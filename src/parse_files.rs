use bio::io::fasta;
use bio_types::genome::AbstractInterval;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::{Cigar as hts_Cigar, CigarString, Record};
use rust_htslib::bam::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::prelude::*;
use std::io::prelude::*;
use std::io::BufReader;
use std::str;
use std::sync::Mutex;

pub struct Params {
    pub output_name: String,
}

pub struct ChainData {
    pub score: usize,
    pub start_ref: i64,
    pub end_ref: i64,
    pub start_query: i64,
    pub end_query: i64,
    pub primary_name: String,
    pub secondary_name: String,
    pub blocks: Vec<(i64, i64, i64)>,
}

pub fn read_chain_file(chain_file: &str) -> ChainData {
    let file = File::open(chain_file).unwrap();
    let reader = BufReader::new(file);
    let mut contig_blocks = vec![];
    let mut chain_score = 0;
    let mut chain_start_ref = 0;
    let mut chain_end_ref = 0;
    let mut chain_start_query = 0;
    let mut chain_end_query = 0;
    let mut secondary_chrom_name = String::new();
    let mut primary_chrom_name = String::new();
    let mut run_ref: i64 = 0;
    let mut run_query: i64 = 0;
    for line in reader.lines() {
        let line = line.unwrap();
        let spl: Vec<&str> = line.split_whitespace().collect();
        if spl.len() == 0 {
            continue;
        }
        if spl[0] == "chain" {
            if secondary_chrom_name != "" {
                break;
            }
            secondary_chrom_name = spl[2].to_owned();
            primary_chrom_name = spl[7].to_owned();
            chain_end_query = spl[11].parse().unwrap();
            chain_start_query = spl[10].parse().unwrap();
            chain_score = spl[1].parse().unwrap();
            chain_start_ref = spl[5].parse().unwrap();
            chain_end_ref = spl[6].parse().unwrap();
            run_ref = chain_start_ref;
            run_query = chain_start_query;
        } else {
            let block_len: i64 = spl[0].parse().unwrap();
            contig_blocks.push((run_ref + block_len, run_query + block_len, block_len));
            if spl.len() > 1 {
                let ref_gap = spl[1].parse::<i64>().unwrap();
                let q_gap = spl[2].parse::<i64>().unwrap();
                run_ref += ref_gap + block_len;
                run_query += q_gap + block_len;
            }
        }
    }

    let cinfo = ChainData {
        score: chain_score,
        start_ref: chain_start_ref,
        end_ref: chain_end_ref,
        start_query: chain_start_query,
        end_query: chain_end_query,
        blocks: contig_blocks,
        primary_name: primary_chrom_name,
        secondary_name: secondary_chrom_name,
    };
    cinfo
}

pub fn realign(
    bam_file: &str,
    chain_file: &str,
    primary_contigs: &str,
    secondary_contigs: &str,
    params: Params,
) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(10)
        .build_global()
        .unwrap();

    let reader = fasta::Reader::from_file(secondary_contigs).unwrap();
    let mut secondary_seqs_dict = HashMap::new();
    for result in reader.records() {
        let rec = result.unwrap();
        let ref_seq = rec.seq().to_vec();
        secondary_seqs_dict.insert(rec.id().to_owned(), ref_seq);
    }

    let reader = fasta::Reader::from_file(primary_contigs).unwrap();
    let mut primary_seqs_dict = HashMap::new();
    for result in reader.records() {
        let rec = result.unwrap();
        let ref_seq = rec.seq().to_vec();
        primary_seqs_dict.insert(rec.id().to_owned(), ref_seq);
    }

    let cinfo = read_chain_file(chain_file);

    for contig_block in cinfo.blocks.iter() {
        let s1 = &secondary_seqs_dict[&cinfo.secondary_name]
            [(contig_block.0 - contig_block.2) as usize..(contig_block.0) as usize];
        let s2 = &primary_seqs_dict[&cinfo.primary_name]
            [(contig_block.1 - contig_block.2) as usize..(contig_block.1) as usize];
        //        println!(
        //            "{}",
        //            str::from_utf8(&s1[0..usize::min(s1.len(), 20)]).unwrap()
        //        );
        //        println!(
        //            "{}",
        //            str::from_utf8(&s2[0..usize::min(s2.len(), 20)]).unwrap()
        //        );
    }

    let mut bam = bam::Reader::from_path(bam_file).unwrap();
    let mut header = bam::Header::from_template(bam.header());
    let mut primary_head_rec = HeaderRecord::new(b"SQ");
    primary_head_rec.push_tag(b"SN", &cinfo.primary_name);
    primary_head_rec.push_tag(b"LN", &primary_seqs_dict[&cinfo.primary_name].len());
    //    header.push_record(&primary_head_rec);

    let mut writer =
        bam::Writer::from_path(&params.output_name, &header, bam::Format::Bam).unwrap();

    for rec in bam.records() {
        let record = rec.unwrap();
        if record.mapq() > 20 && record.tid() >= 0 {
            let contig = record.contig().to_owned();
            if contig != cinfo.secondary_name {
                continue;
            }
            //            dbg!(&contig, &cinfo.secondary_name);
            let aligned_pairs = record.aligned_pairs();
            let seq = record.seq();
            let cigar = record.cigar();
            let curr_ref_pos = record.pos();
            let curr_read_pos = cigar.leading_softclips();
            let mut new_cigar: Vec<hts_Cigar> = vec![];
            new_cigar.push(hts_Cigar::SoftClip(cigar.leading_softclips() as u32));
            let match_cb = cinfo.blocks.binary_search_by(|x| x.0.cmp(&(curr_ref_pos)));
            let mut gap = false;
            let mut gap_len = 0;
            let mut len_til_gap = 0;
            let mut len_til_block = 0;
            let mut new_start_pos = 0;
            let mut block_ind = match match_cb {
                Ok(i) => i,
                Err(i) => i,
            };
            if block_ind == cinfo.blocks.len() {
                continue;
            } else {
                let gap_len = cinfo.blocks[block_ind].0 as i64
                    - cinfo.blocks[block_ind].2 as i64
                    - curr_ref_pos;
                if gap_len > 0 {
                    gap = true;
                    len_til_block = gap_len;
                } else {
                    len_til_gap = cinfo.blocks[block_ind].0 - curr_ref_pos;
                }
            }

            let new_start_pos;
            if gap{
                new_start_pos = cinfo.blocks[block_ind].1 - cinfo.blocks[block_ind].2;
            }
            else{
                new_start_pos = cinfo.blocks[block_ind].1 - len_til_gap;
            }

            //            if gap {
            //                continue;
            //            }

            for c in cigar.iter() {
                if matches!(c, hts_Cigar::Match(u32))
                    || matches!(c, hts_Cigar::Ins(u32))
                    || matches!(c, hts_Cigar::Del(u32))
                {
                    let mut len = c.len() as i64;
                    loop {
                        if matches!(c, hts_Cigar::Ins(u32)) {
                            new_cigar.push(hts_Cigar::Ins(c.len() as u32));
                            break;
                        }
                        if gap {
                            if len > len_til_block {
                                if matches!(c, hts_Cigar::Match(u32)) {
                                    new_cigar.push(hts_Cigar::Ins(len_til_block as u32));
                                }
                                len -= len_til_block;
                                gap = false;
                                len_til_block = 0;
                                len_til_gap = cinfo.blocks[block_ind].2;
                                if block_ind == 0{
                                    new_cigar.push(hts_Cigar::SoftClip(len_til_block as u32));
                                }
                                else{
                                    let primary_gap_len = cinfo.blocks[block_ind].1 as i64
                                        - cinfo.blocks[block_ind].2 as i64
                                        - cinfo.blocks[block_ind-1].1;
                                    if primary_gap_len > 0 {
                                        new_cigar.push(hts_Cigar::Del(primary_gap_len as u32));
                                    }
                                }
                            } else {
                                if matches!(c, hts_Cigar::Match(u32)){
                                    if block_ind != 0{
                                        new_cigar.push(hts_Cigar::Ins(len as u32));
                                    }
                                    else{
                                        new_cigar.push(hts_Cigar::SoftClip(len as u32));
                                    }
                                }
                                len_til_block -= len;
                                break;
                            }
                        } else {
                            if len > len_til_gap {
                                if matches!(c, hts_Cigar::Match(u32)) {
                                    new_cigar.push(hts_Cigar::Match(len_til_gap as u32));
                                }
                                if matches!(c, hts_Cigar::Del(u32)) {
                                    new_cigar.push(hts_Cigar::Del(len_til_gap as u32));
                                }
                                len -= len_til_gap;
                                len_til_gap = 0;
                                block_ind += 1;
                                len_til_block = cinfo.blocks[block_ind].2;
                                let gap_len = cinfo.blocks[block_ind].0 as i64
                                    - cinfo.blocks[block_ind].2 as i64
                                    - cinfo.blocks[block_ind - 1].0;
                                if gap_len > 0 {
                                    gap = true;
                                    len_til_block = gap_len;
                                } else {
                                    //                                    block_ind += 1;
                                    len_til_gap = cinfo.blocks[block_ind].2;
                                    let primary_gap_len = cinfo.blocks[block_ind].1 as i64
                                        - cinfo.blocks[block_ind].2 as i64
                                        - cinfo.blocks[block_ind - 1].1;
                                    new_cigar.push(hts_Cigar::Del(primary_gap_len as u32));
                                }
                            } else {
                                if matches!(c, hts_Cigar::Match(u32)) {
                                    new_cigar.push(hts_Cigar::Match(len as u32));
                                }
                                if matches!(c, hts_Cigar::Del(u32)) {
                                    new_cigar.push(hts_Cigar::Del(len as u32));
                                }
                                len_til_gap -= len;
                                break;
                            }
                        }
                    }
                }
            }
            //            dbg!(&new_cigar, new_start_pos, curr_ref_pos);
            let s1 = &secondary_seqs_dict[&cinfo.secondary_name]
                [curr_ref_pos as usize..curr_ref_pos as usize + 20];
            let s2 = &primary_seqs_dict[&cinfo.primary_name]
                [new_start_pos as usize..new_start_pos as usize + 20];
            //            println!("{}", str::from_utf8(&s1).unwrap());
            //            println!("{}", str::from_utf8(&s2).unwrap());

            let mut seq_len_from_cig = 0;
            for cigar in new_cigar.iter() {
                if matches!(cigar, hts_Cigar::Match(u32)) {
                    seq_len_from_cig += cigar.len();
                } else if matches!(cigar, hts_Cigar::Ins(u32)) {
                    seq_len_from_cig += cigar.len();
                } else if matches!(cigar, hts_Cigar::SoftClip(u32)) {
                    seq_len_from_cig += cigar.len();
                }
            }
            //            dbg!(seq_len_from_cig, seq.len());
            new_cigar.push(hts_Cigar::SoftClip(cigar.trailing_softclips() as u32));
            let mut new_rec = Record::new();
            let quals = record.qual();
            let hts_cigar_view = CigarString(new_cigar);
            let hts_cigar = Some(hts_cigar_view);
            if record.is_reverse() {
                new_rec.set_reverse();
            }
            new_rec.set(record.qname(), hts_cigar.as_ref(), &seq.as_bytes(), quals);
            let hv = HeaderView::from_header(&header);
            let tid = hv.tid(cinfo.primary_name.as_bytes()).unwrap();
            new_rec.set_tid(tid as i32);
            new_rec.set_pos(new_start_pos);
            new_rec.set_mapq(60);
            writer.write(&new_rec).unwrap();
        }
    }
}
