pub mod parse_files;
use crate::parse_files::*;
use clap::{App, AppSettings, Arg, Command, SubCommand};

fn main() {
    let matches = Command::new("rs-liftover")
        .setting(AppSettings::ArgRequiredElseHelp)
        .version("0.1")
        .about("A simple tool for lifting over reads given a chain file with heuristics around indels.")
        .arg(
            Arg::new("chain file")
                .short('c')
                .help("Input UCSC chain file")
                .takes_value(true)
                .multiple_values(true)
                .required(true),
        )
        .arg(
            Arg::new("primary contig")
                .help("Contig being lifted TO")
                .short('p')
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::new("secondary contig")
                .help("Contig being lifted FROM")
                .short('s')
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::new("bam file")
                .help("BAM file to lift")
                .short('b')
                .takes_value(true)
                .required(true),

        )
        .arg(
            Arg::new("output name")
                .help("Output bam file name") 
                .short('o')
                .takes_value(true)
                .required(false)
        )
        .get_matches();

    let bam_file = matches.value_of("bam file").unwrap();
    let chain_files : Vec<&str>  = matches.values_of("chain file").unwrap().collect();
    let primary_file = matches.value_of("primary contig").unwrap();
    let secondary_file = matches.value_of("secondary contig").unwrap();
    let output_name = matches.value_of("output name").unwrap_or("lift.bam");

    let params = Params{output_name: output_name.to_owned()};
    realign(bam_file, &chain_files, primary_file, secondary_file, params);
}
