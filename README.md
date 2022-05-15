# rs-liftover
Simple sequence mapping liftover tool designed to be general and work with even noisy data.

## Introduction

I had trouble working with liftover tools like leviosam and CrossMap, especially in the presence of SVs, so I hacked together a very simple liftover tool. `rs-liftover` lifts over the mapping to the other reference in a very straightforward (i.e. stupid) manner.


Given a

1. BAM file with mapped sequences to liftover 
2. UCSC chain file between two references

Gaps in the chain file are turned into deletions and insertions without much thought and the reads are extended appropriately. Works well if reads map decently onto the initial reference.

Tested on nanopore data with standard CIGAR strings. 
I would reccomend using other liftover tools like leviosam or CrossMap unless they don't work with your data for some reason.

### Requirements 

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.
### Install

```
git clone https://github.com/bluenote-1577/rs-liftover
cd rs-liftover
cargo build --release
./target/release/liftover -h
```

### Running liftover

`./target/release/liftover -b my_output.bam chainfile.chain bamfile.bam lift-destination.fa lift-source.fa` outputs an unsorted bam file called `my_output.bam`. 

rs-liftover only works with chain files with a single long chain right now. Will fix in a bit...