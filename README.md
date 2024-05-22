## Introduction
Bamliftover lifts over SAM/BAM alignments from custom assemblies, such as subspecies genomes, to the reference genome.
It is designed for **haplotype-resolved NGS analysis**, where custom assemblies are built by adding SNPs/indels to the reference genome.
Compared to the UCSC liftOver tool, `bamliftover` is faster and updates CIGAR alignment strings in SAM/BAM format so that SNP/indel sites and spliced alignments can be visualized in downstream analyses. 

## Requirements

* [samtools](https://www.htslib.org/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [mawk](https://invisible-island.net/mawk/) 

## Usage
```
$ ruby bamliftover.rb
ruby bamliftover.rb [options] <chain.bedpe> <bamfile|->
ruby bamliftover.rb [options] -n <RNAMEs_table> <bamfile|->
  Version: 0.5.1
  chain.bedpe is created by chain2bedpe.awk
  Output is in the SAM format with header

        --sam                        The format of the input is SAM with a header
    -h, --header File                The file including a SAM header for new @SQ
    -n, --noliftcoord File           RNAMEs conversion table for RNAME-only liftOver. Each line consists of tab-delimited old_RNAME and new_RNAME
    -z, --zn                         Add ZN tag at second last optional tag to store the original RNAMEs
    -s, --soft                       Use a soft-clipping operation instead of an insertion operation for the end of reads
    -@, --threads NUM                Number of samtools threads [4]
    -m, --memory SIZE                Memory size for sort command. If SIZE is set to 0, the default memory size of the sort command is used. [1G]
        --dry                        Dry run
```

## Limitations
* Bamliftover only supports the case where the two genomic coordinates are linearly one-to-one and does not support genomic rearrangements such as translocations and inversions
* Bamliftover is designed for single-end mapping and removes mate information such as RNEXT, PNEXT, and TLEN 
* Bamliftover internally sorts whole data twice. For large data like Hi-C mapping, splitting the BAM file is recommended
