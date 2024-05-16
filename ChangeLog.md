# Change log

## bamliftover_samtobed6.awk

### v0.1.1, 2015.03.16
* fixed a bug. The samheader must be closed before processing the first SAM alignment. It is needed for piping to bed6tosam
* fixed a bug with duplicate @PG 

### v0.1.2, 2015.04.13
* change the output header to include only @PG for this script

### v0.2.0, 2015.04.16
* add -n/--noconv option for simple liftOver

### v0.3.0, 2024.05.15
* public release
* clarify usage

## bamliftover_bed6liftover.awk

### v0.1.1, 2015.03.16
* add -H/--header option
* fixed a bug where -c/--clipop does not work

### v0.2.0, 2017.10.17
* fixed a bug. Clipling operations in CIGAR should be at the start/end of the CIGAR string. '20M10S1I' should be '20M1I10S'. This is important for bwa mem alignments

### v0.3.0, 2024.05.15
* public release
* clarify usage

## bamliftover_bed6tosam.awk

### v0.1.1, 2015.03.16
* fixed a bug where the SAM header output is empty before retrieving the first SAM record
* accept multiple SAM headers

### v0.1.2, 2015.03.16
* fixed a bug for a very very rare case in gap refinement

### v0.2.0, 2015.04.13
* check consistency between RNAMEs and SAM header
* add -p/--chainbedpe option to create a new SAM header
* add -Q/--sqheader option to inherit the sort order
* change RNAMEs to "*" according to the SAM spec if not mapped
* SAM header is imported from options, -s/--sam File or -b/--bam File
* disable -c/--chrsize option
* output unmapped reads

### v0.3.0, 2015.04.16
* adapt to bamliftover_samtobed6.awk v0.2.0
* add -n/--nobed6 option for simple liftOver
* add a ZN tag to mapped alignments at $(NF-1) to store the original RNAME if -z/--zn is specified
* fixed a bug where the last unmapped records are not processed
* convert the RNAMEs of unmapped reads '*' in any case

### v0.4.0 2023.09.28
* fixed a bug in -z/--zn option. If no optional tags, add the ZN tag to the last optional field

### v0.5.0, 2024.05.15
* public release
* clarify usage
* remove mate information to avoid RNEXT and PNEXT not described in the SAM header.

## bamliftover_simple.awk

### v0.2.0, 2015.03.16
* public release
* clarify usage
* remove mate information to avoid RNEXT and PNEXT not described in the SAM header.

## bamliftover.rb

### v0.2.0, 2015.03.16
* change the required argument, chromsize, to an optional argument
* if the header in the BAM file contains liftovered RNAMAs, no chromsize file is needed and the sort order between the input and output BAM files is consistent.

### v0.2.1, 2015.03.16
* SAM input is allowed

### v0.2.2, 2015.04.13
* fixed a bug where SAM input (change in v0.2.1) does not work
* adapt to bamliftover_samtobed6.awk v0.2.0 and remove the chromsize option 
* disable the -c/--chromsize option and add the -h/--header option

### v0.3.0, 2015.04.16
* adapt to bamliftover_bed6tosam_v0.3.0 and bamlifover_santobed6 v0.2.0
* add a ZN tag to mapped alignments to the second last optional tag to store the original RNAME if the option -z/--zn is specified
* add the -n/--nolift option
* sort command uses -s 1G option
* change awk library names

### v0.3.1, 2015.04.20
* add the -t/--tmpfs option

### v0.3.2, 2015.05.22
* use a temporary BAM as an uncompressed BAM

### v0.3.3 2016.06.06
* add -c/--compress option

### v0.4.0 2020.04.02
* add no chain.bedpe mode

### v0.5.0, 2024.05.16
* public release
* clarify usage
* remove the -t/--tmpfs option
* remove the -c/--compress option and use a temporary BAM as an uncompressed BAM by default
* remove dependency on tmpfile.rb
* add -@/--threads option
