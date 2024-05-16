# samtools view -h [bam_file] | awk -f bamliftover_samtobed6.awk 

## operation
# 1. convert SAM alignments to split BED6 format like bedtools bamtobed -splitD
# 2. bed6tobe12.awk compatible
# 3. CIAGR output
#
# Note about bedtools bamtobed -splitD
# 1. No CIGAR output
# 2. CIAGR like xDxN makes zero length BED record
#
# Test
# $ diff -U 0 <(bedtools bamtobed -splitD -i <bamfile> | mawk '($2!=$3){print}' | cut -f1-4,6 | head -n 100) \
#   <(samtools view -F4 <bamfile> | mawk -f bamliftover_samtobed6.awk -- -q | cut -f1-4,6 | head -n 100)
# $ diff -U 0 <(samtools view -F4 <bamfile> | head -n 100 | cut -f6) \
#   <(samtools view -F4 <bamfile> | head -n 100 | mawk -f bamliftover_samtobed6.awk | \
#    sort -k4,4n -k5,5nr -k2,2n | mawk -f bed6tobed12.awk -- -c | mawk '{gsub(",", "", $9); print $9}')

BEGIN {
  version = "0.3.0";
  program_name = "bamliftover_samtobed6";
  FS="\t";
  OFS="\t";
  use_qname = 0; samheader = ""; show_seq = 0;
  noconv_file = "";
  cmdargs = ""; 

  for (i = 1; i < ARGC; i++) {
    if (ARGV[i]=="-n" || ARGV[i] == "--noconv"){
      delete ARGV[i];
      noconv_file = ARGV[++i];
      cmdargs = cmdargs sprintf("--noconv %s ", noconv_file);
      delete ARGV[i];      
    } else if (ARGV[i]=="-H" || ARGV[i] == "--header"){
      delete ARGV[i];
      samheader = ARGV[++i];
      cmdargs = cmdargs sprintf("--header %s ", samheader);
      delete ARGV[i];      
    } else if (ARGV[i]=="-q" || ARGV[i] == "--qname"){
      use_qname = 1;
      cmdargs = cmdargs "--qname ";
      delete ARGV[i];
    } else if (ARGV[i]=="-s" || ARGV[i] == "--seq"){
      show_seq = 1;
      cmdargs = cmdargs "--seq ";
      delete ARGV[i];
    } else if (ARGV[i]=="-h" || ARGV[i] == "--help"){
      usage();
    } else if (ARGV[i] ~ /^-./){
      raise("unrecognized option " ARGV[i]);
    } else break;
  }

  split("", noconv);
  if (noconv_file != ""){
    while((getline < noconv_file) > 0) noconv[$1] = $2;
    close(noconv_file);
  }
  split("", report); # initialize report
  n_records = 0; 
  if (samheader != ""){
    print sprintf("@PG\tID:%s\tPN:%s\tVN:%s\tCL:%s", program_name, program_name, version, cmdargs) > samheader;
    close(samheader);
  }
}

function usage(){
  print sprintf("mawk -f %s.awk -- [options] <samfile_with_header>", program_name) > "/dev/stderr";
  print sprintf("  version:%s\n", version) > "/dev/stderr";
  print "Convert SAM alignments as BED6 format with CIGAR in the 7th column like bedtools bamtobed -splitD" > "/dev/stderr";
  print "This BED6 can be reverted with bamliftover_bed6tosam.awk" > "/dev/stderr";
  print > "/dev/stderr";
  print "options" > "/dev/stderr";
  print "  -n/--noconv File :       List of RNAMEs to ignore for conversion. Same as the RNAME conversion table for bamliftover_simple.awk" > "/dev/stderr";  
  print "  -h/--header File :       File to save SAM header" > "/dev/stderr";  
  print "  -q/--qname :             Show the QNAME in the name column" > "/dev/stderr";  
  print "  -s/--seq :               Show the sequence in the 8th column" > "/dev/stderr";  
  _exit = 0;
  exit(1);
}
function raise(msg){
  print msg > "/dev/stderr";
  _exit = 1;
  exit(1);
}

function flag_unmapped(flag){ return( (int(flag/4))%2 ) }
function flag_reverse(flag){ return( (int(flag/16))%2 ) }

/^@/{ next; }
##################### main loop
{ n_records += 1; }
(flag_unmapped($2) || $3 in noconv){ 
  report["no_conv"]++;
  next;
}
{ cigar = $6; strand = (flag_reverse($2) ? "-" : "+");
  if (cigar ~ /^[0-9]+M$/){ # simple case. 
    st = $4-1;
    en = st + substr(cigar, 1, length(cigar)-1);
    if (show_seq) cigar = cigar "\t" $10;
    print $3, st, en, (use_qname ? $1 : n_records), 1, strand, cigar ",";
    report["n_bed6"] += 1;
    next;
  }
  idx = 1; st = $4-1; # 0-based 
  while(1){
    if (match(cigar, /([0-9]+[DN])+/)){
      if (RSTART == 1) raise("Unexpected CIGAR " $0);
      cigar1 = substr(cigar, 1, RSTART-1);
      cigar2 = substr(cigar, RSTART, RLENGTH);
      cigar_out = cigar1 "," cigar2;
      cigar = substr(cigar, RSTART+RLENGTH);
    } else {
      cigar1 = cigar; cigar2 = ""; 
      cigar_out = cigar1 ",";
    }
    gsub("[0-9]+[^0-9M]", "", cigar1);
    n = split(cigar1, arr, "M");
    en = st;
    for (i = 1; i < n; i++) en += arr[i];
    if (show_seq) cigar_out = cigar_out "\t" $10;
    print $3, st, en, (use_qname ? $1 : n_records), idx, strand, cigar_out;
    if (cigar2 == "") break;
    n = split(cigar2, arr, "[DN]");
    for (i=1; i < n; i++) en += arr[i];
    st = en; idx += 1;
  }
  report["n_bed6"] += idx;
}
END {
  print sprintf("# %s,  version %s", program_name, version) > "/dev/stderr";
  print sprintf("# SAM records: %d", n_records) > "/dev/stderr";
  print sprintf("# Converted SAM records: %d", n_records - report["no_conv"]) > "/dev/stderr";
  print sprintf("# BED6 records: %d", report["n_bed6"]) > "/dev/stderr";
  
  if (_exit){
    print "Failed to complete" > "/dev/stderr";
    exit(1);
  }
}

