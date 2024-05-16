BEGIN{ 
  version = "0.3.0";
  program_name = "bamliftover_bed6liftover";
  default_endclip_cigar_op = "I";

  endclip_cigar_op = default_endclip_cigar_op;
  samheader = "";
  cmdargs = "";
  
  OFS="\t";
  for (i = 1; i < ARGC; i++) {
    if (ARGV[i]=="-c" || ARGV[i] == "--clipop"){
      delete ARGV[i];
      endclip_cigar_op = ARGV[++i];
      cmdargs = cmdargs sprintf("--clipop %s ", ARGV[i]);
      delete ARGV[i];
    } else if (ARGV[i]=="-H" || ARGV[i] == "--header"){
      delete ARGV[i];
      samheader = ARGV[++i];
      cmdargs = cmdargs sprintf("--header %s ", ARGV[i]);
      delete ARGV[i];
    } else if (ARGV[i]=="-h" || ARGV[i] == "--help"){
      usage();
    } else if (ARGV[i] ~ /^-./){
      raise("unrecognized option " ARGV[i]);
    } else break;
  }
  if (samheader != ""){
    print sprintf("@PG\tID:%s\tPN:%s\tVN:%s\tCL:%s", program_name, program_name, version, cmdargs) > samheader;
    close(samheader);
  }

  name = ""; score = ""; chrom = "";
  split("", stat);
}
  
function usage(){
  print "sort -k1,1 -k2,2n <BED6_with_CIGAR> | \\" > "/dev/stderr";
  print "bedtools intersect -wao -sorted -a stdin -b <CHAIN.bedpe> | \\" > "/dev/stderr";
  print sprintf("mawk -f %s", program_name) > "/dev/stderr";
  print sprintf("  version:%s\n", version) > "/dev/stderr";
  print "Note: Assume BED6_with_CIGAR derived from bamliftover_samtobed6.awk" > "/dev/stderr";
  print "      CHAIN.bedpe is created by chain2bedpe.awk" > "/dev/stderr";
  print "      The output is BED6 format with CIGAR in the 7th column" > "/dev/stderr";
  print > "/dev/stderr";
  print "options" > "/dev/stderr";
  print sprintf("  -c/--clipop CHAR :          Single letter CIGAR operation character for clipping ends [%s]", default_endclip_cigar_op) > "/dev/stderr";
  print         "  -H/--header File :          File to write @PG as SAM header" > "/dev/stderr";
  _exit = 0;
  exit(1);
}
function raise(msg){
  print msg > "/dev/stderr";
  _exit = 1;
  exit(1);
}

function output(   x){
  if (chrom == ".") stat["unmap"] += 1;
  ## print cliplen_begin, cigar, qcigar1, cliplen_end, qcigar2; # for debug
  gsub("[0-9]+M", "", qcigar1); # qcigar1 should not contain 'N' and 'D'
  cigar = cigar qcigar1;
  if (cliplen_begin > 0){
    if (match(cigar, "^([0-9]+[SH])+")){
      x = substr(cigar, 1, RLENGTH);
      cigar = substr(cigar, RLENGTH+1);
    } else x = "";
    cigar = x cliplen_begin endclip_cigar_op cigar;
  }
  if (cliplen_end > 0){
    if (match(cigar, "([0-9]+[SH])+$")){
      x = substr(cigar, RSTART);
      cigar = substr(cigar, 1, RSTART-1);
    } else x = "";
    cigar = cigar cliplen_end endclip_cigar_op x;
  }
  print chrom, st, en, name, score, strand, cigar "," qcigar2;
  chrom = ""; 
  stat["total"] += 1;
}

# take 'len' bases CIGAR string and store it to the global value '_frontside_cigar'
# return the back side CIGAR
# if unset keepall, CIGAR operations moving on the reference side are removed from '_fronside_cigar'
function split_cigar_by_reflen(cigar, len, keepall,   n, op){
  _frontside_cigar = "";
  while((len > 0 && cigar !="") || cigar ~ /^[0-9]+[ISHP]/){
    match(cigar, /^[0-9]+/);
    n = substr(cigar, 1, RLENGTH)+0;
    op = substr(cigar, RLENGTH+1, 1);
    cigar = substr(cigar, RLENGTH+2);
    if (op ~ /[DMXN=]/){
      if (len < n){
        cigar = (n-len) op cigar;
        n = len;
      }
      len -= n;
      if (keepall) _frontside_cigar = _frontside_cigar n op;
    } else _frontside_cigar = _frontside_cigar n op; 
  }
  return(cigar);
}
##################### main loop
## ex. bedtools intersect -a BED6 -b bedpe
# 1       2       3       4       5       6       7       -10     -9      -8      -7      -6      -5      -4      -3      -2      -1      0
# chr10j  5950054 5950153 860     1       +       99M     chr10j  5947717 5950136 chr10   5949466 5951885 10      0       +       +       82
# chr10j  5950054 5950153 860     1       +       99M     chr10j  5950136 5950335 chr10   5951886 5952085 10      0       +       +       17
# chr10j  5985290 5985389 458     1       +       99M     chr10j  5984717 5985530 chr10   5986590 5987403 10      0       +       +       99
#
# query 860 
#   old chr10j  |=== 82 bp ===| |=== 17 bp ===|
#   new chr10   |=== 82 bp ===|-|=== 17 bp ===|
#                              1bp insertion
(name != $4 || score != $5){ # new block
  if (name != "") output();
  qchr = $1; qstart = $2+0; qend = $3+0;
  name = $4; score = $5+0; strand = $6; 
  if ((i = index($7, ",")) == 0) raise("Unexpected CIGAR " $0);
  qcigar1 = substr($7, 1, i-1); qcigar2 = substr($7, i+1);
}
($(NF-7) == "."){ # unmap => CIGAR covered by read (qcigar1) is converted by end-clip CIAGR
  # print "Unmap: " $0 > "/dev/stderr";
  gsub("M", endclip_cigar_op, qcigar1);  # qcigar1 should not contain 'N' and 'D'
  chrom = "."; cigar = qcigar1; qcigar1 = "";
  st = -1; en = 0; 
  cliplen_begin = 0; cliplen_end = 0;
  next; # unmap means no additional matching block
}
{
  if (chrom == ""){ # initial overlap
    chrom = $(NF-7); cigar = ""; 
    st = $(NF-6)+0; 
    cliplen_begin = $(NF-9) - $2; 
    if (cliplen_begin < 0){ # completely covered by chain block
      st -= cliplen_begin;  
      cliplen_begin = 0;
    } else {
      qcigar1 = split_cigar_by_reflen(qcigar1, cliplen_begin, 0);
      qcigar1 = _frontside_cigar qcigar1;
    }
    en = st;
  } else {
    qgapsize = $(NF-9) - prev_qend; tgapsize = $(NF-6) - prev_tend;
    if (qgapsize > 0){
      qcigar1 = split_cigar_by_reflen(qcigar1, qgapsize, 0);
      cigar = cigar qgapsize "I" _frontside_cigar;
    }
    if (tgapsize > 0){
      cigar = cigar tgapsize "D";
      en += tgapsize;
    }
  }
  cliplen_end = $3 - $(NF-8);
  if (cliplen_end < 0){
    cliplen_end = 0; 
    cigar = cigar qcigar1;
    qcigar1 = "";
  } else {
    qcigar1 = split_cigar_by_reflen(qcigar1, $(NF)+0, 1);
    cigar = cigar _frontside_cigar;
  }
  en += $(NF);
  prev_qend = $(NF-8)+0; prev_tend = $(NF-5)+0;
}

END {
  if (chrom != "") output();
  print sprintf("# %s, version %s", program_name, version) > "/dev/stderr";
  print sprintf("# Total records: %d", stat["total"]) > "/dev/stderr"; 
  print sprintf("# Non-liftovered records: %d", stat["unmap"]) > "/dev/stderr"; 
		
  if (_exit){
    print "Failed to complete" > "/dev/stderr";
    exit(1);
  }
}

