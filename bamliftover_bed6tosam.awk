## operation
# 1. revert the lifed BED6 with its source bamfile back to samfile
# 2. refine gaps

# Test
# $ diff <(samtools view <bamfile>) \
#        <(samtools view <bamfile> | mawk -f bamliftover_samtobed6.awk | sort -k4,4n -k5,5nr | \
#          mawk -f bamliftover_bed6tosam.awk -- -b <bamfile> )
BEGIN {
  version = "0.5.1";
  program_name = "bamliftover_bed6tosam";
  FS="\t";
  OFS="\t";
  sampipe = ""; chain_bedpefile = ""; nobed6_file = "";
  split("", samheader); n_samheader = 0; use_zntag = 0;
  cmdargs = ""; out_pg = 0;

  for (i = 1; i < ARGC; i++) {
    if (ARGV[i]=="-b" || ARGV[i] == "--bam"){
      delete ARGV[i];
      sampipe = sprintf("samtools view -h %s", ARGV[++i]);
      cmdargs = cmdargs sprintf("--bam %s ", ARGV[i]);
      delete ARGV[i];
    } else if (ARGV[i]=="-s" || ARGV[i] == "--sam"){
      delete ARGV[i];
      sampipe = sprintf("cat %s", ARGV[++i]);
      cmdargs = cmdargs sprintf("--sam %s ", ARGV[i]);
      delete ARGV[i];
    } else if (ARGV[i]=="-n" || ARGV[i] == "--nobed6"){
      delete ARGV[i];
      nobed6_file = ARGV[++i];
      cmdargs = cmdargs sprintf("--nobed6 %s ", nobed6_file);
      delete ARGV[i];
    } else if (ARGV[i]=="-z" || ARGV[i] == "--zn"){
      use_zntag = 1;
      cmdargs = cmdargs "--zn ";
      delete ARGV[i];
    } else if (ARGV[i]=="-Q" || ARGV[i] == "--sqheader"){
      delete ARGV[i];
      sqheader = ARGV[++i];
      cmdargs = cmdargs sprintf("--sqheader %s ", ARGV[i]);
      delete ARGV[i];
    } else if (ARGV[i]=="-p" || ARGV[i] == "--chainbedpe"){
      delete ARGV[i];
      chain_bedpefile = ARGV[++i];
      cmdargs = cmdargs sprintf("--chain %s ", ARGV[i]);
      delete ARGV[i];
    } else if (ARGV[i]=="-H" || ARGV[i] == "--header"){
      delete ARGV[i];
      n_samheader += 1;
      samheader[n_samheader] = ARGV[++i];
      cmdargs = cmdargs sprintf("--header %s ", ARGV[i]);
      delete ARGV[i];
    } else if (ARGV[i]=="-h" || ARGV[i] == "--help"){
      usage();
    } else if (ARGV[i] ~ /^-./){
      raise("unrecognized option " ARGV[i]);
    } else break;
  }
  if (sampipe == "") raise("Bam/Samfile are required");

  # read SN and SL from header, chrsize, chain_bedpe
  split("", old2new_rname); split("", chrsize); split("", newrname_order); n_newrname = 0; split("", used_rnames);
  if (sqheader != ""){
    while((getline < sqheader) > 0){
      if ($1 != "@SQ") continue;
      sn = substr($2, 4);
      newrname_order[++n_newrname] = sn;
      if (!(sn in chrsize)) chrsize[sn] = substr($3, 4)+0;
    }
    close(sqheader);
  }
  if (chain_bedpefile != ""){
    bedpecmd = sprintf("%s %s", (chain_bedpefile ~ /\.gz$/ ? "gzip -dc" : "cat"), chain_bedpefile);
    while((bedpecmd | getline) > 0){
      if ($1 in old2new_rname && old2new_rname[$1] != $4) raise("RNAMEs between old and new assemblies in the CHAIN.bedpe are not one-to-one: " $0);
      old2new_rname[$1] = $4;
      if (chrsize[$4]+0 < $6+0) chrsize[$4] = $6+0;
    }
    close(bedpecmd);
  } 

  split("", nobed6);
  if (nobed6_file != ""){
    while((getline < nobed6_file) > 0){
      if ($1 in old2new_rname && old2new_rname[$1] != $2) raise("RNAMEs conversion table is conflict with CHAIN.bedpe: " $0);
      nobed6[$1] = $2;
      old2new_rname[$1] = $2;
    }
    close(nobed6_file);
  }
  split("", report); # initialize report
  n_bamrecord = 0; 
}

function usage(){
  print sprintf("mawk -f %s.awk -- -b/s <Bam/Samfile> [options] <sorted_BED6>", program_name) > "/dev/stderr";
  print sprintf("  version:%s\n", version) > "/dev/stderr";
  print "Convert BED6 + original BAM/SAM to SAM for bemliftover" > "/dev/stderr";
  print "Prior to execution, BED6 should be sorted by -k4,4n -k5,5nr" > "/dev/stderr";
  print "Bam/Samfile is required" > "/dev/stderr";
  print > "/dev/stderr";
  print "options" > "/dev/stderr";
  print "  -b/--bam File:        Bamfile used for bamliftover_samtobed6" > "/dev/stderr";  
  print "  -s/--sam File:        Samfile with header used for bamliftover_samtobed6" > "/dev/stderr";  
  print "  -n/--nobed6 File :    RNAMEs conversion table for RNAME-only liftOver. Each line consists of tab-delimited old_RNAME and new_RNAME" > "/dev/stderr";  
  print "  -H/--header File:     Samfile header for @PG. Multiple options are allowed" > "/dev/stderr";  
  print "  -Q/--sqheader File:   Samfile header for @SQ" > "/dev/stderr";  
  print "  -p/--chainbedpe File: CHAIN.bedpe to create a new SAM header" > "/dev/stderr";  
  print "  -z/--zn               Add a ZN tag to the second last optional tag to store the original RNAME" > "/dev/stderr";  
  _exit = 0;
  exit(1);
}
function raise(msg){
  print msg > "/dev/stderr";
  _exit = 1;
  exit(1);
}

function flag_unmapped(flag){ return( (int(flag/4))%2 ) }
function set_unmapped_flag(flag){ return( flag_unmapped(flag) ? flag : flag+4) }

function chop_gap(by, chop_len,  l){
  while(chop_len > 0 && match(qcigar2, "[0-9]+" by)){
    l = substr(qcigar2, RSTART, RLENGTH-1)+0;
    if (l > chop_len){ l -= chop_len; chop_len = 0;  }
    else { chop_len -= l; l = 0; }
    # using lowercase not to reuse CIGAR like 0N
    qcigar2 = substr(qcigar2, 1, RSTART-1) l tolower(by) substr(qcigar2, RSTART+RLENGTH);
  }
  qcigar2 = toupper(qcigar2);
  return(chop_len);
}

# return ReadLength
function cigar2len(cigar,   buf, i, n, len){
  len = 0;
  gsub(/[0-9]+[DNHP]/, "", cigar);
  n = split(cigar, buf, /[MIS=X]/)+0;
  for (i = 1; i < n; i++) len += buf[i];
  return(len);
} 

function create_header(   sn, pg, i){
  if ($1 == "@SQ"){
    if (n_newrname > 0){
      for (i=1; i <= n_newrname; i++){
        sn = newrname_order[i];
        print "@SQ", "SN:" sn, "LN:" chrsize[sn];
        used_rnames[sn] = 1;
      }
      n_newrname = 0;
    }
    sn = substr($2, 4);
    if (sn in old2new_rname) sn = old2new_rname[sn];
    if (!(sn in used_rnames)){
      print "@SQ", "SN:" sn, (sn in chrsize ? "LN:" chrsize[sn] : $3);
      used_rnames[sn] = 1;
    }
  } else if ($1 == "@PG" && !out_pg){
    pg = $0;
    print sprintf("@PG\tID:%s\tPN:%s\tVN:%s\tCL:%s", program_name, program_name, version, cmdargs);
    for (i=1; i <= n_samheader; i++){
      while((getline < samheader[i]) > 0) if ($1 == "@PG") print;
      close(samheader[i]);
    }
    print pg;
    out_pg = 1;
  } else {
    print;
  }
}

function get_samrecord(){
  while((sampipe | getline) > 0){
    if ($1 ~ /^@/) create_header()
    else {
      n_bamrecord++;
      if (use_zntag && $3 !="*"){
        if (NF > 11) $(NF) = "ZN:Z:" $3 OFS $(NF)
        else $(NF) = $(NF) OFS "ZN:Z:" $3 
      } 
      $7 = "*"; $8 = 0; $9 = 0; # remove mate information
      if (flag_unmapped($2)) $3 = "*"
      else if ($3 in nobed6) $3 = nobed6[$3]
      else break;
      print;
    }
  }
}
##################### main loop
{ bedname = $4; blockidx = $5+0; 
  rname = $1;
  if ($7 !~ /,$/) raise("Unexpected CIGAR " $0);
  cigar = substr($7, 1, length($7)-1);
  st = $2+0;
  while($5+0 > 1){
    if ((getline) == 0) raise("Truncated file");
    if (rname == ".") rname = $1; # looking for mapped block
    blockidx -= 1;
    if (bedname != $4 || $5+0 != blockidx) raise("BED6 is not sorted: " $0);
    if ((i = index($7, ",")) == 0) raise("Unexpected CIGAR " $0);

    qcigar1 = substr($7, 1, i-1); qcigar2 = substr($7, i+1);
    gap_len = 0
    n = split(qcigar2, arr, "[DN]");
    for (i=1; i < n; i++) gap_len += arr[i]
    if ($3+0 > 0 && st >= 0){ # both prev. and curr. block are liftovered
      # print (st-$3), gap_len, qcigar2 > "/dev/stderr";
      diff_len = (st-$3) - gap_len;
      if (diff_len > 0){ # gap extend, extend 'N' if exist, otherwise extend 'D'
        if (match(qcigar2, "[0-9]+N")) 
          qcigar2 = substr(qcigar2, 1, RSTART-1) substr(qcigar2, RSTART, RLENGTH-1)+diff_len substr(qcigar2, RSTART+RLENGTH-1)
        else qcigar2 = qcigar2 diff_len "D"
      } else { # gap removal, at first from N, next from D
        diff_len = -diff_len;
        diff_len = chop_gap("N", diff_len);
        # print sprintf("remove %d bases from non-N gap %s", -diff_len, qcigar2) > "/dev/stderr";
        if (diff_len > 0) diff_len = chop_gap("D", diff_len);
        if (diff_len > 0){ # very very rare case
          qcigar2 = qcigar2 cigar; cigar = "";
          diff_len = chop_gap("D", diff_len);
          if (diff_len > 0) diff_len = chop_gap("N", diff_len);
          if (diff_len > 0){ 
            while((sampipe | getline) > 0){
              if ($1 ~ /^@/) continue;
              n_bamrecord += 1;
              if (!flag_unmapped($2)) break;
            }
            if (n_bamrecord == bedname+0) print > "/dev/stderr";
            raise(sprintf("Can not remove %d bases from gap any more", diff_len));
          }
        }
      }
    } else if ($3+0 < 1 && st >= 0){ # curr. block is not liftovered
      $2 = st - gap_len;
    } else if ($3+0 > 0 && st < 0 ){ # prev. block is not liftovered. 
      # this happens only if the last block is not liftovered (due to input BED6 is sorted in reverse)
      # nothing has to do
    }
    cigar = qcigar1 qcigar2 cigar;
    st = $2+0;
  }

  # get SamRecord
  get_samrecord();
  if (n_bamrecord != bedname+0) raise("Conflict between BED6 and Bamfile. Is BED6 sorted?");

  # Output new sam records after checking sequence length == length infering from CIGAR 
  if (cigar2len(cigar) != length($10)) raise(sprintf("Inconsistent CIGAR %s\n%s", cigar, $0));
  $3 = rname; $4 = st+1; $6 = cigar;
  if ($3 == "."){
    $3 = "*";
    $2 = set_unmapped_flag($2);
    report["unmap"] += 1;
  } else if (!($3 in used_rnames)) raise("RNAME: " $3 " is not in header");
  print;
  report["bed6sam"] += 1;
}

END {
  get_samrecord();
  if ((sampipe | getline) > 0) raise("Not all SAM records are processed");
  close(sampipe);
  print sprintf("# %s,  version %s", program_name, version) > "/dev/stderr";
  print sprintf("# SAM records: %d", n_bamrecord) > "/dev/stderr";
  print sprintf("# SAM records from BED6: %d", report["bed6sam"]) > "/dev/stderr";
  print sprintf("# Non-liftovered records: %d", report["unmap"]) > "/dev/stderr"; 

  if (_exit){
    print "Failed to complete" > "/dev/stderr";
    exit(1);
  }
}

