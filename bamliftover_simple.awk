BEGIN {
  version = "0.2.0";
  program_name = "bamliftover_simple";
  FS="\t"; OFS="\t";
  src_type = "bam"; rname_file = ""; use_zntag = 0;
  cmdargs = ""; out_pg = 0;

  for (i = 1; i < ARGC; i++) {
    if (ARGV[i]=="-s" || ARGV[i] == "--sam"){
      delete ARGV[i];
      src_type = "sam";
      cmdargs = cmdargs "--sam ";
    } else if (ARGV[i]=="-n" || ARGV[i] == "--name"){
      delete ARGV[i];
      rname_file = ARGV[++i];
      cmdargs = cmdargs sprintf("--name %s ", rname_file);
      delete ARGV[i];
    } else if (ARGV[i]=="-z" || ARGV[i] == "--zn"){
      delete ARGV[i];
      use_zntag = 1;
      cmdargs = cmdargs "--zn ";
    } else if (ARGV[i]=="-h" || ARGV[i] == "--help"){
      usage();
    } else if (ARGV[i] ~ /^-./){
      raise("unrecognized option " ARGV[i]);
    } else break;
  }
  if (rname_file=="") raise("No RNAME table");
  if (ARGV[i]=="") raise("No souce file");
  sampipe = sprintf("%s %s 2> /dev/null", (src_type=="sam" ? "cat" : "samtools view -h"), ARGV[i]);
  delete ARGV[i];

  # read RNAME table
  split("", old2new_rname); split("", used_rnames);
  while((getline < rname_file) > 0) old2new_rname[$1] = $2;
  close(rname_file);

  ## main loop
  n_bamrecord = 0; 
  while((sampipe | getline) > 0){
    if ($1 ~ /^@/) create_header()
    else {
      n_bamrecord++;
      if (use_zntag && $3 !="*") $(NF) = "ZN:Z:" $3 OFS $(NF);
      $7 = "*"; $8 = 0; $9 = 0; # remove mate information
      if ((int($2/4))%2) $3 = "*"
      else if ($3 in old2new_rname) $3 = old2new_rname[$3]
      else raise("Unknown name: " $3);
      print;
    }
  }
  close(sampipe);
  print sprintf("# %s,  version %s", program_name, version) > "/dev/stderr";
  print sprintf("# SAM records: %d", n_bamrecord) > "/dev/stderr";
  exit;
}

function create_header(   sn){
  if ($1 == "@SQ"){
    sn = substr($2, 4);
    if (!(sn in old2new_rname)) raise(sprintf("Unknown name: %s in header", sn));
    sn = old2new_rname[sn];
    if (sn in used_rnames) return;
    $2 = "SN:" sn;
    used_rnames[sn] = 1;
  } else if ($1 == "@PG" && !out_pg){
    print sprintf("@PG\tID:%s\tPN:%s\tVN:%s\tCL:%s", program_name, program_name, version, cmdargs);
    out_pg = 1;
  }
  print;
}

function usage(){
  print sprintf("mawk -f %s.awk -- -n <RNAMEtable> [options] <SAMfile|->", program_name) > "/dev/stderr";
  print sprintf("  version:%s\n", version) > "/dev/stderr";
  print > "/dev/stderr";
  print "options" > "/dev/stderr";
  print "  -s/--sam File:        Source is SAM with header [default BAM]" > "/dev/stderr";  
  print "  -n/--name File :      Required. RNAMEs conversion table for RNAME-only liftOver. Each line consists of tab-delimited old_RNAME and new_RNAME." > "/dev/stderr";  
  _exit = 0;
  exit(1);
}
function raise(msg){
  print msg > "/dev/stderr";
  _exit = 1;
  exit(1);
}

END {
  if (_exit){
    print "Failed to complete" > "/dev/stderr";
    exit(1);
  }
}

