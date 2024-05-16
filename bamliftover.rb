#!/usr/bin/env ruby

require 'optparse'
require 'tmpdir'
require 'fileutils'

class Tmpfile
  @index = 0
  @tmpdir = nil
  class << self
    def create(suffix=''); "#{tmpdir}/#{@index += 1}_#{suffix}"; end
    def tmpdir; @tmpdir ||= Dir.mktmpdir; end
    def clear;
      FileUtils.remove_entry_secure(@tmpdir) if @tmpdir
      @tmpdir = nil
    end
  end
end

at_exit do
  Tmpfile.clear
end

Version = '0.5.0'
params = {
  :softclip => false,
  :sam => false,
  :nolift => nil,
  :zn => false,
  :header => nil,
  :sort_size => '1G',
  :dry => false,
  :threads => 4
}

opts = OptionParser.new(<<HEADER){ |opt|
ruby #{File.basename($0)} [options] <chain.bedpe> <bamfile|->
ruby #{File.basename($0)} [options] -n <RNAMEs_table> <bamfile|->
  Version: #{Version}
  chain.bedpe is created by chain2bedpe.awk
  Output is in the SAM format with header

HEADER
  opt.on('--sam',
         'The format of the input is SAM with a header'){ params[:sam] = true }
  opt.on('-h', '--header File', String,
         'The file including a SAM header for new @SQ'){ |v| params[:header] = v }
  opt.on('-n', '--noliftcoord File', String,
         'RNAMEs conversion table for RNAME-only liftOver. Each line consists of tab-delimited old_RNAME and new_RNAME'){ |v| params[:nolift] = v }
  opt.on('-z', '--zn',
         'Add ZN tag at second last optional tag to store the original RNAMEs'){ params[:zn] = true }
  opt.on('-s', '--soft',
         "Use a soft-clipping operation instead of an insertion operation for the end of reads"){ params[:softclip] = true }
  opt.on('-@', '--threads NUM', Integer,
         "Number of samtools threads [#{params[:threads]}]"){ |v| params[:threads] = v }
  opt.on('-m', '--memory SIZE', String,
         "Memory size for sort command [#{params[:sort_size]}]"){ |v| params[:sort_size] = v }
  opt.on('--dry', 
         "Dry run"){ params[:dry] = true }
  opt.parse!(ARGV)
}

case ARGV.size
when 1
  raise "The -n/--nolift option is required if the chain.bedpe is not provided" unless params[:nolift]
  cmd = "mawk -f #{File.dirname($0)}/bamliftover_simple.awk -- -n #{params[:nolift]} #{params[:sam] ? '--sam' : ''} #{params[:zn] ? '-z' : ''} #{ARGV[0]}"
when 2
  tmpheader1 = Tmpfile.create
  tmpheader2 = Tmpfile.create
  case 
  when ARGV[1] == '-'
    bamtmp = Tmpfile.create
    cmd = "tee #{bamtmp} | samtools view -h -@#{params[:threads]} - 2> /dev/null | "
    cmd = "samtools view -Sb -@#{params[:threads]} - 2> /dev/null | " + cmd if params[:sam]
    input1 = ''
    input2 = "-b #{bamtmp}"
  when params[:sam]
    cmd = ''
    input1 = ARGV[1]
    input2 = "-s #{ARGV[1]}"
  else
    cmd = "samtools view -h -@#{params[:threads]} #{ARGV[1]} 2> /dev/null | "
    input1 = ''
    input2 = "-b #{ARGV[1]}"
  end

  awkopts = [ [ "-H #{tmpheader1}",
                (params[:nolift] ? "-n #{params[:nolift]}" : nil),
              ],
              [ "-H #{tmpheader2}",
                (params[:softclip] ? '-c S' : nil)
              ],
              [ input2,
                "-p #{ARGV[0]}",
                (params[:header] ? "-Q #{params[:header]}" : nil),
                "-H #{tmpheader2}",
                "-H #{tmpheader1}",
                (params[:nolift] ? "-n #{params[:nolift]}" : nil),
                (params[:zn] ? '-z' : nil)
              ] 
            ].collect{ |opt| opt.compact.join(' ') }
  
  cmd += <<CMD
mawk -f #{File.dirname($0)}/bamliftover_samtobed6.awk -- #{awkopts[0]} #{input1} |
sort -k1,1 -k2,2n -S #{params[:sort_size]} |
bedtools intersect -wao -sorted -a stdin -b #{ARGV[0]} |
mawk -f #{File.dirname($0)}/bamliftover_bed6liftover.awk -- #{awkopts[1]} |
sort -k4,4n -k5,5nr -S #{params[:sort_size]} |
mawk -f #{File.dirname($0)}/bamliftover_bed6tosam.awk -- #{awkopts[2]}
CMD
else
  $stderr.puts opts.to_s
  exit!
end

if params[:dry]
  $stderr.puts cmd
else
  system(cmd) || exit!
end
