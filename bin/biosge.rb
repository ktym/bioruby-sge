#!/usr/bin/env ruby
#
# = Bio::SGE -- Sun Grid Engine array job submitter (Bio::FlatFile query to SGE)
#
# Copyright::	Copyright (C) 2009, 2010 Toshiaki Katayama <mailto:ktym at hgc dot jp>
# License::	Distributes under the same terms as Ruby
# Site::	http://kanehisa.hgc.jp/~k/sge/
# Download::	http://kanehisa.hgc.jp/~k/sge/sge.rb
#
# == USAGE (AS A COMMAND)
#
# Usage:
#     % sge.rb \[options...\] -q input_file -t db_file -c 'command --opts #{query} #{target}'
# 
# Options:
#     -q or --query file
#        Specify a flatfile including multiple entries.
#     -t or --target file
#        Specify a database file to be used.
#     -c or --command 'string'
#        Specify a command line to be executed.
#        The following identifiers can be used in the command line 'string'.
#          '#{query}'       fragmented query file name (== input_file)
#          '#{target}'      target database file name
#          '#{work_dir}'    current working directory
#          '#{task_id}'     SGE_TASK_ID
#          '#{slice}'       -- task_id / @@slice (integer >= 1)
#          '#{input_file}'  -- 'input/#{slice}/#{task_id}'
#          '#{output_file}' -- 'output/#{slice}/#{task_id}'
#          '#{error_file}'  -- 'error/#{slice}/#{task_id}'
#     -o or --sge_opts 'string'
#        Additional options for the qsub command.
#          '-l s_vmem=16G -l mem_req=16' to reserve 16GB RAM for each job
#          '-l cpu_arch=xeon'            to limit to use xeon CPUs only
#        Resource reservation and backfill options:
#          '-R y -l s_rt=12:0:0'         to limit max exec time to 12h (SIGUSER1)
#          '-R y -l h_rt=12:0:0'         to limit max exec time to 12h (SIGKILL)
#          '-R y -pe mpi-fillup 4'       to reserve 4 threads for MPI
#     -m or --task_min integer
#        Start number of tasks (default is 1, increase to start from halfway).
#     -M or --taks_max integer
#        Last value (default is a total number of entries in query).
#     -s or --task_step integer
#        Number of processes per one job (default is 1000). Large value is
#        recommended for short tasks with a large number of queries, and
#        a small value (minimum is 1) can be used for time consuming tasks
#        with a small number of queries.
#     --clear
#        Remove a SGE script and output/error/log directories
#     --clean
#        Remove a count file and the extracted input directory
#     --distclean
#        Exec both of --clear and --clean
#     -h or --help
#        Print this help message.
# 
# Examples:
#     % sge.rb -q data/query.pep -t data/target.pep -c 'blastall -p blastp -i #{query} -d #{target}' -o '-l cpu_arch=xeon'
#     % sge.rb -q data/query.nuc -t /usr/local/db/blast/ncbi/nr -c 'blastall -p blastx -s 10 -i #{query} -d #{target}' -o '-l cpu_arch=xeon -l sjob -l s_vmem=4G,mem_req=4'
#     % sge.rb -q data/dme.nuc -t data/dme.genome -s 1 -c 'exonerate --bestn 1 --model est2genome --showtargetgff 1 --showvulgar yes #{query} #{target}'
#     % sge.rb -q data/hsa.pep -t data/Pfam-A.hmm -m 1000 -M 2000 -s 10 -c 'hmmscan --tblout output/#{slice}/#{task_id}.tbl #{target} #{query}'
#     % sge.rb -q data/refseq.gb -c 'bp_genbank2gff3.pl -out stdout #{query}'
#     % sge.rb --distclean
# 
# See also:
#     http://kanehisa.hgc.jp/~k/sge/
#
# == RESULTS
# 
# The execution results will be stored in the following files and directories.
# 
#   count.txt     # correspondence table of the file numbers and entry IDs
#   input/        # extracted sequence files (one file, one sequence)
#   output/       # outputs of the command (numberd same as the input files)
#   error/        # errors of the command (numberd same as the input files)
#   log/          # log files of the qsub run (stdout and stderr)
# 
# You can confirm whether there were no system errors during the SGE execution
# by sizes and contents of files in the log/ directory.
# 
# Then, check the error/ directory whether there was a problem or not in your
# jobs (some command may utilize the stderr to another purpose).
# 
# Finally, main results can be obtained from files in the output/ directory.
# 

require 'bio-sge'
require 'getoptlong'

def show_usage
  prog  = File.basename($0)
  usage = %Q[
Usage:
    % #{prog} \[options...\] -q input_file -t db_file -c 'command --opts \#{query} \#{target}'

Options:
    -q or --query file
       Specify a flatfile including multiple entries.
    -t or --target file
       Specify a database file to be used.
    -c or --command 'string'
       Specify a command line to be executed.
       The following identifiers can be used in the command line 'string'.
         '\#{query}'       fragmented query file name (== input_file)
         '\#{target}'      target database file name
         '\#{work_dir}'    current working directory
         '\#{task_id}'     SGE_TASK_ID
         '\#{slice}'       -- task_id / @@slice (integer >= 0)
         '\#{input_file}'  -- "input/\#{slice}/\#{task_id}"
         '\#{output_file}' -- "output/\#{slice}/\#{task_id}"
         '\#{error_file}'  -- "error/\#{slice}/\#{task_id}"
    -o or --sge_opts 'string'
       Additional options for the qsub command.
         '-l s_vmem=16G -l mem_req=16' to reserve 16GB RAM for each job
         '-l cpu_arch=xeon'            to limit to use xeon CPUs only
       Resource reservation and backfill options:
         '-R y -l s_rt=12:0:0'         to limit max exec time to 12h (SIGUSER1)
         '-R y -l h_rt=12:0:0'         to limit max exec time to 12h (SIGKILL)
         '-R y -pe mpi-fillup 4'       to reserve 4 threads for MPI
    -m or --task_min integer
       Start number of tasks (default is 1, increase to start from halfway).
    -M or --taks_max integer
       Last value (default is a total number of entries in query).
    -s or --task_step integer
       Number of processes per one job (default is 1000). Large value is
       recommended for short tasks with a large number of queries, and
       a small value (minimum is 1) can be used for time consuming tasks
       with a small number of queries.
    -h or --help
       Print this help message.
    --clear
       Remove a SGE script and output/error/log directories
    --clean
       Remove a count file and the extracted input directory
    --distclean
       Exec both of --clear and --clean

Examples:
    % #{prog} -q data/query.pep -t data/target.pep -c 'blastall -p blastp -i \#{query} -d \#{target}' -o '-l cpu_arch=xeon'
    % #{prog} -q data/query.nuc -t /usr/local/db/blast/ncbi/nr -c 'blastall -p blastx -s 10 -i \#{query} -d \#{target}' -o '-l cpu_arch=xeon -l sjob -l s_vmem=4G,mem_req=4'
    % #{prog} -q data/dme.nuc -t data/dme.genome -s 1 -c 'exonerate --bestn 1 --model est2genome --showtargetgff 1 --showvulgar yes \#{query} \#{target}'
    % #{prog} -q data/hsa.pep -t data/Pfam-A.hmm -m 1000 -M 2000 -s 10 -c 'hmmscan --tblout output/\#{slice}/\#{task_id}.tbl \#{target} \#{query}'
    % #{prog} -q data/refseq.gb -c 'bp_genbank2gff3.pl -out stdout \#{query}'
    % #{prog} --distclean

See also:
    http://kanehisa.hgc.jp/~k/sge/

]
  puts usage
  exit
end

$opts = Hash.new

args = GetoptLong.new(
  [ '--query',     '-q',  GetoptLong::REQUIRED_ARGUMENT ],
  [ '--target',    '-t',  GetoptLong::REQUIRED_ARGUMENT ],
  [ '--command',   '-c',  GetoptLong::REQUIRED_ARGUMENT ],
  [ '--sge_opts',  '-o',  GetoptLong::REQUIRED_ARGUMENT ],
  [ '--task_min',  '-m',  GetoptLong::REQUIRED_ARGUMENT ],
  [ '--task_max',  '-M',  GetoptLong::REQUIRED_ARGUMENT ],
  [ '--task_step', '-s',  GetoptLong::REQUIRED_ARGUMENT ],
  [ '--clear',            GetoptLong::NO_ARGUMENT ],
  [ '--clean',            GetoptLong::NO_ARGUMENT ],
  [ '--distclean',        GetoptLong::NO_ARGUMENT ],
  [ '--help',      '-h',  GetoptLong::NO_ARGUMENT ]
)

args.each_option do |name, value|
  case name
  when /--query/
    $opts[:query] = value
  when /--target/
    $opts[:target] = value
  when /--command/
    $opts[:command] = value
  when /--sge_opts/
    $opts[:sge_opts] = value
  when /--task_min/
    $opts[:task_min] = value.to_i
  when /--task_max/
    $opts[:task_max] = value.to_i
  when /--task_step/
    $opts[:task_step] = value.to_i
  when /--clear/
    $opts[:clear] = true
  when /--clean/
    $opts[:clean] = true
  when /--distclean/
    $opts[:clear] = true
    $opts[:clean] = true
  when /--help/
    $opts[:help] = true
  end
end

if $opts[:clear]
  sge = Bio::SGE.new
  sge.clear
end

if $opts[:clean]
  sge = Bio::SGE.new
  sge.clean
end

show_usage if $opts[:help] or !$opts[:command]

sge = Bio::SGE.new { |opt|
  opt.query     = $opts[:query]     if $opts[:query]
  opt.target    = $opts[:target]    if $opts[:target]
  opt.command   = $opts[:command]   if $opts[:command]
  opt.sge_opts  = $opts[:sge_opts]  if $opts[:sge_opts]
  opt.task_min  = $opts[:task_min]  if $opts[:task_min]
  opt.task_max  = $opts[:task_max]  if $opts[:task_max]
  opt.task_step = $opts[:task_step] if $opts[:task_step]
}
sge.prepare
sge.submit
