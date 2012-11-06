require 'rake/clean'
require 'debugger'
require 'subroutines'
require 'progress'

probe_category = 'comprehensive'

CLOBBER.include(probe_category, '*.CEL', '*.gz', '*_*', 'replicates.txt')

celDir = File.absolute_path('.')
celFiles = File.join(celDir,'*.CEL')
libBase = 'HuEx-1_0-st-v2.r2'
pgf = libBase + '.pgf'
clf = libBase + '.clf'
ps = "HuEx-1_0-st-v2.r2.dt1.hg18.#{probe_category}.ps"
qcc = libBase + '.qcc'
bgp = libBase + '.antigenomic.bgp'
apt = 'apt-probeset-summarize'
action = 'rma-sketch'

zipGSE = '/Volumes/Spin/GSE/GSE24391_RAW.tar'

file 'GSM601318.CEL' do

end

cel_files = %w[
  GSM601318.CEL
  GSM601319.CEL
  GSM601320.CEL
  GSM601321.CEL
  GSM601322.CEL
  GSM601323.CEL
  GSM601324.CEL
  GSM601325.CEL
  GSM601326.CEL
  GSM601327.CEL
  GSM601328.CEL
  GSM601329.CEL
  GSM601330.CEL
  GSM601331.CEL
  GSM601332.CEL
  GSM601333.CEL
  GSM601334.CEL
  GSM601335.CEL
  GSM601336.CEL
  GSM601337.CEL
  GSM601338.CEL
]

conditions = %w[E000 E020 E040 E060 E120 E240 E480]
replicates = conditions.each_with_object([]) do |cond, acc|
  3.times {|i| acc << "#{cond}_#{i+1}"}
end

file_hash = Hash[cel_files.zip(replicates)]

file replicates[0] do
  sh "tar xvf #{zipGSE}"
  sh "gunzip -v *.gz"
  file_hash.each do |orig_name, new_name|
    puts "moving #{orig_name} to #{new_name}"
    mv orig_name, new_name
  end
end

cel_file_txt = 'replicates.txt'
file "#{probe_category}/rma-sketch.summary.txt" => replicates[0] do
  File.open(cel_file_txt,'w') do |io|
    io << "cel_files\n"
    io << replicates.join("\n")
  end
  cmd = "#{apt} -a #{action} -p #{pgf} -c #{clf} -s #{ps} -b #{bgp} -o #{probe_category} --qc-probesets #{qcc} --cel-files #{cel_file_txt}"
  sh cmd
end

file "#{probe_category}/expr.ps" => "#{probe_category}/rma-sketch.summary.txt" do |f|
  sh "removeAffyHeader.rb -f #{f.prerequisites[0]} -o #{f.name}"
end

def parse_header(header)
  # This function expect that the replicates will be deliminated by '_'
  # e.g. E000_1 E000_2 are replicates
  head_base = header.map{|h| h.split('_')[0]}
  head_base.each_with_index.with_object({}) do |(c, i), h|
    h[c] ||= []
    h[c] << i
  end
end

def avg_line_reps(data_row_ar, cond_hash)
  cond_hash.each_value.with_object([]) do |rep_cols, acc|
    acc << rep_cols.map{|col| data_row_ar[col].to_f}.mean
  end
end


def average_affy_replicates(in_file_name, out_file_name)
  File.open(in_file_name) do |f_in|
    header = f_in.gets.split
    cond_hash = parse_header(header[1..-1])
    cond_hash.each {|cond, cols| puts "#{cond} => #{cols}"}
    File.open(out_file_name, 'w') do |f_out|
      f_out << "ps\t" << cond_hash.keys.join("\t") << "\n"
      f_in.each.with_progress('Averaging Replicates') do |line|
        line_ar = line.split
        out_ar = [line_ar.shift]
        out_ar << avg_line_reps(line_ar, cond_hash)
        f_out << out_ar.join("\t") << "\n"
      end
    end
  end
end

file "#{probe_category}/rep_average.ps" => "#{probe_category}/expr.ps" do |f|
  average_affy_replicates(f.prerequisites[0], f.name)
end

def add_bed_header(file_io, name, desc)
  file_io << 'track type=bedGraph visibility=full color=0,0,0 altColor=127,127,127 autoScale=on graphType=bar alwaysZero=off windowingFunction=mean smoothingWindow=8 maxHeightPixels=64:64:8 '
  file_io << %Q[description='#{desc}' ] unless desc.nil?
  file_io << %Q[name='#{name}' \n] unless name.nil?
end

def make_line_hash(line, header)
  Hash[header.zip(line.split)]
end

def form_bed_line(line, header, col, dict)
  ps_col = header[0]
  line_hash = make_line_hash(line, header)
  val = line_hash[col]
  ps = line_hash[ps_col]
  coord = dict[ps]
  coord.nil? ? '' : "#{coord}\t#{val}\n"
end
    
def fill_bed(in_file_name, col, out_file_io)
  puts "Loading Dictionary"
  dict = readKyotoHash('ps2bed')
  File.open(in_file_name) do |f_in|
    header = f_in.gets.split
    progress_msg = "Filling Bed with #{col} from #{in_file_name}"
    f_in.each.with_progress(progress_msg) do |line|
      out_file_io << form_bed_line(line, header, col, dict)
    end
  end
end

def bed_from_ps(in_file_name, col, out_file_name, track_name = nil, desc = nil)
  track_name ||= in_file_name
  File.open(out_file_name, 'w') do |out_file_io|
    add_bed_header(out_file_io, track_name, desc)
    fill_bed(in_file_name, col, out_file_io)
  end
end

def bed_name(cond, probe_category)
  "#{probe_category}/#{cond}.bed"
end

example_bed = "#{probe_category}/#{conditions[0]}.bed"
file example_bed => "#{probe_category}/rep_average.ps" do |f|
  conditions.each do |cond|
    track_name = "raw_#{cond}"
    track_desc = "Raw Expression Values for #{cond}"
    out_file_name = bed_name(cond, probe_category)
    bed_from_ps(f.prerequisites[0], cond, out_file_name, track_name, track_desc)
    sh "gzip -cv #{out_file_name} > #{out_file_name}.gz"
  end
end

multi_bed_file = "#{probe_category}/raw_multi.bed"
file "#{multi_bed_file}.gz" => example_bed do
  conditions.each do |cond|
    bed_file = bed_name(cond, probe_category)
    sh "cat #{bed_file} >> #{multi_bed_file}"
  end
  sh "gzip -v #{multi_bed_file}"
end

task :default => "#{multi_bed_file}.gz"