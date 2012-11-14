require_relative 'microarray_config'

# Replace this if your data are not in simple replicates
# Meaning you did a different number of reps for each condition
cel_files = Conditions.inject([]) do |acc, cond|
  N_reps.times {|i| acc << File.join(Cel_dir, "#{cond}_#{i+1}")}
  acc
end

#######################################
## Below here might not need editing if all goes well.

orig_cel_files = Orig_cel_names.map{|cel| File.join(Orig_cel_dir, cel)}

require 'rake/clean'
require_relative 'lib/mar_subs'

CLOBBER.include('ps','mps')

file_hash = Hash[orig_cel_files.zip(cel_files)]

######################## Get the GSE from NCBI
if defined? GSE
  gse_url = "http://www.ncbi.nlm.nih.gov/geosuppl/?acc=#{GSE}"
  
  directory Orig_cel_dir
  tar_file = File.join(Orig_cel_dir, "#{GSE}.tar")
  file tar_file => Orig_cel_dir do |f|
    sh "curl -o #{f.name} #{gse_url}"
  end
  #CLOBBER.include(Orig_cel_dir, tar_file) uncomment to clobber
  # the the downloaded tar file
  
  
  # the m switch tells tar to give the extracted
  # data a current timestamp.  This will keep the cel_files
  # newer than the tar file.  Otherwise they keep getting downloaded
  orig_cel_files.map{|cel| "#{cel}.gz"}.each do |cel_gz|
    file cel_gz => tar_file do |f|
      sh "tar -mxvf #{f.prerequisites[0]} -C #{Orig_cel_dir}"
    end
  end
  
  orig_cel_files.each do |cel|
    file cel => "#{cel}.gz" do |f|
      sh "gunzip -cv #{f.prerequisites[0]} > #{f.name}"
    end
  end
  
  CLOBBER.include(File.join(Orig_cel_dir, '*.gz'))
  CLOBBER.include(File.join(Orig_cel_dir, '*.CEL'))
else
  # put some code in here to fetch/unpack your data if needed
end

######################## Copy Files
directory Cel_dir
CLOBBER.include(Cel_dir)
file_hash.each do |orig_name, new_name|
  file new_name => [orig_name, Cel_dir] do |f|
    puts "copying #{f.prerequisites[0]} to #{f.name}"
    cp f.prerequisites[0], f.name
  end
  CLOBBER.include(new_name)
end

######################### Create cel_files.txt
cel_file_txt = 'cel_files.txt'
file 'cel_files.txt' => cel_files do
  File.open(cel_file_txt,'w') do |io|
    io << "cel_files\n"
    io << cel_files.join("\n")
  end
end
CLOBBER.include(cel_file_txt)

######### Configuration for APT, requires AFFX_ANALYSIS_FILES_PATH to be set
celDir = File.absolute_path('.')
celFiles = File.join(celDir,'*.CEL')
libBase = 'HuEx-1_0-st-v2.r2'
pgf = libBase + '.pgf'
clf = libBase + '.clf'
qcc = libBase + '.qcc'
bgp = libBase + '.antigenomic.bgp'
apt = 'apt-probeset-summarize'
action = 'rma-sketch'

############################## APT Summarize
summary_proc = ->(ps_or_mps) {"#{ps_or_mps}/rma-sketch.summary.txt"}
ps_mps_task( summary_proc, cel_files + [cel_file_txt]) do |ps_or_mps|
    ps_or_mps_flag = ps_or_mps == 'ps' ? "-s" : "-m"
    ps_or_mps_lib = "#{libBase}.dt1.hg18.#{Probe_category}.#{ps_or_mps}"
    cmd = "#{apt} -a #{action} -p #{pgf} -c #{clf} #{ps_or_mps_flag} #{ps_or_mps_lib} -b #{bgp} -o #{ps_or_mps} --qc-probesets #{qcc} --cel-files #{cel_file_txt}"
    sh cmd
end

############################### RemoveAffyHeader
expr_proc = ->(ps_or_mps) {"#{ps_or_mps}/expr.#{ps_or_mps}"}
ps_mps_task(expr_proc, summary_proc) do |ps_or_mps, f|
  puts "Striping Affy Header"
  sh "removeAffyHeader.rb -f #{f.prerequisites[0]} -o #{f.name}"
end

############################## Average Replicates
rep_average_proc = ->(ps_or_mps) {"#{ps_or_mps}/rep_average.#{ps_or_mps}"}
ps_mps_task(rep_average_proc, expr_proc) do |ps_or_mps, f|
  average_replicates(f.prerequisites[0], f.name)
end


############################## Make condition specific bedfiles
def cond_bed_name(ps_or_mps, cond)
  "#{ps_or_mps}/#{cond}.bed"
end

Conditions.each do |cond|
  # this proc will close over cond
  cond_bed_proc = ->(ps_or_mps) { cond_bed_name(ps_or_mps, cond) }
  ps_mps_task(cond_bed_proc, rep_average_proc) do |ps_or_mps, f|
    track_name = "raw_#{cond}"
    track_desc = "Raw Expression Values for #{cond}, #{Probe_category} probesets"
    bed_from_ps(f.prerequisites[0], cond, f.name, track_name, track_desc)
    sh "gzip -cv #{f.name} > #{f.name}.gz"
  end
end

######################## Combine the condition specific bed files and compress
multi_bed_proc = ->(ps_or_mps) { "#{ps_or_mps}/raw_multi.bed" }
all_cond_beds_proc = ->(ps_or_mps) do
  Conditions.map{|cond| cond_bed_name(ps_or_mps, cond)}
end
ps_mps_task(multi_bed_proc, all_cond_beds_proc) do |ps_or_mps, f|
  f.prerequisites.each do |bed_file|
    sh "cat #{bed_file} >> #{f.name}"
  end
end

multi_bed_gz_proc = ->(ps_or_mps) { "#{multi_bed_proc.(ps_or_mps)}.gz"}
ps_mps_task(multi_bed_gz_proc, multi_bed_proc) do |ps_or_mps, f|
  sh "gzip -cv #{f.prerequisites[0]} > #{f.name}"
end


########################### Make mean bed files across all datapoints
mean_bed_proc = ->(ps_or_mps) {"#{ps_or_mps}/all_tp_mean.bed"}
  # we average the raw data, because there might be a different number
  # of replicates in each condition
ps_mps_task(mean_bed_proc, expr_proc) do |ps_or_mps, f|
    track_name = "raw_tp_mean"
    track_desc = "Raw Expression Values Averaged Across All Time Points, #{Probe_category}.#{ps_or_mps}"
    make_average_bed(f.prerequisites[0], f.name, track_name, track_desc)
end

#################### Normalize data

data_normalizer = make_data_normalizer(Control_name)
#this is a proc/lambda
norm_avg_proc = ->(ps_or_mps) {"#{ps_or_mps}/norm_avg.#{ps_or_mps}"}
ps_mps_task(norm_avg_proc, rep_average_proc) do |ps_or_mps, f|
  data_normalizer.(f.prerequisites[0], f.name)
end


##################### 2-fold up and down





=begin

now we need to get down to business,  Find up and down regulated genes

1) Normalize data vs control
2) Label with gs_nm_ps or gs_nm_mps
3) For each condition, two fold up, two fold down
4) Same thing but with fdr, rely on matlab for this
5) 1 bed files for each condition, signal relative to control ps and mps
    This show up/down regulation at the probe or mps level
6) Ratio bedfiles
    ps (full) / mps (comprehensive) helps to see intron vs exon (pre vs mature)
    ps (comprehensive) / mps (comprehensive) helps to look for isoform variation

=end