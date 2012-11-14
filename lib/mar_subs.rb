
require 'debugger'
require 'subroutines'
require 'progress'

############## ps_mps_task

def call_if_proc(procedure, inputs)
  (procedure.is_a? Proc) ? procedure.(inputs) : procedure
end

def ps_mps_task(target, deps, &block)
  # creates a file task for a block both for mps and ps
  # target can be a filename or a procedure which accepts ps_or_mps and
  # evaluates to a filename
  # same for deps
  ['ps', 'mps'].each do |ps_or_mps|
    target_file = call_if_proc(target, ps_or_mps)
    dep_files = call_if_proc(deps, ps_or_mps)
    file target_file => dep_files do |f|
      puts "Creating #{target_file} from #{dep_files}"
      case block.parameters.count
        when 0 then block.()
        when 1 then block.(ps_or_mps)
        when 2 then block.(ps_or_mps, f)
      end
    end
    task default: target_file
    CLOBBER.include(target_file)
  end
end

####################### average_replicates

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

def average_replicates(in_file_name, out_file_name)
  File.open(in_file_name) do |f_in|
    header = f_in.gets.split
    cond_hash = parse_header(header[1..-1])
    cond_hash.each {|cond, cols| puts "#{cond} => #{cols}"}
    File.open(out_file_name, 'w') do |f_out|
      f_out << "ps_or_tc\t" << cond_hash.keys.join("\t") << "\n"
      f_in.each.with_progress('Averaging Replicates') do |line|
        line_ar = line.split
        out_ar = [line_ar.shift]
        out_ar << avg_line_reps(line_ar, cond_hash)
        f_out << out_ar.join("\t") << "\n"
      end
    end
  end
end


############################# bed_from_ps

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


############### make_average_bed
def form_avg_bed_line(line, dict)
  line_ar = line.split
  ps = line_ar.shift
  mean = line_ar.map{|e| e.to_f}.mean
  coord = dict[ps]
  coord.nil? ? '' : "#{coord}\t#{mean}\n"
end

def fill_avg_bed(in_file_name, out_file_io)
  puts "Loading Dictionary"
  dict = readKyotoHash('ps2bed')
  File.open(in_file_name) do |f_in|
    header = f_in.gets.split
    progress_msg = "Filling Average Bed from #{in_file_name}"
    f_in.each.with_progress(progress_msg) do |line|
      out_file_io << form_avg_bed_line(line, dict)
    end
  end
end

def make_average_bed(in_file_name, out_file_name, track_name, track_desc)
  track_name ||= in_file_name
  File.open(out_file_name, 'w') do |out_file_io|
    add_bed_header(out_file_io, track_name, track_desc)
    fill_avg_bed(in_file_name, out_file_io)
  end
end

############ normalize data
def normalize_data(control_val, treat_data)
  treat_data.map{|e| e ? e.to_f-control_val.to_f : e}
end

def save_ar(f, ar)
  f << ar.join("\t") << "\n"
end

def add_norm_header(fio, id_name, control_name, treatments)
  save_ar(fio, [id_name, control_name] + treatments)
end

def parse_line_hash(line_hash, control_name, treatments)
  control_val = line_hash[control_name]
  treat_data = line_hash.values_at(*treatments)
  return control_val, treat_data
end

def get_treatment_cols(header_ar, id, control)
  #remove id, control, assume the rest are treaments
  header_ar.reject{|e| [id, control].include?(e)}
end

=begin
Data_normalizer expects data that looks like the following:
where the data values are log2, so we can subtract them
ps_or_tc control treat1 treat2 treat3
323149   1.234   2.345  2.3345 5.333
...

control_col_name input tells which column is control

id_col_name tells which column contains the probe id
=end

def make_data_normalizer(control_col_name = 'control', id_col_name = 'ps_or_tc' )
  # data_normalizer will close over id_col_name and control_col_name
  # once created, it can be called on input and output datafiles
  ->(datafile_in, datafile_out) do
    File.open(datafile_in, 'r') do |fio_in|
      header_ar = fio_in.gets.split
      treatment_cols = get_treatment_cols(header_ar, id_col_name, control_col_name)
      File.open(datafile_out, 'w') do |fio_out|
        add_norm_header(fio_out, id_col_name, control_col_name, treatment_cols)
        msg = "Normalizing #{datafile_in} -> #{datafile_out}"
        fio_in.with_progress(msg).inject(fio_out) do |fout, line|
          line_hash = make_line_hash(line, header_ar)
          id = line_hash[id_col_name]
          control_val = line_hash[control_col_name]
          treat_data = line_hash.values_at(*treatment_cols)
          norm_data = normalize_data(control_val, treat_data)
          save_ar(fout, [id] + norm_data)
        end
      end
    end
  end
end


=begin
This function takes data that looks like the following:

ps_or_tc control treat1 treat2 treat3
323149   1.234   2.345  2.3345 5.333
...

control_col_name input tells which column is control

id_col_name tells which column contains the probe id

=end




