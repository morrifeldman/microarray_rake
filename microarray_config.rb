# Edit the following to adjust to your microarray data:
#####################################################

# uncomment and set to gse number to use a gse
# GSE =

Orig_cel_names = %w[
  GFP_0.CEL
  GFP_1.CEL
  HER2_0.CEL
  HER2_1.CEL
]

Probe_category = 'comprehensive'

Orig_cel_dir = if defined? GSE
  '../gse'
  else
  # This points to your input cel files if you not downloading a gse
  '../Raw_HER2'
end

Conditions = %w[control her2]

Control_name = 'control'

N_reps = 2

Cel_dir = 'cel_data'
