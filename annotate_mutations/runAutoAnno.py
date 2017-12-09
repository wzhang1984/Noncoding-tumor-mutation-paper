import os

fn_prefix = "TCGA_snv_mnv"
merge_size = '50'
occurrence_cutoff = 5

print "1. mut4merge_alt.py"
os.system("python ./mut4merge_alt.py {} {} {}".format(fn_prefix,merge_size,occurrence_cutoff))

print "2. annotate_mut_GeneHancer.py"
os.system("python ./annotate_mut_GeneHancer.py {} {} {}".format(fn_prefix, merge_size, occurrence_cutoff))

