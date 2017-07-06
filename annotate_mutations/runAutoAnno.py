import os

fn_prefix="TCGA_snv_mnv"
merge_size='50'
occurrence_cutoff=3

print "0. mut4merge_alt.py"
os.system("python ./mut4merge_alt.py {} {} {}".format(fn_prefix,merge_size,occurrence_cutoff))

print "2. anotate_mut_islands_motif.py"
os.system("python ./annotate_mut_islands_motif.py {} {}".format(fn_prefix,merge_size))

print "3. anno_promoter_gtex_hic.py"
os.system("python ./anno_promoter_gtex_hic.py {} {}".format(fn_prefix,merge_size))

