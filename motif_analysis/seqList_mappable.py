
line_out=''
line_info=''
header=True
for line in open('../autoAnno_GeneHancer/TCGA_snv_mnv_merged_50_anno_promoter_genehancer.txt'):
    if header:
        header=False
        continue
    a=line.split('\n')[0].split('\t')
    for mut in a[4].split(","):
        [pat,Chr_StartEnd,source,alt]=mut.split("__")
        [Chr,StartEnd]=Chr_StartEnd.split(":")
        [Start,End]=[int(i) for i in StartEnd.split("-")]
        temp='{}:{}-{}'.format(Chr,Start-8,End+7)
        line_out+=temp+'\n'
        line_info+='{}_{}_{}\t{}\t{}:{}-{}\n'.format(Chr_StartEnd,alt,pat[:12],a[-1],a[0],a[1],a[2])

open('seqList_mappable_50.txt','wb').write(line_out)
open('info_mappable_50.txt','wb').write(line_info)

