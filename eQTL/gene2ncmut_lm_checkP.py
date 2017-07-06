
import glob
import numpy as np
import os

index2gene={}
for line in open('./index2gene4lm.txt').read().splitlines():
    a=line.split('\t')
    index2gene[a[0]]=a[1]

model2p={}
fns=glob.glob('lm_model_p/*')
for fn in fns:
    index=fn.split('/')[-1].split('.')[0]
    gene=index2gene[index]
    for line in open(fn).read().splitlines():
        model2p[gene]=line

for prefix in ['lm_coef_p/']:
    fns=glob.glob(prefix+'/*')
    line_out=''
    line_out_all=''
    for fn in fns:
        index=fn.split('/')[-1].split('.')[0]
        gene=index2gene[index]
        if gene.split('|')[0]=='?':
            continue
#        p=1
        nominalP=1
        coef=0
        fdr=1
        for line in open(fn).read().splitlines():
            a=line.split('\t')
            nominalP_tmp=float(a[2])
            if nominalP_tmp<nominalP:
                nominalP=nominalP_tmp
                coef=a[1]
#                p=a[-1]
                fdr=a[3]
            line_out_all+='{}\t{}\n'.format(gene,line)
        line_out+='{}\t{}\t{}\t{}\t{}\n'.format(gene,coef,nominalP,fdr,model2p[gene])
#        line_out+='{}\t{}\t{}\t{}\n'.format(gene,coef,nominalP,p)
    open('./gene2ncmut_lm_p.txt','wb').write(line_out)
    open('./gene2ncmut_lm_p_allPairs.txt','wb').write(line_out_all)

os.system('Rscript run_qvalue.R')

