
import os

os.system('mkdir -p data4lm_peer lm_coef_p lm_model_p log')
os.system('rm data4lm_peer/* lm_coef_p/* lm_model_p/* log/*')

index2gene={}
for line in open('./index2gene4lm.txt').read().splitlines():
    a=line.split('\t')
    index2gene[a[0]]=a[1]

gene2expr_zscore={}
for line in open('PEER_zscore/expr_zscore.txt').read().splitlines():
    a=line.split('\t')
    gene2expr_zscore[a[0]]='exp\t{}\n'.format('\t'.join(a[1:]))

line_factors=''
for line in open('PEER_zscore/factors_35.txt').read().splitlines():
    a=line.split('\t')
    if a[0]=='mean':
        continue
    if a[0][0]=='h' and int(a[0][1:])>30:
        continue
    line_factors+=line+'\n'

for i in index2gene:
    fn=i+'.txt'
    gene=index2gene[i]
    line_out=''
    line_CNA=''
    for line in open('data4lm/{}'.format(fn)).read().splitlines():
        a=line.split('\t')
        if a[0]=='CNA':
            line_CNA=line+'\n'
        else:
            line_out+=line+'\n'
    line_out+=line_CNA
    line_out+=line_factors
    line_out+=gene2expr_zscore[gene]
    open('data4lm_peer/{}'.format(fn),'wb').write(line_out)

