
import numpy as np
import os

pats_ncm=set()
pats_12_ncm=set()
for line in open('../autoAnno/TCGA2pat.txt').read().rstrip().split('\n'):
    a=line.split('\t')
    if a[-1]!='0':
        pats_ncm.add(a[0][:15])
        pats_12_ncm.add(a[0][:12])

print 'Reading pat2CNA'
header=True
gene2CNA={}
pats_CNA=set()
for line in open('../CNA/pat2CNA.txt').read().rstrip().split('\n'):
    a=line.split('\t')
    if header:
        header=False
        pats=a[1:]
        pats_CNA=set(pats)
        continue
    gene=a[0].split('|')[0]
    if gene not in gene2CNA:
        gene2CNA[gene]={}
    v=[float(i) for i in a[1:]]
    for i in range(len(v)):
        if v[i]>=1.5:
            gene2CNA[gene][pats[i]]=1
        elif v[i]<=-1.5:
            gene2CNA[gene][pats[i]]=-1

print 'Reading ncmut'
gene2MI={}
MI2ncmut={}
MI2patPerBp={}
MI2concentrate={}
header=True
occur2freq={}
for line in open('../autoAnno/TCGA_snv_mnv_merged_50_anno_promoter_gtex_noHic.txt'):
    if header:
        header=False
        continue
    a=line.split('\n')[0].split('\t')
    MI=a[5]
    pats=set()
    mut2pats={}
    for mut in a[4].split(','):
        [pat,Chr_StartEnd,source,alt,TF_motifs]=mut.split("__")
        pat=pat[:15]
        pats.add(pat)
        if Chr_StartEnd not in mut2pats:
            mut2pats[Chr_StartEnd]=set()
        mut2pats[Chr_StartEnd].add(pat)
    MI2ncmut[MI]=pats
    MI2patPerBp[MI]=len(pats)/(float(a[2])-float(a[1])-1)
    MI2concentrate[MI]=0
    for mut in mut2pats:
        MI2concentrate[MI]=max(MI2concentrate[MI],len(mut2pats[mut])/float(len(pats)))

    for g in a[-1].split(','):
        [entrez,symbol,tp]=g.split('|')
        for gene in [entrez,symbol]:
            if gene not in ['','NA']:
                if gene not in gene2MI:
                    gene2MI[gene]=set()
                gene2MI[gene].add(MI)

    if a[3] not in occur2freq:
        occur2freq[a[3]]=0
    occur2freq[a[3]]+=1

open('occur2freq.txt','wb').write(''.join(['{}\t{}\n'.format(i,occur2freq[i]) for i in occur2freq]))

print 'Reading pat2disease'
disease2pats={}
pat2gender={}
for line in open('../pat2clin4surv_gender.txt').read().splitlines():
    a=line.split('\t')
    if a[0] in pats_12_ncm:
        if a[1] in ['STAD','UCEC']:
            continue
        if a[-1]=='':
            continue
        gender=''
        if a[1] not in disease2pats:
            disease2pats[a[1]]=set()
        disease2pats[a[1]].add(a[0])
        if a[-1]=='FEMALE':
            gender=1
        elif a[-1]=='MALE':
            gender=0
        pat2gender[a[0]]=gender

MI2motif={}
for line in open('../motifAnalysis/summarize_instances_50_compare_MI.txt').read().rstrip().split('\n'):
    a=line.split('\t')
    MI=a[1]
    if int(a[4])<4:
        continue
    motif_gain_loss_nPats='|'.join(a[2:5])
    if MI not in MI2motif:
        MI2motif[MI]=set()
    MI2motif[MI].add(motif_gain_loss_nPats)

print 'Reading pat2exp'
header_flag=True
gene2line_exp={}
gene2line_CNA={}
gene2line_ncmut={}
gene2exp_median={}
ncmut_pat=set()
MIs_all=set()
for line in open('./pat2exp.txt'):
    a=line.split('\n')[0].split('\t')
    if header_flag:
        header_flag=False
        header=a[0]
        pats_orig=a[1:]
        pats=[]
        line_gender='FEMALE'
        for pat in pats_orig:
            if pat in pats_CNA and pat in pats_ncm and pat[:12] in pat2gender:
                header+='\t'+pat
                pats.append(pat)
                line_gender+='\t{}'.format(pat2gender[pat[:12]])
        line_gender+='\n'
        header+='\n'
        print len(pats_orig)
        print len(pats)
        pats_set=set(pats)
        line_diseases=''
        for disease in sorted(disease2pats):
            if disease=='BRCA':
                continue
            line_diseases+=disease
            for pat in pats:
                if pat[:12] in disease2pats[disease]:
                    line_diseases+='\t1'
                else:
                    line_diseases+='\t0'
            line_diseases+='\n'
        continue

    [symbol,entrez]=a[0].split('|')
    v=[float(i) for i in a[1:]]
    v_median=np.median(v)
    if v_median<=1:
        continue
    gene2exp_median[a[0]]=v_median
    line_CNA=''
    line_ncmut=''

    line_exp=a[0]
    exps=[]
    for i in range(len(v)):
        if pats_orig[i] in pats:
            logexp=np.log2(max(v[i],1.0/8.0))
            line_exp+='\t{}'.format(logexp)
            exps.append(logexp)
    if np.median(exps)<=0:
        print a[0]
        continue
    line_exp+='\n'
    gene2line_exp[a[0]]=line_exp

    if symbol!='?' and symbol in gene2CNA:
        line_CNA='CNA'
        for pat in pats:
            if pat in gene2CNA[symbol]:
                line_CNA+='\t'+str(gene2CNA[symbol][pat])
            else:
                line_CNA+='\t0'
        line_CNA+='\n'
    gene2line_CNA[a[0]]=line_CNA

    gene2line_ncmut[a[0]]=line_ncmut
    MIs=set()
    if symbol!='?' and symbol in gene2MI:
        MIs=MIs | gene2MI[symbol]
    if entrez in gene2MI:
        MIs=MIs | gene2MI[entrez]
    if MIs==set():
        continue
    for MI in MIs:
        ncmut_pats=MI2ncmut[MI]
        if len(ncmut_pats&pats_set)<5:
            continue
#        if MI2patPerBp[MI]<0.5:
#            continue        
        if MI2concentrate[MI]<0.35:
            continue
#        if len(ncmut_pats&pats_set)<10 and (MI not in MI2motif):
#            continue
#        if MI not in MI2motif:
#            continue
        line_ncmut+=MI
        for pat in pats:
            if pat in ncmut_pats:
                line_ncmut+='\t1'
                ncmut_pat.add('{}\t{}'.format(MI,pat))
                MIs_all.add(MI)
            else:
                line_ncmut+='\t0'
        line_ncmut+='\n'
    if line_ncmut=='':
        continue
    gene2line_ncmut[a[0]]=line_ncmut

print len(ncmut_pat)
print len(MIs_all)

os.system('mkdir -p data4lm lm_coef_p log')
os.system('rm data4lm/* lm_coef_p/* log/*')

print 'output'
line_index2gene=''
index=1
for gene in sorted(gene2line_ncmut):
    [symbol,entrez]=gene.split('|')
    gene2line=gene2line_CNA[gene]+gene2line_ncmut[gene]
    if gene2line_ncmut[gene]=='':
        continue
    line_out=header+gene2line
    open('./data4lm/{}.txt'.format(index),'wb').write(line_out)
    line_index2gene+='{}\t{}\t{}\n'.format(index,gene,gene2exp_median[gene])
    index+=1

open('index2gene4lm.txt','wb').write(line_index2gene)

open('covariates.txt','wb').write(header+line_gender+line_diseases)

line_out=''
for gene in sorted(gene2line_exp):
    line_out+=gene2line_exp[gene]
open('pat2exp_log.txt','wb').write(header+line_out)
