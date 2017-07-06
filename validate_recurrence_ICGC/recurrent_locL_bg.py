
import random
import numpy as np

binSize=1000000
n_perm=10000

bin2ann={}
print 'Reading annotation file'
for line in open('/cellar/users/wzhang1984/soft/homer/data/genomes/hg19/hg19.basic.annotation').read().splitlines():
    a=line.split('\t')
    Chr=a[1]
    Start=int(a[2])
    End=int(a[3])
    ann=a[5]
    startBin=(Start-200)/binSize
    endBin=(End+200)/binSize
    for Bin in range(startBin,endBin+1):
        chr_bin='{}:{}'.format(Chr,Bin)
        if chr_bin not in bin2ann:
            bin2ann[chr_bin]={}
        if ann not in bin2ann[chr_bin]:
            bin2ann[chr_bin][ann]=[]
        bin2ann[chr_bin][ann].append([Start,End])

eQTL2info={}
eQTL2bg={}
print 'Reading eQTLs'
for line in open('../oncodrivefml/regions_MI5_glmnet.txt').read().splitlines():
    a=line.split('\t')
    Chr='chr'+a[0]
    Start=int(a[1])
    End=int(a[2])
    mid=int(round((Start+End)/2.0))
    half_len=int(round((End-Start)/2.0))
    startBin=(mid-200)/binSize
    endBin=(mid+200)/binSize
    eQTL_ann=''
    for Bin in range(startBin,endBin+1):
        chr_bin='{}:{}'.format(Chr,Bin)
        for ann in bin2ann[chr_bin]:
            for start_end in bin2ann[chr_bin][ann]:
                [ann_start,ann_end]=start_end
                if mid>=ann_start and mid<=ann_end:
                    eQTL_ann=ann
                    eQTL2info[a[-1]]=[Chr,Start,End,mid,half_len,ann]
                    break
            if eQTL_ann!='':
                break
        if eQTL_ann!='':
            break
    eQTL_bg_start=mid-1000000
    eQTL_bg_end=mid+1000000
    startBin=(eQTL_bg_start-200)/binSize
    endBin=(eQTL_bg_end+200)/binSize
    ranges=[]
    for Bin in range(startBin,endBin+1):
        chr_bin='{}:{}'.format(Chr,Bin)
        if eQTL_ann not in bin2ann[chr_bin]:
            continue
        for start_end in bin2ann[chr_bin][eQTL_ann]:
            [ann_start,ann_end]=start_end
            if ann_start>=eQTL_bg_start and ann_end<=eQTL_bg_end:
                ranges.append(range(ann_start,ann_end+1))
            elif ann_start<eQTL_bg_start and ann_end>=eQTL_bg_start:
                ranges.append(range(eQTL_bg_start,ann_end+1))
            elif ann_start<=eQTL_bg_end and ann_end>eQTL_bg_end:
                ranges.append(range(ann_start,eQTL_bg_end+1))
    ranges=set().union(*ranges)
    eQTL2bg[a[-1]]=np.random.choice(list(ranges), n_perm)
    print a[-1],eQTL_ann,len(ranges),eQTL2bg[a[-1]]


bin2mut={}
print 'Reading mutations'
cline=0
for line in open('../oncodrivefml/input_muts.txt'):
    cline+=1
    if cline%1000000==0:
        print cline
    a=line.rstrip().split('\t')
    Chr='chr'+a[0]
    Start=int(a[1])
    pat=a[-1][:12]
    startBin=(Start-200)/binSize
    endBin=(Start+200)/binSize
    for Bin in range(startBin,endBin+1):
        chr_bin='{}:{}'.format(Chr,Bin)
        if chr_bin not in bin2mut:
            bin2mut[chr_bin]=[]
        bin2mut[chr_bin].append([Start,pat])

line_log=''
line_out=''
print 'Calculating p-values'
for eQTL in eQTL2info:
    [Chr,Start,End,mid,half_len,ann]=eQTL2info[eQTL]
    Bin=mid/binSize
    chr_bin='{}:{}'.format(Chr,Bin)
    pats=set()
    for mut_pat in bin2mut[chr_bin]:
        [mut,pat]=mut_pat
        if mut>=Start and mut<=End:
            pats.add(pat)
    n_pats=len(pats)
    line_log+='{}\t{}\t{}\t{}\t{}\n'.format(eQTL,Chr,Start,End,n_pats)
    if n_pats==0:
        line_out+='{}\t{}\t{}\t{}\t{}\t{}\n'.format(eQTL,Chr,Start,End,n_pats,1)
        print eQTL,n_pats,0
        continue
    p=0
    n=0.0
    for eQTL_rand in eQTL2bg[eQTL]:
        n+=1
        eQTL_rand_start=eQTL_rand-half_len
        eQTL_rand_end=eQTL_rand+half_len
        Bin=eQTL_rand/binSize
        chr_bin='{}:{}'.format(Chr,Bin)
        pats_rand=set()
        if chr_bin not in bin2mut:
            line_log+='{}\t{}\t{}\t{}\t{}\n'.format(eQTL,Chr,eQTL_rand_start,eQTL_rand_end,0)
            continue
        for mut_pat in bin2mut[chr_bin]:
            [mut,pat]=mut_pat
            if mut>=eQTL_rand_start and mut<=eQTL_rand_end:
                pats_rand.add(pat)
        n_pats_rand=len(pats_rand)
        line_log+='{}\t{}\t{}\t{}\t{}\n'.format(eQTL,Chr,eQTL_rand_start,eQTL_rand_end,n_pats_rand)
        if n_pats_rand>=n_pats:
            p+=1
        if n>=1000 and p>=15:
            break
    print eQTL,n_pats,p/n
    line_out+='{}\t{}\t{}\t{}\t{}\t{}\n'.format(eQTL,Chr,Start,End,n_pats,p/n)

open('eQTL2recurrent_log.txt','w').write(line_log)
open('eQTL2recurrent_p.txt','w').write(line_out)
