import os
import sys
import glob
from multiprocessing import Pool

fn_out_prefix=sys.argv[1]
merge_size=sys.argv[2]
occurrence_cutoff=int(sys.argv[3])


binSize=10000
bin2CDS={}
for line in open("/cellar/data/users/wzhang1984/bcbio/genomes/Hsapiens/GRCh37/rnaseq-2014-07-14/ref-transcripts.gtf"):
    a=line.split("\t")
    if a[2]!="CDS":
        continue
    Chr="chr"+a[0]
    Start=int(a[3])-2
    End=int(a[4])+2
    startBin=Start/binSize
    endBin=End/binSize
    Bin=""
    for i in range(startBin,endBin+1):
        Bin=Chr+":"+str(i)
        if not Bin in bin2CDS:
            bin2CDS[Bin]=[]
        bin2CDS[Bin].append([Start,End])

def check_cds(Chr,Start,End):
    startBin=Start/binSize
    endBin=End/binSize
    cds_flag=False
    for i in range(startBin,endBin+1):
        Bin=Chr+":"+str(i)
        if Bin in bin2CDS:
            for CDS in bin2CDS[Bin]:
                [cdsStart,cdsEnd]=CDS
                if (Start>=cdsStart and Start<=cdsEnd) or (End>=cdsStart and End<=cdsEnd) or (cdsStart>=Start and cdsStart<=End) or (cdsEnd>=Start and cdsEnd<=End):
                    cds_flag=True
                    break
            if cds_flag:
                break
    return cds_flag

TCGA2line_out={}
TCGA2nMuts={}

def parse_mp(TCGA,fns,func,n_processes):
    pool = Pool(processes=n_processes)
    args = zip([TCGA]*len(fns),fns)
    rtn = pool.map(func,args)
    pool.close()
    pool.join()
    line_out_union=set.union(*rtn)
    return line_out_union

###############################################################################################

TCGA2pat={}
pat2TCGA={}

fn2size={}
for line in open('/cellar/data/users/wzhang1984/WGS/reanalyze/gdc_manifest.2016-10-05T21_21_43.100096.tsv').read().splitlines()[1:]:
    a=line.split('\t')
    fn2size[a[1]]=int(a[-2])

case_id2tcga={}
header=True
columns={}
pat2size={}
for line in open('/cellar/data/users/wzhang1984/WGS/reanalyze/files.txt').read().splitlines():
    row=line.split('\t')
    if header:
        header=False
        for i in range(len(row)):
            columns[row[i]]=i
        continue
    tcga=row[columns['cases_0_submitter_id']]
    case_id2tcga[row[columns['cases_0_case_id']]]=tcga
    fn=row[columns['file_name']]
    tcga_full=row[columns['cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id']]
    tissue_code=int(tcga_full[13:15])
    if tissue_code<10 and fn in fn2size:
        if not tcga in pat2TCGA:
            pat2TCGA[tcga]=tcga_full
            TCGA2pat[tcga_full]=tcga
            pat2size[tcga]=fn2size[fn]
        else:
            if fn2size[fn]>pat2size[tcga]:
                pat2TCGA[tcga]=tcga_full
                TCGA2pat[tcga_full]=tcga

for line in open('/cellar/data/users/wzhang1984/WGS/reanalyze/tnPairs.txt').read().splitlines():
    row=line.split('\t')
    TCGA2pat[row[1]]=row[0]
    pat2TCGA[row[0]]=row[1]

def parse_Snyder(args):
    [TCGA,fn]=args[:2]
    line_out_local=set()
    for line in open(fn).read().splitlines()[3:]:
        row=line.split('\t')
        if row[-1]=='Somatic':
            Chr='chr'+row[0]
            Start=int(row[1])
            End=Start+len(row[3])-1
            cds_flag=False
            cds_flag=check_cds(Chr,Start,End)
            if cds_flag:
                continue
            Start=str(Start)
            End=str(End)
            line_out_local.add('\t'.join([Chr,Start,End,TCGA+'__'+Chr+':'+Start+'-'+End+'__snv_mnv__'+row[3]+'>'+row[4],'.','.'])+'\n')
    return line_out_local

def parse_Snyder_wrapper(case_id,fns,TCGA2line_out,TCGA2nMuts):
    pat=case_id2tcga[case_id]
    TCGA=pat2TCGA[pat]
    if TCGA not in TCGA2line_out:
        TCGA2line_out[TCGA]=set()
    for fn in fns:
        line_out_tmp=parse_Snyder([TCGA,fn])
        TCGA2line_out[TCGA]=TCGA2line_out[TCGA]|line_out_tmp
    TCGA2nMuts[TCGA]=len(TCGA2line_out[TCGA])
    print pat,len(TCGA2line_out),TCGA2nMuts[TCGA]

for Dir in glob.glob('/cellar/data/users/wzhang1984/WGS/reanalyze/Snyder/MergedMutationFiles/*'):
    case_id=Dir.split('id')[1][:-1]
    fns=glob.glob(Dir+'/*')
    parse_Snyder_wrapper(case_id,fns,TCGA2line_out,TCGA2nMuts)

for Dir in glob.glob('/cellar/data/users/wzhang1984/WGS/reanalyze/Snyder/FromCluster229MergedMutationFiles/*'):
    case_id=Dir.split('/')[-1]
    fns=glob.glob(Dir+'/*')
    parse_Snyder_wrapper(case_id,fns,TCGA2line_out,TCGA2nMuts)

for Dir in glob.glob('/cellar/data/users/wzhang1984/WGS/reanalyze/Snyder/FromCluster211MergedMutationFiles/*'):
    if Dir.split('/')[-1]=='1231AnnotationProblems':
        continue
    for fn in glob.glob(Dir+'/*'):
        case_id=fn.split('/')[-1].split('.')[0]
        fns=[fn]
        parse_Snyder_wrapper(case_id,fns,TCGA2line_out,TCGA2nMuts)

for Dir in glob.glob('/cellar/data/users/wzhang1984/WGS/reanalyze/Snyder/FromCluster26MergedMutationFiles/*'):
    case_id=Dir.split('/')[-1]
    fns=glob.glob(Dir+'/*')
    parse_Snyder_wrapper(case_id,fns,TCGA2line_out,TCGA2nMuts)

###############################################################################################

tumor_fn2TCGA={}
for line in open('/cellar/data/users/wzhang1984/WGS/reanalyze/tnPairs_360.txt').read().splitlines()[1:]:
    row=line.split('\t')
    tumor_fn2TCGA[row[3]]=row[1]

groupID2TCGA={}
for line in open('/cellar/data/users/wzhang1984/WGS_reanalyze/groupID.txt').read().splitlines():
    row=line.split('\t')
    groupID2TCGA[row[1]]=tumor_fn2TCGA[row[5]+'.bam']

def parse_mutect(args):
    [TCGA,fn]=args[:2]
    line_out_local=set()
    for line in os.popen('gunzip -c {}'.format(fn)).read().splitlines():
        if line=='':
            continue
        if line[0]=='#':
            continue
        row=line.split('\t')
        if row[6]=='PASS' and len(row[7].split('SOMATIC'))>1:
            Chr='chr'+row[0]
            Start=int(row[1])
            End=Start+len(row[3])-1
            cds_flag=False
            cds_flag=check_cds(Chr,Start,End)
            if cds_flag:
                continue
            Start=str(Start)
            End=str(End)
            line_out_local.add('\t'.join([Chr,Start,End,TCGA+'__'+Chr+':'+Start+'-'+End+'__snv_mnv__'+row[3]+'>'+row[4],'.','.'])+'\n')
    return line_out_local

for Dir in glob.glob('/cellar/data/users/wzhang1984/WGS_reanalyze/vcf_files_358/*'):
    groupID=Dir.split('/')[-1]
    TCGA=groupID2TCGA[groupID]
    if TCGA not in TCGA2line_out:
        TCGA2line_out[TCGA]=set()
    fns=glob.glob(Dir+'/*')
    line_out_tmp=parse_mp(TCGA,fns,parse_mutect,13)
    TCGA2line_out[TCGA]=TCGA2line_out[TCGA]|line_out_tmp
    TCGA2nMuts[TCGA]=len(TCGA2line_out[TCGA])
    print TCGA[:12],len(TCGA2line_out),TCGA2nMuts[TCGA]

###############################################################################################

line_out=''
for TCGA in TCGA2line_out:
    line_out+=''.join(TCGA2line_out[TCGA])

open(fn_out_prefix+'.bed','w').write(line_out)

line_out=''
for TCGA in sorted(TCGA2nMuts):
    line_out+='{}\t{}\t{}\n'.format(TCGA,TCGA2pat[TCGA],TCGA2nMuts[TCGA])
open('TCGA2pat.txt','w').write(line_out)

print 'sort'
os.system('sort -k1,1 -k2,2n '+fn_out_prefix+'.bed > '+fn_out_prefix+'.sorted.bed')
os.system('rm '+fn_out_prefix+'.bed')

print 'bedtools merge'
os.system('bedtools merge -d '+merge_size+' -c 1,4 -o count,collapse -i '+fn_out_prefix+'.sorted.bed >'+fn_out_prefix+'_merged_'+merge_size+'.txt')
os.system('rm '+fn_out_prefix+'.sorted.bed')

if occurrence_cutoff>0:
    fn=fn_out_prefix+'_merged_'+merge_size+'.txt'
    line_out=''
    for line in open(fn).read().splitlines():
        a=line.split('\t')
        if int(a[3])>=occurrence_cutoff:
            line_out+=line+'\n'
    open(fn,'w').write(line_out)

print 'annotatePeaks.pl'
os.system('annotatePeaks.pl '+fn_out_prefix+'_merged_'+merge_size+'.txt hg19 -annStats annStats.txt >'+fn_out_prefix+'_merged_'+merge_size+'_homerAnno.txt')
# os.system('annotatePeaks.pl '+fn_out_prefix+'_merged_'+merge_size+'.txt hg19 >'+fn_out_prefix+'_merged_'+merge_size+'_homerAnno.txt')
