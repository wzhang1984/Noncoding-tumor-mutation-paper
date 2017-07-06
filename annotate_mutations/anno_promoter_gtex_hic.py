
import sys

binSize=10000

fn_prefix=sys.argv[1]
merge_size=sys.argv[2]

ENSG2Entrez={}
for line in open("/cellar/data/users/wzhang1984/TCGA_enhancer_annotations/HiC_GTEx/ENSG2Entrez.txt").read().rstrip().split("\n")[1:]:
    a=line.split("\t")
    if len(a)==2 and a[1]:
        ENSG2Entrez[a[0]]=a[1]

mask={}
for line in open("/cellar/data/users/wzhang1984/bcbio/genomes/Hsapiens/GRCh37/rnaseq-2014-07-14/ref-transcripts-mask.gtf"):
    a=line.split("\t")
    mask[a[-1].split('gene_id "')[1].split('"')[0]]=1

bin2TSS={}
coding_genes={}
for line in open("/cellar/data/users/wzhang1984/bcbio/genomes/Hsapiens/GRCh37/rnaseq-2014-07-14/ref-transcripts.gtf"):
    a=line.split("\t")
    if a[1]!="protein_coding" or a[2]!="transcript":
    # if a[2]!="transcript":
        continue
    ENSG=a[-1].split('gene_id "')[1].split('"')[0]
    gene_name=a[-1].split('gene_name "')[1].split('"')[0]
    if ENSG in mask:
        continue
    entrez=""
    if ENSG in ENSG2Entrez:
        entrez=ENSG2Entrez[ENSG]
    Chr="chr"+a[0]
    if a[6]=="+":
        tss=int(a[3])
    elif a[6]=="-":
        tss=int(a[4])
    Bin=Chr+"\t"+str(tss/binSize)
    if not Bin in bin2TSS:
        bin2TSS[Bin]=[]
    bin2TSS[Bin].append([Chr,tss,entrez+"|"+gene_name])
    if gene_name:
        coding_genes[gene_name]=1
    if entrez:
        coding_genes[entrez]=1

mut_island2promoter={}
header=True
for line in open(fn_prefix+"_merged_"+merge_size+"_homerAnno.txt"):
    if header:
        header=False
        continue
    a=line.split("\t")
    if a[9]=='NA':
        continue
    dist2tss=int(a[9])
    if dist2tss>=-1000 and dist2tss<=1000 and a[-1].rstrip()=="protein-coding":
    # if a[7].split(" ")[0]=="promoter-TSS" and a[-1].rstrip()=="protein-coding":
    # if a[7].split(" ")[0]=="promoter-TSS":
        mut_island2promoter["\t".join([a[1],str(int(a[2])-1),a[3]])]=a[11]+"|"+a[15]
        if a[11]:
            coding_genes[a[11]]=1
        if a[15]:
            coding_genes[a[15]]=1

# binSize=10000
# bin2CDS={}
# for line in open("/cellar/data/users/wzhang1984/bcbio/genomes/Hsapiens/GRCh37/rnaseq-2014-07-14/ref-transcripts.gtf"):
    # a=line.split("\t")
    # if a[2]!="CDS":
        # continue
    # Chr="chr"+a[0]
    # Start=int(a[3])-2
    # End=int(a[4])+2
    # startBin=Start/binSize
    # endBin=End/binSize
    # Bin=""
    # for i in range(startBin,endBin+1):
        # Bin=Chr+":"+str(i)
        # if not Bin in bin2CDS:
            # bin2CDS[Bin]=[]
        # bin2CDS[Bin].append([Start,End])

line_out=""
for line in open(fn_prefix+"_merged_"+merge_size+"_anno.txt").read().rstrip().split("\n"):
    a=line.split("\t")
    if a[0]=="Chr":
        line_out+="\t".join(["\t".join(a[:-4]),a[-2],a[-1],"genes"])+"\n"
    else:
        Chr=a[0]
        Start=int(a[1])
        End=int(a[2])
        # startBin=Start/binSize
        # endBin=End/binSize
        # cds_flag=False
        # for i in range(startBin,endBin+1):
            # Bin=Chr+":"+str(i)
            # if Bin in bin2CDS:
                # for CDS in bin2CDS[Bin]:
                    # [cdsStart,cdsEnd]=CDS
                    # if (Start>=cdsStart and Start<=cdsEnd) or (End>=cdsStart and End<=cdsEnd) or (cdsStart>=Start and cdsStart<=End) or (cdsEnd>=Start and cdsEnd<=End):
                        # cds_flag=True
                        # break
                # if cds_flag:
                    # break
        # if cds_flag:
            # continue
        genes={}
        isPromoter=False
        mut_island="\t".join(a[:3])
        mid=(Start+End)/2.0
        if a[6] and mut_island in mut_island2promoter:
        # if mut_island in mut_island2promoter:
            isPromoter=True
            genes[mut_island2promoter[mut_island]+"|TSS"]=1
        if a[-3]:
            for assignment in a[-3].split(","):
                if assignment[:3]=="chr":
                    continue
                    # if not a[6]:
                        # continue
                    # if isPromoter:
                        # continue
                    # [hicChr,start_end]=assignment.split(":")
                    # [hicStart,hicEnd]=[int(i) for i in start_end.split("-")]
                    # hicMid=(hicStart+hicEnd)/2.0
                    # if (mid>=hicStart-0.5 and mid<=hicEnd+0.5) or (hicMid>=Start-0.5 and hicMid<=End+0.5):
                        # local_remote="_l"
                    # else:
                        # local_remote="_r"
                    # startBin=hicStart/binSize
                    # endBin=hicEnd/binSize
                    # for binPos in range(startBin,endBin+1):
                        # Bin=hicChr+"\t"+str(binPos)
                        # if Bin in bin2TSS:
                            # for tss in bin2TSS[Bin]:
                                # [tssChr,position,gene]=tss
                                # if position>=hicStart-0.5 and position<=hicEnd+0.5:
                                    # genes[gene+"|HiC"+local_remote]=1
                # elif (assignment.split("|")[0] in coding_genes or assignment.split("|")[1] in coding_genes) and (a[6] or len(a[-3].split("epig"))>1):
                if (assignment.split("|")[0] in coding_genes or assignment.split("|")[1] in coding_genes) and (a[6] or len(a[-3].split("epig"))>1):
                # elif (assignment.split("|")[0] in coding_genes or assignment.split("|")[1] in coding_genes):
                # else:
                    genes[assignment]=1
        if genes!={}:
            line_out+="\t".join(["\t".join(a[:-4]),a[-2],a[-1],",".join(sorted(genes.keys()))])+"\n"

open(fn_prefix+"_merged_"+merge_size+"_anno_promoter_gtex_noHic.txt","w").write(line_out)


