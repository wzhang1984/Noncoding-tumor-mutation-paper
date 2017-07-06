
import glob
import numpy as np
import scipy.spatial.distance
import scipy.stats
import sys

fn_prefix=sys.argv[1]
merge_size=sys.argv[2]

def euclidean2states(scores):
    euclidean_P=scipy.spatial.distance.euclidean(scores,[1,0,0])
    euclidean_E=scipy.spatial.distance.euclidean(scores,[0,1,0])
    euclidean_I=scipy.spatial.distance.euclidean(scores,[0,0,1])
    euclidean_PE=scipy.spatial.distance.euclidean(scores,[0.5,0.5,0])
    euclidean_PI=scipy.spatial.distance.euclidean(scores,[0.5,0,0.5])
    euclidean_EI=scipy.spatial.distance.euclidean(scores,[0,0.5,0.5])
    euclidean_PEI=scipy.spatial.distance.euclidean(scores,[1.0/3,1.0/3,1.0/3])
    euclidean_min=min([euclidean_P,euclidean_E,euclidean_I,euclidean_PE,euclidean_PI,euclidean_EI,euclidean_PEI])
    if euclidean_PEI==euclidean_min:
        states="PEI"
    elif euclidean_PE==euclidean_min:
        states="PE"
    elif euclidean_EI==euclidean_min:
        states="EI"
    elif euclidean_PI==euclidean_min:
        states="PI"
    elif euclidean_E==euclidean_min:
        states="E"
    elif euclidean_P==euclidean_min:
        states="P"
    elif euclidean_I==euclidean_min:
        states="I"
    return states


complement={}
complement["A"]="T"
complement["C"]="G"
complement["G"]="C"
complement["T"]="A"

ENSG2Entrez={}
for line in open("/cellar/data/users/wzhang1984/TCGA_enhancer_annotations/HiC_GTEx/ENSG2Entrez.txt").read().rstrip().split("\n")[1:]:
    a=line.split("\t")
    if len(a)==2 and a[1]:
        ENSG2Entrez[a[0]]=a[1]

binSize=10000

fns=glob.glob("/cellar/users/wzhang1984/projects/20150511_TCGA_enhancer/HiC_GTEx/GTEx_v6/*")

bin2roadmap_assignments={}
print "Reading roadmap_stringent_enhancers_with_Ensembl_ID_and_name_and_cell_types.bed.txt"
for line in open("/cellar/data/users/wzhang1984/TCGA_enhancer_annotations/roadmap_assignments/roadmap_stringent_enhancers_with_Ensembl_ID_and_name_and_cell_types.bed.txt").read().rstrip().split("\n")[1:]:
    a=line.split("\t")
    if a[3]==".":
        continue
    Chr=a[0]
    Start=int(a[1])
    End=int(a[2])
    position="\t".join(a[1:3])
    startBin=Start/binSize
    endBin=End/binSize
    for binPos in range(startBin,endBin+1):
        Bin=Chr+"\t"+str(binPos)
        if not Bin in bin2roadmap_assignments:
            bin2roadmap_assignments[Bin]={}
        if not position in bin2roadmap_assignments[Bin]:
            bin2roadmap_assignments[Bin][position]={}
        for g in a[3].split(" "):
            [ENSG,name]=g.split("_")
            entrez=""
            if ENSG in ENSG2Entrez:
                entrez=ENSG2Entrez[ENSG]
            bin2roadmap_assignments[Bin][position][entrez+"|"+name]=1

bin2motif={}
print "Reading motifs"
for line in open("/cellar/users/wzhang1984/projects/20150511_TCGA_enhancer/ActiveGenomicRegions/factorbookMotifPos.txt").read().rstrip().split("\n")[1:]:
    a=line.split("\t")
    Chr=a[1]
    Start=int(a[2])
    End=int(a[3])
    motif=a[4]
    strand=a[-1]
    startBin=Start/binSize
    endBin=End/binSize
    for binPos in range(startBin,endBin+1):
        Bin=Chr+"\t"+str(binPos)
        if not Bin in bin2motif:
            bin2motif[Bin]=[]
        bin2motif[Bin].append([Start,End,motif,strand])

motif2TF_ENCODE={}
print "Reading motifs"
for line in open("/cellar/users/wzhang1984/projects/20150511_TCGA_enhancer/ActiveGenomicRegions/ENCODE.tf.bound.union.bed").read().rstrip().split("\n"):
    a=line.split("\t")
    Chr=a[0]
    Start=int(a[1])
    End=int(a[2])
    motif=a[3]
    strand=a[-2]
    TF=a[-1]
    if motif not in motif2TF_ENCODE:
        motif2TF_ENCODE[motif]=set()
    motif2TF_ENCODE[motif].add(TF)
    startBin=Start/binSize
    endBin=End/binSize
    for binPos in range(startBin,endBin+1):
        Bin=Chr+"\t"+str(binPos)
        if not Bin in bin2motif:
            bin2motif[Bin]=[]
        bin2motif[Bin].append([Start,End,motif,strand])

motif2TFs={}
print "Readin motifCanonical"
for line in open("/cellar/users/wzhang1984/projects/20150511_TCGA_enhancer/ActiveGenomicRegions/factorbookMotifCanonical.txt").read().rstrip().split("\n")[1:]:
    a=line.split("\t")
    TF=a[0]
    if len(a)>1 and a[1]:
        for motif in a[1].split(","):
            if not motif in motif2TFs:
                motif2TFs[motif]={}
            motif2TFs[motif][TF]=1

bin2GTEx={}
for fn in fns:
    print "Reading "+fn.split("/")[-1]
    for line in open(fn).read().rstrip().split("\n")[1:]:
        a=line.split("\t")
        Chr="chr"+a[13]
        position=int(a[14])
        Bin=Chr+"\t"+str(position/binSize)
        ENSG=a[1].split(".")[0]
        entrez=""
        if ENSG in ENSG2Entrez:
            entrez=ENSG2Entrez[ENSG]
        if not Bin in bin2GTEx:
            bin2GTEx[Bin]={}
        if not position in bin2GTEx[Bin]:
            bin2GTEx[Bin][position]={}
        bin2GTEx[Bin][position][entrez+"|"+a[26]]=1

fns=glob.glob("/cellar/users/wzhang1984/projects/20150511_TCGA_enhancer/HiC_GTEx/HiC/*")

bin2HiC={}
for fn in fns:
    print "Reading "+fn.split("/")[-1]
    for line in open(fn).read().rstrip().split("\n")[1:]:
        a=line.split("\t")
        chr1="chr"+a[0]
        start1=int(a[1])
        end1=int(a[2])
        chr2="chr"+a[3]
        start2=int(a[4])
        end2=int(a[5])
        chrPosition1="chr"+"\t".join(a[:3])
        chrPosition2="chr"+"\t".join(a[3:6])
        startBin=start1/binSize
        endBin=end1/binSize
        for binPos in range(startBin,endBin+1):
            Bin=chr1+"\t"+str(binPos)
            if not Bin in bin2HiC:
                bin2HiC[Bin]={}
            if not chrPosition1 in bin2HiC[Bin]:
                bin2HiC[Bin][chrPosition1]={}
            bin2HiC[Bin][chrPosition1][chrPosition2]=1
        startBin=start2/binSize
        endBin=end2/binSize
        for binPos in range(startBin,endBin+1):
            Bin=chr2+"\t"+str(binPos)
            if not Bin in bin2HiC:
                bin2HiC[Bin]={}
            if not chrPosition2 in bin2HiC[Bin]:
                bin2HiC[Bin][chrPosition2]={}
            bin2HiC[Bin][chrPosition2][chrPosition1]=1

bin2segment={}
print "Reading ActiveGenomicRegions.bed"
for line in open("/cellar/users/wzhang1984/projects/20150511_TCGA_enhancer/ActiveGenomicRegions/ActiveGenomicRegions.bed").read().rstrip().split("\n"):
    a=line.split("\t")
    scores=np.array([float(i) for i in a[4].split(",")])
    Chr=a[0]
    Start=int(a[1])
    End=int(a[2])
    startBin=Start/binSize
    endBin=End/binSize
    for binPos in range(startBin,endBin+1):
        Bin=Chr+"\t"+str(binPos)
        if not Bin in bin2segment:
            bin2segment[Bin]={}
        bin2segment[Bin]["\t".join([str(Start),str(End)])]=scores

bin2RegulomeDB={}
print "Reading RegulomeDB.dbSNP141.txt"
for line in open("/cellar/users/wzhang1984/projects/20150511_TCGA_enhancer/RegulomeDB/RegulomeDB.dbSNP141.txt"):
    a=line.rstrip().split("\n")[0].split("\t")
    Chr=a[0]
    Start=int(a[1])
    End=int(a[2])
    startBin=Start/binSize
    endBin=End/binSize
    for binPos in range(startBin,endBin+1):
        Bin=Chr+"\t"+str(binPos)
        if not Bin in bin2RegulomeDB:
            bin2RegulomeDB[Bin]={}
        bin2RegulomeDB[Bin]["\t".join([str(Start),str(End)])]=a[-1]

peak2nearestGene={}
print "Reading "+fn_prefix+"_merged_"+merge_size+"_homerAnno.txt"
for line in open(fn_prefix+"_merged_"+merge_size+"_homerAnno.txt").read().rstrip().split("\n")[1:]:
    a=line.split("\t")
    peak2nearestGene["\t".join(a[1:4])]=a[15]

bin2TF={}
print "Reading TFs"
for line in open("/cellar/users/wzhang1984/projects/20150511_TCGA_enhancer/ActiveGenomicRegions/wgEncodeRegTfbsClusteredV3.bed").read().rstrip().split("\n"):
    a=line.split("\t")
    Chr=a[0]
    Start=int(a[1])
    End=int(a[2])
    TF=a[3]
    startBin=Start/binSize
    endBin=End/binSize
    for binPos in range(startBin,endBin+1):
        Bin=Chr+"\t"+str(binPos)
        if not Bin in bin2TF:
            bin2TF[Bin]=[]
        bin2TF[Bin].append([Start,End,TF])

line_out="Chr\tStart\tEnd\tcount\tTF_motif\tID\tstates\tscores\tRegulomeDB_score\tnearest_gene\tinteractions\tnTFs\tnMotifs\n"

print "Reading "+fn_prefix+"_merged_"+merge_size+".txt"
cline=0
c1=0
c2=0
for line in open(fn_prefix+"_merged_"+merge_size+".txt").read().rstrip().split("\n"):
    cline+=1
    if cline/1000==cline/1000.0:
        print cline
    a=line.split("\t")
    Chr=a[0]
    Start=int(a[1])-100
    End=int(a[2])+100
    mid=(Start+End)/2.0
    # if int(a[3])<3:
        # continue
    startBin=Start/binSize
    endBin=End/binSize
    mapped_segments={}
    scores=np.array([0.0,0.0,0.0])
    interactions={}
    RegulomeDB_score=""
    region2TFs={}
    region2motifs={}
    region="\t".join(a[:3])
    for binPos in range(startBin,endBin+1):
        Bin=Chr+"\t"+str(binPos)
        if Bin in bin2segment:
            for segment in bin2segment[Bin]:
                if segment in mapped_segments:
                    continue
                [segmentStart,segmentEnd]=[int(i) for i in segment.split("\t")]
                segmentMid=(segmentStart+segmentEnd)/2.0
                if (mid>=segmentStart-0.5 and mid<=segmentEnd+0.5) or (segmentMid>=Start-0.5 and segmentMid<=End+0.5):
                    mapped_segments[segment]=1
                    segmentScores=bin2segment[Bin][segment]
                    overlap=min(End,segmentEnd)-max(Start,segmentStart)
                    scores+=segmentScores/np.array([134.0,134.0,6.0])*overlap/(End-Start)
    # if max(scores)==0:
        # continue
    if max(scores)==0:
        states=""
        c2+=1
    else:
        states=euclidean2states(scores/sum(scores))
        c1+=1

    TF_motif={}
    for patInfo in a[4].split(","):
        [pat,position,Type,ref_alt]=patInfo.split("__")
        [ref,alt]=ref_alt.split(">")
        [mutChr,start_end]=position.split(":")
        [mutStart,mutEnd]=[int(i) for i in start_end.split("-")]
        mutMid=(mutStart+mutEnd)/2.0
        [mutStartBin,mutEndBin]=[mutStart/binSize,mutEnd/binSize]
        for binPos in range(mutStartBin,mutEndBin+1):
            Bin=mutChr+"\t"+str(binPos)
            if Bin in bin2motif:
                for motif in bin2motif[Bin]:
                    [motifStart,motifEnd,motifName,motifStrand]=motif
                    motifMid=(motifStart+motifEnd)/2.0
                    if (mutEnd>=motifStart-2.0 and mutStart<=motifEnd+2.0):
                        TFs_covered=set()
                        if motifName in motif2TF_ENCODE:
                            TFs_covered=motif2TF_ENCODE[motifName]
                        else:
                            if motifName not in motif2TFs:
                                continue
                            TFs=motif2TFs[motifName]
                            if Bin in bin2TF:
                                for TF_i in bin2TF[Bin]:
                                    [TFStart,TFEnd,TF]=TF_i
                                    if TF not in TFs:
                                        continue
                                    TFMid=(TFStart+TFEnd)/2.0
                                    if (motifEnd>=TFStart and motifStart<=TFEnd):
                                        TFs_covered.add(TF)
                        # if TFs_covered==set():
                            # continue
                        if motifStrand=="+":
                            mut_rel_motif=[mutStart-motifStart,mutEnd-motifStart]
                            ref_motif=ref
                            alt_motif=alt
                        elif motifStrand=="-":
                            mut_rel_motif=[motifEnd-mutEnd,motifEnd-mutStart]
                            ref_motif="".join([complement[i] for i in ref[::-1]])
                            alt_motif="".join([complement[i] for i in alt[::-1]])
                        if not patInfo in TF_motif:
                            TF_motif[patInfo]={}
                        for TF in TFs_covered:
                            TF_motif[patInfo]["{}>{}|{}|{}|{}|{}>{}".format(TF,motifName,mut_rel_motif[0],mut_rel_motif[1],motifStrand,ref_motif,alt_motif)]=1
    # if TF_motif=={}:
        # continue

    for binPos in range(startBin,endBin+1):
        Bin=Chr+"\t"+str(binPos)
        if Bin in bin2GTEx:
            for position in bin2GTEx[Bin]:
                if position>=Start and position<=End:
                    for gene in bin2GTEx[Bin][position]:
                        interactions[gene+"|eQTL"]=1
        if Bin in bin2HiC:
            for chrPosition1 in bin2HiC[Bin]:
                [Chr1,hicStart1,hicEnd1]=chrPosition1.split("\t")
                hicStart1=int(hicStart1)
                hicEnd1=int(hicEnd1)
                hicMid1=(hicStart1+hicEnd1)/2.0
                if (mid>=hicStart1-0.5 and mid<=hicEnd1+0.5) or (hicMid1>=Start-0.5 and hicMid1<=End+0.5):
                    interactions[Chr1+":"+str(hicStart1)+"-"+str(hicEnd1)]=1
                    for chrPosition2 in bin2HiC[Bin][chrPosition1]:
                        [Chr2,hicStart2,hicEnd2]=chrPosition2.split("\t")
                        interactions[Chr2+":"+hicStart2+"-"+hicEnd2]=1
        if Bin in bin2roadmap_assignments:
            for position in bin2roadmap_assignments[Bin]:
                [roadmapStart,roadmapEnd]=[int(i) for i in position.split("\t")]
                roadmapMid=(roadmapStart+roadmapEnd)/2.0
                if (roadmapMid>=Start-0.5 and roadmapMid<=End+0.5) or (mid>=roadmapStart-0.5 and mid<=roadmapEnd+0.5):
                    for gene in bin2roadmap_assignments[Bin][position]:
                        interactions[gene+"|epig"]=1
        if Bin in bin2RegulomeDB:
            for snp in bin2RegulomeDB[Bin]:
                [snpStart,snpEnd]=[int(i) for i in snp.split("\t")]
                snpMid=(snpStart+snpEnd)/2.0
                if (mid>=snpStart-0.5 and mid<=snpEnd+0.5) or (snpMid>=Start-0.5 and snpMid<=End+0.5):
                    if RegulomeDB_score!="":
                        if bin2RegulomeDB[Bin][snp]<RegulomeDB_score:
                            RegulomeDB_score=bin2RegulomeDB[Bin][snp]
                    else:
                        RegulomeDB_score=bin2RegulomeDB[Bin][snp]
        if Bin in bin2TF:
            for TF_i in bin2TF[Bin]:
                [TFStart,TFEnd,TF]=TF_i
                TFMid=(TFStart+TFEnd)/2.0
                if (mid>=TFStart-0.5 and mid<=TFEnd+0.5) or (TFMid>=Start-0.5 and TFMid<=End+0.5):
                    region2TFs[TF]=1
        if Bin in bin2motif:
            for motif_i in bin2motif[Bin]:
                [motifStart,motifEnd,motif,motifStrand]=motif_i
                motifMid=(motifStart+motifEnd)/2.0
                if (mid>=motifStart-0.5 and mid<=motifEnd+0.5) or (motifMid>=Start-0.5 and motifMid<=End+0.5):
                    region2motifs[motif]=1
    line_out+="\t".join([region,a[3],",".join([i+"__"+"#".join(TF_motif[i].keys()) if i in TF_motif else i+'__' for i in a[4].split(",")]),Chr+":"+a[1]+"-"+a[2],states,",".join([str(round(i,4)) for i in scores]),RegulomeDB_score,peak2nearestGene["\t".join([Chr,str(int(a[1])+1),a[2]])],",".join(interactions),str(len(region2TFs)),str(len(region2motifs))])+"\n"

open(fn_prefix+"_merged_"+merge_size+"_anno.txt","w").write(line_out)

print c1
print c2

