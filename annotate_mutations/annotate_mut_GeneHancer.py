
# Annotate noncoding mutations using GeneHancer
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5467550/

import sys

fn_prefix = sys.argv[1]
merge_size = sys.argv[2]
occurrence_cutoff = int(sys.argv[3])

binsize = 10000

print 'Reading GeneHancer_hg19'
bin2GeneHancer = {}
fn = '/cellar/data/users/wzhang1984/GeneHancer/GeneHancer_hg19.bed'
for line in open(fn).read().rstrip().splitlines():
    row = line.split('\t')
    Chr = row[0]
    Start = int(row[1])
    End = int(row[2])
    Pos = '\t'.join(row[1:3])
    genes = []
    for entry in row[3].split(';'):
        entry_split = entry.split('connected_gene=')
        if len(entry_split) == 2:
            genes.append(entry_split[1])
    binstart = Start / binsize
    binend = End / binsize
    for binpos in range(binstart, binend+1):
        Bin = '{}:{}'.format(Chr, binpos)
        if Bin not in bin2GeneHancer:
            bin2GeneHancer[Bin] = {}
        if Pos not in bin2GeneHancer[Bin]:
            bin2GeneHancer[Bin][Pos] = genes

print 'Reading Homer TSS annotation'
locus2promoter = {}
header = True
for line in open("./" + fn_prefix + "_merged_" + merge_size + "_homerAnno.txt"):
    if header:
        header = False
        continue
    row = line.split("\t")
    if row[9] == 'NA':
        continue
    dist2tss = int(row[9])
    if dist2tss>=-1000 and dist2tss<=1000:
        locus2promoter["\t".join([row[1], str(int(row[2])-1), row[3]])] = row[15]


print 'Reading noncoding mutation loci'
nline = 0
line_out = ""
for line in open("./" + fn_prefix + "_merged_" + merge_size + ".txt").read().rstrip().splitlines():
    row = line.split("\t")
    nline += 1
    if nline%1000 == 0:
        print nline
    if int(row[3]) < occurrence_cutoff:
        continue
    row = line.split("\t")
    Chr = row[0]
    Start = int(row[1]) - 100
    End = int(row[2]) + 100
    mid = (Start+End) / 2.0
    locus = "\t".join(row[:3])
    genes = set()
    if locus in locus2promoter:
        genes.add(locus2promoter[locus] + '|TSS')

    binstart = Start / binsize
    binend = End / binsize
    for binpos in range(binstart, binend+1):
        Bin = '{}:{}'.format(Chr, binpos)
        if Bin in bin2GeneHancer:
            for Pos in bin2GeneHancer[Bin]:
                [genehancerstart,genehancerend] = [int(i) for i in Pos.split("\t")]
                genehancermid = (genehancerstart+genehancerend) / 2.0
                if (genehancermid>=Start-0.5 and genehancermid<=End+0.5) or (mid>=genehancerstart-0.5 and mid<=genehancerend+0.5):
                    for gene in bin2GeneHancer[Bin][Pos]:
                        genes.add(gene + '|enhancer')

    if genes != set():
        line_out += line + '\t' + ','.join(genes) + '\n'

open(fn_prefix+"_merged_"+merge_size+"_anno_promoter_genehancer.txt","w").write(line_out)

