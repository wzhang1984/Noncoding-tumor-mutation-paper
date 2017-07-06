
ref_fn='./motif_instances_mappable_50_ref.txt'
alt_fn='./motif_instances_mappable_50_alt.txt'

mut2motif={}
mut2motif={}
for fn in [ref_fn,alt_fn]:
    header=True
    c2i={}
    for line in open(fn):
        a=line.split('\n')[0].split('\t')
        if header:
            if a[0]!='FASTA ID':
                continue
            header=False
            for i in range(len(a)):
                c2i[a[i]]=i
            continue
        mut=a[c2i['FASTA ID']][:-4]
        motif=a[c2i['Motif Name']]
        if mut not in mut2motif:
            mut2motif[mut]=[set(),set()]
        if fn==ref_fn:
            mut2motif[mut][0].add(motif)
        elif fn==alt_fn:
            mut2motif[mut][1].add(motif)

line_out=''
motif2count={}
for mut in sorted(mut2motif):
    motifs_gain=mut2motif[mut][1]-mut2motif[mut][0]
    motifs_lost=mut2motif[mut][0]-mut2motif[mut][1]
    for motif in sorted(motifs_gain):
        line_out+='{}\t{}\tgain\n'.format(mut,motif)
        if motif not in motif2count:
            motif2count[motif]=[0,0]
        motif2count[motif][0]+=1
    for motif in sorted(motifs_lost):
        line_out+='{}\t{}\tlost\n'.format(mut,motif)
        if motif not in motif2count:
            motif2count[motif]=[0,0]
        motif2count[motif][1]+=1
open('./motif_instances_mappable_50_compare.txt','wb').write(line_out)

line_out=''
for motif in sorted(motif2count):
    line_out+='{}\t{}\t{}\n'.format(motif,motif2count[motif][0],motif2count[motif][1])
open('./motif_instances_mappable_50_count.txt','wb').write(line_out)

