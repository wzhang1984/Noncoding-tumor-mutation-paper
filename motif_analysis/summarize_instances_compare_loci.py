
mut2gene={}
mut2locus={}
for line in open('./info_mappable_50.txt').read().rstrip().split('\n'):
    a=line.split('\t')
    if a[0] not in mut2gene:
        mut2gene[a[0]]=set()
    for gene in a[1].split(','):
        mut2gene[a[0]].add(gene.split('|')[0])
    mut2locus[a[0]]=a[2]


gene2motif2info={}
for line in open('./motif_instances_mappable_50_compare.txt').read().rstrip().split('\n'):
    a=line.split('\t')
    mut=a[0]
    locus=mut2locus[mut]
    [mut_loc,alt,pat]=mut.split('_')
    motif=a[1]
    for gene in mut2gene[mut]:
        key='{}\t{}\t{}\t{}'.format(gene,locus,motif,a[2])
        if key not in gene2motif2info:
            gene2motif2info[key]={'pats':set(),'muts':set(),'locuss':set()}
        gene2motif2info[key]['pats'].add(pat)
        gene2motif2info[key]['muts'].add(mut_loc)

line_out=''
for key in gene2motif2info:
    line_out+='{}\t{}\t{}\t{}\t{}\n'.format(key,len(gene2motif2info[key]['pats']),len(gene2motif2info[key]['muts']),','.join(gene2motif2info[key]['pats']),','.join(gene2motif2info[key]['muts']))

open('summarize_instances_50_compare_loci.txt','wb').write(line_out)
