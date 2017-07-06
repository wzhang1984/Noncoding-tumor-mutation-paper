
mut2gene={}
mut2MI={}
for line in open('./info_mappable_50.txt').read().rstrip().split('\n'):
    a=line.split('\t')
    if a[0] not in mut2gene:
        mut2gene[a[0]]=set()
    for gene in a[1].split(','):
        mut2gene[a[0]].add('|'.join(gene.split('|')[:2]))
    mut2MI[a[0]]=a[2]


gene2motif2info={}
for line in open('./motif_instances_mappable_50_compare.txt').read().rstrip().split('\n'):
    a=line.split('\t')
    mut=a[0]
    MI=mut2MI[mut]
    [mut_loc,alt,pat]=mut.split('_')
    motif=a[1]
    for gene in mut2gene[mut]:
        key='{}\t{}\t{}\t{}'.format(gene,MI,motif,a[2])
        if key not in gene2motif2info:
            gene2motif2info[key]={'pats':set(),'muts':set(),'MIs':set()}
        gene2motif2info[key]['pats'].add(pat)
        gene2motif2info[key]['muts'].add(mut_loc)
        # gene2motif2info[key]['MIs'].add(MI)

line_out=''
for key in gene2motif2info:
    line_out+='{}\t{}\t{}\t{}\t{}\n'.format(key,len(gene2motif2info[key]['pats']),len(gene2motif2info[key]['muts']),','.join(gene2motif2info[key]['pats']),','.join(gene2motif2info[key]['muts']))

open('summarize_instances_50_compare_MI.txt','wb').write(line_out)
