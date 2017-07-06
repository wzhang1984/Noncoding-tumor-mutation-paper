
info=[]
for line in open('./info_mappable_50.txt').read().rstrip().split('\n'):
    a=line.split('\t')
    info.append(a[0])

line_out=''
line_out2=''
seq=''
index=0
for line in open('./seqList_mappable_50.fa'):
    if line[0]=='>':
        if seq:
            replace=seq[7:-7]
            if replace!=ref:
                print header
            seq_alt=seq[:7]+alt+seq[-7:]
            line_out+='>'+header+'_ref'+'\n'+seq+'\n'
            line_out2+='>'+header+'_alt'+'\n'+seq_alt+'\n'
        header=info[index]
        ref_alt=header.split('_')[1]
        [ref,alt]=ref_alt.split('>')
        index+=1
        seq=''
        if index/1000000==index/1000000.0:
            print index
    else:
        seq+=line.split('\n')[0]

if seq:
    replace=seq[7:-7]
    if replace!=ref:
        print header
    seq_alt=seq[:7]+alt+seq[-7:]
    line_out+='>'+header+'_ref'+'\n'+seq+'\n'
    line_out2+='>'+header+'_alt'+'\n'+seq_alt+'\n'

open('./seqList_mappable_50_ref.fa','wb').write(line_out)
open('./seqList_mappable_50_alt.fa','wb').write(line_out2)

