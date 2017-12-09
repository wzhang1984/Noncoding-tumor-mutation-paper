
library(qvalue)

p=read.delim('gene2ncmut_lm_p.txt',header=F,row.names=1)
qobj=qvalue(p[,4])
p$q=qobj$qvalues
p$fdr=p.adjust(p[,4],method='BH')
write.table(p,'gene2ncmut_lm_p_q.txt',quote=F,sep='\t',col.names=F)
