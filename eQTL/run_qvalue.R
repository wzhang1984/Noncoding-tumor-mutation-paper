
library(qvalue)

p=read.delim('gene2ncmut_lm_p.txt',header=F,row.names=1)
qobj=qvalue(p[,4], lambda=seq(0.00,0.99,0.01))
$qobj=qvalue(p[,4], lambda=seq(0.00,0.99,0.01), pi0.meth='bootstrap')
p$q=qobj$qvalues
p$fdr=p.adjust(p[,4],method='BH')
write.table(p,'gene2ncmut_lm_p_q.txt',quote=F,sep='\t',col.names=F)



#library(qvalue)
#
#p=read.delim('gene2ncmut_lm_p_allPairs.txt',header=F)
#qobj=qvalue(p[,4], lambda=seq(0.2,0.8,0.01))
#p$q=qobj$qvalues
#p$fdr=p.adjust(p[,4],method='BH')
#write.table(p,'gene2ncmut_lm_p_allPairs_q.txt',quote=F,sep='\t',col.names=F,row.names=F)
#
