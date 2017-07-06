library(peer)

expr = read.delim('../pat2exp_log.txt', header=T, row.names=1)
expr_zscore = apply(expr, 1, function(y) (y-mean(y)) / sd(y))

model = PEER()

PEER_setPhenoMean(model,as.matrix(expr_zscore))

PEER_setNk(model,35)

PEER_setAdd_mean(model, TRUE)

covs = read.delim('../covariates.txt',header=T,row.names=1)
covs_matrix=t(matrix(as.numeric(unlist(covs)),nrow=nrow(covs)))
PEER_setCovariates(model, covs_matrix)

PEER_update(model)

factors = PEER_getX(model)
rownames(factors)=rownames(expr_zscore)
colnames(factors)=c(rownames(covs),'mean',paste("h", seq(1:35), sep = ""))

residuals = PEER_getResiduals(model)
rownames(residuals)=rownames(expr_zscore)
colnames(residuals)=colnames(expr_zscore)

write.table(t(factors),'factors_35.txt',sep='\t',quote=F,col.names=F)
write.table(t(residuals),'residuals_35.txt',sep='\t',quote=F,col.names=F)
write.table(t(expr_zscore),'expr_zscore.txt',sep='\t',quote=F,col.names=F)

alpha = PEER_getAlpha(model)
rownames(alpha)=colnames(factors)
write.table(alpha,'precision_35.txt',quote=F,sep='\t',row.names=T,col.names=F)
alpha=alpha[22:length(alpha)]
pdf('PEER_plotModel_35.pdf')
plot(alpha,xlab="Factors",ylab="Inverse variance of factor weights", type="b", col="blue", lwd=4, xaxp=c(1,length(alpha), length(alpha)-1))
dev.off()
