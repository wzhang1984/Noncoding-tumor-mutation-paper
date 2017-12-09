
library(survival)


sink("survival/coxph.txt")

d=read.delim("survival/pat2surv2labels.txt",header=T,row.names=1)

d=d[which(d$K10!="NaN"),]
d=d[which(d$DFS_MONTH!="NA"),]
d=d[which(d$DFS_MONTH>0),]

d$DISEASE=droplevels(d$DISEASE)

H0_OS=tryCatch(coxph(Surv(OS_MONTHS,OS_STATUS)~DISEASE,data=d), error = function(e) e)
H0_DFS=tryCatch(coxph(Surv(DFS_MONTHS,DFS_STATUS)~DISEASE,data=d), error = function(e) e)

# output file format: K, OS_simple, DFS_simple, OS_strata, DFS_strata, OS_complete, DFS_complete
out=as.data.frame(matrix(ncol=7))
v=colnames(d)[c(-1,-2,-3,-4,-5,-6)]
for (i in 1:length(v)) {
    d2=d[,c(1,2,3,4,5,6,i+6)]
    d2$cluster=factor(d2[[v[i]]])
    d2$cluster_relevel <- relevel(d2$cluster, ref = names(which.max(summary(d2$cluster))))
    cr_OS_simple=tryCatch(coxph(Surv(OS_MONTHS,OS_STATUS)~cluster_relevel,data=d2), error = function(e) e)
    cr_DFS_simple=tryCatch(coxph(Surv(DFS_MONTHS,DFS_STATUS)~cluster_relevel,data=d2), error = function(e) e)
    cr_OS=tryCatch(coxph(Surv(OS_MONTHS,OS_STATUS)~DISEASE+cluster_relevel,data=d2), error = function(e) e)
    cr_DFS=tryCatch(coxph(Surv(DFS_MONTHS,DFS_STATUS)~DISEASE+cluster_relevel,data=d2), error = function(e) e)
    out[i,1]=v[i]
    out[i,2]=summary(cr_OS_simple)$logtest[3]
    out[i,3]=summary(cr_DFS_simple)$logtest[3]
    out[i,4]=summary(cr_OS)$logtest[3]
    out[i,5]=summary(cr_DFS)$logtest[3]
    a=anova(H0_OS,cr_OS)
    df=a$Df[2]
    p=pchisq(-2*(H0_OS$loglik[2]-cr_OS$loglik[2]),df,lower.tail=FALSE)
    out[i,6]=p
    a=anova(H0_DFS,cr_DFS)
    df=a$Df[2]
    p=pchisq(-2*(H0_DFS$loglik[2]-cr_DFS$loglik[2]),df,lower.tail=FALSE)
    out[i,7]=p
    pdf(paste("survival/",v[i],"_OS.pdf"))
    fit <- survfit(Surv(OS_MONTHS,OS_STATUS)~cluster, data=d2)
    print(fit)
    print("#############################################################################################################################################################################################")
    plot(fit,col=rainbow(i+1),xlab="Overall survival (months)",ylab="Survival probability",xlim=c(0,200),mark.time=TRUE)
    legend("topright",levels(factor(d2[,8])),col=rainbow(i+1),lty=rep(1,i+1))
    dev.off()
    pdf(paste("survival/",v[i],"_DFS.pdf"))
    fit <- survfit(Surv(DFS_MONTHS,DFS_STATUS)~cluster, data=d2)
    print(fit)
    print("#############################################################################################################################################################################################")
    plot(fit,col=rainbow(i+1),xlab="Disease free survival (months)",ylab="Survival probability",xlim=c(0,200),mark.time=TRUE)
    legend("topright",levels(factor(d2[,8])),col=rainbow(i+1),lty=rep(1,i+1))
    dev.off()
}
write.table(out,"survival/coxph_FTest.txt",sep="\t",row.name=F, quote =F,col.name=F)

sink()
