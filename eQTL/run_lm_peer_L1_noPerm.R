

args <- commandArgs(trailingOnly = TRUE)

fn=args[1]

#fn='data4lm_peer/1255.txt'

library(glmnet)
library(methods)

d=read.delim(fn,header=T,row.names=1)
d=t(d)
x=subset(d, select=-c(exp))
y=d[,'exp']

factors=c()
RMEs=c()
for (i in colnames(x)) {
    if (grepl('chr', i)) {
        RMEs=c(RMEs,i)
    } else {
        factors=c(factors,i)
    }
}

out=data.frame(coef=numeric(),p=numeric())

cvfit = cv.glmnet(x, y)

#mat = coef(cvfit, s = cvfit$lambda.min)
#summ <- summary(mat)
#df_coef = data.frame(feature = rownames(mat)[summ$i], coef = summ$x)
#
#factors_L1 = c()
#RMEs_L1 = c()
#for (i in df_coef$feature) {
#    if (i=="(Intercept)") {
#        next
#    }
#    if (grepl('chr', i)) {
#        RMEs_L1 = c(RMEs_L1,i)
#    } else {
#        factors_L1 = c(factors_L1,i)
#    }
#}
#
#if (length(RMEs_L1)==0) {
#    fn_out = paste('./lm_coef_p/',strsplit(fn,'/')[[1]][2],sep='')
#    write.table('0\t1\t1', fn_out, quote = F, sep = '\t', col.names=F)
#    fn_out = paste('./lm_model_p/',strsplit(fn,'/')[[1]][2],sep='')
#    write.table('1', fn_out, quote = F, sep = '\t', col.names=F, row.names=F)
#    stop('0 RMEs')
#}

for (lambda in cvfit$lambda) {
    if (lambda>cvfit$lambda.min) {next}
    mat = coef(cvfit, s = lambda)
    summ <- summary(mat)
    df_coef = data.frame(feature = rownames(mat)[summ$i], coef = summ$x)

    factors_L1 = c()
    RMEs_L1 = c()
    for (i in df_coef$feature) {
        if (i=="(Intercept)") {
            next
        }
        if (grepl('chr', i)) {
            RMEs_L1 = c(RMEs_L1,i)
        } else {
            factors_L1 = c(factors_L1,i)
        }
    }
    if (length(RMEs_L1)>0) {break}
}

if (length(RMEs_L1)==0) {
    RMEs_L1=RMEs
}


x_L1 = x[,c(RMEs_L1,factors_L1)]

df=data.frame(x_L1,y)
fit=lm(y~.,data=df)
s=summary(fit)
coef=data.frame(s$coefficients)
out=coef[2:(length(RMEs_L1)+1),c(1,4)]
colnames(out)=c('coef','p')
rownames(out)=RMEs_L1

out$fdr=p.adjust(out$p,method='BH')

fn_out = paste('./lm_coef_p/',strsplit(fn,'/')[[1]][2],sep='')
write.table(out, fn_out, quote = F, sep = '\t', col.names=F)



x_factors_L1 = x[,factors_L1]
df=data.frame(x_factors_L1,y)
fit0=lm(y~.,data=df)
a=anova(fit0,fit)
fn_out = paste('./lm_model_p/',strsplit(fn,'/')[[1]][2],sep='')
write.table(a$'Pr(>F)'[2], fn_out, quote = F, sep = '\t', col.names=F, row.names=F)

