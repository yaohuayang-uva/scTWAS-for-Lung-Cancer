rm(list=ls())
library(glmnet)

exp <- read.delim("gene1.exp", header=T, stringsAsFactors=F, sep="\t")
gt <- read.delim("dosage.traw", header=T, row.names=2, stringsAsFactors=F, sep="\t", check.names=F)
snp_info <- gt[, c("COUNTED","ALT")]
gt <- as.data.frame(t(gt[,6:ncol(gt)]), stringsAsFactors=F)
gt <- gt[as.character(exp$Sample),]
table(rownames(gt)==exp$Sample)
gt <- apply(gt, 2, function(x) ifelse(is.na(x),mean(x,na.rm=T),x))

dnam_val <- as.numeric(exp$Exp)
if (shapiro.test(dnam_val)$p.value < 0.05) {dnam_val <- qnorm(rank(dnam_val, ties.method="average")/(length(dnam_val)+1))}
set.seed(123)
fit <- cv.glmnet(gt,dnam_val,nfolds=5,alpha=0.5,keep=T,parallel=F)
fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm))
best.lam <- fit.df[which.min(fit.df[,1]),]; cvm.best = best.lam[,1]; lambda.best = best.lam[,2]; nrow.best = best.lam[,3]
ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best])
ret[ret == 0.0] <- NA; bestbetas = as.vector(ret[which(!is.na(ret)),]); names(bestbetas) = rownames(ret)[which(!is.na(ret))] # vector of non-zero betas
pred.mat <- fit$fit.preval[,nrow.best]

ext_head <- c("Gene", "cvm", "lambda.iter", "lambda.min", "n.snps", "R", "R2", "Pval")
wt_head <- c("Gene", "SNP", "EffectA", "RefA", "beta")
extra_file <- 'extra.txt'; write(ext_head,file=extra_file,ncolumns=8,sep="\t")
weights_file <- 'weights.txt'; write(wt_head,file=weights_file,ncolumns=5,sep="\t")

res <- cor.test(dnam_val, pred.mat); rval <- res$estimate[[1]]; pval <- res$p.value
sum_res <- c("gene1", cvm.best, nrow.best, lambda.best, length(bestbetas), rval, rval**2, pval)
write(sum_res,file=extra_file,ncolumns=8,append=T,sep="\t")

bestbetalist <- names(bestbetas); bestbetainfo <- snp_info[bestbetalist,]; betatable<-as.matrix(cbind(bestbetainfo,bestbetas))
betafile <- cbind("gene1",rownames(betatable),betatable[,1],betatable[,2],betatable[,3]) 
write(t(betafile),file=weights_file,ncolumns=5,append=T,sep="\t")
