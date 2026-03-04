rm(list=ls())
pacman::p_load(glmnet, foreach, doParallel, readr)
args = commandArgs(trailingOnly=TRUE)
#mycores <- detectCores() - 2
mycores <- 40

gt_fd = paste0('/project/yanglab/scTWAS/Geno/', args[1])
exp_f = paste0('/project/yanglab/scTWAS/Residuals/', args[2])
res_folder <- paste0("/project/yanglab/scTWAS/Output/", args[3])
cis_reg <- NULL
if (args[3] == "100K") {cis_reg <- 100000}
if (args[3] == "500K") {cis_reg <- 500000}
if (args[3] == "Mega") {cis_reg <- 1000000}

gt_f <- paste(gt_fd, '/geno.traw', sep='')
gt <- read.delim(gt_f, header=T, stringsAsFactors=F, row.names=2, check.names=F, sep="\t")
snp_info <- gt[, c("CHR","POS","COUNTED","ALT")]
gt <- data.matrix(t(gt[,6:ncol(gt)])); rownames(gt) <- as.character(t(data.frame(strsplit(rownames(gt),"_"), stringsAsFactors=F))[,2])

peak <- read.delim(exp_f, header=T, stringsAsFactors=F, row.names=1, check.names=F, sep="\t")
info <- readRDS("/project/yanglab/scTWAS/Residuals/info.rds")
op_gene <- intersect(rownames(peak), rownames(info)); peak <- peak[op_gene,]; info <- info[op_gene,]; table(rownames(peak)==rownames(info))

peak <- peak[,rownames(gt)]; table(rownames(gt)==colnames(peak))

gt <- apply(gt, 2, function(x) ifelse(is.na(x),mean(x,na.rm=T),x))

ext_head <- c("CpG", "cvm", "lambda.iter", "lambda.min", "n.snps", "R", "R2", "Pval")
wt_head <- c("CpG", "SNP", "EffectA", "RefA", "beta"); cov_head <- c('CpG', 'RSID2','RSID1','VALUE')

func <- function (x, m=peak, g=gt, sinfo=snp_info, mc=mycores, cis=cis_reg) {
  every <- ceiling(nrow(m)/mc); st <- (x-1) * every + 1; en <- st + every - 1
  if (en > nrow(m)) { en <- nrow(m) }
  extra_file <- paste(res_folder, "/", args[1], '__part_', x, '.extra.txt', sep=''); write(ext_head,file=extra_file,ncolumns=8,sep="\t")
  weights_file <- paste(res_folder, "/", args[1], '__part_', x, '.weights.txt', sep=''); write(wt_head,file=weights_file,ncolumns=5,sep="\t")
  cov_file <- paste(res_folder, "/", args[1], '__part_', x, '.cov.txt', sep=''); write(cov_head,file=cov_file,ncolumns=4,sep="\t")

  for (j in st:en) {
    dnam_val <- as.numeric(m[j,])
    if (shapiro.test(dnam_val)$p.value < 0.05) {dnam_val <- qnorm(rank(dnam_val, ties.method="average")/(length(dnam_val)+1))}
    dnam_val <- scale(dnam_val, center=T, scale=T); dnam_val[is.na(dnam_val)] <- 0; rownames(dnam_val) <- colnames(m)
    peak_info <- info[rownames(m)[j],]; peak_chr <- as.numeric(peak_info$V4); left <- as.numeric(peak_info$V5) - cis; right <- as.numeric(peak_info$V6) + cis
    snps <- rownames(subset(sinfo, CHR == peak_chr & POS >= left & POS <= right)); geno <- g[,snps]
    if(is.null(dim(geno)) | (!is.null(dim(geno)) && dim(geno)[2] == 0)) {next}
    set.seed(as.numeric(sub("ENSG","",rownames(peak_info))))
    isbad <- FALSE
    tryCatch({fit <- cv.glmnet(geno,dnam_val,nfolds=5,alpha=0.5,keep=T,parallel=F)}, error=function(e){isbad <<- TRUE })
    if (isbad) {next}
    fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm))
    best.lam <- fit.df[which.min(fit.df[,1]),]; cvm.best = best.lam[,1]; lambda.best = best.lam[,2]; nrow.best = best.lam[,3]
    ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best]) # get betas from best lambda
    ret[ret == 0.0] <- NA; bestbetas = as.vector(ret[which(!is.na(ret)),]); names(bestbetas) = rownames(ret)[which(!is.na(ret))] # vector of non-zero betas
    pred.mat <- fit$fit.preval[,nrow.best]
    if(length(bestbetas) > 0) {
      res <- cor.test(dnam_val, pred.mat); rval <- res$estimate[[1]]; pval <- res$p.value
      #if (rval>0.1 & pval<0.05) {
        sum_res <- c(rownames(m)[j], cvm.best, nrow.best, lambda.best, length(bestbetas), rval, rval**2, pval)
        write(sum_res,file=extra_file,ncolumns=8,append=T,sep="\t") # write extra results

        bestbetalist <- names(bestbetas); bestbetainfo <- sinfo[bestbetalist,]; betatable<-as.matrix(cbind(bestbetainfo,bestbetas))
        betafile <- cbind(rownames(m)[j],rownames(betatable),betatable[,3],betatable[,4],betatable[,5]) ##output "CpG","SNP","effectAllele","refAllele","beta"
        write(t(betafile),file=weights_file,ncolumns=5,append=T,sep="\t") # t() necessary for correct output from write() function

        if (length(bestbetalist)>1) {
          dsg <- g[,bestbetalist]; cov = cov(as.matrix(dsg)); cov[upper.tri(cov)] <- NA; cov = cbind(expand.grid(dimnames(cov)), value = as.vector(cov))
          colnames(cov) <- c('RSID2','RSID1','VALUE'); cov = cov[!is.na(cov$VALUE),]; cov$GENE <- rownames(m)[j]; cov = cov[,c('GENE','RSID1','RSID2','VALUE')]
          write(t(cov),file=cov_file,ncolumns=4,append=T,sep="\t") # write cov results
        }     
      #}
    }
  }
}

cl <- makeCluster(mycores); registerDoParallel(cl)
x <- foreach(x=1:mycores,.packages='glmnet') %dopar% func(x)
stopCluster(cl)
