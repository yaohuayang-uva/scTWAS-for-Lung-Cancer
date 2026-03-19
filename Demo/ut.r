rm(list=ls())
library(glmnet)
source('glasso.r')

exp <- readRDS("gene1.allcelltype.exp.rds")
Yt <- names(exp); ntune = 5; fold = 5; P = length(Yt); T_num = length(Yt)
gtfile <- read.delim("dosage.traw", header=T, row.names=2, stringsAsFactors=F, check.names=F, sep="\t")
dose <- as.data.frame(t(gtfile[,6:ncol(gtfile)]), stringsAsFactors=F)
dose$sampleid <- rownames(dose); rownames(dose) <- NULL; dose <- dose[, c("sampleid", colnames(dose)[1:ncol(dose)-1])]
dose$id<-dose[,1]; dose<-dose[,c(1,ncol(dose),2:(ncol(dose)-1))]
for(j in 3:ncol(dose)) {
  dose[,j][is.na(dose[,j])] <- mean(dose[,j], na.rm=T)
  dose[,j] = dose[,j] - mean(dose[,j])
}
N = nrow(dose)

tmp = as.matrix(dose[,-(1:2)]); XX = t(tmp)%*%as.matrix(tmp)/N; Xnorm = diag(XX); remove(tmp); remove(XX)
sub_id = dose[,1]; M = ncol(dose) - 2; M = ncol(dose) - 2 #no of snps

sub_id_map = list()  #subject id map for each tissue
for(t in 1:T_num){
  tmp = rep(0, nrow(exp[[t]]))
  for(j in 1:length(tmp)){
    tmp[j] = which(sub_id==exp[[t]][j,1])
  }
  sub_id_map[[t]] = tmp
}

cv_config = cv_helper(N, fold); cv_perm = cv_config$perm; cv_idx = cv_config$idx
single_res_test = list(); single_lam = matrix(0,fold,P); single_theta_est = list()
multi_res_test = list(); multi_lam = matrix(0,fold,2); multi_theta_est = list()
multi_res_test2 = list(); multi_lam2 = array(0, dim=c(fold, P, 2)); multi_theta_est2 = list()
res_tune = list(); rec_lamv = matrix(0, fold, ntune)
avg_tune_res<-list() 
single_initial_est<-list()

for(f in 1:fold){
  test_index = cv_perm[cv_idx[f,1]:cv_idx[f,2]]
  test_id = sub_id[test_index] 
  tuning_index = cv_perm[cv_idx[f%%fold+1,1]:cv_idx[f%%fold+1,2]]
  tuning_id = sub_id[tuning_index] 

  X_test = list(); Y_test = list(); X_tune = list(); Y_tune = list(); X_train = list(); Y_train = list(); X_all = list(); Y_all = list()

  for(t in 1:T_num){
    X_train_tmp = sub_id_map[[t]][!(sub_id_map[[t]]%in%c(test_index,tuning_index))]
    Y_train_tmp = !(sub_id_map[[t]]%in%c(test_index,tuning_index))
    X_tuning_tmp = sub_id_map[[t]][(sub_id_map[[t]]%in%tuning_index)]
    Y_tuning_tmp = (sub_id_map[[t]]%in%tuning_index)
    X_test_tmp = sub_id_map[[t]][(sub_id_map[[t]]%in%test_index)]
    Y_test_tmp = (sub_id_map[[t]]%in%test_index)
    X_train[[t]] = apply(as.matrix(dose[X_train_tmp,-c(1,2)]),2,as.numeric)
    Y_train[[t]] = exp[[t]][Y_train_tmp, 2]
    X_tune[[t]] = apply(as.matrix(dose[X_tuning_tmp,-c(1,2)]),2,as.numeric)
    Y_tune[[t]] = exp[[t]][Y_tuning_tmp, 2]
    X_test[[t]] = apply(as.matrix(dose[X_test_tmp,-c(1,2)]),2,as.numeric)
    Y_test[[t]] = exp[[t]][Y_test_tmp, 2]
    X_all_tmp = sub_id_map[[t]]
    X_all[[t]] = apply(as.matrix(dose[X_all_tmp,-c(1,2)]),2,as.numeric)
    Y_all[[t]] = exp[[t]][,2]
    
    X_train[[t]]=rbind(X_train[[t]],X_test[[t]])
    Y_train[[t]]=c(Y_train[[t]],Y_test[[t]])
  }
    if (f==1){
      single_initial_est_all = matrix(0, ncol(X_train[[1]]), T_num) 
      single_summary_all = list()
      for(t in 1:T_num){
        set.seed(123)
        tt = cv.glmnet(X_all[[t]], Y_all[[t]], alpha = 0.5, nfolds = 5) 
        single_summary_all[[t]] = tt
        single_initial_est_all[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)] 
      }
    }

    single_initial_est[[f]] = matrix(0, ncol(X_train[[1]]), T_num) 
    for(t in 1:T_num){
    set.seed(123)
    tt = cv.glmnet(X_train[[t]], Y_train[[t]], alpha = 0.5, nfolds = 5) 
    single_initial_est[[f]][,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)] 
  }
 lam_range = minmax_lambda(single_summary_all) #lambda range
  sig_norm = apply(single_initial_est[[f]], 1, function(x){sqrt(sum(x^2))}) 
  sig_norm[sig_norm==0] = rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2 
  sig_norm = sig_norm/sum(sig_norm) 
  weights2 = 1/sig_norm; weights2 = weights2/sum(weights2);

  tis_norm = apply(single_initial_est[[f]], 2, function(x){sum(abs(x))})
  tis_norm[tis_norm==0] = rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
  tis_norm = tis_norm/sum(tis_norm)
  weights1 = 1/tis_norm; weights1 = weights1/sum(weights1);
  lam_V = seq(lam_range[1], lam_range[2], length.out = ntune)

  initial_numeric = as.numeric(single_initial_est[[f]])

  XY = grad_prep(X_train, Y_train)
  XX_train = lapply(X_train, function(x){t(x)%*%x/nrow(x)})
  spsz = unlist(lapply(X_train,nrow)) 
  res_tune[[f]] = array(NA, dim=c(ntune, ntune, P)) 

  rec_lamv[f,] = lam_V  

  for(lam1 in 1:ntune){
    for(lam2 in 1:ntune){
      single_est = matrix(initial_numeric, M, P) 
      ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[lam1]/spsz, lambda2=lam_V[lam2], theta=single_est)

      if (lam1==1 & lam2==1){
        multi_res_test[[f]]<-list()
      }
      multi_res_test[[f]][[paste0(lam1,'_',lam2)]] = multi_mse(ans$est, X_tune, Y_tune)
      if(sum(ans$est!=0)>0){
        res_tune[[f]][lam1,lam2, ] = ans$tune_err
        remove(single_est); remove(ans);
      }else{
        res_tune[[f]][lam1,lam2, ] = ans$tune_err
        remove(single_est); remove(ans);
      }
    }
  }
  avg_tune_res[[f]] = apply(res_tune[[f]], c(1,2), mean) 
}
avg_tune_cross_fold<-matrix(rep(1,ntune*ntune),nrow = ntune,ncol=ntune)
  for (i_tune in 1:ntune){
    for (j_tune in 1:ntune){
      avg_tune_cross_fold[i_tune,j_tune]=mean(c(avg_tune_res[[1]][i_tune,j_tune],avg_tune_res[[2]][i_tune,j_tune],avg_tune_res[[3]][i_tune,j_tune],avg_tune_res[[4]][i_tune,j_tune],avg_tune_res[[5]][i_tune,j_tune]),na.rm = T)
      avg_tune_cross_fold[i_tune,j_tune]<-ifelse(is.na(avg_tune_cross_fold[i_tune,j_tune]),1,avg_tune_cross_fold[i_tune,j_tune])
    }
  }

  best.lam = which(avg_tune_cross_fold == min(avg_tune_cross_fold), arr.ind = TRUE)[1,]

  multi_res<-list()
  for(f in 1:ntune){
    multi_res[[f]]<-multi_res_test[[f]][[paste0(as.numeric(best.lam)[1],'_',as.numeric(best.lam)[2])]]
  }

cv_df<-list()
  cv_r<-cv_p<-c()
  for (tissue_i in 1:T_num){
    cv_df[[tissue_i]]<-multi_res[[1]][tissue_i][[1]]

    for(fold_i in 2:fold){
      cv_df[[tissue_i]]<-rbind(cv_df[[tissue_i]],multi_res[[fold_i]][tissue_i][[1]])
    }
    fit<-cor.test(cv_df[[tissue_i]][,1],cv_df[[tissue_i]][,2])
    cv_r[tissue_i]<-as.numeric(fit$estimate)
    cv_p[tissue_i]<-as.numeric(fit$p.value)

    cv_r[tissue_i]<-ifelse(is.na(cv_r[tissue_i]),0,cv_r[tissue_i])
    cv_p[tissue_i]<-ifelse(is.na(cv_p[tissue_i]),1,cv_p[tissue_i])
  }

  XY = grad_prep(X_all, Y_all)
  XX_all = lapply(X_all, function(x){t(x)%*%x/nrow(x)})
  ans = glasso_no_early_stopping(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=lam_V[best.lam[1]]/spsz, lambda2=lam_V[best.lam[2]], theta=single_initial_est_all)

info <- gtfile[, c("ALT","COUNTED")]; colnames(info) <- c("ref_allele","counted_allele")
  info$chr_bp <- paste(gtfile$CHR, gtfile$POS, sep="_"); info$rsid <- rownames(info)
  rownames(info) <- NULL; info <- info[, c("rsid","chr_bp","ref_allele","counted_allele")]
  combined_res <- rep(list(NULL), length(Yt)) 
  for (tissue_i in 1:length(Yt)){
    tissue_name=Yt[tissue_i]
    weight_df<-info
    weight_df$gene="gene1"
    weight_df$weight=ans$est[,tissue_i]
    #weight_df$r2=round((cv_r[tissue_i])^2,3)
    weight_df$r=cv_r[tissue_i]; weight_df$r2=(cv_r[tissue_i])^2
    weight_df$p=cv_p[tissue_i]
    weight_df$lambda=paste0(round(lam_V[best.lam[1]],2),';',round(lam_V[best.lam[2]],2))
    weight_df$iter=ans$iter
    weight_df<-weight_df[,c('gene','rsid','chr_bp','ref_allele','counted_allele','weight','r','r2','p','lambda','iter')]
    weight_df<-weight_df[weight_df$weight!=0,]
    if(nrow(weight_df)>0){
      combined_res[[tissue_i]] <- weight_df 
    }
  }
names(combined_res) <- Yt
saveRDS(combined_res, file="ut.res.rds")

