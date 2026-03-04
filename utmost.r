args=as.numeric(commandArgs(TRUE)) #subjob id
options(stringsAsFactors=F)
#pacman::p_load(glmnet,foreach,stringr, doParallel, readr)
library(glmnet); library(foreach); library(stringr); library(doParallel); library(readr)


allm <- readRDS(paste0("/scratch/vta8we/scTWAS/UTMOST/Exp/chr",args[1],"_rn.rds")) # No NA are allowed in exp values!
chromposition <- readRDS("/scratch/vta8we/scTWAS/UTMOST/Exp/info.rds")
mycores <- 20

func <- function (x, allmet=allm, chrpos=chromposition, mc=mycores) {
source('/scratch/vta8we/CPTAC/GLASSO/glasso.r')
gploc <- "/scratch/vta8we/scTWAS/UTMOST/Geno/chr"
every <- ceiling(nrow(allmet[[1]])/mc); st <- (x-1) * every + 1; en <- st + every - 1
if (st <= nrow(allmet[[1]])) {
if (en > nrow(allmet[[1]])) { en <- nrow(allmet[[1]]) }

cpgs <- c(); for (ii in 1:length(allmet)) {cpgs <- c(cpgs, rownames(allmet[[ii]])[st:en])}; cpgs <- unique(cpgs)
#scss <- read.table("/nobackup/sbcs/yangy23/Genetics_LungCA/DNAm/Model/UTMOST/success.cpgs", header=F, stringsAsFactors=F); scss <- as.character(scss$V1)

for (gene_i in cpgs){
#  if (gene_i %in% scss) {next}
  chrom <- as.numeric(chrpos[gene_i,]$V4); mapinfo_st <- as.numeric(chrpos[gene_i,]$V5); mapinfo_en <- as.numeric(chrpos[gene_i,]$V6)
  st_pos <- mapinfo_st-100000; en_pos <- mapinfo_en+100000
  dir.create(paste0('/scratch/vta8we/scTWAS/UTMOST/Model/100K/Processed_Met/', gene_i))
  plink_cmd <- paste("plink2 --bfile ", gploc, chrom, " --chr ", chrom,  " --from-bp ", st_pos, " --to-bp ", en_pos, " --export A-transpose --out /scratch/vta8we/scTWAS/UTMOST/Model/100K/Processed_Met/", gene_i, "/geno ", "--threads 1 --memory 800 ", sep="")
  system(plink_cmd,ignore.stdout=T,ignore.stderr=T,wait=T)
  tmp_files <- system(paste("ls /scratch/vta8we/scTWAS/UTMOST/Model/100K/Processed_Met/", gene_i, sep=""), intern=T)
  geno_file <- paste("/scratch/vta8we/scTWAS/UTMOST/Model/100K/Processed_Met/", gene_i, "/geno.traw", sep="")
  clean_cmd <- paste("rm -rf /scratch/vta8we/scTWAS/UTMOST/Model/100K/Processed_Met/", gene_i, sep="")
  
  if (!("geno.traw" %in% tmp_files)) {system(clean_cmd, ignore.stdout=T,ignore.stderr=T,wait=T); next}

  tissue_list <- c()
  for (jj in 1:length(allmet)) {
    if (gene_i %in% rownames(allmet[[jj]])) {tissue_list <- c(tissue_list, names(allmet)[jj])}
  } 
  Yt <- tissue_list; ntune = 5; fold = 5
  
  P = length(Yt)
  if(P==0){system(clean_cmd, ignore.stdout=T,ignore.stderr=T,wait=T); next}
  
  ## expr files ##
  Y = list()
  for(t in Yt){
    tmp_exp <- as.data.frame(t(allmet[[t]][gene_i,]), stringsAsFactors=F); tmp_exp$tissue <- rep(t, nrow(tmp_exp))
    tmp_exp$sampleid <- rownames(tmp_exp); colnames(tmp_exp)[1] <- "exp"; tmp_exp <- tmp_exp[, c("sampleid","exp","tissue")]
    rownames(tmp_exp) <- NULL; tmp_exp$exp <- as.numeric(tmp_exp$exp); Y[[t]] <- tmp_exp
  }
  ssize = unlist(lapply(Y, nrow))  #N of samples for each tissue
  T_num = length(Yt)  #N of tissues
  
  ##Load genotype files
  gtfile <- read.delim(geno_file, header=T, row.names=2, stringsAsFactors=F, check.names=F, sep="\t")
  system(clean_cmd, ignore.stdout=T,ignore.stderr=T,wait=T)
  dose <- as.data.frame(t(gtfile[,6:ncol(gtfile)]), stringsAsFactors=F)
  tmp_rnames <- c()
  for (name in rownames(dose)) {tmp_rnames <- c(tmp_rnames, strsplit(name,"_")[[1]][2])} # NEED TO CHECK!
  dose$sampleid <- tmp_rnames; rownames(dose) <- NULL; dose <- dose[, c("sampleid", colnames(dose)[1:ncol(dose)-1])]

  if (ncol(dose) <= 2) {system(clean_cmd, ignore.stdout=T,ignore.stderr=T,wait=T); next}
  dose$id<-dose[,1]; dose<-dose[,c(1,ncol(dose),2:(ncol(dose)-1))]
  # Impute NA; center dosage to 0
  for(j in 3:ncol(dose)) {
    dose[,j][is.na(dose[,j])] <- mean(dose[,j], na.rm=T)
    dose[,j] = dose[,j] - mean(dose[,j])
  }
  N = nrow(dose)
  
  ## covariance matrix ##
  tmp = as.matrix(dose[,-(1:2)]); XX = t(tmp)%*%as.matrix(tmp)/N; Xnorm = diag(XX); remove(tmp); remove(XX)
  sub_id = dose[,1]; M = ncol(dose) - 2; M = ncol(dose) - 2 #no of snps
  
  sub_id_map = list()  #subject id map for each tissue
  for(t in 1:T_num){
    tmp = rep(0, nrow(Y[[t]]))
    for(j in 1:length(tmp)){
      tmp[j] = which(sub_id==Y[[t]][j,1])
    }
    sub_id_map[[t]] = tmp
  }
  
  cv_config = cv_helper(N, fold); cv_perm = cv_config$perm; cv_idx = cv_config$idx
  
  single_res_test = list(); single_lam = matrix(0,fold,P); single_theta_est = list()
  
  multi_res_test = list(); multi_lam = matrix(0,fold,2); multi_theta_est = list()
  
  multi_res_test2 = list(); multi_lam2 = array(0, dim=c(fold, P, 2)); multi_theta_est2 = list()
  
  res_tune = list(); rec_lamv = matrix(0, fold, ntune)
  
  avg_tune_res<-list() #average tuning error
  single_initial_est<-list()
  
  #----------------------
  
  for(f in 1:fold){
    print(f)
    #bgt = Sys.time()
    test_index = cv_perm[cv_idx[f,1]:cv_idx[f,2]]
    test_id = sub_id[test_index] #sample id in the testing set
    tuning_index = cv_perm[cv_idx[f%%fold+1,1]:cv_idx[f%%fold+1,2]]
    tuning_id = sub_id[tuning_index] #sample in the 'early stop' set
    
    #dfs for each set
    X_test = list()
    Y_test = list()
    X_tune = list()
    Y_tune = list()
    X_train = list()
    Y_train = list()
    X_all = list()
    Y_all = list()
    
    for(t in 1:T_num){
      X_train_tmp = sub_id_map[[t]][!(sub_id_map[[t]]%in%c(test_index,tuning_index))]
      Y_train_tmp = !(sub_id_map[[t]]%in%c(test_index,tuning_index))
      X_tuning_tmp = sub_id_map[[t]][(sub_id_map[[t]]%in%tuning_index)]
      Y_tuning_tmp = (sub_id_map[[t]]%in%tuning_index)
      X_test_tmp = sub_id_map[[t]][(sub_id_map[[t]]%in%test_index)]
      Y_test_tmp = (sub_id_map[[t]]%in%test_index)
      X_train[[t]] = apply(as.matrix(dose[X_train_tmp,-c(1,2)]),2,as.numeric)
      Y_train[[t]] = Y[[t]][Y_train_tmp, 2]
      X_tune[[t]] = apply(as.matrix(dose[X_tuning_tmp,-c(1,2)]),2,as.numeric)
      Y_tune[[t]] = Y[[t]][Y_tuning_tmp, 2]
      X_test[[t]] = apply(as.matrix(dose[X_test_tmp,-c(1,2)]),2,as.numeric)
      Y_test[[t]] = Y[[t]][Y_test_tmp, 2]
      X_all_tmp = sub_id_map[[t]]
      X_all[[t]] = apply(as.matrix(dose[X_all_tmp,-c(1,2)]),2,as.numeric)
      Y_all[[t]] = Y[[t]][,2]
      
      #rbind training and test set
      X_train[[t]]=rbind(X_train[[t]],X_test[[t]])
      Y_train[[t]]=c(Y_train[[t]],Y_test[[t]])
      
    }
    
    
    ## model training using all of the samples to get lambda range and initial est ##
    if (f==1){
      single_initial_est_all = matrix(0, ncol(X_train[[1]]), T_num) #N of snp by N tissue matrix
      single_summary_all = list()
      for(t in 1:T_num){
        set.seed(05162021)
        #tt = cv.glmnet(X_all[[t]], Y_all[[t]], alpha = 0.5, nfolds = 5,nlambda=50,pmax=200) #cv in 100% of the data, do not find the lambda here #
        tt = cv.glmnet(X_all[[t]], Y_all[[t]], alpha = 0.5, nfolds = 5) # as.numeric(Y_all[[t]]) if esg "Error in y - predmat : non-numeric argument to binary operator"
        single_summary_all[[t]] = tt
        single_initial_est_all[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)] #single tissue EN, beta and lambda with best cv error
      }
    }
    
    ## for each fold
    single_initial_est[[f]] = matrix(0, ncol(X_train[[1]]), T_num) #N of snp by N tissue matrix
    #single_summary = list()
    for(t in 1:T_num){
      set.seed(f)
      #print(length(Y_tune[[t]]))
      #tt = cv.glmnet(X_train[[t]], Y_train[[t]], alpha = 0.5, nfolds = 5,nlambda=50,pmax=200) #cv in 100% of the data, do not find the lambda here #,pmax=200
      tt = cv.glmnet(X_train[[t]], Y_train[[t]], alpha = 0.5, nfolds = 5) # as.numeric(Y_all[[t]]) if esg "Error in y - predmat : non-numeric argument to binary operator"
      #single_summary[[t]] = tt
      single_initial_est[[f]][,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)] #single tissue EN, beta and lambda with best cv error
    }
    
    
    ## use elastic net ests row norm as weights ## 
    lam_range = minmax_lambda(single_summary_all) #lambda range
    sig_norm = apply(single_initial_est[[f]], 1, function(x){sqrt(sum(x^2))}) #norm for each snp
    sig_norm[sig_norm==0] = rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2  #replace =0 with min?
    sig_norm = sig_norm/sum(sig_norm) #??? length=N of snps
    weights2 = 1/sig_norm; weights2 = weights2/sum(weights2);
    
    tis_norm = apply(single_initial_est[[f]], 2, function(x){sum(abs(x))})
    tis_norm[tis_norm==0] = rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
    tis_norm = tis_norm/sum(tis_norm)
    weights1 = 1/tis_norm; weights1 = weights1/sum(weights1);
    lam_V = seq(lam_range[1], lam_range[2], length.out = ntune)
    
    initial_numeric = as.numeric(single_initial_est[[f]])
    
    ## preparation
    XY = grad_prep(X_train, Y_train)
    XX_train = lapply(X_train, function(x){t(x)%*%x/nrow(x)})
    spsz = unlist(lapply(X_train,nrow)) #N of samples for each training set
    res_tune[[f]] = array(NA, dim=c(ntune, ntune, P)) #5*5 matrix *49 tissue, to find the best comb of lam1 and lam2
    
    rec_lamv[f,] = lam_V  #should be the same across 5 folds
    
    for(lam1 in 1:ntune){
      for(lam2 in 1:ntune){
        print(paste0('lam1= ',lam1,' lam2= ',lam2))
        single_est = matrix(initial_numeric, M, P) #snps by tissue matrix. initial iteration with single tissue estimates, save compute time
        ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[lam1]/spsz, lambda2=lam_V[lam2], theta=single_est)
        
        
        #apply to tuning set
        if (lam1==1 & lam2==1){
          multi_res_test[[f]]<-list()
        }
        multi_res_test[[f]][[paste0(lam1,'_',lam2)]] = multi_mse(ans$est, X_tune, Y_tune)
        #print(dim(multi_mse(ans$est, X_tune, Y_tune)))
        
        #ans$est: snps by tissue matrix
        #ans$tune_err: tuning error for each tissue
        #ans$avg_tune_err: average tuning error among 49 tissues
        
        if(sum(ans$est!=0)>0){
          res_tune[[f]][lam1,lam2, ] = ans$tune_err #tuning err for each lam pair
          remove(single_est); remove(ans);
        }else{
          res_tune[[f]][lam1,lam2, ] = ans$tune_err
          remove(single_est); remove(ans);
          #break
        }			
      }
    }
    
    avg_tune_res[[f]] = apply(res_tune[[f]], c(1,2), mean) #average tuning error across the tissues for each lambda pair
    #edt = Sys.time()
    #print(edt-bgt)
  }
  #print(multi_lam)
  print(rec_lamv)  #should be the same across 5 folds
  
  
  
  #-----------------------------
  #find the best comb of lambda1 and lambda2 across 5 folds
  avg_tune_cross_fold<-matrix(rep(1,ntune*ntune),nrow = ntune,ncol=ntune)
  for (i_tune in 1:ntune){
    for (j_tune in 1:ntune){
      avg_tune_cross_fold[i_tune,j_tune]=mean(c(avg_tune_res[[1]][i_tune,j_tune],avg_tune_res[[2]][i_tune,j_tune],avg_tune_res[[3]][i_tune,j_tune],avg_tune_res[[4]][i_tune,j_tune],avg_tune_res[[5]][i_tune,j_tune]),na.rm = T)
      avg_tune_cross_fold[i_tune,j_tune]<-ifelse(is.na(avg_tune_cross_fold[i_tune,j_tune]),1,avg_tune_cross_fold[i_tune,j_tune])
    }
  }
  
  #find the best comb for each fold
  best.lam = which(avg_tune_cross_fold == min(avg_tune_cross_fold), arr.ind = TRUE)[1,] 
  
  #predicted expression in tuning set under best lambda comb
  multi_res<-list()
  for(f in 1:ntune){
    multi_res[[f]]<-multi_res_test[[f]][[paste0(as.numeric(best.lam)[1],'_',as.numeric(best.lam)[2])]]
  }
  
  
  #combine 5 fold results
  cv_df<-list()
  cv_r<-cv_p<-c()
  for (tissue_i in 1:T_num){
    cv_df[[tissue_i]]<-multi_res[[1]][tissue_i][[1]]
    
    for(fold_i in 2:fold){
      cv_df[[tissue_i]]<-rbind(cv_df[[tissue_i]],multi_res[[fold_i]][tissue_i][[1]])
    }
    print(nrow(cv_df[[tissue_i]]))
    fit<-cor.test(cv_df[[tissue_i]][,1],cv_df[[tissue_i]][,2])
    cv_r[tissue_i]<-as.numeric(fit$estimate)
    cv_p[tissue_i]<-as.numeric(fit$p.value)
    
    cv_r[tissue_i]<-ifelse(is.na(cv_r[tissue_i]),0,cv_r[tissue_i])
    cv_p[tissue_i]<-ifelse(is.na(cv_p[tissue_i]),1,cv_p[tissue_i])
  }
  
  
  ## generate an estimate with whole data ##
  XY = grad_prep(X_all, Y_all)
  XX_all = lapply(X_all, function(x){t(x)%*%x/nrow(x)})
  ans = glasso_no_early_stopping(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=lam_V[best.lam[1]]/spsz, lambda2=lam_V[best.lam[2]], theta=single_initial_est_all)

  #-------------------------------------------------------------

  #output 
  info <- gtfile[, c("ALT","COUNTED")]; colnames(info) <- c("ref_allele","counted_allele")
  info$chr_bp <- paste(gtfile$CHR, gtfile$POS, sep="_"); info$rsid <- rownames(info)
  rownames(info) <- NULL; info <- info[, c("rsid","chr_bp","ref_allele","counted_allele")]
  combined_res <- rep(list(NULL), length(Yt)) ## Try combine
  for (tissue_i in 1:length(Yt)){
    tissue_name=Yt[tissue_i]
    weight_df<-info
    weight_df$gene=gene_i
    weight_df$weight=ans$est[,tissue_i]
    #weight_df$r2=round((cv_r[tissue_i])^2,3)
    weight_df$r=cv_r[tissue_i]; weight_df$r2=(cv_r[tissue_i])^2
    weight_df$p=cv_p[tissue_i]
    weight_df$lambda=paste0(round(lam_V[best.lam[1]],2),';',round(lam_V[best.lam[2]],2))
    weight_df$iter=ans$iter
    weight_df<-weight_df[,c('gene','rsid','chr_bp','ref_allele','counted_allele','weight','r','r2','p','lambda','iter')]
    
    #rm 0
    weight_df<-weight_df[weight_df$weight!=0,]
    
    #if(nrow(weight_df)>0 & cv_r[tissue_i]>0.1 & cv_p[tissue_i]<0.05){
    if(nrow(weight_df)>0){
      #saveRDS(weight_df,paste0('/project/yanglab/DNAm/Model/', tissue_name, '/', gene_i, '.rds')) # 1 file per CpG per tissue
      combined_res[[tissue_i]] <- weight_df # Try combine
    }
  }
  names(combined_res) <- Yt
  if (length(unique(combined_res)) > 1) {saveRDS(combined_res, paste0('/scratch/vta8we/scTWAS/UTMOST/Model/100K/Model/', gene_i, '.rds'))} # Try combine
}
}
}

cl <- makeCluster(mycores)
registerDoParallel(cl)
x <- foreach(x=1:mycores,.packages='glmnet') %dopar% func(x)
stopCluster(cl)
