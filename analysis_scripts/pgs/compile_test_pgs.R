source(file.path(Sys.getenv("PGSC_HOME"), "source_file.R"))
library("dplyr")

pheno_table <- read.table(PHENO_FILE, header=TRUE)
cat("pheno:", PHENO_LOWER, "\n")

for(context in CONTEXT_LIST){
  cat("context:", context, "\n")
  dir.create(R2, showWarnings = FALSE)
  dir.create(paste0(R2,"all_thresh/"), showWarnings = FALSE)
  out_file <- paste0(R2,"processed_",PHENO_LOWER,"_",TEST_POP,"_",context,"_PGS.txt")
  all_thresh <- paste0(R2,"all_thresh/all_thresh_",PHENO_LOWER,"_",TEST_POP,"_",context,"_PGS.txt")
  if(file.exists(out_file)) next

  ### build phenos+covariates for PGS individuals
  covars <- read.table(paste0(COVARS,"covariates_sa40PC_",context,".pheno"),header=TRUE)
  context_col <- CONTEXT_MAP[[context]]$col
  covs_table  <- covars[, CONTEXT_MAP[[context]]$keep]
  rm(covars)

  PGS_file <- paste0(PGS, PHENO_LOWER,"_",context,"_pgs.",THRESH[1],".sscore")  
  PGS_IDs <- read.table(PGS_file,header=FALSE, col.names=PGS_COLS)[,1:2]
  pheno_cov <- merge(merge(PGS_IDs, covs_table, by=c("FID","IID")),pheno_table,by=c("FID", "IID"))
  if(context %in% c("sex","age")){
    colnames(pheno_cov)[15] <- "pheno_code"
  }else{
    colnames(pheno_cov)[16] <- "pheno_code"
  }
  y <- pheno_cov$pheno_code
  X <- pheno_cov[,!names(pheno_cov)%in%c("FID","IID","pheno_code")]
  Z <- pheno_cov[[context_col]]
  rm(PGS_IDs,PGS_file,covs_table)
  
  ### Step 1: Score all of the tau0
  # test pgs
  prs.resultA <- NULL
  for(i in THRESH)try({ 
    prs_tableA <- read.table(paste0(PGS, PHENO_LOWER,"_",context,"_pgs.",i,".sscore"), header=FALSE,col.names=PGS_COLS)[,c(1,6)]
    prsA <- merge(pheno_cov,prs_tableA,by="FID")$SCORE1_SUM
    prs.R2A <- R2fxn(x=prsA,y=y,X=X)
    prs.resultA <- rbind(prs.resultA, data.frame(Threshold=i,R2=prs.R2A$r.squared,w0=prs.R2A$coef["x","Estimate"],type="pgs"))
  })
  best_std_threshA <- prs.resultA[which.max(prs.resultA$R2),]
  best_prs_tableA <- read.table(paste0(PGS, PHENO_LOWER,"_",context,"_pgs.",best_std_threshA[[1]],".sscore"), header=FALSE, col.names=PGS_COLS)[,c(1,6)]
  best_prsA <- merge(pheno_cov,best_prs_tableA,by="FID")$SCORE1_SUM
  rm(prs_tableA,best_prs_tableA)

  ### Step 2: ampPGS
  prs.R2amp <- R2fxn(x=cbind( best_prsA, best_prsA * as.numeric(Z) ),y=y,X=X)
  ampprs.result <- data.frame(Threshold=best_std_threshA[[1]],R2=prs.R2amp$r.squared,w0a=prs.R2amp$coef["xbest_prsA","Estimate"],w0b=prs.R2amp$coef["x","Estimate"], type="ampPGS")
  
  ### Step 3a: PGSC model variations vs PGS
  prsc.resultA <- NULL

  for(i in THRESH)try({ 
    PGxCS_table <- read.table(paste0(PGS_GXE, PHENO_LOWER,"_",context,"_pgxcs.",i,".sscore"),header=FALSE,col.names=PGS_COLS)[,c(1,6)] 
    PGxCS <- merge(pheno_cov,PGxCS_table,by="FID")$SCORE1_SUM
    context_max <- PGxCS * as.numeric(Z==max(Z))
    context_min <- PGxCS * as.numeric(Z==min(Z))
    rm(PGxCS_table)

    # generalized PGSC: w0 PGS + w1a PGxCS * max(Ca) + w1b PGxCS * min(Cb)
    prs.R2A <- R2fxn(x=cbind( best_prsA, context_max , context_min ),y=y,X=X)
    prsc.resultA <- rbind(prsc.resultA,data.frame(Threshold=i,R2=prs.R2A$r.squared,w0=prs.R2A$coef["xbest_prsA","Estimate"],w1a=prs.R2A$coef["xcontext_max","Estimate"],w1b=prs.R2A$coef["xcontext_min","Estimate"],type="PGSC"))
  })
  
  best_prscA <- prsc.resultA[which.max(prsc.resultA$R2),]

  prs_full <- dplyr::bind_rows(prs.resultA, ampprs.result, prsc.resultA)
  prs.update <- dplyr::bind_rows(best_std_threshA, ampprs.result, best_prscA)
  
  write.table(prs_full, file = all_thresh)
  write.table(prs.update, file = out_file)
} 
