source(file.path(Sys.getenv("PGSC_HOME"), "source_file.R"))
source("bootstrap_prs.R")
library("dplyr")
pheno_table <- read.table(PHENO_FILE, header=TRUE)

cat("pheno:", PHENO_LOWER, "\n")
for(valid_pop in VALID_POPS){
  cat("pop:",valid_pop, "\n")
  for(context in CONTEXT_LIST){
    cat("context:",context, "\n")
    dir.create(paste0(R2,"valid/"), showWarnings = FALSE)
    out_file <- paste0(R2,"valid/processed_",PHENO_LOWER,"_",valid_pop,"_",context,"_PGS.txt")
    misc_out <- paste0(R2,"valid/misc_",PHENO_LOWER,"_",valid_pop,"_",context,".txt")
    boots_out <- paste0(R2,"valid/r2boots_",PHENO_LOWER,"_",valid_pop,"_",context,".txt")
    if(file.exists(out_file)) next 

    ### build phenos+covariates for PGS individuals
    covars <- read.table(paste0(COVARS,"covariates_sa40PC_",context,".pheno"),header=TRUE)
    context_col <- CONTEXT_MAP[[context]]$col
    covs_table  <- covars[, CONTEXT_MAP[[context]]$keep]
    rm(covars)

    best_threshes <- read.table(paste0(R2,"processed_",PHENO_LOWER,"_",TEST_POP,"_",context,"_PGS.txt"),header=TRUE) 
    thresh_list <- format(as.numeric(best_threshes$Threshold),scientific=F)
    PGS_fileA <- paste0(PGS,"valid/",PHENO_LOWER,"_",context,"_",valid_pop,"_pgs_valid.",format(as.numeric(thresh_list[1]), scientific=F),".sscore")
    PGS_IDs <- read.table(PGS_fileA,header=FALSE, col.names=PGS_COLS)[,1:2]
    pheno_cov <- merge(merge(PGS_IDs, covs_table, by=c("FID","IID")),pheno_table,by=c("FID", "IID"))
    if(context %in% c("sex","age")){
      colnames(pheno_cov)[15] <- "pheno_code"
    }else{
      colnames(pheno_cov)[16] <- "pheno_code"
    } 
    rm( PGS_IDs,covs_table)
    
    y <- pheno_cov$pheno_code
    X <- pheno_cov[,!colnames(pheno_cov)%in%c("FID","IID","pheno_code")]
    Z <- pheno_cov[,colnames(pheno_cov)%in%c(context_col)]
    
    ### Step 1: Score all of the tau0
    #std pgs
    prs_tableA <- read.table(PGS_fileA, header=FALSE,col.names=PGS_COLS)[,c(1,6)]
    best_prsA <- merge(pheno_cov,prs_tableA,by="FID")$SCORE1_SUM
    prs.R2A <- r2_boot(x=best_prsA,y=y,X=X)#,B=2)
    prs.result <- data.frame(Threshold=format(as.numeric(thresh_list[1]), scientific=F),R2=prs.R2A$r2,w0=prs.R2A$pgs_w0,sd=prs.R2A$r2_sd,se=prs.R2A$r2_se,pval=prs.R2A$r2_pval,type="pgs")
    pgs_bootsA <- prs.R2A$r2_boots
    pgs_pv <- summary( lm( y ~ 1 + best_prsA + as.matrix(X) )  )$coef['best_prsA',4]

    ### Step 2: ampPGS
    ampprs.result <- NULL
    amp_rho <- best_threshes[3,6]/best_threshes[3,5] # w0b/w0a
    ampPGS <- best_prsA + amp_rho * (best_prsA * as.numeric(Z))
    prs.R2amp <- r2_diff_boot(x0=best_prsA,x1=ampPGS,y=y,X=X,vali=TRUE)#,B=2)
    ampprs.result <- data.frame(Threshold=format(as.numeric(thresh_list[3]), scientific=F),R2=prs.R2amp$r2,r2_delta=prs.R2amp$delta,w0=prs.R2amp$pgsC_w0,w1=prs.R2amp$pgsC_w1,sd=prs.R2amp$delta_sd,CI_25=prs.R2amp$delta_25,CI_75=prs.R2amp$delta_75,pval=prs.R2amp$delta_pval,type="ampPGS")
    amp_boots <- prs.R2amp$delta_boots
    
    ### Step 3: PGSC
    PGxCS_file <- paste0(PGS_GXE,"valid/",PHENO_LOWER,"_",context,"_",valid_pop,"_pgxcs_valid.",format(as.numeric(thresh_list[4]), scientific=F),".sscore")
    PGxCS_table <- read.table(PGxCS_file,header=FALSE,col.names=PGS_COLS)[,c(1,6)]

    # generalized PGSC: PGS + w1a PGxCS * Ca + w1b PGxCS * Cb
    PGxCS <- merge(pheno_cov,PGxCS_table,by="FID")$SCORE1_SUM
    rm(PGxCS_table)
    pgs_gxc_cor <- cor(best_prsA, PGxCS) # misc summ stat
    context_max <- PGxCS * as.numeric(Z==max(Z))
    context_min <- PGxCS * as.numeric(Z==min(Z))
    gen_rho_maxA <- best_threshes[4,7]/best_threshes[4,3] # w1a/w0
    gen_rho_minA <- best_threshes[4,8]/best_threshes[4,3] # w1b/w0
    PGSC <- best_prsA + gen_rho_minA * context_min + gen_rho_maxA * context_max
  
    ### Step 4: Test for significance
    prs.update <- dplyr::bind_rows(prs.result,ampprs.result)
    pgsc_list <- c("ampPGS","PGSC")
    threshold_map <- c(PGSC = thresh_list[4])
    prsc.result <- NULL
    # PGSCs vs pgs 
    for(pgsc_name in pgsc_list[-1]){
      thresh <- threshold_map[pgsc_name]
      pgsc.R2A <- r2_diff_boot(x0=best_prsA,x1=get(pgsc_name),y=y,X=X,vali=TRUE)#, B=2)
      prsc.result <- rbind(prsc.result,data.frame(Threshold=format(as.numeric(thresh), scientific=F),R2=pgsc.R2A$r2,r2_delta=pgsc.R2A$delta,w0=pgsc.R2A$pgsC_w0,w1=pgsc.R2A$pgsC_w1,sd=pgsc.R2A$delta_sd,CI_25=pgsc.R2A$delta_25,CI_75=pgsc.R2A$delta_75,pval=pgsc.R2A$delta_pval,type=paste0(pgsc_name,"_v_pgs")))
    }

    # gen PGSC vs. everything but itself
    for(pgsc_name in pgsc_list[-2]){
      pgsc.R2C <- r2_diff_boot(x0=get(pgsc_name),x1=PGSC,y=y,X=X,vali=TRUE)#,B=2)
      prsc.result <- rbind(prsc.result,data.frame(Threshold=format(as.numeric(thresh_list[4]), scientific=F),R2=pgsc.R2C$r2,r2_delta=pgsc.R2C$delta,w0=pgsc.R2C$pgsC_w0,w1=pgsc.R2C$pgsC_w1,sd=pgsc.R2C$delta_sd,CI_25=pgsc.R2C$delta_25,CI_75=pgsc.R2C$delta_75,pval=pgsc.R2C$delta_pval,type=paste0("PGSC_v_",pgsc_name)))
    }
    
    prs.update <- dplyr::bind_rows(prs.update, prsc.result)
    all_boots <- cbind(pgs_bootsA, amp_boots)
    write.table(prs.update, file = out_file)
    write.table(all_boots, file = boots_out)
    write.table(c(pgs_pv,pgs_gxc_cor), file = misc_out,row.names=c("pgs_pv","pgs_gxc_cor"))
    rm(PGxCS,context_max,context_min, prs.update, all_boots)
  }
}



