## example script
source("bootstrap_prs.R")
library("dplyr")

phenotype <- commandArgs(trailingOnly=T)[1]
context <- "sex"

pheno_table <- read.table(paste0(phenotype,".pheno"), header=TRUE)

dir.create("valid_r2/", showWarnings = FALSE)
out_file <- paste0("valid_r2/",phenotype,"_",context,"_R2.txt")
if(file.exists(out_file)) quit(save = "no")
  
### Build Phenotype & covariate table
# covar file should contain sex, age, PCs 1-10, and the context
covars <- read.table(paste0("covars/covariates_sa40PC_",context,".pheno"),header=TRUE)
context_col <- 4 # context col index. replace w appropriate col number
covs_table <- covars[,1:14] # FID, IID, age, sex, 10 PCs. replace w appropriate col numbers
rm(covars)
  
# Load in optimized tau threshes & rhos
tuning_vals <- read.table(paste0("rhos/processed_",phenotype,"_whitebrit_sex_PGS.txt"),header=TRUE)
PGS_thresh <- format(tuning_vals[1,1],scientific=F) # should import the pgs threshold
PGS_file <- paste0("valid_pgs/",phenotype,"_pgs.",PGS_thresh,".sscore") # import your local pgs
PGS_IDs <- read.table(PGS_file)[,1:2]

covs <- covs_table[covs_table$FID %in% PGS_IDs$FID, ] # only keep the people in the pgs
pheno <- pheno_table[pheno_table$FID %in% PGS_IDs$FID, ] # only keep the people in the pgs
pheno_cov <- merge(pheno, covs, by=c("FID", "IID"))
colnames(pheno_cov)[3] <- "pheno_code"
rm( PGS_IDs,covs,covs_table,pheno )
  
y <- pheno_cov$pheno_code
X <- pheno_cov[,!colnames(pheno_cov)%in%c("FID","IID","pheno_code")]
Z <- pheno_cov[,context_col]  

### Step 1: Standard PGS
pgs.result <- NULL	
pgs_table <- read.table(PGS_file)[,c(1,6)]
best_pgs <- merge(pheno_cov,pgs_table,by="FID")$SCORE1_SUM
pgs.R2 <- r2_boot(x=best_pgs,y=y,X=X)
pgs_pv <- summary(lm(y ~ 1 + best_pgs + as.matrix(X)))$coef['best_pgs',4]
pgs.result <- rbind(pgs.result, data.frame(Threshold=PGS_thresh,R2=pgs.R2$r2,
                                           w0=pgs.R2$pgs_w0,sd=pgs.R2$r2_sd,
                                           se=pgs.R2$r2_se,pval=pgs_pv,type="pgs"))

### Step 2: ampPGS
PGS_thresh <- format(tuning_vals[3,1],scientific=F) # should import the ampPGS threshold
amp_rho <- tuning_vals[3,6]/tuning_vals[3,5] # w0b/w0a
ampPGS <- best_pgs + amp_rho * (best_pgs * as.numeric(Z))
pgs.R2amp <- r2_diff_boot(x0=best_pgs,x1=ampPGS,y=y,X=X,vali=TRUE)
amppgs.result <- data.frame(Threshold=PGS_thresh,R2=pgs.R2amp$r2,r2_delta=pgs.R2amp$delta,
                            w0=pgs.R2amp$pgsC_w0,w1=pgs.R2amp$pgsC_w1,sd=pgs.R2amp$delta_sd,
                            CI_25=pgs.R2amp$delta_25,CI_75=pgs.R2amp$delta_75,
                            pval=pgs.R2amp$delta_pval,type="ampPGS")

### Step 3: PGSC
pgsc.result <- NULL	
  
PGSC_thresh <- format(tuning_vals[4,1],scientific=F) # should import the PGSC threshold
PGxCS_file <- paste0("valid_pgxcs/",phenotype,"_pgxcs.",PGSC_thresh,".sscore")
PGxCS_table <- read.table(PGxCS_file)[,c(1,6)] 
PGxCS <- merge(pheno_cov,PGxCS_table,by="FID")$SCORE1_SUM
rm(PGxCS_table)

context_max <- PGxCS * as.numeric(Z==max(Z))
context_min <- PGxCS * as.numeric(Z==min(Z))
gen_rho_max <- tuning_vals[4,7]/tuning_vals[4,3] # w1a/w0
gen_rho_min <- tuning_vals[4,8]/tuning_vals[4,3] # w1b/w0
PGSC <- best_pgs + gen_rho_min * context_min + gen_rho_max * context_max

# Compare PGSC to pgs
pgsc.R2 <- r2_diff_boot(x0=best_pgs,x1=PGSC,y=y,X=X,vali=TRUE)
pgsc.result <- rbind(pgsc.result,data.frame(Threshold=PGSC_thresh,R2=pgsc.R2$r2,
                                            r2_delta=pgsc.R2$delta,w0=pgsc.R2$pgsC_w0,
                                            w1=pgsc.R2$pgsC_w1,sd=pgsc.R2$delta_sd,
                                            CI_25=pgsc.R2$delta_25,CI_75=pgsc.R2$delta_75,
                                            pval=pgsc.R2$delta_pval,type="PGSC_v_pgs"))

# Compare PGSC to ampPGS 
pgsc.R2_v_amp <- r2_diff_boot(x0=ampPGS,x1=PGSC,y=y,X=X,vali=TRUE)
pgsc.result <- rbind(pgsc.result,data.frame(Threshold=PGSC_thresh,R2=pgsc.R2_v_amp$r2,
                                            r2_delta=pgsc.R2_v_amp$delta,w0=pgsc.R2_v_amp$pgsC_w0,
                                            w1=pgsc.R2_v_amp$pgsC_w1,sd=pgsc.R2_v_amp$delta_sd,
                                            CI_25=pgsc.R2_v_amp$delta_25,
                                            CI_75=pgsc.R2_v_amp$delta_75,
                                            pval=pgsc.R2_v_amp$delta_pval,type="PGSC_v_amp"))
rm(PGxCS,context_max,context_min)

### Output
pgs.update <- dplyr::bind_rows(pgs.result, amppgs.result,pgsc.result)
write.table(pgs.update, file = out_file)



