# calculate R2
r2fxn <- function(x,y,X){
  X <- as.matrix(X)
  y <- resid( lm( y ~ X ,na.action = na.exclude) )
  x <- resid( lm( x ~ X ,na.action = na.exclude) )
  summary( lm( y ~ x )  )
} 

# bootstrap R2
r2_boot <- function(x,y,B=1e4,X){
  ## x=pgs, y=pheno, X=covars
  lm_out <- r2fxn(x,y,X)
  coef <- ifelse("x" %in% rownames(lm_out$coefficients), lm_out$coefficients["x", "Estimate"], NA)
  if (is.na(coef)) {cat("Error: Coefficient 'x' not found in lm_out, pgs might be 0s\n")}
  r2 <- lm_out$r.squared
  r2_boots <- NULL
  for( b in 1:B ){
    is <- sample( 1:length(x), replace=T )
    x1 <- x[is]
    y1 <- y[is]
    X1 <- X[is,]
    r2_boots[b] <- r2fxn( x1, y1, X1 )$r.squared
  }
  list( r2=r2, pgs_w0=coef, r2_boots=r2_boots, r2_sd=sd(r2_boots), r2_se=lm_out$coefficients["x",2], r2_pval=lm_out$coefficients["x",4])
}

# bootrap R2 diff
r2_diff_boot <- function(x0,x1,y,B=1e4,X,vali=FALSE){
  if(vali==TRUE){
    lm_diff_out <- r2fxn(x1,y,X)
    w0_coef <- NA
    w1_coef <- ifelse("x" %in% rownames(lm_diff_out$coefficients), lm_diff_out$coefficients["x", "Estimate"], NA)
    if (is.na(w1_coef)) {cat("Error: Coefficient 'x' not found in lm_diff_out, GxC pgs might be 0s\n")}
    coef_update <- NA
  }else{
    lm_diff_out <- r2fxn(cbind(x0,x1),y,X)
    w0_coef <-lm_diff_out$coefficients["xx0","Estimate"]
    w1_coef <- ifelse("xx1" %in% rownames(lm_diff_out$coefficients), lm_diff_out$coefficients["xx1", "Estimate"], NA)
    if (is.na(w1_coef)) {cat("Error: Coefficient 'xx1' not found in lm_diff_out, GxC pgs might be 0s\n")}
    coef_update <- w1_coef/w0_coef
  }
  r2 <- lm_diff_out$r.squared
  delta <- r2 - r2fxn(x0,y,X)$r.squared
  delta_boots <- NULL
  for( b in 1:B ){
    is <- sample( 1:length(x0), replace=T )
    x0b <- x0[is]
    x1b <- x1[is]
    yb  <- y[is]
    Xb  <- X[is,]
    if(vali==TRUE){
	    delta_boots[b] <- r2fxn(x1b,yb,Xb)$r.squared - r2fxn(x0b,yb,Xb)$r.squared
    }else{
	    delta_boots[b] <- r2fxn(cbind(x0b,x1b),yb,Xb)$r.squared - r2fxn(x0b,yb,Xb)$r.squared
    }
  }
  delta_lower_pval <- (1 + sum(delta_boots <= 0))/ (1+B)
  delta_2pval <- 2*min(delta_lower_pval,1 - delta_lower_pval)
  list( r2=r2, pgsC_w0=w0_coef, pgsC_w1=w1_coef,coef_diff=coef_update, delta=delta, delta_boots=delta_boots, delta_sd=sd(delta_boots), delta_25=quantile(delta_boots,0.025)[[1]], delta_75=quantile(delta_boots,0.975)[[1]], delta_pval=delta_2pval)
}




