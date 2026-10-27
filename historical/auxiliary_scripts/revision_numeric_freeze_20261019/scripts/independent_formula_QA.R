# Public-package copy: see metadata/code_changes.csv for provenance.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else getwd(), winslash='/', mustWork=TRUE)
out <- file.path(root,'revision_followup_20261019','bridge43')
resdir <- file.path(out,'results'); qadir <- file.path(out,'qa')
dir.create(qadir,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(limma);library(jsonlite)})
E <- readRDS(file.path(root,'revision_scope_lock_20261018/bridge43/inputs_local/expression43.rds'))
Xs <- readRDS(file.path(root,'revision_scope_lock_20261018/bridge43/designs_local/design_matrices.rds'))
sets <- readRDS(file.path(root,'revision_scope_lock_20261018/bridge43/inputs_local/GO_indices.rds'))
all_checks <- list(); add <- function(name,value,detail='') {all_checks[[length(all_checks)+1L]] <<- list(check=name,passed=isTRUE(value),detail=detail)}
for (mid in names(Xs)) {
  D <- Xs[[mid]]; p <- ncol(D); n <- nrow(D); G <- nrow(E); gcol <- which(colnames(D)=='groupSCI')
  Xr <- D[,c(setdiff(seq_len(p),gcol),gcol),drop=FALSE]
  qr0 <- qr(Xr); eff <- qr.qty(qr0,t(E)); ut <- eff[p,]
  if(qr0$qr[p,p]<0) ut <- -ut
  U <- eff[-seq_len(p),,drop=FALSE]; sig2 <- colMeans(U^2)
  sv <- squeezeVar(sig2,df=n-p,covariate=NULL)
  stdev <- sqrt(diag(solve(crossprod(Xr))))[p]
  b <- as.numeric(solve(crossprod(D),crossprod(D,t(E)))[which(colnames(D)=="groupSCI"),])
  modt <- ut/sqrt(sv$var.post)
  tab <- fread(file.path(resdir,paste0('gene_models_',mid,'.csv')))
  add(paste0(mid,'_base_group_effect'),max(abs(b-tab$logFC))<=1e-10,paste0('max_abs=',format(max(abs(b-tab$logFC)),digits=16)))
  add(paste0(mid,'_moderated_t_reconstruction'),max(abs(modt-tab$t))<=1e-10,paste0('max_abs=',format(max(abs(modt-tab$t)),digits=16)))
  dft <- pmin((n-p)+sv$df.prior,G*(n-p))
  ses <- sqrt(sv$var.post)*sqrt(diag(solve(crossprod(D))))[gcol]
  ci_low <- b-qt(.975,df=dft)*ses
  ci_high <- b+qt(.975,df=dft)*ses
  pp <- 2*pt(-abs(b/ses),df=dft)
  add(paste0(mid,'_moderated_CI_reconstruction'),max(abs(c(ci_low-tab$CI.L,ci_high-tab$CI.R)))<=1e-10,paste0('max_abs=',format(max(abs(c(ci_low-tab$CI.L,ci_high-tab$CI.R))),digits=16)))
  add(paste0(mid,'_moderated_P_reconstruction'),max(abs(pp-tab$P.Value))<=1e-10,paste0('max_abs=',format(max(abs(pp-tab$P.Value)),digits=16)))
  add(paste0(mid,'_BH_reconstruction'),max(abs(p.adjust(tab$P.Value,'BH')-tab$adj.P.Val))<=5e-14,paste0('max_abs=',format(max(abs(p.adjust(tab$P.Value,'BH')-tab$adj.P.Val)),digits=16)))
  ## CAMERA aggregation reconstruction (fixed rho=.01)
  df.total <- min((n-p)+sv$df.prior,G*(n-p))
  Stat <- limma:::zscoreT(modt,df=df.total,approx=TRUE,method='hill')
  meanStat <- mean(Stat); varStat <- var(Stat); G2 <- G
  pvals <- numeric(length(sets)); dirs <- character(length(sets)); nset <- integer(length(sets))
  for(i in seq_along(sets)) {
    ii <- sets[[i]]; m <- length(ii); m2 <- G2-m; vif <- max(1,1+(m-1)*.01)
    delta <- G2/m2*(mean(Stat[ii])-meanStat)
    varP <- ((G2-1)*varStat-delta^2*m*m2/G2)/(G2-2)
    st <- delta/sqrt(varP*(vif/m+1/m2))
    dn <- pt(st,df=G2-2); up <- pt(st,df=G2-2,lower.tail=FALSE)
    pvals[i] <- 2*min(dn,up); dirs[i] <- ifelse(dn<up,'Down','Up'); nset[i] <- m
  }
  cam <- fread(file.path(resdir,paste0('camera_models_',mid,'.csv')))
  add(paste0(mid,'_CAMERA_P_reconstruction'),max(abs(pvals-cam$PValue))<=1e-10,paste0('max_abs=',format(max(abs(pvals-cam$PValue)),digits=16)))
  add(paste0(mid,'_CAMERA_direction_reconstruction'),all(dirs==cam$Direction),paste0('mismatches=',sum(dirs!=cam$Direction)))
  add(paste0(mid,'_CAMERA_BH_reconstruction'),max(abs(p.adjust(pvals,'BH')-cam$FDR))<=5e-14,paste0('max_abs=',format(max(abs(p.adjust(pvals,'BH')-cam$FDR)),digits=16)))
  fwrite(data.table(model=mid,gene_effect_max_abs=max(abs(b-tab$logFC)),moderated_t_max_abs=max(abs(modt-tab$t)),gene_q_max_abs=max(abs(p.adjust(tab$P.Value,'BH')-tab$adj.P.Val)),camera_p_max_abs=max(abs(pvals-cam$PValue)),camera_q_max_abs=max(abs(p.adjust(pvals,'BH')-cam$FDR)),camera_direction_mismatches=sum(dirs!=cam$Direction)),file.path(qadir,paste0('independent_QA_',mid,'.csv')))
}
qa <- list(passed=all(vapply(all_checks,`[[`,logical(1),'passed')),checks=all_checks,scope='Formula reconstruction with base QR; limma squeezeVar/zscoreT retained for the locked moderated-variance and CAMERA transform; not a second biological dataset or independent calibration.',software=list(R=R.version.string,limma=as.character(packageVersion('limma'))))
write_json(qa,file.path(qadir,'independent_formula_QA.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat('independent QA passed:',qa$passed,'checks:',length(all_checks),'\n')

