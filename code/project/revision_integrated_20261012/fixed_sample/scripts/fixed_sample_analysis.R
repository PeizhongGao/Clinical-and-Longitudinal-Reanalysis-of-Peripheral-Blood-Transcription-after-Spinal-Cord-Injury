# Public-package copy: see metadata/code_changes.csv for provenance.
# Fixed-sample revision-stage diagnostics. Scores remain on the frozen TC scale.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) normalizePath(args[1], winslash='/', mustWork=TRUE) else normalizePath(getwd(),winslash='/',mustWork=TRUE)
base <- file.path(root,'revision_R2_20261008')
out <- file.path(base,'fixed_sample_20261012')
for (d in c('inputs','tables','logs','figures','scripts')) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
wcsv <- function(x,n) fwrite(x,file.path(out,'tables',paste0(n,'.csv')),na='NA')
logfile <- file.path(out,'logs','analysis_run.txt')
sink(logfile,split=TRUE)
cat('Started UTC:',format(Sys.time(),tz='UTC'), '\n')
files <- c(clinical_scores=file.path(base,'tables','ODC405_locked_scores_with_clinical_and_cell_estimates.csv'),
           score_reference=file.path(base,'supplement','discovery_locked_scores_descriptive.csv'),
           prior_models=file.path(base,'supplement','ODC405_locked_score_models.csv'),
           frozen_genes=file.path(base,'supplement','locked_signature_reused_S14.csv'),
           standardization=file.path(base,'supplement','standardization_parameters.csv'),
           protocol=file.path(out,'analysis-protocol.md'))
stopifnot(all(file.exists(files)))
input_hashes <- data.table(role=names(files),source=unname(files),md5=unname(tools::md5sum(files)))
wcsv(input_hashes,'input_hashes_before')
for (i in seq_along(files)) if (names(files)[i]!='protocol') stopifnot(file.copy(files[i],file.path(out,'inputs',basename(files[i])),overwrite=TRUE))
d <- fread(files['clinical_scores'],na.strings=c('NA',''))
ref <- fread(files['score_reference'],na.strings=c('NA',''))
prior <- fread(files['prior_models'])
outcomes <- c('innate_score','adaptive_score','balance_score')
stopifnot(nrow(d)==48L,uniqueN(d$sample)==48L,uniqueN(d$subject)==48L,uniqueN(d$gsm)==48L,
          sum(d$group=='SCI')==38L,sum(d$group=='TC')==10L)
reference_join <- merge(d[,c('sample',outcomes),with=FALSE],ref[,c('sample',outcomes),with=FALSE],by='sample',suffixes=c('_input','_reference'))
stopifnot(nrow(reference_join)==48L)
audit <- rbindlist(lapply(outcomes,function(y) data.table(check=paste0('frozen_',y),max_abs_difference=max(abs(reference_join[[paste0(y,'_input')]]-reference_join[[paste0(y,'_reference')]])))))
audit <- rbind(audit,data.table(check='balance_identity',max_abs_difference=max(abs(d$balance_score-d$innate_score+d$adaptive_score))))
stopifnot(all(audit$max_abs_difference<1e-11), all(is.finite(as.matrix(d[,..outcomes]))))
wcsv(audit,'frozen_score_audit')
d[,group:=factor(group,levels=c('TC','SCI'))]
d[,sex:=factor(sex,levels=c('F','M'))]
for(v in c('tbi','prior_cns')) set(d,j=v,value=factor(d[[v]],levels=c('No','Yes','Unknown')))
centers <- data.table(variable=c('age','iss','draw_hours'),center=vapply(c('age','iss','draw_hours'),function(v)mean(d[[v]],na.rm=TRUE),numeric(1)))
for(i in seq_len(nrow(centers))) set(d,j=c('age_c','iss_c','draw_hours_c')[i],value=d[[centers$variable[i]]]-centers$center[i])
wcsv(centers,'covariate_center_constants')
clinical_vars <- c('age','sex','iss','tbi','prior_cns')
d[,complete_clinical:=complete.cases(.SD),.SDcols=clinical_vars]
d[,complete_time:=complete_clinical & is.finite(draw_hours)]
d[,missing_reason_clinical:=vapply(seq_len(.N),function(i)paste(clinical_vars[vapply(clinical_vars,function(v)is.na(d[[v]][i]),logical(1))],collapse=';'),character(1))]
d[,missing_reason_time:=ifelse(complete_time,'',ifelse(!complete_clinical,paste0(missing_reason_clinical,ifelse(is.na(draw_hours),';draw_hours','')),'draw_hours'))]
wcsv(d[,.(sample,gsm,subject,group,complete_clinical,complete_time,missing_reason_clinical,missing_reason_time)],'sample_membership')
sets <- list(full48=copy(d),clinical43=copy(d[complete_clinical==TRUE]),time41=copy(d[complete_time==TRUE]))
stopifnot(nrow(sets$clinical43)==43L,nrow(sets$time41)==41L,
          sum(sets$clinical43$group=='TC')==7L,sum(sets$time41$group=='TC')==5L,
          identical(sets$clinical43[group=='SCI',sample],sets$time41[group=='SCI',sample]))
for(s in names(sets)) wcsv(sets[[s]],paste0('fixed_input_',s))
rhs <- c(unadjusted='group',age_sex='group + age_c + sex',clinical='group + age_c + sex + iss_c + tbi + prior_cns',clinical_time='group + age_c + sex + iss_c + tbi + prior_cns + draw_hours_c')
specs <- list(full48=c('unadjusted','age_sex'),clinical43=c('unadjusted','age_sex','clinical'),time41=c('unadjusted','age_sex','clinical','clinical_time'))

hc3_fit <- function(dat,y,spec){
  fit <- lm(as.formula(paste(y,'~',rhs[spec])),data=dat,na.action=na.fail,x=TRUE,y=TRUE)
  X <- model.matrix(fit); j <- match('groupSCI',colnames(X)); p <- ncol(X)
  if(fit$rank!=p || is.na(j) || nobs(fit)!=nrow(dat)) stop('Invalid model design.')
  B <- solve(crossprod(X)); h <- rowSums((X%*%B)*X); e <- residuals(fit)
  stopifnot(max(h)<1-1e-8)
  weights <- e^2/(1-h)^2
  V <- B %*% crossprod(X,X*weights) %*% B
  se <- sqrt(V[j,j]); est <- unname(coef(fit)[j]); df <- df.residual(fit); crit <- qt(.975,df)
  group_resid <- if(p>2) residuals(lm.fit(X[,-j,drop=FALSE],X[,j])) else X[,j]-mean(X[,j])
  total_group_ss <- sum((X[,j]-mean(X[,j]))^2)
  group_resid_ss <- sum(group_resid^2)
  Xscaled <- cbind(1,scale(X[,-1,drop=FALSE]))
  group_vif <- total_group_ss/group_resid_ss
  result <- data.table(outcome=y,spec=spec,estimate=est,robust_SE_HC3=se,
    lower_CL=est-crit*se,upper_CL=est+crit*se,CI_width=2*crit*se,statistic=est/se,
    df=df,p_value=2*pt(-abs(est/se),df),OLS_SE=sqrt(vcov(fit)[j,j]),n_samples=nrow(dat),
    n_participants=uniqueN(dat$subject),SCI_n=sum(dat$group=='SCI'),TC_n=sum(dat$group=='TC'),
    hat_max=max(h),design_rank=fit$rank,design_columns=p,condition_number_scaled=kappa(Xscaled,exact=TRUE),
    group_VIF=group_vif,group_R2_from_covariates=1-group_resid_ss/total_group_ss,
    group_residual_SS=group_resid_ss,model=paste(y,'~',rhs[spec]),postselection=TRUE)
  inf <- data.table(sample=dat$sample,gsm=dat$gsm,group=as.character(dat$group),hat=h,
    residual=e,fitted=as.numeric(fitted(fit)),studentized_residual=as.numeric(rstudent(fit)),
    Cook_distance=as.numeric(cooks.distance(fit)),DFBETAS_group=as.numeric(dfbetas(fit)[,j]),
    group_beta_change_if_deleted=-as.numeric(X%*%B[,j])*e/(1-h),
    leverage_screen_2p_over_n=2*p/nrow(dat),Cook_screen_4_over_n=4/nrow(dat),DFBETAS_screen_2_over_sqrt_n=2/sqrt(nrow(dat)))
  list(fit=fit,X=X,result=result,influence=inf,V=V)
}

results <- list(); fits <- list(); infl <- list(); assumptions <- list(); design_columns <- list()
for(s in names(sets)) for(spec in specs[[s]]) for(y in outcomes){
  k <- paste(s,spec,y,sep='__'); z <- hc3_fit(sets[[s]],y,spec)
  z$result[,subset:=s]; results[[k]] <- z$result; fits[[k]] <- z
  z$influence[,`:=`(subset=s,spec=spec,outcome=y)]
  infl[[k]] <- z$influence
  sh <- shapiro.test(residuals(z$fit)); bp <- lmtest::bptest(z$fit,studentize=TRUE)
  assumptions[[k]] <- data.table(subset=s,spec=spec,outcome=y,n=nrow(sets[[s]]),
    Shapiro_W=unname(sh$statistic),Shapiro_P=sh$p.value,BP_studentized=unname(bp$statistic),BP_df=unname(bp$parameter),BP_P=bp$p.value,
    residual_RMSE=sqrt(mean(residuals(z$fit)^2)),max_abs_studentized_residual=max(abs(z$influence$studentized_residual)),
    n_leverage_screen=sum(z$influence$hat>z$influence$leverage_screen_2p_over_n),
    n_Cook_screen=sum(z$influence$Cook_distance>z$influence$Cook_screen_4_over_n),
    n_DFBETAS_screen=sum(abs(z$influence$DFBETAS_group)>z$influence$DFBETAS_screen_2_over_sqrt_n))
  if(y==outcomes[1]){
    X <- z$X
    vifs <- diag(solve(cor(X[,-1,drop=FALSE])))
    design_columns[[paste(s,spec)]] <- data.table(subset=s,spec=spec,column=colnames(X)[-1],VIF=as.numeric(vifs))
    wcsv(data.table(sample=sets[[s]]$sample,X),paste0('design_',s,'_',spec))
  }
}
res <- rbindlist(results)
res[,q_BH_3_within_model:=p.adjust(p_value,'BH'),by=.(subset,spec)]
res[,q_BH_18_primary_models:=NA_real_]
res[outcome!='balance_score',q_BH_18_primary_models:=p.adjust(p_value,'BH')]
wcsv(res,'fixed_sample_models'); wcsv(rbindlist(infl),'all_participant_influence')
wcsv(rbindlist(assumptions),'model_assumption_diagnostics'); wcsv(rbindlist(design_columns),'design_column_VIF')
saveRDS(fits,file.path(out,'logs','fitted_models.rds'))

# Exact source reproduction, independent of re-centering continuous covariates.
mapping <- data.table(subset=c('full48','full48','clinical43','time41'),spec=c('unadjusted','age_sex','clinical','clinical_time'),
                      old_spec=c('unadjusted','age_sex','clinical_complete_case','clinical_plus_draw_time'))
checks <- list()
for(i in seq_len(nrow(mapping))) for(y in outcomes){
  a <- res[subset==mapping$subset[i] & spec==mapping$spec[i] & outcome==y]
  b <- prior[spec==mapping$old_spec[i] & outcome==y]
  stopifnot(nrow(a)==1L,nrow(b)==1L,a$n_samples==b$n_samples)
  fields <- c('estimate','robust_SE_HC3','lower_CL','upper_CL','p_value','hat_max')
  delta <- max(abs(as.numeric(a[,..fields])-as.numeric(b[,..fields])))
  q_delta <- abs(a$q_BH_3_within_model-b$q_BH)
  checks[[paste(i,y)]] <- data.table(subset=a$subset,spec=a$spec,outcome=y,max_abs_numeric_difference=delta,q_difference=q_delta,passed=delta<1e-9 & q_delta<1e-9)
}
repro <- rbindlist(checks); stopifnot(all(repro$passed)); wcsv(repro,'prior_model_reproduction')

transitions <- list(
 list(label='48_to_43_selection_unadjusted',from=c('full48','unadjusted'),to=c('clinical43','unadjusted'),type='sample_selection'),
 list(label='43_age_sex_added',from=c('clinical43','unadjusted'),to=c('clinical43','age_sex'),type='fixed_sample_adjustment'),
 list(label='43_full_clinical_added',from=c('clinical43','age_sex'),to=c('clinical43','clinical'),type='fixed_sample_adjustment'),
 list(label='43_to_41_selection_clinical',from=c('clinical43','clinical'),to=c('time41','clinical'),type='sample_selection'),
 list(label='41_draw_time_added',from=c('time41','clinical'),to=c('time41','clinical_time'),type='fixed_sample_adjustment'),
 list(label='48_to_41_selection_unadjusted',from=c('full48','unadjusted'),to=c('time41','unadjusted'),type='sample_selection'),
 list(label='41_age_sex_added',from=c('time41','unadjusted'),to=c('time41','age_sex'),type='fixed_sample_adjustment'),
 list(label='41_full_clinical_added',from=c('time41','age_sex'),to=c('time41','clinical'),type='fixed_sample_adjustment'))
changes <- rbindlist(lapply(transitions,function(t) rbindlist(lapply(outcomes,function(y){
  a <- res[subset==t$from[1]&spec==t$from[2]&outcome==y]; b <- res[subset==t$to[1]&spec==t$to[2]&outcome==y]
  data.table(transition=t$label,type=t$type,outcome=y,from_subset=a$subset,from_spec=a$spec,to_subset=b$subset,to_spec=b$spec,
    from_n=a$n_samples,to_n=b$n_samples,estimate_before=a$estimate,estimate_after=b$estimate,
    estimate_change=b$estimate-a$estimate,abs_estimate_ratio=abs(b$estimate)/abs(a$estimate),
    HC3_SE_before=a$robust_SE_HC3,HC3_SE_after=b$robust_SE_HC3,HC3_SE_ratio=b$robust_SE_HC3/a$robust_SE_HC3,
    CI_width_before=a$CI_width,CI_width_after=b$CI_width,CI_width_ratio=b$CI_width/a$CI_width,
    interpretation='Descriptive, path-dependent change; not a causal attribution or test between independent estimates')
}))))
wcsv(changes,'model_stage_changes')

# Refit every TC deletion in the fuller models; never choose a favorable deletion.
loo_specs <- list(c('clinical43','clinical'),c('time41','clinical'),c('time41','clinical_time'))
loo <- list(); loo_checks <- list()
for(pair in loo_specs){
  s <- pair[1]; spec <- pair[2]; dat <- sets[[s]]
  for(id in dat[group=='TC',sample]) for(y in outcomes){
    k <- paste(s,spec,y,sep='__'); full <- fits[[k]]
    z <- hc3_fit(dat[sample!=id],y,spec); row <- z$result
    row[,`:=`(subset=s,deleted_TC=id,deleted_gsm=dat[sample==id,gsm],full_estimate=full$result$estimate,
      estimate_change=estimate-full$result$estimate,full_lower_CL=full$result$lower_CL,full_upper_CL=full$result$upper_CL,
      direction_changed=sign(estimate)!=sign(full$result$estimate),diagnostic_only=TRUE)]
    loo[[paste(k,id)]] <- row
    expected <- full$influence[sample==id,group_beta_change_if_deleted]
    loo_checks[[paste(k,id)]] <- data.table(subset=s,spec=spec,outcome=y,deleted_TC=id,
      analytic_delta=expected,refit_delta=row$estimate_change,max_abs_difference=abs(expected-row$estimate_change))
  }
}
loo <- rbindlist(loo)
loo[,q_BH_3_within_deletion:=p.adjust(p_value,'BH'),by=.(subset,spec,deleted_TC)]
loo[,multiplicity_caveat:='Within-deletion BH only; does not adjust across deletions, models or discovery selection']
wcsv(loo,'TC_leave_one_out_models')
loochk <- rbindlist(loo_checks); stopifnot(max(loochk$max_abs_difference)<1e-9); wcsv(loochk,'TC_leave_one_out_identity_check')
loosum <- loo[,.(n_deletions=.N,estimate_min=min(estimate),estimate_max=max(estimate),
  max_abs_change=max(abs(estimate_change)),most_influential_TC=deleted_TC[which.max(abs(estimate_change))],
  any_direction_changed=any(direction_changed),all_HC3_CI_include_zero=all(lower_CL<=0 & upper_CL>=0),
  min_lower_CL=min(lower_CL),max_upper_CL=max(upper_CL),min_HC3_SE=min(robust_SE_HC3),max_HC3_SE=max(robust_SE_HC3)),by=.(subset,spec,outcome)]
wcsv(loosum,'TC_leave_one_out_summary')

# Descriptive support and selection audit; no hypothesis tests on tiny omitted groups.
descriptive <- list(); categorical <- list(); overlap <- list()
for(s in names(sets)){
  dat <- sets[[s]]
  for(g in c('TC','SCI')) for(v in c(outcomes,'age','iss','draw_hours')){
    x <- dat[group==g][[v]]; x <- x[is.finite(x)]
    descriptive[[paste(s,g,v)]] <- data.table(subset=s,group=g,variable=v,n_group=sum(dat$group==g),n_available=length(x),n_missing=sum(dat$group==g)-length(x),
      mean=mean(x),SD=sd(x),median=median(x),Q1=as.numeric(quantile(x,.25)),Q3=as.numeric(quantile(x,.75)),min=min(x),max=max(x))
  }
  for(v in c('sex','tbi','prior_cns')) for(g in c('TC','SCI')){
    categorical[[paste(s,v,g)]] <- data.table(subset=s,group=g,variable=v,level=levels(dat[[v]]),
      n=vapply(levels(dat[[v]]),function(l)sum(dat$group==g & dat[[v]]==l),integer(1)))
  }
  for(v in c('age','iss','draw_hours')){
    tc <- dat[group=='TC'][[v]]; tc <- tc[is.finite(tc)]
    sci <- dat[group=='SCI'][[v]]; sci <- sci[is.finite(sci)]
    low <- max(min(tc),min(sci)); high <- min(max(tc),max(sci))
    overlap[[paste(s,v)]] <- data.table(subset=s,variable=v,TC_n=length(tc),SCI_n=length(sci),
      TC_min=min(tc),TC_max=max(tc),SCI_min=min(sci),SCI_max=max(sci),overlap_low=low,overlap_high=high,
      TC_n_in_range_overlap=sum(tc>=low & tc<=high),SCI_n_in_range_overlap=sum(sci>=low & sci<=high),
      caveat='Univariate observed-range intersection; not sufficient multivariable common support')
  }
}
wcsv(rbindlist(descriptive),'descriptive_statistics'); wcsv(rbindlist(categorical),'categorical_counts'); wcsv(rbindlist(overlap),'covariate_range_overlap')
selection <- list()
for(s in c('clinical43','time41')) for(g in c('TC','SCI')) for(kept in c(TRUE,FALSE)) for(v in c(outcomes,'age','iss','draw_hours')){
  ids <- sets[[s]]$sample; z <- d[group==g & (sample%in%ids)==kept]; x <- z[[v]]; x <- x[is.finite(x)]
  selection[[paste(s,g,kept,v)]] <- data.table(subset=s,group=g,retained=kept,variable=v,n_group=nrow(z),n_available=length(x),
    mean=if(length(x))mean(x) else NA_real_,SD=if(length(x)>1)sd(x) else NA_real_,
    median=if(length(x))median(x) else NA_real_,min=if(length(x))min(x) else NA_real_,max=if(length(x))max(x) else NA_real_)
}
wcsv(rbindlist(selection),'selection_descriptives')
fullinf <- rbindlist(infl)[spec%in%c('clinical','clinical_time')]
infsum <- fullinf[,.(max_hat=max(hat),max_hat_sample=sample[which.max(hat)],max_hat_group=group[which.max(hat)],
  max_abs_group_change_all=max(abs(group_beta_change_if_deleted)),most_influential_all_sample=sample[which.max(abs(group_beta_change_if_deleted))],
  max_abs_group_change_TC=max(abs(group_beta_change_if_deleted[group=='TC'])),max_Cook=max(Cook_distance)),by=.(subset,spec,outcome)]
wcsv(infsum,'fuller_model_influence_summary')

after <- unname(tools::md5sum(files)); stopifnot(identical(after,input_hashes$md5))
wcsv(data.table(role=names(files),source=unname(files),md5_before=input_hashes$md5,md5_after=after,unchanged=after==input_hashes$md5),'input_integrity_after')
qa <- list(status='PASS',analysis_unit='one participant per row',n_models=nrow(res),n_unique_designs=uniqueN(res[,.(subset,spec)]),
  fixed_sets=list(full48=48,clinical43=43,time41=41),no_score_restandardization=TRUE,
  score_max_difference=max(audit$max_abs_difference),prior_max_difference=max(repro$max_abs_numeric_difference),
  all_designs_full_rank=all(res$design_rank==res$design_columns),n_TC_deletion_fits=nrow(loo),
  deletion_identity_max_difference=max(loochk$max_abs_difference),inputs_unchanged=TRUE,
  exploratory_postselection=TRUE,independent_validation=FALSE,protocol_MD5=input_hashes[role=='protocol',md5])
write_json(qa,file.path(out,'logs','analysis_QA.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'logs','R_sessionInfo.txt'))
print(res[outcome!='balance_score',.(subset,spec,outcome,n_samples,SCI_n,TC_n,estimate,lower_CL,upper_CL,robust_SE_HC3,hat_max,group_VIF)])
print(loosum); print(rbindlist(overlap)); print(infsum)
cat('Completed UTC:',format(Sys.time(),tz='UTC'),'\nPASS\n')
sink()
