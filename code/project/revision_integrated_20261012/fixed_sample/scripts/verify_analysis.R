# Public-package copy: see metadata/code_changes.csv for provenance.
# Numerical QA using direct participant-deletion refits, separate from HC3 matrix code.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args))args[1] else getwd(),winslash='/',mustWork=TRUE)
out <- file.path(root,'revision_R2_20261008','fixed_sample_20261012')
fits <- readRDS(file.path(out,'logs','fitted_models.rds'))
checks <- list()
for(k in names(fits)){
  z <- fits[[k]]; X <- z$X; y <- z$fit$y; j <- match('groupSCI',colnames(X))
  deleted <- vapply(seq_len(nrow(X)),function(i){
    fit <- lm.fit(X[-i,,drop=FALSE],y[-i]); stopifnot(fit$rank==ncol(X));fit$coefficients[j]
  },numeric(1))
  delta <- deleted-z$result$estimate
  # HC3 variance equals the uncentered sum of squared deletion coefficient changes.
  direct_SE <- sqrt(sum(delta^2))
  analytic_delta <- z$influence$group_beta_change_if_deleted
  checks[[k]] <- data.table(model=k,n_subject_refits=length(deleted),
    deletion_SE=direct_SE,matrix_HC3_SE=z$result$robust_SE_HC3,
    SE_difference=abs(direct_SE-z$result$robust_SE_HC3),max_delta_difference=max(abs(delta-analytic_delta)))
}
checks <- rbindlist(checks)
stopifnot(max(checks$SE_difference)<1e-9,max(checks$max_delta_difference)<1e-9)
fwrite(checks,file.path(out,'tables','HC3_direct_refit_verification.csv'))
models <- fread(file.path(out,'tables','fixed_sample_models.csv'))
scores <- fread(file.path(out,'inputs','ODC405_locked_scores_with_clinical_and_cell_estimates.csv'))
sets <- list(full48=fread(file.path(out,'tables','fixed_input_full48.csv')),
             clinical43=fread(file.path(out,'tables','fixed_input_clinical43.csv')),
             time41=fread(file.path(out,'tables','fixed_input_time41.csv')))
for(s in names(sets)){
  dat <- sets[[s]]
  for(o in c('innate_score','adaptive_score','balance_score')){
    raw <- dat[group=='SCI',mean(get(o))]-dat[group=='TC',mean(get(o))]
    beta <- models[subset==s & spec=='unadjusted' & outcome==o,estimate]
    stopifnot(abs(raw-beta)<1e-10)
    original <- scores[[o]][match(dat$sample,scores$sample)]
    stopifnot(identical(original,dat[[o]]))
  }
}
for(s in names(sets)) for(spec_name in unique(models[subset==s,spec])){
  z <- models[subset==s & spec==spec_name]
  stopifnot(abs(z[outcome=='balance_score',estimate]-(z[outcome=='innate_score',estimate]-z[outcome=='adaptive_score',estimate]))<1e-10)
}
qa <- list(status='PASS',n_subject_deletion_refits=sum(checks$n_subject_refits),
  HC3_SE_max_difference=max(checks$SE_difference),deletion_effect_max_difference=max(checks$max_delta_difference),
  raw_mean_difference_matches_unadjusted=TRUE,frozen_score_values_bit_identical=TRUE,balance_coefficient_identity=TRUE)
write_json(qa,file.path(out,'logs','independent_numeric_QA.json'),pretty=TRUE,auto_unbox=TRUE)
print(qa)
