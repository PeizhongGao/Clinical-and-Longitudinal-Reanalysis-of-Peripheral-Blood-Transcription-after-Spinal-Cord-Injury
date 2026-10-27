# Public-package copy: see metadata/code_changes.csv for provenance.
args<-commandArgs(trailingOnly=TRUE)
root<-normalizePath(if(length(args))args[1] else getwd(),winslash='/',mustWork=TRUE)
out<-file.path(root,'revision_numeric_freeze_20261019')
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
specs<-list(
 list('BTM_discovery_all88.csv','SE','lower_CL','upper_CL','p_value'),
 list('BTM_longitudinal_all240.csv','SE','lower_CL','upper_CL','p_value'),
 list('original15plus15_fixed_sample27.csv','robust_SE_HC3','lower_CL','upper_CL','p_value'),
 list('original15plus15_longitudinal36.csv','SE','lower.CL','upper.CL','p.value'),
 list('original15plus15_longitudinal_sensitivity42.csv','SE','lower.CL','upper.CL','p.value'),
 list('original15plus15_clinical_composition15.csv','robust_SE_HC3','lower_CL','upper_CL','p_value'))
checks<-lapply(specs,function(s){
 d<-fread(file.path(out,'tables',s[[1]])); se<-d[[s[[2]]]];df<-d$df
 plo<-2*pt(-abs(d$estimate/se),df);lo<-d$estimate-qt(.975,df)*se;hi<-d$estimate+qt(.975,df)*se
 dp<-max(abs(plo-d[[s[[5]]]]));dc<-max(abs(c(lo-d[[s[[3]]]],hi-d[[s[[4]]]])))
 data.table(file=s[[1]],n=nrow(d),max_abs_P=dp,max_abs_CI=dc,passed=dp<1e-10&&dc<1e-10)
})
z<-rbindlist(checks);fwrite(z,file.path(out,'qa','consolidated_inference_QA.csv'))
write_json(list(passed=all(z$passed),n_rows=sum(z$n),checks=z,scope='Recompute two-sided t P and pointwise 95% CI from saved estimates, SE and df. Existing mixed-model fits additionally have previous cross-engine QA. This check is arithmetic, not independent re-estimation of all fits.'),file.path(out,'qa','consolidated_inference_QA.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
print(z);stopifnot(all(z$passed))
