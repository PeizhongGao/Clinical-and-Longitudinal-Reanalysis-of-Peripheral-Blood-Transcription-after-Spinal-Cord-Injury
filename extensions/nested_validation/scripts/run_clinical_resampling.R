args<-commandArgs(trailingOnly=TRUE);ROOT<-normalizePath(args[1],winslash='/',mustWork=TRUE)
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(lmtest)})
all<-fread(file.path(ROOT,'inputs/btm_model_inputs.csv'));frozen<-fread(file.path(ROOT,'inputs/btm_frozen_models.csv'))
mods<-c('M37.1','M163','M7.0','M223');BPERM<-9999L
d0<-all[subset=='clinical43' & spec=='clinical' & variant=='published' & module_id==mods[1]][order(sample)]
stopifnot(nrow(d0)==43,uniqueN(d0$subject)==43,sum(d0$group=='SCI')==36,sum(d0$group=='TC')==7)
for(v in c('tbi','prior_cns'))d0[,(v):=factor(get(v),levels=c('No','Yes','Unknown'))]
d0[,group:=factor(group,levels=c('TC','SCI'))];d0[,sex:=factor(sex,levels=c('F','M'))]
X<-model.matrix(~group+age_c+sex+iss_c+tbi+prior_cns,d0);j<-match('groupSCI',colnames(X));Z<-X[,-j,drop=FALSE]
stopifnot(qr(X)$rank==9,ncol(X)==9,qr(Z)$rank==8)
inv<-solve(crossprod(X));A<-inv%*%t(X);h<-rowSums(X*t(A));az<-solve(crossprod(Z),t(Z));h0<-rowSums(Z*t(az))
stopifnot(all(h<1-1e-8),all(h0<1-1e-8))
studentized<-function(Y) {
  Y<-as.matrix(Y);bb<-A%*%Y;ee<-Y-X%*%bb
  se<-sqrt(colSums((ee*(A[j,]/(1-h)))^2))
  stopifnot(all(is.finite(se)),all(se>0))
  list(t=bb[j,]/se,estimate=bb[j,],se=se,residual=ee)
}
set.seed(20261030);perms<-matrix(seq_len(nrow(d0)),nrow(d0),BPERM)
for(b in seq_len(BPERM))for(s in levels(d0$sex)) {ix<-which(d0$sex==s);perms[ix,b]<-sample(ix)}
set.seed(20261031);weights<-matrix(sample(c(-1,1),nrow(d0)*BPERM,replace=TRUE),nrow(d0),BPERM)
saveRDS(list(samples=d0$sample,sex=as.character(d0$sex),permutations=perms,wild_weights=weights),file.path(ROOT,'tables/clinical_resampling_assignments.rds'))
results<-nulls<-diagnostics<-inputs<-list();verification<-list()
for(id in mods) {
  d<-all[subset=='clinical43' & spec=='clinical' & variant=='published' & module_id==id][match(d0$sample,sample)]
  stopifnot(identical(d$sample,d0$sample));y<-d$score;obs<-studentized(y)
  ref<-frozen[subset=='clinical43' & spec=='clinical' & variant=='published' & module_id==id]
  stopifnot(nrow(ref)==1,abs(obs$estimate-ref$estimate)<1e-10,abs(obs$se-ref$SE)<1e-10)
  f0<-as.vector(Z%*%(az%*%y));e0<-y-f0
  yperm<-matrix(f0,nrow(d0),BPERM)+matrix(e0[perms],nrow(d0),BPERM)
  ywild<-matrix(f0,nrow(d0),BPERM)+weights*as.vector(e0/sqrt(1-h0))
  tp<-studentized(yperm)$t;tw<-studentized(ywild)$t
  bp<-lmtest::bptest(lm.fit<-lm(y~X-1))$p.value
  for(method in c('Freedman_Lane_within_sex','null_wild_Rademacher_HC2')) {
    tt<-if(method=='Freedman_Lane_within_sex')tp else tw
    exceed<-sum(abs(tt)>=abs(obs$t)-1e-12);p<-(exceed+1)/(BPERM+1);ci<-binom.test(exceed,BPERM)$conf.int
    key<-paste(id,method)
    results[[key]]<-data.table(module_id=id,method=method,n=43,SCI_n=36,TC_n=7,estimate=obs$estimate,HC3_SE=obs$se,df=34,lower_CL=ref$lower_CL,upper_CL=ref$upper_CL,observed_t=obs$t,original_p=ref$p_value,original_q_BH16=ref$q_BH_primary16,resamples=BPERM,exceedances=exceed,p_resampling=p,mc_binomial_lower=ci[1],mc_binomial_upper=ci[2])
    nulls[[key]]<-data.table(module_id=id,method=method,resample_id=seq_len(BPERM),t=tt)
  }
  diagnostics[[id]]<-data.table(module_id=id,full_rank=ncol(X),nuisance_rank=ncol(Z),max_leverage=max(h),max_null_leverage=max(h0),group_VIF=sum((X[,j]-mean(X[,j]))^2)/sum(qr.resid(qr(Z),X[,j])^2),full_BP_p=bp,full_shapiro_p=shapiro.test(obs$residual[,1])$p.value,null_shapiro_p=shapiro.test(e0)$p.value,null_residual_variance_SCI=var(e0[d$group=='SCI']),null_residual_variance_TC=var(e0[d$group=='TC']))
  inputs[[id]]<-data.table(module_id=id,sample=d$sample,group=d$group,sex=as.character(d0$sex),score=y,null_fitted=f0,null_residual=e0,full_residual=obs$residual[,1],full_hat=h,null_hat=h0)
  verification[[id]]<-list(estimate_difference=unname(obs$estimate-ref$estimate),SE_difference=unname(obs$se-ref$SE),passed=TRUE)
}
r<-rbindlist(results);r[,p_Holm4:=p.adjust(p_resampling,'holm'),by=method]
fwrite(r,file.path(ROOT,'tables/clinical_resampling_results.csv'));fwrite(rbindlist(nulls),file.path(ROOT,'tables/clinical_resampling_null.csv'))
fwrite(rbindlist(diagnostics),file.path(ROOT,'tables/clinical_resampling_diagnostics.csv'));fwrite(rbindlist(inputs),file.path(ROOT,'tables/clinical_resampling_inputs.csv'))
fwrite(cbind(data.table(sample=d0$sample),as.data.table(X)),file.path(ROOT,'tables/clinical_design_matrix.csv'))
write_json(list(status='PASS',frozen_estimates_reproduced=verification,full_rank=ncol(X),nuisance_rank=ncol(Z),permutation_count=BPERM,wild_count=BPERM,shared_resampling_across_modules=TRUE,all_null_statistics_finite=TRUE),file.path(ROOT,'audit/clinical_resampling_verification.json'),pretty=TRUE,auto_unbox=TRUE)
print(r[,.(module_id,method,estimate,original_q_BH16,p_resampling,p_Holm4)])
cat('CLINICAL RESAMPLING COMPLETE\n')
