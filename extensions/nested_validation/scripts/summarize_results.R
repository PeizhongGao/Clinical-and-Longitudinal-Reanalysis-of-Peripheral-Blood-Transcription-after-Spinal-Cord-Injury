args<-commandArgs(trailingOnly=TRUE);ROOT<-normalizePath(args[1],winslash='/',mustWork=TRUE)
source(file.path(ROOT,'scripts/nested_core.R'))
m<-fread(file.path(ROOT,'tables/nested_repeat_metrics.csv'));pr<-fread(file.path(ROOT,'tables/nested_oof_predictions.csv'));au<-fread(file.path(ROOT,'tables/nested_outer_fit_audit.csv'))
stopifnot(nrow(m)==80,nrow(pr)==3840,nrow(au)==400)
summaries<-cirows<-freqrows<-jacrows<-list();set.seed(20261029)
for(src in c('raw','processed'))for(cn in c('SCI_HC','SCI_TC')) {
  d<-m[source==src & contrast==cn];a<-au[source==src & contrast==cn];p<-pr[source==src & contrast==cn,.(y=unique(y),prediction=mean(prediction),baseline=mean(baseline)),by=.(sample,subject)]
  stopifnot(nrow(p)==48,sum(p$y)==38,all(p$sample==sort(p$sample)))
  agg<-metrics(p$y,p$prediction);draws<-matrix(NA_real_,5000,4,dimnames=list(NULL,names(agg)))
  pos<-which(p$y==1);neg<-which(p$y==0)
  for(b in 1:5000) {ix<-c(sample(pos,length(pos),TRUE),sample(neg,length(neg),TRUE));draws[b,]<-metrics(p$y[ix],p$prediction[ix])}
  for(k in names(agg)) {
    ci<-quantile(draws[,k],c(.025,.975),names=FALSE)
    cirows[[paste(src,cn,k)]]<-data.table(source=src,contrast=cn,metric=k,mean_oof_prediction_metric=agg[k],conditional_interval_low=ci[1],conditional_interval_high=ci[2],bootstrap_resamples=5000,interval_scope='stratified participant bootstrap of fixed mean OOF predictions; excludes model-training uncertainty')
  }
  s<-data.table(source=src,contrast=cn,n_participants=48,SCI_n=38,control_n=10,repeats=20,outer_folds=5,inner_folds=3,mean_repeat_AUC=mean(d$AUC),sd_repeat_AUC=sd(d$AUC),min_repeat_AUC=min(d$AUC),max_repeat_AUC=max(d$AUC),mean_repeat_AP=mean(d$average_precision),mean_repeat_Brier=mean(d$Brier),mean_repeat_log_loss=mean(d$log_loss),mean_baseline_Brier=mean(d$baseline_Brier),mean_baseline_log_loss=mean(d$baseline_log_loss),mean_selected_genes=mean(a$selected_n),min_selected_genes=min(a$selected_n),max_selected_genes=max(a$selected_n),intercept_only_fits=sum(a$selected_n==0),outer_fits=nrow(a),inner_small_class_warnings=sum(a$inner_small_class_warnings),mean_oof_AUC=agg['AUC'])
  summaries[[paste(src,cn)]]<-s
  lists<-strsplit(ifelse(is.na(a$selected),'',a$selected),';',fixed=TRUE);lists<-lapply(lists,function(x)x[nzchar(x)])
  genes<-rownames(read_dataset(ROOT,src,cn)$E);cnt<-table(unlist(lists));fr<-data.table(source=src,contrast=cn,gene=genes,selected_count=as.integer(cnt[match(genes,names(cnt))]));fr[is.na(selected_count),selected_count:=0L];fr[,`:=`(outer_training_fits=100,selection_frequency=selected_count/100)];freqrows[[paste(src,cn)]]<-fr
  comb<-combn(1:100,2);jacs<-apply(comb,2,function(ix){u<-union(lists[[ix[1]]],lists[[ix[2]]]);if(!length(u))1 else length(intersect(lists[[ix[1]]],lists[[ix[2]]]))/length(u)})
  jacrows[[paste(src,cn)]]<-data.table(source=src,contrast=cn,mean_jaccard=mean(jacs),median_jaccard=median(jacs),min_jaccard=min(jacs),max_jaccard=max(jacs),dependent_fit_pairs=length(jacs),scope='descriptive overlap; fits and pairs are dependent')
}
summary<-rbindlist(summaries);fwrite(summary,file.path(ROOT,'tables/nested_performance_summary.csv'));fwrite(rbindlist(cirows),file.path(ROOT,'tables/nested_conditional_intervals.csv'));fwrite(rbindlist(freqrows),file.path(ROOT,'tables/nested_gene_selection_frequencies.csv'));fwrite(rbindlist(jacrows),file.path(ROOT,'tables/nested_selection_overlap.csv'))
fwrite(pr[,.(y=unique(y),mean_prediction=mean(prediction),mean_baseline=mean(baseline)),by=.(source,contrast,sample,subject)],file.path(ROOT,'tables/nested_mean_oof_predictions.csv'))
if(!file.exists(file.path(ROOT,'tables/full_pipeline_permutation_null.csv'))) {print(summary);cat('OBSERVED SUMMARIES COMPLETE; PERMUTATION SUMMARY PENDING\n');quit(status=0)}
nu<-fread(file.path(ROOT,'tables/full_pipeline_permutation_null.csv'));stopifnot(nrow(nu)==1998,all(nu[, .N,by=contrast]$N==999))
pres<-rbindlist(lapply(c('SCI_HC','SCI_TC'),function(cn){
  obs<-m[source=='raw' & contrast==cn & repeat_id==1,AUC];v<-nu[contrast==cn,AUC];ex<-sum(v>=obs-1e-12);ci<-binom.test(ex,length(v))$conf.int
  data.table(contrast=cn,observed_first_repeat_AUC=obs,null_mean_AUC=mean(v),null_sd_AUC=sd(v),permutations=length(v),exceedances=ex,p_permutation=(ex+1)/(length(v)+1),mc_binomial_lower=ci[1],mc_binomial_upper=ci[2],exchangeability='labels within sex; marginal association conditional on sex, not clinically adjusted')
}));pres[,p_Holm2:=p.adjust(p_permutation,'holm')];fwrite(pres,file.path(ROOT,'tables/full_pipeline_permutation_results.csv'))
jsonlite::write_json(list(nested_summary=summary,conditional_intervals=rbindlist(cirows),selection_overlap=rbindlist(jacrows),full_pipeline_permutation=pres,clinical=fread(file.path(ROOT,'tables/clinical_resampling_results.csv')),clinical_diagnostics=fread(file.path(ROOT,'tables/clinical_resampling_diagnostics.csv'))),file.path(ROOT,'audit/result_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
fwrite(data.table(family=c('Nested_raw_classifier_Holm2','Clinical_BTM_FreedmanLane_Holm4','Clinical_BTM_wild_Holm4'),tests=c(2L,4L,4L),method=c('within-sex complete-pipeline label permutation','within-sex nuisance-residual permutation','HC2-adjusted Rademacher null wild bootstrap'),resamples_per_test=c(999L,9999L,9999L),p_field=c('p_permutation','p_resampling','p_resampling'),adjusted_field=c('p_Holm2','p_Holm4','p_Holm4'),original_BH_families_unchanged=TRUE),file.path(ROOT,'tables/new_resampling_family_registry.csv'))
print(summary);print(pres);cat('SUMMARY COMPLETE\n')
