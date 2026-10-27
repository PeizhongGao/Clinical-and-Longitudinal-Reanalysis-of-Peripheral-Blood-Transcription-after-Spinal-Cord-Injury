args<-commandArgs(trailingOnly=TRUE);ROOT<-normalizePath(args[1],winslash='/',mustWork=TRUE);mode<-if(length(args)>1)args[2] else 'pilot'
source(file.path(ROOT,'scripts/nested_core.R'))
dir.create(file.path(ROOT,'logs/permutations'),recursive=TRUE,showWarnings=FALSE)
if(mode=='pilot') {
  ds<-read_dataset(ROOT);tr<-1:38;te<-39:48
  # Deterministic explicit stratification produces a valid test of held-out invariance.
  f<-stratified_folds(ds$y,5,20261027);tr<-which(f!=1);te<-which(f==1)
  a<-training_state(ds,tr,ds$y);altered<-ds;altered$E[,te]<-altered$E[,te]*8+100;altered$CPM[,te]<-altered$CPM[,te]*50
  yp<-ds$y;yp[te]<-1L-yp[te];b<-training_state(altered,tr,yp)
  stopifnot(identical(a,b));fit1<-fit_path(a,ds$y[tr]);fit2<-fit_path(b,yp[tr]);stopifnot(identical(fit1,fit2))
  set.seed(1);xr<-rnorm(100);stopifnot(abs(auc_rank(rep(0:1,50),xr)-as.numeric(pROC::auc(rep(0:1,50),xr,direction='<',quiet=TRUE)))<1e-12)
  tm<-system.time(z<-nested_once(ds))
  jsonlite::write_json(list(status='PASS',heldout_expression_label_invariance=TRUE,training_state_fit_invariance=TRUE,auc_formula_check=TRUE,fold_disjointness_assertions=TRUE,seconds=unname(tm['elapsed']),n=length(ds$y),source_gene_count=nrow(ds$E)),file.path(ROOT,'audit/implementation_tests.json'),pretty=TRUE,auto_unbox=TRUE)
  saveRDS(z,file.path(ROOT,'audit/pilot_result.rds'));cat('PILOT PASSED; nested repeat seconds:',tm['elapsed'],'\n')
}
if(mode=='observed') {
  allmetrics<-allpred<-allaudit<-allsplit<-allinner<-list();counter<-0L
  for(source in c('raw','processed')) for(contrast in c('SCI_HC','SCI_TC')) {
    ds<-read_dataset(ROOT,source,contrast)
    # Comparison-specific seeds; processed sensitivity uses matching splits for the same contrast.
    offset<-if(contrast=='SCI_HC')0 else 100000
    for(rep in 1:20) {
      z<-nested_once(ds,seed=20261027+offset+(rep-1)*10000)
      tag<-list(source=source,contrast=contrast,repeat_id=rep);counter<-counter+1L
      allmetrics[[counter]]<-cbind(as.data.table(tag),as.data.table(as.list(z$metrics)),as.data.table(setNames(as.list(z$baseline),paste0('baseline_',names(z$baseline)))))
      allpred[[counter]]<-cbind(as.data.table(tag),z$predictions);allaudit[[counter]]<-cbind(as.data.table(tag),z$audit)
      allsplit[[counter]]<-cbind(as.data.table(tag),z$splits);allinner[[counter]]<-cbind(as.data.table(tag),z$inner_audit)
      cat(source,contrast,'repeat',rep,'of 20 completed\n');flush.console()
    }
  }
  fwrite(rbindlist(allmetrics),file.path(ROOT,'tables/nested_repeat_metrics.csv'))
  write_predictions(rbindlist(allpred),file.path(ROOT,'tables/nested_oof_predictions.csv'))
  fwrite(rbindlist(allaudit),file.path(ROOT,'tables/nested_outer_fit_audit.csv'))
  fwrite(rbindlist(allsplit),file.path(ROOT,'tables/nested_split_assignments.csv'))
  fwrite(rbindlist(allinner),file.path(ROOT,'tables/nested_inner_screen_audit.csv'))
  writeLines(capture.output(sessionInfo()),file.path(ROOT,'logs/R_sessionInfo.txt'))
  cat('OBSERVED COMPLETE\n')
}
run_perm_one <- function(job) {
  fname<-file.path(ROOT,'logs/permutations',paste0(job$contrast,'_',sprintf('%04d',job$id),'.rds'))
  if(file.exists(fname))return(invisible(TRUE))
  ds<-DATASETS[[job$contrast]];offset<-if(job$contrast=='SCI_HC')0 else 100000
  yp<-permute_sex(ds$y,ds$meta$sex,20261028+offset+job$id)
  z<-nested_once(ds,y=yp,seed=20261027+offset,detail=FALSE)
  saveRDS(list(id=job$id,contrast=job$contrast,seed=20261028+offset+job$id,labels=yp,metrics=z$metrics,audit=z$audit),fname)
  TRUE
}
if(mode=='permutations') {
  DATASETS<-setNames(lapply(c('SCI_HC','SCI_TC'),function(v)read_dataset(ROOT,'raw',v)),c('SCI_HC','SCI_TC'))
  jobs<-list();for(cn in names(DATASETS))for(i in 1:999)jobs[[length(jobs)+1L]]<-list(contrast=cn,id=i)
  jobs<-Filter(function(j)!file.exists(file.path(ROOT,'logs/permutations',paste0(j$contrast,'_',sprintf('%04d',j$id),'.rds'))),jobs)
  cl<-parallel::makePSOCKcluster(6,outfile=file.path(ROOT,'logs/permutation_workers.log'))
  parallel::clusterExport(cl,c('ROOT','DATASETS','run_perm_one'),envir=environment())
  parallel::clusterEvalQ(cl,{source(file.path(ROOT,'scripts/nested_core.R'));data.table::setDTthreads(1);NULL})
  for(start in seq.int(1,max(1,length(jobs)),by=40)) {
    if(!length(jobs))break
    batch<-jobs[start:min(start+39,length(jobs))]
    ok<-parallel::parLapplyLB(cl,batch,run_perm_one);stopifnot(all(unlist(ok)))
    cat('PERMUTATIONS completed files:',length(list.files(file.path(ROOT,'logs/permutations'),pattern='rds$')),'/ 1998\n');flush.console()
  }
  parallel::stopCluster(cl)
  res<-lapply(list.files(file.path(ROOT,'logs/permutations'),pattern='rds$',full.names=TRUE),readRDS)
  stopifnot(length(res)==1998)
  null<-rbindlist(lapply(res,function(z)cbind(data.table(contrast=z$contrast,permutation_id=z$id,seed=z$seed),as.data.table(as.list(z$metrics)))))
  fwrite(null,file.path(ROOT,'tables/full_pipeline_permutation_null.csv'))
  fwrite(rbindlist(lapply(res,function(z)cbind(data.table(contrast=z$contrast,permutation_id=z$id),z$audit))),file.path(ROOT,'tables/permutation_outer_fit_audit.csv'))
  fwrite(rbindlist(lapply(res,function(z)data.table(contrast=z$contrast,permutation_id=z$id,sample=DATASETS[[z$contrast]]$meta$sample,permuted_label=z$labels))),file.path(ROOT,'tables/permutation_label_assignments.csv'))
  cat('FULL PIPELINE PERMUTATIONS COMPLETE\n')
}
