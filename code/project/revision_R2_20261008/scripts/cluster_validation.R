# Public-package copy: see metadata/code_changes.csv for provenance.
args<-commandArgs(TRUE);root<-normalizePath(if(length(args))args[1] else '.',winslash='/');out<-file.path(root,'revision_R2_20261008');r1<-file.path(root,'revision_R1_20261007')
suppressPackageStartupMessages({library(data.table);library(cluster)})
set.seed(20261008);sink(file.path(out,'logs/cluster_validation.txt'),split=TRUE)
s<-readRDS(file.path(r1,'tables/analysis_state.rds'));cl<-fread(file.path(r1,'tables/original_S3_Cluster_Assignment.csv'));cl$Gene[cl$Gene=='SEPT15']<-'SEP15'
stopifnot(!anyDuplicated(cl$Gene),all(cl$Gene %in% rownames(s$discovery)))
E<-s$discovery[cl$Gene,,drop=FALSE];Z<-t(scale(t(E)));stopifnot(all(is.finite(Z)))
write<-function(d,n)fwrite(d,file.path(out,'tables',paste0(n,'.csv')))
ari<-function(a,b){tt<-table(a,b);f<-function(x)sum(x*(x-1)/2);n<-sum(tt);obs<-f(tt);expected<-f(rowSums(tt))*f(colSums(tt))/(n*(n-1)/2);(obs-expected)/((f(rowSums(tt))+f(colSums(tt)))/2-expected)}
fits<-list();metrics<-list();init<-list();D<-dist(Z)
for(k in 2:6){
  fit<-kmeans(Z,k,nstart=50,iter.max=200);fits[[as.character(k)]]<-fit
  metrics[[as.character(k)]]<-data.frame(k=k,silhouette=mean(silhouette(fit$cluster,D)[,3]),withinss=fit$tot.withinss,min_size=min(fit$size),max_size=max(fit$size),n_genes=nrow(Z),n_participants=ncol(Z))
  for(b in 1:20){fb<-kmeans(Z,k,nstart=1,iter.max=200);init[[paste(k,b)]]<-data.frame(k=k,iteration=b,ARI=ari(fit$cluster,fb$cluster),withinss=fb$tot.withinss)}
  cat('reference k',k,'completed\n');flush.console()
}
write(rbindlist(metrics),'cluster_metrics');write(rbindlist(init),'cluster_initializations')
write(data.frame(gene=rownames(Z),original=cl$Cluster,k2=fits[['2']]$cluster,k3=fits[['3']]$cluster,k4=fits[['4']]$cluster,k5=fits[['5']]$cluster,k6=fits[['6']]$cluster),'cluster_assignments')
cat('gap B100 started\n');flush.console()
gap<-clusGap(Z,FUNcluster=function(x,k)list(cluster=kmeans(x,k,nstart=10,iter.max=200)$cluster),K.max=6,B=100,d.power=2,spaceH0='scaledPCA',verbose=FALSE)
write(data.frame(k=1:6,gap$Tab),'cluster_gap')
gap_k<-maxSE(gap$Tab[,'gap'],gap$Tab[,'SE.sim'],method='Tibs2001SEmax');cat('gap selected k',gap_k,'\n');flush.console()
boots<-list();overall<-list();groups<-split(seq_len(ncol(E)),s$s151$group)
for(b in 1:200){
  ii<-unlist(lapply(groups,function(ix)sample(ix,length(ix),replace=TRUE)),use.names=FALSE);zb<-t(scale(t(E[,ii,drop=FALSE])));stopifnot(all(is.finite(zb)))
  for(k in 2:6){
    fb<-kmeans(zb,k,nstart=25,iter.max=200);a<-fits[[as.character(k)]]$cluster
    overall[[paste(b,k)]]<-data.frame(iteration=b,k=k,ARI=ari(a,fb$cluster))
    for(j in 1:k){ref<-which(a==j);jac<-vapply(1:k,function(l){bb<-which(fb$cluster==l);length(intersect(ref,bb))/length(union(ref,bb))},numeric(1));boots[[paste(b,k,j)]]<-data.frame(iteration=b,k=k,cluster=j,n_genes=length(ref),max_jaccard=max(jac))}
  }
  if(b%%20==0){cat('bootstrap',b,'/200\n');flush.console()}
}
boot<-rbindlist(boots);write(boot,'cluster_bootstrap');write(rbindlist(overall),'cluster_bootstrap_ARI')
summ<-boot[,.(n_genes=first(n_genes),mean_jaccard=mean(max_jaccard),median_jaccard=median(max_jaccard),q025=quantile(max_jaccard,.025),q975=quantile(max_jaccard,.975),fraction_below_half=mean(max_jaccard<.5)),by=.(k,cluster)];write(summ,'cluster_bootstrap_summary')
write(data.frame(silhouette_k=which.max(vapply(metrics,function(x)x$silhouette,numeric(1)))+1,gap_k=gap_k,original_vs_refit_k3_ARI=ari(cl$Cluster,fits[['3']]$cluster),n_genes=nrow(E),n_samples=ncol(E),bootstrap_iterations=200,gap_references=100),'cluster_decision')
saveRDS(list(fits=fits,gap=gap),file.path(out,'tables/cluster_state.rds'));capture.output(sessionInfo(),file=file.path(out,'logs/cluster_sessionInfo.txt'));cat('COMPLETE\n');sink()
