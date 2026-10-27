# Public-package copy: see metadata/code_changes.csv for provenance.
# Locked revision-stage external module analysis; R1-10 remains frozen.
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(if(length(args))args[1] else '.',winslash='/')
out<-file.path(root,'revision_R1_20261007')
suppressPackageStartupMessages({library(data.table);library(AnnotationDbi);library(org.Hs.eg.db);library(edgeR);library(lmerTest);library(emmeans);library(jsonlite)})
stopifnot(file.exists(file.path(out,'tables/external_module_protocol_lock.json')))
sink(file.path(out,'logs/external_BTM_analysis.txt'),split=TRUE)
cat('Analysis started UTC:',format(Sys.time(),tz='UTC',usetz=TRUE),'\n')
write_table<-function(x,n)fwrite(as.data.frame(x),file.path(out,'tables',paste0(n,'.csv')),na='NA')
members<-as.data.frame(fread(file.path(out,'tables/external_BTM_source_members.csv')))
locked<-as.data.frame(fread(file.path(out,'tables/candidate_gene_pools.csv')))
source(file.path(out,'scripts/candidate_pools.R'))
valid_symbols<-keys(org.Hs.eg.db,keytype='SYMBOL');valid_alias<-keys(org.Hs.eg.db,keytype='ALIAS')
map_symbols<-function(x){
  x<-unique(x);direct<-x%in%valid_symbols;lookup<-x[!direct & !grepl('///',x,fixed=TRUE) & x%in%valid_alias]
  aliases<-if(length(lookup))suppressMessages(AnnotationDbi::select(org.Hs.eg.db,keys=lookup,keytype='ALIAS',columns='SYMBOL')) else data.frame(ALIAS=character(),SYMBOL=character())
  split_alias<-split(aliases$SYMBOL,aliases$ALIAS)
  do.call(rbind,lapply(seq_along(x),function(i){
    g<-x[i];opts<-unique(na.omit(split_alias[[g]]));reason<-if(grepl('///',g,fixed=TRUE))'ambiguous_source_token' else if(direct[i])'exact_SYMBOL' else if(length(opts)==1)'unique_ALIAS' else if(length(opts)>1)'ambiguous_ALIAS' else 'unmapped'
    data.frame(source_token=g,symbol=if(reason=='exact_SYMBOL')g else if(reason=='unique_ALIAS')opts else NA_character_,mapping=reason)
  }))
}
mapping<-map_symbols(unique(c(members$source_token,locked$gene)));write_table(mapping,'external_BTM_gene_mapping')
members<-merge(members,mapping,by='source_token',all.x=TRUE,sort=FALSE)
pool<-unique(na.omit(mapping$symbol[match(locked$gene,mapping$source_token)]))
members$overlaps_original_106<-!is.na(members$symbol)&members$symbol%in%pool
write_table(members,'external_BTM_members_annotated')

meta<-as.data.frame(fread(file.path(out,'tables/sample_manifest.csv')))
d<-fread(file.path(root,'explore/Report/GSE151371_normalized_expression_matrix_comBat_de-ID.csv'))
E<-as.matrix(d[,-1,with=FALSE]);rownames(E)<-d[[1]];storage.mode(E)<-'double';E<-E[!duplicated(rownames(E)),,drop=FALSE]
tested<-fread(file.path(out,'tables/SCI_vs_TC_processed_limma.csv'))$gene;stopifnot(all(tested%in%rownames(E)));E<-E[tested,,drop=FALSE]
emap<-map_symbols(rownames(E));write_table(emap,'external_BTM_discovery_symbol_mapping')
sy<-emap$symbol[match(rownames(E),emap$source_token)];ok<-!is.na(sy)
E151<-rowsum(E[ok,,drop=FALSE],sy[ok],reorder=FALSE);E151<-E151/as.numeric(table(sy[ok])[rownames(E151)])
m151<-meta[meta$dataset=='GSE151371',];m151<-m151[match(colnames(E151),m151$sample),];stopifnot(!anyNA(m151$sample))

raw<-read.delim(gzfile(file.path(root,'data/GSE226238/GSE226238_raw_counts_GRCh38.p13_NCBI.tsv.gz')),check.names=FALSE,colClasses=c('character',rep('numeric',46)))
ann<-read.delim(gzfile(file.path(root,'data/GSE226238/Human.GRCh38.p13.annot.tsv.gz')),check.names=FALSE,colClasses='character')
m226<-as.data.frame(fread(file.path(out,'tables/NCBI_GSE226238_sample_map.csv')))
counts<-as.matrix(raw[,-1]);rownames(counts)<-raw$GeneID
stopifnot(identical(colnames(counts),m226$gsm))
symbols<-ann$Symbol[match(rownames(counts),ann$GeneID)]
amap<-map_symbols(symbols);write_table(amap,'external_BTM_NCBI_symbol_mapping');symbols<-amap$symbol[match(symbols,amap$source_token)]
requested<-unique(na.omit(members$symbol));take<-!is.na(symbols)&symbols%in%requested
mc<-rowsum(counts[take,,drop=FALSE],symbols[take],reorder=FALSE)
norm<-as.data.frame(fread(file.path(out,'tables/NCBI_GSE226238_normalization.csv')));stopifnot(identical(norm$gsm,m226$gsm))
lib<-norm$lib.size*norm$norm.factors
E226<-cpm(mc,lib.size=lib,log=TRUE,prior.count=.5);colnames(E226)<-m226$sample
write_table(data.frame(gene=rownames(E226),mc,check.names=FALSE),'external_BTM_NCBI_aggregated_counts')
write_table(data.frame(gene=rownames(E226),E226,check.names=FALSE),'external_BTM_NCBI_logCPM')
write_table(data.frame(gene=rownames(E151)[rownames(E151)%in%requested],E151[rownames(E151)%in%requested,,drop=FALSE],check.names=FALSE),'external_BTM_discovery_log_expression')
coverage<-list();scores<-list();usedrows<-list()
for(dataset in c('GSE151371','GSE226238'))for(variant in c('published','exclude_original106'))for(id in unique(members$module_id)){
  mm<-members[members$module_id==id,];if(variant=='exclude_original106')mm<-mm[!mm$overlaps_original_106,]
  gs<-unique(na.omit(mm$symbol));X<-if(dataset=='GSE151371')E151 else E226;md<-if(dataset=='GSE151371')m151 else m226
  g<-intersect(gs,rownames(X));g<-g[vapply(g,function(z)all(is.finite(X[z,])),logical(1))]
  denominator<-nrow(mm);cov<-if(denominator)length(g)/denominator else 0;eligible<-length(g)>=5 & cov>=.8
  ky<-paste(dataset,variant,id)
  coverage[[ky]]<-data.frame(dataset=dataset,variant=variant,module_id=id,lineage=members$lineage[match(id,members$module_id)],n_source_tokens=denominator,n_unique_mapped=length(gs),n_used=length(g),coverage=cov,eligible=eligible,missing_tokens=paste(mm$source_token[is.na(mm$symbol)|!mm$symbol%in%g],collapse=';'),genes_used=paste(g,collapse=';'))
  if(length(g))usedrows[[ky]]<-data.frame(dataset=dataset,variant=variant,module_id=id,gene=g)
  s<-md;s$score<-if(eligible)colMeans(X[g,,drop=FALSE]) else NA_real_;s$variant<-variant;s$module_id<-id;s$lineage<-members$lineage[match(id,members$module_id)];scores[[ky]]<-s
}
covtab<-rbindlist(coverage,fill=TRUE);sc<-as.data.frame(rbindlist(scores,fill=TRUE));write_table(covtab,'external_BTM_coverage');write_table(sc,'external_BTM_scores');write_table(rbindlist(usedrows),'external_BTM_used_genes')
visits<-c('acute','3m','6m','12m');ids<-unique(m226$subject[m226$group=='SCI']);complete<-ids[vapply(ids,function(id)all(visits%in%m226$time[m226$subject==id]),logical(1))]
res<-list();diag<-list();resids<-list();fits<-list()
empty<-function(v,id,dataset,contrast,scope='all_available')data.frame(variant=v,module_id=id,dataset=dataset,scope=scope,contrast=contrast,estimate=NA_real_,lower.CL=NA_real_,upper.CL=NA_real_,df=NA_real_,p.value=NA_real_,n_participants=NA_integer_,n_observations=NA_integer_,status='coverage_failed')
for(v in c('published','exclude_original106'))for(id in unique(members$module_id)){
  dd<-sc[sc$dataset=='GSE151371' & sc$variant==v & sc$module_id==id & sc$group%in%c('SCI','TC'),]
  row<-empty(v,id,'GSE151371','SCI - TC');row$n_participants<-nrow(dd);row$n_observations<-nrow(dd)
  if(all(is.finite(dd$score))){
    tt<-t.test(dd$score[dd$group=='SCI'],dd$score[dd$group=='TC'],var.equal=FALSE)
    row$estimate<-unname(diff(rev(tt$estimate)));row$lower.CL<-tt$conf.int[1];row$upper.CL<-tt$conf.int[2];row$df<-unname(tt$parameter);row$p.value<-tt$p.value;row$status<-'estimated'
  }
  res[[paste(v,id,'discovery')]]<-row
  for(scope in c('all_available',if(v=='published')'complete_cases' else NULL)){
    dd<-sc[sc$dataset=='GSE226238' & sc$variant==v & sc$module_id==id & sc$group=='SCI',]
    if(scope=='complete_cases')dd<-dd[dd$subject%in%complete,]
    dd$time<-factor(dd$time,levels=visits);dd$subject<-factor(dd$subject)
    ky<-paste(v,id,scope)
    rows<-do.call(rbind,lapply(paste(visits[-1],'- acute'),function(con)empty(v,id,'GSE226238',con,scope)))
    rows$n_observations<-nrow(dd);rows$n_participants<-nlevels(dd$subject)
    if(all(is.finite(dd$score))){
      fit<-tryCatch(lmerTest::lmer(score~time+(1|subject),data=dd,REML=TRUE),error=function(e)e)
      if(inherits(fit,'error')){rows$status<-'model_failed';diag[[ky]]<-data.frame(variant=v,module_id=id,scope=scope,error=conditionMessage(fit))}else{
        fits[[ky]]<-fit
        tab<-as.data.frame(summary(contrast(emmeans(fit,~time,lmer.df='satterthwaite'),'trt.vs.ctrl',ref=1),infer=c(TRUE,TRUE),adjust='none'))
        rows[,c('estimate','lower.CL','upper.CL','df','p.value')]<-tab[,c('estimate','lower.CL','upper.CL','df','p.value')];rows$status<-if(isSingular(fit))'estimated_singular' else 'estimated'
        diag[[ky]]<-data.frame(variant=v,module_id=id,scope=scope,singular=isSingular(fit),fixed_design_rank=qr(model.matrix(fit))$rank,fixed_design_columns=ncol(model.matrix(fit)),residual_sd=sigma(fit),random_intercept_variance=as.data.frame(VarCorr(fit))$vcov[1],residual_shapiro_p=shapiro.test(residuals(fit))$p.value,convergence_message=paste(fit@optinfo$conv$lme4$messages,collapse=';'))
        resids[[ky]]<-data.frame(variant=v,module_id=id,scope=scope,sample=dd$sample,subject=dd$subject,time=dd$time,fitted=fitted(fit),residual=residuals(fit))
      }
    }
    res[[ky]]<-rows
  }
}
results<-as.data.frame(rbindlist(res,fill=TRUE));results$family<-ifelse(results$scope=='complete_cases','published_complete_cases_12',paste0(results$variant,'_16'))
results$q_BH<-NA_real_
for(fam in unique(results$family)){
  idx<-which(results$family==fam);nn<-if(grepl('complete',fam))12 else 16;stopifnot(length(idx)==nn)
  results$q_BH[idx]<-p.adjust(results$p.value[idx],method='BH',n=nn)
}
results$lineage<-members$lineage[match(results$module_id,members$module_id)]
results$direction_consistent<-with(results,ifelse(dataset=='GSE151371',ifelse(lineage=='neutrophil',estimate>0,estimate<0),ifelse(lineage=='neutrophil',estimate<0,estimate>0)))
write_table(results,'external_BTM_effects');write_table(rbindlist(diag,fill=TRUE),'external_BTM_model_diagnostics');write_table(rbindlist(resids),'external_BTM_model_residuals')
# Avoid serializing lmer formula environments (recursive state under R 4.6).
# Full coefficients, diagnostics, residuals and executable models are retained separately.
saveRDS(list(scores=sc,coverage=as.data.frame(covtab),effects=results,members=members),file.path(out,'tables/external_BTM_state.rds'))
capture.output(sessionInfo(),file=file.path(out,'logs/external_BTM_sessionInfo.txt'))
print(covtab[,c('dataset','variant','module_id','n_used','n_source_tokens','coverage','eligible')]);print(results[,c('family','module_id','contrast','estimate','lower.CL','upper.CL','p.value','q_BH','status')])
cat('Analysis completed UTC:',format(Sys.time(),tz='UTC',usetz=TRUE),'\n');sink()
