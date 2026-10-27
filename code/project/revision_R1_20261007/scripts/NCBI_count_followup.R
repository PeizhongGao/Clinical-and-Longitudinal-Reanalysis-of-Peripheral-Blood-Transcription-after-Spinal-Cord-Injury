# Public-package copy: see metadata/code_changes.csv for provenance.
# Revision-stage full-coverage GSE226238 analysis; protocol recorded before fitting.
# Distinct NCBI-generated counts, not a reproduction/replacement of discovery preprocessing.
args<-commandArgs(trailingOnly=TRUE)
root<-normalizePath(if(length(args))args[1] else '.',winslash='/')
out<-file.path(root,'revision_R1_20261007')
suppressPackageStartupMessages({library(data.table);library(edgeR);library(lmerTest);library(emmeans)})
sink(file.path(out,'logs/NCBI_count_followup_R.txt'),split=TRUE)
set.seed(20261007)
write_table<-function(x,name)fwrite(as.data.frame(x),file.path(out,'tables',paste0(name,'.csv')),na='NA')
read_gz<-function(path)read.delim(gzfile(path),check.names=FALSE,stringsAsFactors=FALSE,colClasses=c('character',rep('numeric',46)))
raw<-read_gz(file.path(root,'data/GSE226238/GSE226238_raw_counts_GRCh38.p13_NCBI.tsv.gz'))
ann<-read.delim(gzfile(file.path(root,'data/GSE226238/Human.GRCh38.p13.annot.tsv.gz')),check.names=FALSE,stringsAsFactors=FALSE,colClasses='character')
counts<-as.matrix(raw[,-1]);storage.mode(counts)<-'double';rownames(counts)<-raw$GeneID
stopifnot(nrow(counts)==39376,ncol(counts)==46,!anyDuplicated(rownames(counts)),all(is.finite(counts)),all(counts>=0),all(counts==round(counts)))
ai<-match(rownames(counts),ann$GeneID);stopifnot(!anyNA(ai));symbols<-ann$Symbol[ai]
meta<-as.data.frame(fread(file.path(out,'tables/sample_manifest.csv')))
m<-meta[meta$dataset=='GSE226238',];m<-m[match(colnames(counts),m$gsm),]
stopifnot(identical(colnames(counts),m$gsm),nrow(m)==46,sum(m$group=='CTL')==9,sum(m$group=='SCI')==37)
m$time<-ifelse(m$group=='SCI',paste0(sub('^.*[vV]','',m$sample),'m'),'Control');m$time[m$time=='0m']<-'acute'
write_table(m,'NCBI_GSE226238_sample_map')
dge<-DGEList(counts=counts);dge<-calcNormFactors(dge,method='TMM')
Eid<-cpm(dge,log=TRUE,prior.count=.5)
write_table(data.frame(gsm=m$gsm,sample=m$sample,group=m$group,time=m$time,dge$samples,all_zero_genes=sum(rowSums(counts)==0),prior_count=.5),'NCBI_GSE226238_normalization')
locked<-as.data.frame(fread(file.path(out,'tables/locked_signature_reused_S14.csv')));genes<-split(locked$gene,locked$module)
idx<-match(locked$gene,symbols);stopifnot(!anyNA(idx),all(table(symbols)[locked$gene]==1))
E<-Eid[idx,,drop=FALSE];rownames(E)<-locked$gene;colnames(E)<-m$sample
write_table(data.frame(gene=rownames(E),GeneID=rownames(Eid)[idx],E,check.names=FALSE),'NCBI_GSE226238_signature_logCPM')
coverage<-list();parameters<-list()
score_module<-function(gs,source){
  s<-m
  for(mod in names(gs)){
    X<-E[gs[[mod]],,drop=FALSE];ref<-which(m$group=='CTL');mu<-rowMeans(X[,ref,drop=FALSE]);sdr<-apply(X[,ref,drop=FALSE],1,sd)
    stopifnot(all(is.finite(sdr)),all(sdr>1e-6))
    s[[paste0(mod,'_score')]]<-colMeans(sweep(sweep(X,1,mu,'-'),1,sdr,'/'))
    coverage[[paste(source,mod)]]<<-data.frame(source=source,module=mod,n_requested=length(genes[[mod]]),n_present=length(gs[[mod]]),n_used=length(gs[[mod]]),coverage=length(gs[[mod]])/length(genes[[mod]]),reference='9 CTL',missing=paste(setdiff(genes[[mod]],gs[[mod]]),collapse=';'))
    parameters[[paste(source,mod)]]<<-data.frame(source=source,module=mod,gene=gs[[mod]],reference='CTL',ref_mean=mu,ref_sd=sdr)
  }
  s$balance_score<-s$innate_score-s$adaptive_score;s$source<-source;s
}
sfull<-score_module(genes,'NCBI_full_15plus15')
oldmap<-as.data.frame(fread(file.path(out,'tables/official_GSE226238_coverage.csv')))
# The source coverage CSV is gene-level; present genes are derived from original standardization audit.
oldpar<-as.data.frame(fread(file.path(out,'tables/standardization_parameters.csv')))
partial<-lapply(names(genes),function(mod)intersect(genes[[mod]],oldpar$gene[oldpar$dataset=='GSE226238' & oldpar$module==mod]));names(partial)<-names(genes)
stopifnot(length(partial$innate)==8,length(partial$adaptive)==13)
spart<-score_module(partial,'NCBI_restricted_8plus13')
write_table(sfull,'NCBI_longitudinal_full_scores');write_table(spart,'NCBI_longitudinal_restricted_scores');write_table(do.call(rbind,coverage),'NCBI_gene_coverage');write_table(do.call(rbind,parameters),'NCBI_standardization_parameters')
visits<-c('acute','3m','6m','12m');scores<-c('innate_score','adaptive_score','balance_score')
subjects<-unique(m$subject[m$group=='SCI']);complete<-subjects[vapply(subjects,function(id)all(visits%in%m$time[m$subject==id]),logical(1))]
stopifnot(length(subjects)==10,length(complete)==7)
results<-list();diags<-list();models<-list()
fit_scores<-function(s,scope='all_available',markers=NULL){
  df<-s[s$group=='SCI',];if(scope=='complete_cases')df<-df[df$subject%in%complete,]
  df$time<-factor(df$time,levels=visits);df$subject<-factor(df$subject)
  if(!is.null(markers)){
    ii<-match(df$sample,rownames(markers));df$neutrophil_marker<-as.numeric(scale(markers[ii,'Neutrophils']));df$tcell_marker<-as.numeric(scale(markers[ii,'T cells']))
  }
  for(sc in scores){
    f<-as.formula(paste(sc,'~ time',if(is.null(markers))'' else '+ neutrophil_marker + tcell_marker','+ (1|subject)'))
    fit<-lmerTest::lmer(f,data=df,REML=TRUE)
    tab<-as.data.frame(summary(contrast(emmeans(fit,~time,lmer.df='satterthwaite'),'trt.vs.ctrl',ref=1),infer=c(TRUE,TRUE),adjust='none'))
    spec<-if(is.null(markers))'unadjusted' else 'marker_adjusted';key<-paste(unique(s$source),scope,spec,sc)
    tab$score<-sc;tab$source<-unique(s$source);tab$scope<-scope;tab$spec<-spec;tab$n_samples<-nrow(df);tab$n_participants<-length(unique(df$subject));tab$model<-paste(deparse(f),collapse=' ')
    results[[key]]<<-tab
    diags[[key]]<<-data.frame(source=unique(s$source),scope=scope,spec=spec,score=sc,singular=isSingular(fit),fixed_design_rank=qr(model.matrix(fit))$rank,fixed_design_columns=ncol(model.matrix(fit)),kappa=kappa(model.matrix(fit)),random_intercept_variance=as.data.frame(VarCorr(fit))$vcov[1],residual_sd=sigma(fit),residual_shapiro_p=shapiro.test(residuals(fit))$p.value,convergence_message=paste(fit@optinfo$conv$lme4$messages,collapse=';'))
    models[[key]]<<-fit
    capture.output(summary(fit),file=file.path(out,'logs',paste0(gsub(' ','_',key),'.txt')))
  }
}
fit_scores(sfull);fit_scores(sfull,'complete_cases');fit_scores(spart)
marker<-as.data.frame(fread(file.path(out,'inputs/MCPcounter_genes.txt')));names(marker)[1:2]<-c('gene','cell')
sets<-split(marker$gene,marker$cell);sets<-lapply(sets,function(g)setdiff(unique(g),locked$gene))
ms<-list();markercoverage<-list()
for(cell in c('Neutrophils','T cells')){
  g<-intersect(sets[[cell]],symbols);stopifnot(length(g)>=3)
  mc<-rowsum(counts[symbols%in%g,,drop=FALSE],symbols[symbols%in%g],reorder=FALSE)
  me<-cpm(mc,lib.size=dge$samples$lib.size*dge$samples$norm.factors,log=TRUE,prior.count=.5)
  ms[[cell]]<-colMeans(me);markercoverage[[cell]]<-data.frame(cell=cell,genes=paste(g,collapse=';'),n_used=length(g),all_signature_genes_excluded=TRUE)
}
ms<-as.data.frame(ms,check.names=FALSE);rownames(ms)<-m$sample;fit_scores(sfull,markers=ms)
write_table(data.frame(sample=m$sample,ms,check.names=FALSE),'NCBI_MCP_marker_abundances');write_table(do.call(rbind,markercoverage),'NCBI_MCP_marker_coverage')
mixed<-do.call(rbind,results);mixed$q_BH<-ave(mixed$p.value,interaction(mixed$source,mixed$scope,mixed$spec),FUN=function(p)p.adjust(p,'BH'))
write_table(mixed,'NCBI_longitudinal_mixed_contrasts');write_table(do.call(rbind,diags),'NCBI_mixed_model_diagnostics')
paired<-list()
for(sc in scores)for(t in visits[-1]){
  x<-sfull[sfull$time=='acute',c('subject',sc)];y<-sfull[sfull$time==t,c('subject',sc)];z<-merge(x,y,by='subject');delta<-z[[paste0(sc,'.y')]]-z[[paste0(sc,'.x')]]
  tt<-t.test(delta);w<-wilcox.test(delta,exact=FALSE)
  paired[[paste(sc,t)]]<-data.frame(score=sc,time=t,n_pairs=nrow(z),mean_change=mean(delta),ci_low=tt$conf.int[1],ci_high=tt$conf.int[2],median_change=median(delta),wilcox_p=w$p.value)
}
paired<-do.call(rbind,paired);paired$q_BH<-p.adjust(paired$wilcox_p,'BH');write_table(paired,'NCBI_longitudinal_paired_sensitivity')
old<-as.data.frame(fread(file.path(out,'tables/longitudinal_locked_scores.csv')));oi<-match(sfull$gsm,old$gsm);stopifnot(!anyNA(oi));old<-old[oi,]
comparisons<-do.call(rbind,lapply(scores,function(sc)data.frame(score=sc,n=46,NCBI_full_vs_author_partial_spearman=cor(sfull[[sc]],old[[sc]],method='spearman'),NCBI_restricted_vs_author_partial_spearman=cor(spart[[sc]],old[[sc]],method='spearman'),NCBI_full_vs_NCBI_restricted_spearman=cor(sfull[[sc]],spart[[sc]],method='spearman'))))
write_table(comparisons,'NCBI_author_score_comparison')
# Conditioned-on-transform attrition summaries; these cannot rule out informative dropout.
acute<-sfull[sfull$time=='acute',];acute$retained12<-acute$subject%in%complete
attrition<-do.call(rbind,lapply(scores,function(sc)data.frame(score=sc,n_retained=sum(acute$retained12),n_missing12=sum(!acute$retained12),mean_retained=mean(acute[[sc]][acute$retained12]),mean_missing12=mean(acute[[sc]][!acute$retained12]),wilcox_p=wilcox.test(acute[[sc]]~acute$retained12,exact=FALSE)$p.value)))
write_table(attrition,'NCBI_baseline_attrition_descriptive')
saveRDS(list(E=E,Eid=Eid,sfull=sfull,spart=spart,mixed=mixed,models=models,diags=diags,ms=ms,genes=genes,normalization=dge$samples),file.path(out,'tables/NCBI_analysis_state.rds'))
capture.output(sessionInfo(),file=file.path(out,'logs/NCBI_sessionInfo.txt'))
print(comparisons);print(mixed[mixed$source=='NCBI_full_15plus15' & mixed$scope=='all_available' & mixed$spec=='unadjusted',c('contrast','score','estimate','lower.CL','upper.CL','p.value','q_BH')]);print(do.call(rbind,diags))
sink()
