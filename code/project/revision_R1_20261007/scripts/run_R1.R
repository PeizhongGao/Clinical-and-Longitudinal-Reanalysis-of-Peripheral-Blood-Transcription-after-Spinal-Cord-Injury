# Public-package copy: see metadata/code_changes.csv for provenance.
# Reviewer 1 revision. Uses existing processed supplementary data/matrices.
# No read realignment, raw-count filtering, voom or replacement of original preprocessing.
# Exploratory reanalysis decisions recorded before executing this script, 2026-10-07.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1],winslash='/') else normalizePath('..',winslash='/')
out <- file.path(root,'revision_R1_20261007')
for(d in c('tables','figures','logs')) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
sink(file.path(out,'logs/R_run.txt'),split=TRUE)
on.exit <- NULL
suppressPackageStartupMessages({library(data.table);library(limma);library(readxl);library(jsonlite);library(ggplot2);library(patchwork);library(lmerTest);library(emmeans);library(AnnotationDbi);library(org.Hs.eg.db)})
set.seed(20261007)
theme_set(theme_classic(base_size=9,base_family='Arial')+theme(strip.text=element_text(face='bold'),legend.position='bottom',plot.title=element_text(size=10,face='bold')))
palette <- c(HC='#4F6474',TC='#518AA0',SCI='#B9564F',CTL='#4F6474',innate='#B9564F',adaptive='#518AA0',balance='#756783')
write_table <- function(x,name) fwrite(as.data.frame(x),file.path(out,'tables',paste0(name,'.csv')),na='NA')
save_plot <- function(p,name,w=183,h=120){
  stem<-file.path(out,'figures',name)
  ggsave(paste0(stem,'.pdf'),p,width=w/25.4,height=h/25.4,device=cairo_pdf,family='Arial')
  ggsave(paste0(stem,'.svg'),p,width=w/25.4,height=h/25.4,device=svglite::svglite)
  ggsave(paste0(stem,'.png'),p,width=w/25.4,height=h/25.4,dpi=180,device=ragg::agg_png)
  ggsave(paste0(stem,'.tiff'),p,width=w/25.4,height=h/25.4,dpi=600,device=ragg::agg_tiff,compression='lzw')
}
source(file.path(out,'scripts/candidate_pools.R'))
cat('SOURCE POLICY: original Supplementary_Data1-4 and processed expression matrices; no replacement of preprocessing\n')
supp <- list()
suppfiles<-file.path(root,'data/Supplementary_Data',paste0('Supplementary_Data',1:4,'.xlsx'))
for(f in suppfiles) for(s in excel_sheets(f)) {
  tab<-as.data.frame(read_excel(f,sheet=s,.name_repair='minimal'));supp[[s]]<-tab
  write_table(tab,paste0('original_',gsub(' ','_',s)))
}
locked<-supp[['Table S14']]
genes<-split(locked$gene,locked$module)
stopifnot(length(genes$innate)==15,length(genes$adaptive)==15)
allgenes<-unlist(genes,use.names=FALSE)
write_table(locked,'locked_signature_reused_S14')
meta<-as.data.frame(fread(file.path(out,'tables/sample_manifest.csv'),na.strings=c('','NA')))
read_gene_matrix<-function(file){
  d<-fread(file); m<-as.matrix(d[,-1,with=FALSE]);storage.mode(m)<-'double';rownames(m)<-d[[1]]
  # Reproduce original scripts: first duplicate gene retained, no additional filtering.
  m[!duplicated(rownames(m)),,drop=FALSE]
}
E<-read_gene_matrix(file.path(root,'explore/Report/GSE151371_normalized_expression_matrix_comBat_de-ID.csv'))
supp[['S1_DEGs']]$Gene[supp[['S1_DEGs']]$Gene=='SEPT15']<-'SEP15'
write_table(data.frame(supplement_gene='SEPT15',processed_matrix_gene='SEP15',reason='S1 legacy date-conversion repair label; expression values unchanged'),'S1_gene_alias_map')
# Intact gene symbols in this processed copy. Preserve the published S1 tested-gene universe.
stopifnot(all(supp[['S1_DEGs']]$Gene %in% rownames(E)))
E<-E[supp[['S1_DEGs']]$Gene,,drop=FALSE]
m151<-meta[meta$dataset=='GSE151371',];m151<-m151[match(colnames(E),m151$sample),]
stopifnot(!anyNA(m151$sample),identical(colnames(E),m151$sample),all(is.finite(E)))
m151$group<-factor(m151$group,levels=c('HC','TC','SCI'))
fit_pair<-function(mat,metadata,case,ref){
  use<-metadata$group %in% c(ref,case);g<-factor(as.character(metadata$group[use]),levels=c(ref,case));X<-model.matrix(~g)
  f<-eBayes(lmFit(mat[,use,drop=FALSE],X));tt<-topTable(f,coef=2,number=Inf,sort.by='none',confint=TRUE);tt$gene<-rownames(tt)
  list(tt=tt,fit=f,X=X,mat=mat[,use,drop=FALSE])
}
hcfit<-fit_pair(E,m151,'SCI','HC');tcfit<-fit_pair(E,m151,'SCI','TC');thfit<-fit_pair(E,m151,'TC','HC')
orig<-supp[['S1_DEGs']];ix<-match(orig$Gene,hcfit$tt$gene);stopifnot(!anyNA(ix))
audit<-data.frame(check=c('HC_vs_SCI_max_abs_logFC_delta','HC_vs_SCI_max_abs_t_delta','HC_vs_SCI_max_abs_FDR_delta','original_genes','processed_matrix_genes'),value=c(max(abs(orig$logFC-hcfit$tt$logFC[ix])),max(abs(orig$t-hcfit$tt$t[ix])),max(abs(orig$adj.P.Val-hcfit$tt$adj.P.Val[ix])),nrow(orig),nrow(E)))
write_table(audit,'S1_reproduction_audit');print(audit)
write_table(tcfit$tt,'SCI_vs_TC_processed_limma');write_table(hcfit$tt,'SCI_vs_HC_processed_limma');write_table(thfit$tt,'TC_vs_HC_processed_limma')
counts<-do.call(rbind,lapply(list(SCI_vs_TC=tcfit,SCI_vs_HC=hcfit,TC_vs_HC=thfit),function(f)data.frame(up=sum(f$tt$adj.P.Val<.05 & f$tt$logFC>1),down=sum(f$tt$adj.P.Val<.05 & f$tt$logFC< -1))))
counts$contrast<-rownames(counts);write_table(counts,'DEG_counts_processed')
coverage<-list();standardization<-list()
score_module<-function(mat,metadata,dataset,ref,gs=genes,write_audit=TRUE,allow_low_coverage=FALSE){
  metadata<-metadata[match(colnames(mat),metadata$sample),];stopifnot(!anyNA(metadata$sample))
  result<-metadata
  for(mod in names(gs)){
    present<-intersect(gs[[mod]],rownames(mat));X<-mat[present,,drop=FALSE];refs<-which(as.character(metadata$group)==ref)
    if(length(refs)<2)stop('Insufficient reference samples')
    mu<-rowMeans(X[,refs,drop=FALSE]);sdref<-apply(X[,refs,drop=FALSE],1,sd);valid<-is.finite(sdref)&sdref>1e-6
    X<-X[valid,,drop=FALSE];mu<-mu[valid];sdref<-sdref[valid]
    used<-rownames(X)
    if(write_audit){
      coverage[[paste(dataset,mod)]]<<-data.frame(dataset=dataset,module=mod,n_requested=length(gs[[mod]]),n_present=length(present),n_used=length(used),coverage=length(used)/length(gs[[mod]]),missing=paste(setdiff(gs[[mod]],used),collapse=';'),reference=ref,sd_min=1e-6,minimum_coverage=.8)
      standardization[[paste(dataset,mod)]]<<-data.frame(dataset=dataset,module=mod,gene=used,reference=ref,ref_mean=mu,ref_sd=sdref)
    }
    if(length(used)/length(gs[[mod]])<.8 && !allow_low_coverage)stop('Signature coverage below .8: ',dataset,' ',mod)
    if(length(used)/length(gs[[mod]])<.8)warning('LOW COVERAGE PROJECTION, not full signature validation: ',dataset,' ',mod)
    result[[paste0(mod,'_score')]]<-colMeans(sweep(sweep(X,1,mu,'-'),1,sdref,'/'))
  }
  result$balance_score<-result$innate_score-result$adaptive_score
  result
}
s151<-score_module(E,m151,'GSE151371','TC');write_table(s151,'discovery_locked_scores_descriptive')
raw226<-fread(file.path(root,'explore/GSE226238.csv'));E226<-t(as.matrix(raw226[,-1,with=FALSE]));storage.mode(E226)<-'double';colnames(E226)<-raw226[[1]]
transcripts<-rownames(E226);symbols<-sub('^.*>','',transcripts);n_per<-table(symbols)
E226<-rowsum(E226,symbols,reorder=FALSE)/as.numeric(n_per[unique(symbols)])
write_table(data.frame(transcript=transcripts,gene=symbols,collapse='arithmetic mean on author log2/TMM matrix'),'GSE226238_transcript_gene_map')
m226<-meta[meta$dataset=='GSE226238',];m226<-m226[match(tolower(colnames(E226)),tolower(m226$sample)),];m226$sample<-colnames(E226)
stopifnot(!anyNA(m226$gsm),all(is.finite(E226)))
m226$time<-ifelse(m226$group=='SCI',paste0(sub('^.*[vV]','',m226$sample),'m'),'Control');m226$time[m226$time=='0m']<-'acute'
s226<-score_module(E226,m226,'GSE226238','CTL',allow_low_coverage=TRUE)
s226$projection_scope<-'reduced coverage: 8/15 innate and 13/15 adaptive; availability only, no outcome-driven reselection'
write_table(s226,'longitudinal_locked_scores')
write_table(do.call(rbind,coverage),'gene_coverage');write_table(do.call(rbind,standardization),'standardization_parameters')
visit<-do.call(rbind,lapply(c('acute','3m','6m','12m'),function(t){d<-s226[s226$time==t,];data.frame(time=t,n_samples=nrow(d),n_participants=length(unique(d$subject)))}));write_table(visit,'visit_counts')
subjects<-unique(s226$subject[s226$group=='SCI']);visits<-c('acute','3m','6m','12m')
grid<-expand.grid(subject=subjects,time=visits);grid$observed<-paste(grid$subject,grid$time) %in% paste(s226$subject,s226$time);grid$missing_reason<-ifelse(grid$observed,'not_missing','not_available_in_processed_matrix; reason not provided')
write_table(grid,'longitudinal_missingness')
complete<-subjects[vapply(subjects,function(s)all(visits %in% s226$time[s226$subject==s]),logical(1))]
scores<-c('innate_score','adaptive_score','balance_score')
mixed_results<-list();model_diagnostics<-list();models<-list()
for(scope in c('all_available','complete_cases'))for(sc in scores){
  d<-s226[s226$group=='SCI',];if(scope=='complete_cases')d<-d[d$subject %in% complete,]
  d$time<-factor(d$time,levels=visits);d$subject<-factor(d$subject)
  f<-as.formula(paste(sc,'~ time + (1|subject)'))
  fit<-lmerTest::lmer(f,data=d,REML=TRUE)
  em<-emmeans(fit,~time,lmer.df='satterthwaite');cc<-contrast(em,'trt.vs.ctrl',ref=1)
  tab<-as.data.frame(summary(cc,infer=c(TRUE,TRUE),adjust='none'))
  tab$score<-sc;tab$scope<-scope;tab$n_samples<-nrow(d);tab$n_participants<-length(unique(d$subject));tab$model<-deparse(f);tab$df_method<-'Satterthwaite';mixed_results[[paste(scope,sc)]]<-tab
  rr<-residuals(fit);model_diagnostics[[paste(scope,sc)]]<-data.frame(scope=scope,score=sc,singular=lme4::isSingular(fit),random_intercept_variance=as.data.frame(VarCorr(fit))$vcov[1],residual_sd=sigma(fit),residual_shapiro_p=shapiro.test(rr)$p.value,n_samples=nrow(d),n_participants=length(unique(d$subject)))
  capture.output(summary(fit),file=file.path(out,'logs',paste0('mixed_',scope,'_',sc,'.txt')))
  models[[paste(scope,sc)]]<-fit
}
mixed<-do.call(rbind,mixed_results);mixed$q_BH<-ave(mixed$p.value,mixed$scope,FUN=function(p)p.adjust(p,'BH'));write_table(mixed,'longitudinal_mixed_contrasts');write_table(do.call(rbind,model_diagnostics),'mixed_model_diagnostics')
paired<-list()
for(sc in scores)for(t in visits[-1]){
  acute<-s226[s226$time=='acute',c('subject',sc)];later<-s226[s226$time==t,c('subject',sc)];d<-merge(acute,later,by='subject');delta<-d[[paste0(sc,'.y')]]-d[[paste0(sc,'.x')]]
  tt<-t.test(delta);w<-wilcox.test(delta,exact=FALSE)
  paired[[paste(sc,t)]]<-data.frame(score=sc,time=t,n_pairs=nrow(d),mean_change=mean(delta),ci_low=tt$conf.int[1],ci_high=tt$conf.int[2],median_change=median(delta),wilcox_p=w$p.value)
}
paired<-do.call(rbind,paired);paired$q_BH<-p.adjust(paired$wilcox_p,'BH');write_table(paired,'longitudinal_paired_sensitivity')
acute<-s226[s226$time=='acute',];acute$retained12<-acute$subject %in% complete
attrition<-do.call(rbind,lapply(scores,function(sc)data.frame(score=sc,n_retained=sum(acute$retained12),n_missing12=sum(!acute$retained12),mean_retained=mean(acute[[sc]][acute$retained12]),mean_missing12=mean(acute[[sc]][!acute$retained12]),wilcox_p=wilcox.test(acute[[sc]]~acute$retained12,exact=FALSE)$p.value)))
write_table(attrition,'baseline_attrition_descriptive')
# Use existing S16 scores and exact reuse provenance, not the faulty legacy 0-case CSV.
s293<-supp[['Table S16']];s293$overlap<-ifelse(s293$reused_from_gse226238=='yes','reused','new')
gfun<-function(x,y){n<-length(x);m<-length(y);sp<-sqrt(((n-1)*var(x)+(m-1)*var(y))/(n+m-2));(mean(x)-mean(y))/sp*(1-3/(4*(n+m)-9))}
cross<-list()
for(scope in c('all_available','exclude_reused_SCI'))for(sc in scores){
  d<-s293;if(scope=='exclude_reused_SCI')d<-d[d$group=='CTL'|d$overlap=='new',]
  x<-d[[sc]][d$group=='SCI'];y<-d[[sc]][d$group=='CTL'];b<-replicate(2000,gfun(sample(x,replace=TRUE),sample(y,replace=TRUE)));ci<-quantile(b,c(.025,.975),na.rm=TRUE)
  cross[[paste(scope,sc)]]<-data.frame(scope=scope,score=sc,n_SCI=length(x),n_CTL=length(y),mean_SCI=mean(x),mean_CTL=mean(y),hedges_g=gfun(x,y),ci_low=ci[1],ci_high=ci[2],p=wilcox.test(x,y,exact=FALSE)$p.value,reference='9 CTL reused from GSE226238',uncertainty='stratified participant bootstrap; score transformation held fixed')
}
cross<-do.call(rbind,cross);cross$q_BH<-ave(cross$p,cross$scope,FUN=function(p)p.adjust(p,'BH'));write_table(cross,'GSE293559_reconciled_statistics')
# Covariate availability: public sample-level fields only. No imputation of absent clinical data.
covars<-c('age','sex','ais','nli','injury_time_hours','infection','treatment','batch','rin','platform','instrument')
av<-do.call(rbind,lapply(unique(meta$dataset),function(ds)do.call(rbind,lapply(covars,function(c){d<-meta[meta$dataset==ds,];x<-d[[c]];data.frame(dataset=ds,covariate=c,n_samples=nrow(d),n_observed=sum(!is.na(x)&x!=''),n_unique=length(unique(na.omit(x))),limitation=if(sum(!is.na(x)&x!='')==0)'unavailable; no adjustment possible' else 'see actual values and group/time confounding')}))))
write_table(av,'covariate_availability')
# External marker abundances, with ALL 30 signature genes removed to limit shared-gene coupling.
marker<-fread(file.path(out,'inputs/MCPcounter_genes.txt'));names(marker)[1:2]<-c('gene','cell')
marksets<-split(marker$gene,marker$cell);marksets<-lapply(marksets,function(g)setdiff(unique(g),allgenes))
marker_scores<-function(mat){ans<-list();for(cell in names(marksets)){g<-intersect(marksets[[cell]],rownames(mat));if(length(g)>=3)ans[[cell]]<-colMeans(mat[g,,drop=FALSE])};as.data.frame(ans,check.names=FALSE)}
ms151<-marker_scores(E);ms226<-marker_scores(E226)
write_table(data.frame(sample=colnames(E),ms151,check.names=FALSE),'MCP_marker_abundances_discovery');write_table(data.frame(sample=colnames(E226),ms226,check.names=FALSE),'MCP_marker_abundances_longitudinal')
write_table(do.call(rbind,lapply(names(marksets),function(cell)data.frame(cell=cell,genes_after_signature_exclusion=paste(marksets[[cell]],collapse=';'),n_in_discovery=sum(marksets[[cell]]%in%rownames(E)),n_in_longitudinal=sum(marksets[[cell]]%in%rownames(E226))))),'MCP_marker_coverage')
adj<-list()
for(sc in scores){
  d<-s151[s151$group %in% c('SCI','TC'),];idx<-match(d$sample,colnames(E));d$case<-as.integer(d$group=='SCI')
  stopifnot(all(c('Neutrophils','T cells')%in%colnames(ms151)))
  d$neutrophil_marker<-as.numeric(scale(ms151[idx,'Neutrophils']));d$tcell_marker<-as.numeric(scale(ms151[idx,'T cells']))
  for(spec in c('unadjusted','marker_adjusted')){
    f<-as.formula(paste(sc,'~ case',if(spec=='marker_adjusted')'+ neutrophil_marker + tcell_marker' else ''));fit<-lm(f,data=d);co<-summary(fit)$coefficients;ci<-confint(fit,'case')
    adj[[paste(sc,spec)]]<-data.frame(score=sc,spec=spec,estimate=co['case','Estimate'],se=co['case','Std. Error'],ci_low=ci[1],ci_high=ci[2],p=co['case','Pr(>|t|)'],n=nrow(d),design=deparse(f),rank=fit$rank,kappa=kappa(model.matrix(fit)),inference='discovery-derived score; exploratory association, not independent validation')
  }
}
adj<-do.call(rbind,adj);adj$q_BH<-ave(adj$p,adj$spec,FUN=function(p)p.adjust(p,'BH'));write_table(adj,'composition_marker_sensitivity')
# Longitudinal marker sensitivity, one parsimonious 2-marker extension, no random slopes with n=10.
ladj<-list()
for(sc in scores){d<-s226[s226$group=='SCI',];idx<-match(d$sample,colnames(E226));d$time<-factor(d$time,levels=visits);d$neutrophil_marker<-as.numeric(scale(ms226[idx,'Neutrophils']));d$tcell_marker<-as.numeric(scale(ms226[idx,'T cells']));f<-as.formula(paste(sc,'~ time + neutrophil_marker + tcell_marker + (1|subject)'));fit<-lmer(f,data=d,REML=TRUE);tab<-as.data.frame(summary(contrast(emmeans(fit,~time,lmer.df='satterthwaite'),'trt.vs.ctrl',ref=1),infer=c(TRUE,TRUE),adjust='none'));tab$score<-sc;tab$singular<-lme4::isSingular(fit);ladj[[sc]]<-tab}
ladj<-do.call(rbind,ladj);ladj$q_BH<-p.adjust(ladj$p.value,'BH');write_table(ladj,'longitudinal_marker_adjusted_sensitivity')
# Reuse original ALL-candidate frequencies for sensitivity to retention cutoffs.
stab<-fread(file.path(root,'code/01_discovery_signature/GSE151371_stability_selection_all_candidates.csv'))
thresh<-list()
for(freq in c(.3,.4,.5))for(dir in c(.7,.8,.9))for(mod in names(genes)){
  d<-stab[stab$module==mod,];sign<-if(mod=='innate')1 else -1
  g<-d$gene[d$selection_frequency>=freq & d$direction_consistency>=dir & sign*d$full_logFC_SCI_vs_TC>0 & sign*d$full_logFC_SCI_vs_HC>0]
  thresh[[paste(freq,dir,mod)]]<-data.frame(module=mod,min_frequency=freq,min_direction=dir,n_retained_before_size_cap=length(g),genes=paste(g,collapse=';'),jaccard_locked=length(intersect(g,genes[[mod]]))/length(union(g,genes[[mod]])),source='original all-candidate stability table; no raw data reprocessing')
}
write_table(do.call(rbind,thresh),'original_frequency_threshold_sensitivity')
# New top-N sensitivity is an exploratory perturbation of the processed matrix.
# It does not retroactively redefine the frozen S14 signature or reproduce raw-count voom selection.
B<-500;Ns<-c(10,15,20);pools<-list(innate=intersect(innate_candidate_pool,rownames(E)),adaptive=intersect(adaptive_candidate_pool,rownames(E)))
record<-expand.grid(gene=character(),module=character())
record<-do.call(rbind,lapply(names(pools),function(mod)data.frame(gene=pools[[mod]],module=mod)))
select_count<-array(0,dim=c(nrow(record),length(Ns),3),dimnames=list(NULL,as.character(Ns),c('B100','B250','B500')));correct<-rep(0,nrow(record))
for(b in seq_len(B)){
  ids<-unlist(lapply(levels(m151$group),function(g)sample(which(m151$group==g),floor(sum(m151$group==g)*.8))))
  gd<-factor(m151$group[ids],levels=c('HC','TC','SCI'));X<-model.matrix(~0+gd);colnames(X)<-levels(gd)
  f<-eBayes(contrasts.fit(lmFit(E[,ids,drop=FALSE],X),makeContrasts(SCI-TC,SCI-HC,levels=X)))
  for(mod in names(pools)){
    ii<-which(record$module==mod);jj<-match(record$gene[ii],rownames(E));sign<-if(mod=='innate')1 else -1
    ok<-sign*f$coefficients[jj,1]>0 & sign*f$coefficients[jj,2]>0;correct[ii]<-correct[ii]+ok
    ranked<-ii[order(f$p.value[jj,1])];ranked<-ranked[ok[match(ranked,ii)]]
    for(n in Ns){picked<-head(ranked,n);for(k in which(b<=c(100,250,500)))select_count[picked,as.character(n),k]<-select_count[picked,as.character(n),k]+1}
  }
  if(b%%50==0)cat('Processed-matrix resampling',b,'/',B,'\n')
}
stab_new<-do.call(rbind,lapply(seq_along(Ns),function(ni)do.call(rbind,lapply(1:3,function(k)data.frame(record,top_n=Ns[ni],iterations=c(100,250,500)[k],selection_frequency=select_count[,ni,k]/c(100,250,500)[k],direction_consistency_500=correct/B)))))
write_table(stab_new,'processed_matrix_topN_B_sensitivity')
# Biological labels: gene names support annotation; they are not cell-specific proofs.
ann<-AnnotationDbi::select(org.Hs.eg.db,keys=allgenes,keytype='SYMBOL',columns=c('GENENAME','ENTREZID'))
ann$module<-locked$module[match(ann$SYMBOL,locked$gene)];ann$interpretation<-ifelse(ann$module=='innate','myeloid/innate-associated; not cell-exclusive','lymphocyte/T-cell-associated; not cell-exclusive')
write_table(ann,'locked_gene_biological_annotation')
# Threshold-free competitive gene-set test, a valid alternative to cutoff-dependent ORA.
cat('Building GO BP membership; camera on existing processed matrix\n')
gokeys<-intersect(rownames(E),keys(org.Hs.eg.db,keytype='SYMBOL'))
gomap<-AnnotationDbi::select(org.Hs.eg.db,keys=gokeys,keytype='SYMBOL',columns=c('GOALL','ONTOLOGYALL'))
gomap<-unique(gomap[!is.na(gomap$GOALL)&gomap$ONTOLOGYALL=='BP',c('SYMBOL','GOALL')]);bp<-as.data.frame(fread(file.path(out,'inputs/go_bp_terms.csv')))
sets<-split(gomap$SYMBOL,gomap$GOALL);sets<-lapply(sets,unique);sets<-sets[lengths(sets)>=15 & lengths(sets)<=500];sets<-sets[names(sets)%in%bp$go_id]
idxsets<-lapply(sets,function(g)match(g,rownames(E)));cam<-camera(tcfit$mat,index=idxsets,design=tcfit$X,contrast=2,inter.gene.cor=.01,sort=FALSE)
cam$go_id<-rownames(cam);cam$term<-bp$term[match(cam$go_id,bp$go_id)];write_table(cam,'SCI_vs_TC_GO_camera')
write_table(data.frame(go_id=names(sets),n_genes=lengths(sets),genes=vapply(sets,paste,collapse=';',FUN.VALUE=character(1))),'GO_camera_membership')
write_table(data.frame(gene=rownames(E)),'GO_tested_gene_universe')
write_table(as.data.frame(AnnotationDbi::metadata(org.Hs.eg.db)),'annotation_package_metadata')
# Figure contract: quantitative grid; no mechanistic claim; source tables above; all participants shown.
pca<-prcomp(t(E[apply(E,1,sd)>0,]),scale.=TRUE);pdf<-data.frame(pca$x[,1:2],group=m151$group,sample=m151$sample)
write_table(pdf,'PCA_source_data')
p1a<-ggplot(pdf,aes(PC1,PC2,color=group))+geom_point(size=2)+scale_color_manual(values=palette)+labs(title='Processed transcriptomes',subtitle='HC n=10; TC n=10; SCI n=38; descriptive PCA')
vd<-tcfit$tt;vd$direction<-ifelse(vd$adj.P.Val<.05 & abs(vd$logFC)>1,ifelse(vd$logFC>0,'Up','Down'),'Other')
p1b<-ggplot(vd,aes(logFC,-log10(pmax(adj.P.Val,1e-300)),color=direction))+geom_point(size=.6,alpha=.6)+scale_color_manual(values=c(Up=unname(palette['SCI']),Down=unname(palette['TC']),Other='#C3C6C9'))+labs(x='SCI minus TC log2 fold change',y='-log10 BH FDR',title='Direct SCI versus trauma-control comparison')
save_plot(p1a+p1b+plot_annotation(tag_levels='a'),'Figure1_revised',183,90)
focusids<-c('GO:0030098','GO:0030217','GO:0050851','GO:0002376','GO:0006954','GO:0042742','GO:0002446','GO:0032496','GO:0042110','GO:0045087','GO:0002523','GO:0030595')
focus<-cam[cam$go_id %in% focusids,];focus$label<-vapply(focus$term,function(x)paste(strwrap(x,45),collapse='\n'),character(1));focus$label<-factor(focus$label,levels=focus$label[order(focus$FDR,decreasing=TRUE)])
write_table(focus,'GO_camera_figure_panel')
p2<-ggplot(focus,aes(-log10(FDR),label,color=Direction,size=NGenes))+geom_point()+scale_color_manual(values=c(Up=unname(palette['SCI']),Down=unname(palette['TC'])))+labs(x='-log10 BH FDR across all tested BP sets',y=NULL,title='SCI versus trauma controls: threshold-free GO BP test',subtitle='Selected immune terms shown regardless of significance; CAMERA correlation = 0.01')
save_plot(p2,'Figure2_revised',183,135)
heat<-E[allgenes,,drop=FALSE];zheat<-t(scale(t(heat)));hm<-as.data.frame(as.table(zheat));names(hm)<-c('gene','sample','z');hm$gene<-factor(hm$gene,levels=rev(allgenes));hm$sample<-factor(hm$sample,levels=m151$sample);hm$group<-factor(m151$group[match(hm$sample,m151$sample)],levels=c('HC','TC','SCI'))
write_table(hm,'Figure3_source_data')
p3<-ggplot(hm,aes(sample,gene,fill=z))+geom_tile()+facet_grid(~group,scales='free_x',space='free_x')+scale_fill_gradient2(low='#518AA0',mid='white',high='#B9564F',midpoint=0,limits=c(-3,3),oob=scales::squish)+labs(x=NULL,y=NULL,title='Frozen 30-gene signature in the discovery cohort',subtitle='Genes selected using discovery labels; descriptive illustration only',fill='Row z')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_text(size=7),legend.position='right')
save_plot(p3,'Figure3_revised',183,165)
long<-melt(as.data.table(s226[s226$group=='SCI',]),id.vars=c('sample','subject','time'),measure.vars=scores,variable.name='score',value.name='value');long$time<-factor(long$time,levels=visits);long$score<-factor(long$score,levels=scores,labels=c('Innate-associated','Lymphocyte/T-cell-associated','Balance (secondary)'))
p5a<-ggplot(long,aes(time,value,group=subject))+geom_line(color='#A7ACB3',linewidth=.35,alpha=.8)+geom_point(aes(color=score),size=1.6)+facet_wrap(~score,scales='free_y',nrow=1)+scale_color_manual(values=unname(palette[c('innate','adaptive','balance')]))+labs(x='Visit: acute n=10; 3m n=10; 6m n=10; 12m n=7',y='Control-referenced signature score',title='Within-participant trajectories')+theme(legend.position='none')
forest<-mixed[mixed$scope=='all_available',];forest$score<-factor(forest$score,levels=scores,labels=levels(long$score));forest$contrast<-factor(forest$contrast,levels=rev(unique(forest$contrast)))
p5b<-ggplot(forest,aes(estimate,contrast,color=score))+geom_vline(xintercept=0,linetype=2,color='#888888')+geom_errorbar(aes(xmin=lower.CL,xmax=upper.CL),orientation='y',width=.1)+geom_point(size=2)+facet_wrap(~score,scales='free_x',nrow=1)+scale_color_manual(values=unname(palette[c('innate','adaptive','balance')]))+labs(x='Later minus acute estimate (95% CI)',y=NULL,title='Random-intercept model; 10 participants / 37 samples')+theme(legend.position='none')
save_plot(p5a/p5b+plot_annotation(tag_levels='a',subtitle='Reduced-coverage projection: innate 8/15; lymphocyte/T-cell 13/15 genes'),'Figure5_revised',183,160)
p6<-ggplot(s151,aes(innate_score,adaptive_score,color=group))+geom_point(size=1.8,alpha=.85)+scale_color_manual(values=palette)+labs(x='Innate-associated score (TC reference)',y='Lymphocyte/T-cell-associated score (TC reference)',title='Discovery-cohort signature coordinates',subtitle='Descriptive representation; no permutation P value or independent validation claim')
save_plot(p6,'Figure6_revised',183,100)
plot293<-melt(as.data.table(s293),id.vars=c('sample','group','overlap'),measure.vars=scores,variable.name='score',value.name='value');plot293$group<-factor(plot293$group,levels=c('CTL','SCI'));plot293$score<-factor(plot293$score,levels=scores,labels=levels(long$score))
p7<-ggplot(plot293,aes(group,value,color=group))+geom_boxplot(outlier.shape=NA,width=.45)+geom_point(aes(shape=overlap),position=position_jitter(width=.12,seed=20261007),size=1.6)+facet_wrap(~score,scales='free_y',nrow=1)+scale_color_manual(values=palette)+labs(x='CTL n=9; SCI n=36',y='Control-referenced score',title='Additional-cohort reproduction with partial sample reuse',subtitle='10 SCI and all 9 CTL reused from GSE226238; estimates/CI in revision tables')
save_plot(p7,'Figure7_revised',183,100)
pstab<-ggplot(stab_new[stab_new$iterations==500,],aes(selection_frequency,color=factor(top_n)))+stat_ecdf()+facet_wrap(~module)+labs(x='Selection frequency',y='Candidate cumulative fraction',color='Top-N rule',title='Processed-matrix top-N sensitivity; 500 subsamples')
save_plot(pstab,'Supp_topN_sensitivity',183,95)
padj<-ggplot(adj,aes(estimate,spec,color=score))+geom_vline(xintercept=0,linetype=2)+geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation='y',width=.1)+geom_point()+facet_wrap(~score,scales='free_x')+labs(x='SCI minus TC estimate (95% CI)',y=NULL,title='Exploratory marker-adjusted sensitivity',subtitle='Disjoint MCP-counter markers; abundance proxies, not measured cell fractions')+theme(legend.position='none')
save_plot(padj,'Supp_marker_sensitivity',183,95)
# Legacy S8 model for audit only, never used as the unified signature result.
legacy<-supp[['S8_Longitudinal_Scores']];legacy$time<-factor(legacy$time,levels=c('0-3d','3m','6m','12m'))
capture.output(summary(lmer(refined_innate_score~time+(1|subject),data=legacy)),file=file.path(out,'logs/legacy_S8_innate_model_audit_only.txt'))
capture.output(sessionInfo(),file=file.path(out,'logs/sessionInfo.txt'))
saveRDS(list(discovery=E,longitudinal=E226,meta=meta,locked=locked,s151=s151,s226=s226,s293=s293,mixed=mixed,cross=cross,adj=adj,ladj=ladj,cam=cam,genes=genes),file.path(out,'tables/analysis_state.rds'))
cat('COMPLETED. Figures and full numerical tables written.\n')
sink()
