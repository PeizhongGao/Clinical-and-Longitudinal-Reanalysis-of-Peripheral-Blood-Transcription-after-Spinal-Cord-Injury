# Public-package copy: see metadata/code_changes.csv for provenance.
args <- commandArgs(TRUE)
root <- normalizePath(if (length(args)) args[1] else '.', winslash='/', mustWork=TRUE)
r1 <- file.path(root,'revision_R1_20261007')
out <- file.path(root,'revision_R2_20261008')
dir.create(file.path(out,'tables'),showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(out,'supplement'),showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(out,'figures'),showWarnings=FALSE,recursive=TRUE)
suppressPackageStartupMessages({library(data.table);library(limma);library(ggplot2)})
set.seed(20261012)
sink(file.path(out,'logs','odc405_adjusted_analysis.log'),split=TRUE)
write_csv <- function(x,name,where='supplement')
  fwrite(as.data.table(x),file.path(out,where,paste0(name,'.csv')),na='NA')

cat('ODC-SCI 405 revision sensitivity analysis; see ODC405_reanalysis_protocol.md\n')
cat('R:',R.version.string,'\n')
odc_file <- file.path(root,'data','odc-sci_405','TRACK-SCI transcriptomics - final version.csv')
dict_file <- file.path(root,'data','odc-sci_405','odc-sci_405_dict.csv')
matrix_file <- file.path(root,'explore','Report','GSE151371_normalized_expression_matrix_comBat_de-ID.csv')
odc_zip <- file.path(root,'data','odc-sci_405.zip')
stopifnot(file.exists(odc_file),file.exists(dict_file),file.exists(matrix_file),file.exists(odc_zip))

cell_cols <- c('B_cells_naive','B_cells_memory','Plasma_cells','T_cells_CD8',
  'T_cells_CD4_naive','T_cells_CD4_memory_resting','T_cells_CD4_memory_activated',
  'T_cells_follicular_helper','T_cells_regulatory_Tregs','T_cells_gamma_delta',
  'NK_cells_resting','NK_cells_activated','Monocytes','Macrophages_M0','Macrophages_M1',
  'Macrophages_M2','Dendritic_cells_resting','Dendritic_cells_activated',
  'Mast_cells_resting','Mast_cells_activated','Eosinophils','Neutrophils')
tcell_cols <- c('T_cells_CD8','T_cells_CD4_naive','T_cells_CD4_memory_resting',
  'T_cells_CD4_memory_activated','T_cells_follicular_helper','T_cells_regulatory_Tregs',
  'T_cells_gamma_delta')
clinical_cols <- c('Arbitrary_Patient_ID','Group','Sex','Age_at_time_of_Injury','Race',
  'Prior_CNS_Pathology','Concurrent_TBI','ISS','blood_draw_time',
  'Neurological_Level_of_Injury','NLI_grouped','ASIA_Impairment_Scale',
  'AIS_actual_or_estimate','Time_of_exam_in_days',cell_cols)
odc <- fread(odc_file,select=c(clinical_cols,'Gene','expression_value'),na.strings=c('NA',''))
id_col <- 'Arbitrary_Patient_ID'
meta_fields <- setdiff(clinical_cols,id_col)
meta_rows <- odc[is.na(Gene), ..clinical_cols]
if (nrow(meta_rows)!=58L || uniqueN(meta_rows[[id_col]])!=58L)
  stop('Expected one metadata-only row for each ODC participant.')
metadata_consistency <- meta_rows[,lapply(.SD,function(x)uniqueN(x,na.rm=FALSE)),
  by=id_col,.SDcols=meta_fields]
if (max(as.matrix(metadata_consistency[,..meta_fields])) > 1)
  stop('ODC metadata-only rows contain inconsistent values.')
meta <- unique(meta_rows,by=id_col)
setnames(meta,c(id_col,'Group','Sex','Age_at_time_of_Injury','Prior_CNS_Pathology',
  'Concurrent_TBI','ISS','blood_draw_time','Race','Neurological_Level_of_Injury',
  'NLI_grouped','ASIA_Impairment_Scale','AIS_actual_or_estimate','Time_of_exam_in_days'),
  c('sample','odc_group','sex','age','prior_cns','tbi','iss','draw_hours','race',
  'nli','nli_group','ais','ais_source','exam_days'))
meta[,`:=`(group=as.character(odc_group),sex=factor(as.character(sex)),age=as.numeric(age),
  prior_cns=factor(ifelse(is.na(prior_cns)|trimws(as.character(prior_cns))=='','Unknown',trimws(as.character(prior_cns))),levels=c('No','Yes','Unknown')),
  tbi=factor(ifelse(is.na(tbi)|trimws(as.character(tbi))=='','Unknown',trimws(as.character(tbi))),levels=c('No','Yes','Unknown')),
  race=factor(as.character(race)),iss=as.numeric(iss),draw_hours=as.numeric(draw_hours),
  exam_days=as.numeric(exam_days),ais=factor(as.character(ais)),nli=factor(as.character(nli)),
  nli_group=factor(as.character(nli_group)))]
meta[,tcell_fraction:=rowSums(.SD),.SDcols=tcell_cols]
meta[,neutrophil_fraction:=as.numeric(Neutrophils)]

# Existing GEO manifest supplies GSMs for the exact HC01/TC01/SCI01 sample keys.
manifest <- fread(file.path(out,'tables','sample_manifest_R2.csv'),na.strings=c('','NA'))
geo_manifest <- manifest[dataset=='GSE151371',.(dataset,sample,gsm,subject,manifest_group=group)]
stopifnot(nrow(meta)==58L,nrow(geo_manifest)==58L,uniqueN(meta$sample)==58L,
  uniqueN(geo_manifest$sample)==58L)
crosswalk <- merge(geo_manifest,meta,by='sample',all=TRUE,sort=FALSE,suffixes=c('_GEO','_ODC'))
if (nrow(crosswalk)!=58L || anyNA(crosswalk$gsm) || anyNA(crosswalk$odc_group) ||
    any(crosswalk$group!=crosswalk$odc_group) || any(crosswalk$group!=crosswalk$manifest_group))
  stop('ODC-to-GEO sample crosswalk is not a complete group-consistent bijection.')
crosswalk[,mapping_basis:='ODC study-specific group+number ID exactly matches the GEO sample key; GSM from existing GEO manifest']
crosswalk <- crosswalk[match(geo_manifest$sample,sample)]
meta[,gsm:=geo_manifest$gsm[match(sample,geo_manifest$sample)]]
meta[,subject:=geo_manifest$subject[match(sample,geo_manifest$sample)]]
stopifnot(!anyNA(meta$gsm),uniqueN(meta$gsm)==58L)
write_csv(crosswalk[,.(dataset,sample,gsm,subject,group,odc_group,age,sex,race,prior_cns,
  tbi,iss,draw_hours,ais,ais_source,nli,nli_group,exam_days,neutrophil_fraction,
  tcell_fraction,mapping_basis)],'ODC405_sample_crosswalk','tables')
fwrite(crosswalk[,.(n_participants=uniqueN(sample),n_samples=uniqueN(gsm)),by=group],
  file.path(out,'tables','ODC405_group_counts.csv'))

# Keep the original GEO ComBat matrix and exactly the 17,498 genes in the audited R1 fit.
orig_dt <- fread(matrix_file,check.names=FALSE)
Eall <- as.matrix(orig_dt[,-1,with=FALSE]);storage.mode(Eall)<-'double';rownames(Eall)<-orig_dt[[1]]
Eall <- Eall[!duplicated(rownames(Eall)),,drop=FALSE]
prior_fit <- fread(file.path(r1,'tables','SCI_vs_TC_processed_limma.csv'))
genes <- unique(prior_fit$gene)
if (length(genes)!=17498L || !all(genes%in%rownames(Eall))) stop('Original 17,498-gene tested set unavailable or changed.')
E <- Eall[genes,,drop=FALSE]
if (!setequal(colnames(E),meta$sample)) stop('Expression and ODC sample keys do not match.')
meta <- meta[match(colnames(E),meta$sample)]
stopifnot(identical(meta$sample,colnames(E)))

# Audit ODC's long expression field against GEO; do not substitute it for the GEO matrix.
odc_valid <- odc[!is.na(expression_value)&Gene%chin%rownames(Eall),
  .(ODC_value=mean(expression_value)),by=.(Gene,sample=Arbitrary_Patient_ID)]
odc_dups <- odc[!is.na(expression_value)&Gene%chin%rownames(Eall),
  .N,by=.(Gene,sample=Arbitrary_Patient_ID)][N>1]
odc_wide <- dcast(odc_valid,Gene~sample,value.var='ODC_value')
common_genes <- intersect(rownames(Eall),odc_wide$Gene)
odc_mat <- as.matrix(odc_wide[match(common_genes,Gene),colnames(Eall),with=FALSE])
rownames(odc_mat) <- common_genes
geo_common <- Eall[common_genes,,drop=FALSE]
matrix_match <- rbindlist(lapply(colnames(Eall),function(s){
  delta <- odc_mat[,s]-geo_common[,s]
  data.table(sample=s,n_common_genes=sum(is.finite(delta)),
    Pearson_r=cor(odc_mat[,s],geo_common[,s],use='complete.obs'),
    median_abs_difference=median(abs(delta),na.rm=TRUE),max_abs_difference=max(abs(delta),na.rm=TRUE))
}))
write_csv(matrix_match,'ODC405_expression_source_match','tables')
match_summary <- data.table(n_original_tested_genes=length(genes),n_common_exact_gene_symbols=length(common_genes),
  n_duplicate_gene_sample_pairs_in_ODC=nrow(odc_dups),median_sample_correlation=median(matrix_match$Pearson_r),
  minimum_sample_correlation=min(matrix_match$Pearson_r),maximum_abs_difference=max(matrix_match$max_abs_difference),
  interpretation='ODC expression is a concordance check only; models use the original GEO matrix and frozen tested-gene set.')
write_csv(match_summary,'ODC405_expression_source_match_summary','tables')

# Subject-level availability. Preserve explicit Unknown categories; do not impute numeric values.
availability_vars <- c('age','sex','race','prior_cns','tbi','iss','draw_hours','ais','nli',
  'neutrophil_fraction','tcell_fraction')
availability <- rbindlist(lapply(availability_vars,function(v) rbindlist(lapply(c('HC','TC','SCI'),function(g){
  z <- meta[group==g][[v]]
  valid <- if (is.factor(z)||is.character(z)) !is.na(z)&as.character(z)!='' else is.finite(z)
  data.table(variable=v,group=g,n=length(z),observed=sum(valid),missing=sum(!valid),n_unique=uniqueN(z[valid]),
    note=if(v%chin%c('ais','nli'))'SCI-specific; not a common case-control adjustment variable' else '')
}))))
write_csv(availability,'ODC405_covariate_availability','tables')
technical_audit <- data.table(variable=c('platform','instrument','library_strategy','library_selection','library_source',
  'batch','RIN','processing_date','infection','treatment'),available=c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE),
  disposition=c('constant in GEO metadata','constant in GEO metadata','constant in GEO metadata','constant in GEO metadata',
  'constant in GEO metadata','not supplied at participant level','not supplied at participant level',
  'not supplied; GEO submission date is not processing date','not supplied in ODC-405 dictionary/table',
  'not supplied in ODC-405 dictionary/table'),
  note='GEO describes the expression matrix as edgeR-normalized and ComBat batch-corrected; batch labels are unavailable for independent adjustment.')
write_csv(technical_audit,'ODC405_technical_variable_availability','tables')

# Limma models, with rank/df checks and BH handled separately inside each fit.
fit_de <- function(case,ref,spec,rhs,required=character()){
  d <- copy(meta[group%chin%c(case,ref)])
  if(length(required)) d <- d[complete.cases(d[,..required])]
  d[,group:=factor(group,levels=c(ref,case))]
  d[,age_c:=age-mean(age,na.rm=TRUE)];d[,iss_c:=iss-mean(iss,na.rm=TRUE)]
  d[,draw_hours_c:=draw_hours-mean(draw_hours,na.rm=TRUE)]
  d[,neutrophil_z:=as.numeric(scale(neutrophil_fraction))]
  d[,tcell_z:=as.numeric(scale(tcell_fraction))]
  design <- model.matrix(as.formula(paste('~',rhs)),data=d);coef_name<-paste0('group',case)
  if(!coef_name%in%colnames(design)||qr(design)$rank!=ncol(design))
    stop('Rank-deficient or non-estimable design: ',case,' vs ',ref,' ',spec)
  fit<-eBayes(lmFit(E[,d$sample,drop=FALSE],design))
  tt<-topTable(fit,coef=coef_name,number=Inf,sort.by='none',confint=TRUE);tt$gene<-rownames(tt)
  tt$model<-paste('expression ~',rhs);tt$spec<-spec;tt$contrast<-paste(case,'minus',ref)
  tt$n_samples<-nrow(d);tt$n_participants<-uniqueN(d$subject);tt$design_rank<-qr(design)$rank
  tt$design_columns<-ncol(design);tt$residual_df<-fit$df.residual[1]
  tt$case_n<-sum(d$group==case);tt$reference_n<-sum(d$group==ref)
  tt$FDR_lt_0_05_abs_logFC_gt_1<-tt$adj.P.Val<0.05&abs(tt$logFC)>1
  tt$FDR_lt_0_05<-tt$adj.P.Val<0.05
  name<-paste0('ODC405_DE_',case,'_vs_',ref,'_',spec)
  write_csv(tt,name)
  data.table(contrast=paste(case,'minus',ref),spec=spec,model=paste('expression ~',rhs),
    n_samples=nrow(d),n_participants=uniqueN(d$subject),case_n=sum(d$group==case),reference_n=sum(d$group==ref),
    design_rank=qr(design)$rank,design_columns=ncol(design),residual_df=fit$df.residual[1],n_tested_genes=nrow(tt),
    n_FDR_lt_0_05=sum(tt$FDR_lt_0_05),n_FDR_and_abs_logFC_gt_1=sum(tt$FDR_lt_0_05_abs_logFC_gt_1),
    median_logFC=median(tt$logFC))
}
de <- list()
for(pair in list(c('SCI','TC'),c('SCI','HC'),c('TC','HC')))
  de[[paste0(pair[1],'_',pair[2],'_unadjusted')]] <- fit_de(pair[1],pair[2],'unadjusted','group')
de[['SCI_TC_age_sex']] <- fit_de('SCI','TC','age_sex','group + age_c + sex',c('age','sex'))
de[['SCI_TC_clinical']] <- fit_de('SCI','TC','clinical_complete_case',
  'group + age_c + sex + iss_c + tbi + prior_cns',c('age','sex','iss','tbi','prior_cns'))
de[['SCI_TC_clinical_time']] <- fit_de('SCI','TC','clinical_plus_draw_time',
  'group + age_c + sex + iss_c + tbi + prior_cns + draw_hours_c',c('age','sex','iss','tbi','prior_cns','draw_hours'))
de[['SCI_TC_clinical_cells']] <- fit_de('SCI','TC','clinical_plus_CIBERSORTx',
  'group + age_c + sex + iss_c + tbi + prior_cns + neutrophil_z + tcell_z',
  c('age','sex','iss','tbi','prior_cns','neutrophil_fraction','tcell_fraction'))
de[['SCI_HC_adjusted']] <- fit_de('SCI','HC','age_sex_prior_CNS','group + age_c + sex + prior_cns',c('age','sex','prior_cns'))
de[['TC_HC_adjusted']] <- fit_de('TC','HC','age_sex','group + age_c + sex',c('age','sex'))
de_summary<-rbindlist(de,fill=TRUE)
write_csv(de_summary,'ODC405_DE_model_summary','tables')

# Baseline exact-reproduction check for the primary unadjusted SCI-TC fit.
new_unadj<-fread(file.path(out,'supplement','ODC405_DE_SCI_vs_TC_unadjusted.csv'))
old_unadj<-fread(file.path(r1,'tables','SCI_vs_TC_processed_limma.csv'))
cmp<-merge(old_unadj[,.(gene,old_logFC=logFC,old_t=t,old_FDR=adj.P.Val)],
  new_unadj[,.(gene,new_logFC=logFC,new_t=t,new_FDR=adj.P.Val)],by='gene')
reproduction<-data.table(n_genes=nrow(cmp),max_abs_logFC_difference=max(abs(cmp$old_logFC-cmp$new_logFC)),
  max_abs_t_difference=max(abs(cmp$old_t-cmp$new_t)),max_abs_FDR_difference=max(abs(cmp$old_FDR-cmp$new_FDR)))
write_csv(reproduction,'ODC405_unadjusted_reproduction_check','tables')

# Reuse the frozen TC-referenced score file; no genes are selected or restandardized.
scores<-fread(file.path(out,'supplement','discovery_locked_scores_descriptive.csv'))
score_data<-merge(meta,scores[,.(sample,innate_score,adaptive_score,balance_score)],by='sample',all.x=TRUE,sort=FALSE)
score_data<-score_data[group%chin%c('SCI','TC')];score_data[,group:=factor(group,levels=c('TC','SCI'))]
score_specs<-list(
  unadjusted=list(rhs='group',req=character()),
  age_sex=list(rhs='group + age_c + sex',req=c('age','sex')),
  clinical_complete_case=list(rhs='group + age_c + sex + iss_c + tbi + prior_cns',req=c('age','sex','iss','tbi','prior_cns')),
  clinical_plus_draw_time=list(rhs='group + age_c + sex + iss_c + tbi + prior_cns + draw_hours_c',req=c('age','sex','iss','tbi','prior_cns','draw_hours')),
  clinical_plus_CIBERSORTx=list(rhs='group + age_c + sex + iss_c + tbi + prior_cns + neutrophil_z + tcell_z',req=c('age','sex','iss','tbi','prior_cns','neutrophil_fraction','tcell_fraction')))
hc3<-function(fit,term){
  X<-model.matrix(fit);j<-match(term,colnames(X));stopifnot(!is.na(j));bread<-solve(crossprod(X))
  h<-rowSums((X%*%bread)*X);a<-as.vector(X%*%bread[,j]);se<-sqrt(sum((a*residuals(fit)/pmax(1-h,1e-8))^2))
  df<-df.residual(fit);est<-coef(fit)[j];stat<-est/se
  data.table(estimate=est,robust_SE_HC3=se,statistic=stat,df=df,lower_CL=est-qt(.975,df)*se,
    upper_CL=est+qt(.975,df)*se,p_value=2*pt(-abs(stat),df),hat_max=max(h),design_rank=fit$rank,design_columns=ncol(X))
}
score_res<-list()
for(spec in names(score_specs)) for(outcome in c('innate_score','adaptive_score','balance_score')){
  d<-copy(score_data);d[,age_c:=age-mean(age,na.rm=TRUE)];d[,iss_c:=iss-mean(iss,na.rm=TRUE)]
  d[,draw_hours_c:=draw_hours-mean(draw_hours,na.rm=TRUE)]
  d[,neutrophil_z:=as.numeric(scale(neutrophil_fraction))];d[,tcell_z:=as.numeric(scale(tcell_fraction))]
  req<-unique(c(score_specs[[spec]]$req,outcome));d<-d[complete.cases(d[,..req])]
  fit<-lm(as.formula(paste(outcome,'~',score_specs[[spec]]$rhs)),data=d);res<-hc3(fit,'groupSCI')
  res[,`:=`(spec=spec,outcome=outcome,model=paste(deparse(formula(fit)),collapse=' '),n_samples=nobs(fit),
    n_participants=uniqueN(d$subject),SCI_n=sum(d$group=='SCI'),TC_n=sum(d$group=='TC'),postselection=TRUE)]
  score_res[[paste(spec,outcome)]]<-res
}
score_results<-rbindlist(score_res,fill=TRUE);score_results[,q_BH:=p.adjust(p_value,'BH'),by=spec]
write_csv(score_results,'ODC405_locked_score_models')
write_csv(score_data[,.(sample,gsm,subject,group,innate_score,adaptive_score,balance_score,age,sex,race,
  prior_cns,tbi,iss,draw_hours,ais,nli,neutrophil_fraction,tcell_fraction)],
  'ODC405_locked_scores_with_clinical_and_cell_estimates','tables')

# CIBERSORTx group summaries; fractions are relative algorithmic estimates, not measured counts.
cell_summary<-rbindlist(lapply(cell_cols,function(v){
  z<-meta[,.(group,value=as.numeric(get(v)))][is.finite(value)];case<-z[group=='SCI',value];ref<-z[group=='TC',value]
  wt<-suppressWarnings(wilcox.test(case,ref,exact=FALSE));tt<-t.test(case,ref)
  data.table(cell=v,n_SCI=length(case),n_TC=length(ref),mean_SCI=mean(case),sd_SCI=sd(case),
    mean_TC=mean(ref),sd_TC=sd(ref),median_SCI=median(case),median_TC=median(ref),
    mean_diff_SCI_minus_TC=mean(case)-mean(ref),Welch_CI_low=tt$conf.int[1],Welch_CI_high=tt$conf.int[2],wilcoxon_P=wt$p.value)
}))
cell_summary[,q_BH_22:=p.adjust(wilcoxon_P,'BH')]
write_csv(cell_summary,'ODC405_CIBERSORTx_cell_fraction_SCI_vs_TC')
cell_long<-melt(meta,id.vars=c('sample','group'),measure.vars='Neutrophils',variable.name='cell',value.name='fraction')
cell_long<-rbind(cell_long,meta[,.(sample,group,cell='Aggregate T-cell subsets',fraction=tcell_fraction)])
cell_long[,percent:=100*fraction]
write_csv(cell_long,'ODC405_CIBERSORTx_focused_cell_estimates','tables')
focused_cells <- rbindlist(lapply(c('Neutrophils','Aggregate T-cell subsets'),function(cell){
  z <- if(cell=='Neutrophils') meta[,.(group,value=as.numeric(Neutrophils))] else
    meta[,.(group,value=as.numeric(tcell_fraction))]
  z <- z[is.finite(value) & group%chin%c('SCI','TC')]
  case <- z[group=='SCI',value]; ref <- z[group=='TC',value]
  wt <- suppressWarnings(wilcox.test(case,ref,exact=FALSE)); tt <- t.test(case,ref)
  data.table(cell=cell,n_SCI=length(case),n_TC=length(ref),
    mean_SCI=mean(case),sd_SCI=sd(case),median_SCI=median(case),
    mean_TC=mean(ref),sd_TC=sd(ref),median_TC=median(ref),
    mean_difference_SCI_minus_TC=mean(case)-mean(ref),
    wilcoxon_P=wt$p.value,Welch_CI_low=tt$conf.int[1],Welch_CI_high=tt$conf.int[2])
}))
focused_cells[,q_BH_2:=p.adjust(wilcoxon_P,'BH')]
write_csv(focused_cells,'ODC405_CIBERSORTx_focused_group_summary','tables')
fraction_sum<-rowSums(meta[,..cell_cols])
fraction_audit<-data.table(n_samples=length(fraction_sum),sum_min=min(fraction_sum),sum_median=median(fraction_sum),
  sum_max=max(fraction_sum),interpretation='Observed components sum to approximately 1. The dictionary labels them percent; models retain the supplied relative scale.')
write_csv(fraction_audit,'ODC405_CIBERSORTx_fraction_sum_audit','tables')

# Exploratory associations of first five expression PCs with available clinical/cell variables.
pca<-prcomp(t(E),center=TRUE,scale.=FALSE)
pc_scores<-as.data.table(pca$x[,1:5,drop=FALSE],keep.rownames='sample')
pc_data<-merge(pc_scores,meta,by='sample',all.x=TRUE,sort=FALSE)
pc_assoc<-list();pc_vars<-c('group','age','sex','race','prior_cns','tbi','iss','draw_hours','ais','nli','neutrophil_fraction','tcell_fraction')
for(pcn in paste0('PC',1:5)) for(v in pc_vars){
  z<-pc_data[[v]];y<-pc_data[[pcn]];good<-if(is.numeric(z))is.finite(z)&is.finite(y) else !is.na(z)&is.finite(y)
  d<-data.frame(y=y[good],x=z[good]);if(!nrow(d)||length(unique(d$x))<2)next
  d$x<-if(is.numeric(d$x))d$x else droplevels(factor(d$x));fit<-lm(y~x,d);an<-anova(fit)
  pc_assoc[[paste(pcn,v)]]<-data.table(PC=pcn,variable=v,n=nrow(d),df_num=an$Df[1],df_den=df.residual(fit),
    F=an$`F value`[1],p_value=an$`Pr(>F)`[1],R_squared=summary(fit)$r.squared,
    status='exploratory univariate association; not a batch test')
}
pc_assoc<-rbindlist(pc_assoc,fill=TRUE);pc_assoc[,q_BH:=p.adjust(p_value,'BH'),by=PC]
write_csv(pc_assoc,'ODC405_PC_covariate_associations','tables')
pc_var<-data.table(PC=paste0('PC',1:5),variance_fraction=pca$sdev[1:5]^2/sum(pca$sdev^2))
write_csv(pc_var,'ODC405_PCA_variance','tables')

# Supplementary diagnostics: adjusted frozen scores, inferred cells, and PCA.
forest<-copy(score_results)
forest[,outcome:=factor(outcome,levels=c('innate_score','adaptive_score','balance_score'),
  labels=c('Innate-associated score','Adaptive-associated score','Secondary balance score'))]
forest[,spec:=factor(spec,levels=names(score_specs))]
forest[,ypos:=as.numeric(outcome)+(as.numeric(spec)-3)*.12]
p1<-ggplot(forest,aes(y=ypos,x=estimate,colour=spec))+geom_vline(xintercept=0,linetype=2,colour='grey55')+
  geom_segment(aes(x=lower_CL,xend=upper_CL,y=ypos,yend=ypos),linewidth=.55)+geom_point(size=2)+
  scale_y_continuous(breaks=1:3,labels=c('Innate-associated score','Adaptive-associated score','Secondary balance score'))+
  labs(x='Adjusted SCI minus TC score difference (HC3 95% CI)',y=NULL,colour='Model',
    title='Frozen discovery scores under clinical and cell-abundance sensitivities',
    subtitle='Post-selection estimates; CIBERSORTx values are expression-derived proxies')+
  theme_bw(base_size=10)+theme(legend.position='bottom',panel.grid.minor=element_blank())
ggsave(file.path(out,'figures','ODC405_adjusted_score_sensitivities.pdf'),p1,width=9,height=5.6)
ggsave(file.path(out,'figures','ODC405_adjusted_score_sensitivities.png'),p1,width=9,height=5.6,dpi=300)
p2<-ggplot(cell_long,aes(x=group,y=percent,colour=group,fill=group))+geom_boxplot(width=.55,alpha=.16,outlier.shape=NA)+
  geom_jitter(width=.08,height=0,size=2,alpha=.8)+facet_wrap(~cell,scales='free_y')+
  labs(x=NULL,y='CIBERSORTx relative estimate (%)',title='Inferred neutrophil and aggregate T-cell abundance by group',
    subtitle='Algorithmic estimates from bulk expression; not measured blood counts')+
  theme_bw(base_size=10)+theme(legend.position='none',panel.grid.minor=element_blank())
ggsave(file.path(out,'figures','ODC405_CIBERSORTx_focus_by_group.pdf'),p2,width=7.2,height=4.5)
ggsave(file.path(out,'figures','ODC405_CIBERSORTx_focus_by_group.png'),p2,width=7.2,height=4.5,dpi=300)
pca_plot<-ggplot(pc_data,aes(x=PC1,y=PC2,colour=group,shape=sex,size=100*neutrophil_fraction))+geom_point(alpha=.85)+
  labs(x=sprintf('PC1 (%.1f%%)',100*pc_var$variance_fraction[1]),y=sprintf('PC2 (%.1f%%)',100*pc_var$variance_fraction[2]),
    size='Inferred neutrophils (%)',title='GSE151371 expression PCA with linked sample metadata',
    subtitle='Original tested-gene matrix; technical batch labels and RIN are unavailable')+
  theme_bw(base_size=10)+theme(panel.grid.minor=element_blank())
ggsave(file.path(out,'figures','ODC405_PCA_clinical_cell_covariates.pdf'),pca_plot,width=7.2,height=5.5)
ggsave(file.path(out,'figures','ODC405_PCA_clinical_cell_covariates.png'),pca_plot,width=7.2,height=5.5,dpi=300)

write_csv(data.table(file=c(odc_file,dict_file,matrix_file,odc_zip),
  md5=unname(tools::md5sum(c(odc_file,dict_file,matrix_file,odc_zip)))),
  'ODC405_input_hashes','tables')
capture.output(sessionInfo(),file=file.path(out,'logs','odc405_R_sessionInfo.txt'))
cat('Mapping:',nrow(crosswalk),'samples; common expression genes:',length(common_genes),
  '; median per-sample correlation:',median(matrix_match$Pearson_r),
  '; minimum:',min(matrix_match$Pearson_r),'\n')
print(de_summary);print(score_results)
print(cell_summary[cell%chin%c('Neutrophils',tcell_cols)])
cat('COMPLETE\n')
sink()
