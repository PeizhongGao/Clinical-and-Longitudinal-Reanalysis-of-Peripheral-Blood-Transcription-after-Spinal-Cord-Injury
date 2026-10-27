# Public-package copy: see metadata/code_changes.csv for provenance.
args<-commandArgs(TRUE);root<-normalizePath(if(length(args))args[1] else '.',winslash='/');o<-file.path(root,'revision_R2_20261008');r1<-file.path(root,'revision_R1_20261007')
suppressPackageStartupMessages({library(data.table);library(lmerTest);library(emmeans);library(edgeR)})
set.seed(20261008);sink(file.path(o,'logs/sensitivity_analysis.txt'),split=TRUE)
w<-function(d,n)fwrite(as.data.frame(d),file.path(o,'tables',paste0(n,'.csv')),na='NA')
rd<-function(n)as.data.frame(fread(file.path(r1,'tables',paste0(n,'.csv'))))
meta<-as.data.frame(fread(file.path(o,'tables/sample_manifest_R2.csv'),na.strings=c('','NA')))
s<-readRDS(file.path(r1,'tables/analysis_state.rds'));scs<-c('innate_score','adaptive_score','balance_score');visits<-c('acute','3m','6m','12m')
ds<-list(GSE151371=rd('discovery_locked_scores_descriptive'),GSE226238=rd('NCBI_longitudinal_full_scores'),GSE293559=rd('original_Table_S16'))
for(n in names(ds)){d<-ds[[n]];mm<-meta[meta$dataset==n,];idx<-match(d$sample,mm$sample);stopifnot(!anyNA(idx));for(v in c('age','sex','baseline_ais','platform','overlap','subject'))d[[v]]<-mm[[v]][idx];d$sex<-factor(d$sex);d$platform<-factor(d$platform);ds[[n]]<-d}
desc<-list();cors<-list();balances<-list();balancecor<-list();cross<-list();diagnostics<-list();mixed<-list();designed<-list()
gfun<-function(x,y){nn<-length(x);mm<-length(y);(mean(x)-mean(y))/sqrt(((nn-1)*var(x)+(mm-1)*var(y))/(nn+mm-2))*(1-3/(4*(nn+mm)-9))}
for(n in names(ds)){
  d<-ds[[n]];d$stratum<-if(n=='GSE226238')d$time else d$group
  for(gr in unique(d$stratum)){
    dd<-d[d$stratum==gr,]
    for(sc in scs)desc[[paste(n,gr,sc)]]<-data.frame(dataset=n,stratum=gr,score=sc,n_samples=nrow(dd),n_participants=length(unique(dd$subject)),mean=mean(dd[[sc]]),sd=sd(dd[[sc]]),median=median(dd[[sc]]),q25=quantile(dd[[sc]],.25),q75=quantile(dd[[sc]],.75))
    cors[[paste(n,gr)]]<-data.frame(dataset=n,stratum=gr,n_samples=nrow(dd),n_participants=length(unique(dd$subject)),pearson=cor(dd$innate_score,dd$adaptive_score),spearman=cor(dd$innate_score,dd$adaptive_score,method='spearman'),interpretation='within-stratum association; not causal coupling')
  }
  ref<-if(n=='GSE151371')'TC' else 'CTL';ii<-d$group==ref
  d$balance_component_z<-(d$innate_score-mean(d$innate_score[ii]))/sd(d$innate_score[ii])-(d$adaptive_score-mean(d$adaptive_score[ii]))/sd(d$adaptive_score[ii])
  d$balance_rank<-(rank(d$innate_score,ties.method='average')-.5)/nrow(d)-(rank(d$adaptive_score,ties.method='average')-.5)/nrow(d)
  for(v in c('balance_component_z','balance_rank'))balancecor[[paste(n,v)]]<-data.frame(dataset=n,variant=v,spearman_original=cor(d$balance_score,d[[v]],method='spearman'),reference=ref,n_reference=sum(ii))
  balances[[n]]<-data.frame(dataset=n,d[,c('sample','subject','group','stratum','innate_score','adaptive_score','balance_score','balance_component_z','balance_rank')]);ds[[n]]<-d
}
w(rbindlist(desc),'score_descriptives');w(rbindlist(cors),'component_correlations');w(rbindlist(balances),'balance_sensitivity_scores');w(rbindlist(balancecor),'balance_scale_correlations')
fitlong<-function(d,sc,spec,extra=''){
  d<-d[d$group=='SCI',];d$time<-factor(d$time,levels=visits);d$subject<-factor(d$subject);d$age_c<-d$age-mean(d$age);d$severity_AB<-factor(ifelse(d$baseline_ais%in%c('A','B'),'AB','CD'))
  f<-as.formula(paste(sc,'~ time',extra,'+ (1|subject)'));fit<-lmer(f,d,REML=TRUE)
  tab<-as.data.frame(summary(contrast(emmeans(fit,~time,lmer.df='satterthwaite'),'trt.vs.ctrl',ref=1),infer=c(TRUE,TRUE),adjust='none'));tab$score<-sc;tab$spec<-spec;tab$n_samples<-nrow(d);tab$n_participants<-length(unique(d$subject));tab$model<-paste(deparse(f),collapse=' ')
  mixed[[paste(spec,sc)]]<<-tab;vc<-as.data.frame(VarCorr(fit));diagnostics[[paste(spec,sc)]]<<-data.frame(spec=spec,score=sc,rank=qr(model.matrix(fit))$rank,columns=ncol(model.matrix(fit)),singular=isSingular(fit),random_intercept_variance=vc$vcov[1],residual_variance=sigma(fit)^2,ICC=vc$vcov[1]/sum(vc$vcov),residual_Shapiro_P=shapiro.test(residuals(fit))$p.value,convergence=paste(fit@optinfo$conv$lme4$messages,collapse=';'))
}
for(sc in scs){fitlong(ds$GSE226238,sc,'age','+ age_c');fitlong(ds$GSE226238,sc,'sex','+ sex');fitlong(ds$GSE226238,sc,'severity_AB_CD','+ severity_AB')}
for(sc in c('balance_component_z','balance_rank'))fitlong(ds$GSE226238,sc,'balance_scale')
lmfit<-function(d,sc,spec,extra='',dataset='GSE293559'){
  d$case<-as.integer(d$group=='SCI');d$age_c<-d$age-mean(d$age,na.rm=TRUE)
  f<-as.formula(paste(sc,'~ case',extra));fit<-lm(f,d);co<-summary(fit)$coefficients;ci<-confint(fit,'case')
  cross[[paste(dataset,spec,sc)]]<<-data.frame(dataset=dataset,spec=spec,score=sc,n_samples=nobs(fit),n_participants=length(unique(d$subject)),estimate=co['case',1],lower.CL=ci[1],upper.CL=ci[2],p.value=co['case',4],model=paste(deparse(f),collapse=' '),rank=fit$rank,columns=ncol(model.matrix(fit)))
}
for(sc in scs){lmfit(ds$GSE293559,sc,'platform','+ platform');lmfit(ds$GSE293559,sc,'age_sex_platform','+ age_c + sex + platform')}
same<-list()
for(scope in c('same_platform_all','same_platform_new_SCI'))for(sc in scs){
 d<-ds$GSE293559;d<-d[d$platform=='GPL18573',];if(scope=='same_platform_new_SCI')d<-d[d$group=='CTL'|d$overlap=='new_GSE293559',]
 x<-d[[sc]][d$group=='SCI'];y<-d[[sc]][d$group=='CTL'];bs<-replicate(2000,gfun(sample(x,replace=TRUE),sample(y,replace=TRUE)));ci<-quantile(bs,c(.025,.975));tt<-t.test(x,y)
 same[[paste(scope,sc)]]<-data.frame(scope=scope,score=sc,n_SCI=length(x),n_CTL=length(y),g=gfun(x,y),g_ci_low=ci[1],g_ci_high=ci[2],p.value=wilcox.test(x,y,exact=FALSE)$p.value,mean_difference=mean(x)-mean(y),mean_ci_low=tt$conf.int[1],mean_ci_high=tt$conf.int[2])
}
same<-rbindlist(same);same[,q_BH:=p.adjust(p.value,'BH'),by=scope];w(same,'same_platform_sensitivity')
for(sc in c('balance_component_z','balance_rank'))lmfit(ds$GSE293559,sc,'balance_scale')
# All three contrasts are descriptive post-selection checks, not independent signature validation.
pairs<-list(c('TC','HC'),c('SCI','TC'),c('SCI','HC'));dd<-ds$GSE151371;comparisons<-list()
for(pair in pairs)for(sc in scs){x<-dd[[sc]][dd$group==pair[1]];y<-dd[[sc]][dd$group==pair[2]];tt<-t.test(x,y);bs<-replicate(2000,gfun(sample(x,replace=TRUE),sample(y,replace=TRUE)));ci<-quantile(bs,c(.025,.975));comparisons[[paste(pair[1],pair[2],sc)]]<-data.frame(contrast=paste(pair,collapse=' - '),score=sc,n_case=length(x),n_ref=length(y),mean_difference=mean(x)-mean(y),lower.CL=tt$conf.int[1],upper.CL=tt$conf.int[2],hedges_g=gfun(x,y),g_ci_low=ci[1],g_ci_high=ci[2],p.value=tt$p.value,scope='post-selection descriptive; intervals and P values do not account for gene selection')}
cc<-rbindlist(comparisons);cc$q_BH<-p.adjust(cc$p.value,'BH');w(cc,'discovery_three_group_descriptive')
# Official MCP-counter, using the source function and the entire available published marker list.
source(file.path(o,'inputs/MCPcounter.R'));markers<-read.delim(file.path(r1,'inputs/MCPcounter_genes.txt'),check.names=FALSE,stringsAsFactors=FALSE)
raw<-read.delim(gzfile(file.path(root,'data/GSE226238/GSE226238_raw_counts_GRCh38.p13_NCBI.tsv.gz')),check.names=FALSE,colClasses=c('character',rep('numeric',46)))
ann<-read.delim(gzfile(file.path(root,'data/GSE226238/Human.GRCh38.p13.annot.tsv.gz')),check.names=FALSE,colClasses='character');sym<-ann$Symbol[match(raw$GeneID,ann$GeneID)];ok<-!is.na(sym)&sym!='';count<-rowsum(as.matrix(raw[ok,-1]),sym[ok],reorder=FALSE);norm<-rd('NCBI_GSE226238_normalization');stopifnot(identical(colnames(count),norm$gsm));E226<-cpm(count,lib.size=norm$lib.size*norm$norm.factors,log=TRUE,prior.count=.5);colnames(E226)<-norm$sample
coverage<-list()
for(n in c('GSE151371','GSE226238')){
 E<-if(n=='GSE151371')s$discovery else E226;ms<-t(MCPcounter.estimate(E,featuresType='HUGO_symbols',genes=markers));w(data.frame(sample=rownames(ms),ms,check.names=FALSE),paste0(n,'_official_MCPcounter'))
 for(cell in unique(markers[['Cell population']])){gg<-unique(markers[['HUGO symbols']][markers[['Cell population']]==cell]);coverage[[paste(n,cell)]]<-data.frame(dataset=n,cell=cell,n_source=length(gg),n_used=sum(gg%in%rownames(E)),shared_frozen_genes=paste(intersect(gg,unlist(s$genes)),collapse=';'))}
 d<-ds[[n]];ix<-match(d$sample,rownames(ms));stopifnot(!anyNA(ix));d$neutrophil<-as.numeric(scale(ms[ix,'Neutrophils']));d$tcell<-as.numeric(scale(ms[ix,'T cells']))
 if(n=='GSE151371'){d<-d[d$group%in%c('SCI','TC'),];for(sc in scs)lmfit(d,sc,'official_MCPcounter','+ neutrophil + tcell',n)} else for(sc in scs)fitlong(d,sc,'official_MCPcounter','+ neutrophil + tcell')
}
w(rbindlist(coverage),'official_MCPcounter_coverage');mm<-rbindlist(mixed);mm[,q_BH:=p.adjust(p.value,'BH'),by=spec];w(mm,'R2_longitudinal_sensitivities');w(rbindlist(diagnostics),'R2_mixed_model_diagnostics');cr<-rbindlist(cross);cr[,q_BH:=p.adjust(p.value,'BH'),by=.(dataset,spec)];w(cr,'R2_cross_sectional_sensitivities')
# Participant PCA values and technical association audit. No submission date is treated as batch.
pca<-prcomp(t(s$discovery),center=TRUE,scale.=FALSE);pc<-data.frame(sample=rownames(pca$x),pca$x[,1:5]);pc$dataset<-'GSE151371'
d293<-fread(file.path(root,'code/GSE293559_Normalized_Morrison_2025_processed_data.txt'));X<-as.matrix(d293[,-1,with=FALSE]);storage.mode(X)<-'double';rownames(X)<-d293[[1]];stopifnot(nrow(X)==45,all(is.finite(X)));p293<-prcomp(X,center=TRUE,scale.=FALSE);p2<-data.frame(sample=rownames(p293$x),p293$x[,1:5]);p2$dataset<-'GSE293559'
pcs<-rbind(pc,p2);pcs<-merge(pcs,meta[,c('dataset','sample','group','platform','sex','age','overlap')],by=c('dataset','sample'),all.x=TRUE,sort=FALSE);w(pcs,'technical_PCA_scores');assoc<-list()
for(n in unique(pcs$dataset))for(v in c('group','platform','sex','age','overlap'))for(j in 1:5){d<-pcs[pcs$dataset==n,];x<-d[[v]];good<-!is.na(x);d<-d[good,];x<-x[good];if(length(unique(x))<2){assoc[[paste(n,v,j)]]<-data.frame(dataset=n,variable=v,PC=j,n=length(x),p.value=NA_real_,R2=NA_real_,status='unavailable_or_constant');next};d$y<-d[[paste0('PC',j)]];d$x<-if(v=='age')as.numeric(x) else factor(x);f<-lm(y~x,d);an<-anova(f);assoc[[paste(n,v,j)]]<-data.frame(dataset=n,variable=v,PC=j,n=nrow(d),p.value=an$'Pr(>F)'[1],R2=summary(f)$r.squared,status='exploratory_univariate_association')}
aa<-rbindlist(assoc);aa[,q_BH:=p.adjust(p.value,'BH'),by=dataset];w(aa,'technical_PC_associations');w(data.frame(dataset=rep(c('GSE151371','GSE293559'),each=5),PC=rep(1:5,2),variance_fraction=c(pca$sdev[1:5]^2/sum(pca$sdev^2),p293$sdev[1:5]^2/sum(p293$sdev^2))),'technical_PCA_variance')
capture.output(sessionInfo(),file=file.path(o,'logs/sensitivity_sessionInfo.txt'));cat('COMPLETE\n');sink()
