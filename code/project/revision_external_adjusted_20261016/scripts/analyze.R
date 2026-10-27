# Public-package copy: see metadata/code_changes.csv for provenance.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(lmerTest);library(emmeans);library(edgeR)})
root<-normalizePath(commandArgs(TRUE)[1],winslash='/'); base<-file.path(root,'revision_integrated_20261012'); out<-file.path(root,'revision_external_adjusted_20261016')
sink(file.path(out,'logs/analysis_run.txt'),split=TRUE)
cat('Started UTC',format(Sys.time(),tz='UTC'),'\n')
w<-function(x,n)fwrite(as.data.table(x),file.path(out,'tables',paste0(n,'.csv')),na='NA')
rd<-function(p)fread(p,na.strings=c('NA',''))
inputs<-c(scores=file.path(base,'supplement/external_BTM_scores.csv'),members=file.path(base,'supplement/external_BTM_members_annotated.csv'),used=file.path(base,'supplement/external_BTM_used_genes.csv'),coverage=file.path(base,'supplement/external_BTM_coverage.csv'),
 clinical=file.path(base,'fixed_sample/inputs/ODC405_locked_scores_with_clinical_and_cell_estimates.csv'),
 full48=file.path(base,'fixed_sample/tables/fixed_input_full48.csv'),clinical43=file.path(base,'fixed_sample/tables/fixed_input_clinical43.csv'),time41=file.path(base,'fixed_sample/tables/fixed_input_time41.csv'),
 official=file.path(base,'tables/GSE226238_official_MCPcounter.csv'),metadata=file.path(root,'revision_R2_20261008/tables/sample_manifest_R2.csv'),
 markers=file.path(root,'revision_R1_20261007/inputs/MCPcounter_genes.txt'),frozen=file.path(base,'supplement/locked_signature_reused_S14.csv'),
 state=file.path(root,'revision_R1_20261007/tables/analysis_state.rds'),raw=file.path(root,'data/GSE226238/GSE226238_raw_counts_GRCh38.p13_NCBI.tsv.gz'),annotation=file.path(root,'data/GSE226238/Human.GRCh38.p13.annot.tsv.gz'),norm=file.path(root,'revision_R1_20261007/tables/NCBI_GSE226238_normalization.csv'),discovery_map=file.path(base,'supplement/external_BTM_discovery_symbol_mapping.csv'),count_map=file.path(base,'supplement/external_BTM_NCBI_symbol_mapping.csv'))
stopifnot(all(file.exists(inputs)))
w(data.table(role=names(inputs),path=unname(inputs),md5=unname(tools::md5sum(inputs))),'input_hashes_before')
for(k in names(inputs)[!names(inputs)%in%c('state','raw','annotation')])file.copy(inputs[k],file.path(out,'inputs',basename(inputs[k])),overwrite=TRUE)
sc<-rd(inputs['scores']); used<-rd(inputs['used']); members<-rd(inputs['members']); clinical<-rd(inputs['clinical']); meta<-rd(inputs['metadata']); markers<-read.delim(inputs['markers'],check.names=FALSE,stringsAsFactors=FALSE)
mods<-c('M37.1','M163','M7.0','M223'); variants<-c('published','exclude_original106')
stopifnot(setequal(unique(sc$module_id),mods),!anyDuplicated(sc[,.(dataset,sample,variant,module_id)]))
frozen<-rd(inputs['frozen']); cat('Frozen columns',names(frozen),'\n')
frozen_col<-intersect(c('gene','symbol'),names(frozen))[1];stopifnot(!is.na(frozen_col))
excluded<-union(na.omit(members$symbol),frozen[[frozen_col]])
marker_cells<-c('Neutrophils','T cells')
mk<-rbindlist(lapply(marker_cells,function(cellname){genes<-unique(markers[['HUGO symbols']][markers[['Cell population' ]]==cellname]);data.table(cell=cellname,gene=genes,excluded_overlap=genes%in%excluded)}))
s<-readRDS(inputs['state']); E<-s$discovery;emap<-rd(inputs['discovery_map']); sy<-emap$symbol[match(rownames(E),emap$source_token)];ok<-!is.na(sy)
E151<-rowsum(E[ok,,drop=FALSE],sy[ok],reorder=FALSE);E151<-E151/as.numeric(table(sy[ok])[rownames(E151)])
raw<-read.delim(gzfile(inputs['raw']),check.names=FALSE,colClasses=c('character',rep('numeric',46)))
ann<-read.delim(gzfile(inputs['annotation']),check.names=FALSE,colClasses='character');sym<-ann$Symbol[match(raw$GeneID,ann$GeneID)];amap<-rd(inputs['count_map']);sym<-amap$symbol[match(sym,amap$source_token)];ok<-!is.na(sym)&sym!=''
counts<-rowsum(as.matrix(raw[ok,-1]),sym[ok],reorder=FALSE);norm<-rd(inputs['norm']);stopifnot(identical(colnames(counts),norm$gsm))
E226<-cpm(counts,lib.size=norm$lib.size*norm$norm.factors,log=TRUE,prior.count=.5);colnames(E226)<-norm$sample
proxy<-list(); marker_audit<-list();scorecheck<-list()
for(dsid in c('GSE151371','GSE226238')){
 E<-if(dsid=='GSE151371')E151 else E226
 for(cellname in marker_cells){
  allg<-mk[cell==cellname & excluded_overlap==FALSE,gene]
  avail<-intersect(allg,rownames(E));avail<-avail[vapply(avail,function(g)all(is.finite(E[g,])),logical(1))]
  eligible<-length(avail)>=3 & length(avail)/length(allg)>=.5
  marker_audit[[paste(dsid,cellname)]]<-data.table(dataset=dsid,cell=cellname,gene=mk[cell==cellname,gene],excluded_overlap=mk[cell==cellname,excluded_overlap],available=mk[cell==cellname,gene]%in%avail,remaining_denominator=length(allg),n_used=length(avail),coverage=length(avail)/length(allg),eligible=eligible)
  proxy[[paste(dsid,cellname)]]<-data.table(dataset=dsid,sample=colnames(E),cell=cellname,value=if(eligible)colMeans(E[avail,,drop=FALSE])else NA_real_)
 }
 for(v in variants)for(id in mods){
  genes<-used[dataset==dsid & variant==v & module_id==id,gene]
  dd<-sc[dataset==dsid & variant==v & module_id==id];j<-match(dd$sample,colnames(E));stopifnot(!anyNA(j),all(genes%in%rownames(E)))
  scorecheck[[paste(dsid,v,id)]]<-data.table(dataset=dsid,variant=v,module_id=id,n=nrow(dd),max_abs_difference=max(abs(colMeans(E[genes,j,drop=FALSE])-dd$score)))
 }
}
ma<-rbindlist(marker_audit); px<-rbindlist(proxy);ck<-rbindlist(scorecheck);stopifnot(all(ck$max_abs_difference<1e-10),all(ma$eligible))
w(ma,'marker_overlap_coverage');w(px,'disjoint_marker_scores');w(ck,'prior_score_reproduction')
req<-unique(ma[available==TRUE,gene]);w(data.table(gene=intersect(req,rownames(E151)),E151[intersect(req,rownames(E151)),,drop=FALSE]),'disjoint_marker_discovery_expression');w(data.table(gene=intersect(req,rownames(E226)),E226[intersect(req,rownames(E226)),,drop=FALSE]),'disjoint_marker_longitudinal_logCPM')
# Fixed sets: exact previous membership and no new clinical imputations.
sets<-lapply(inputs[c('full48','clinical43','time41')],rd);names(sets)<-c('full48','clinical43','time41')
rhs<-c(unadjusted='group',age_sex='group + age_c + sex',clinical='group + age_c + sex + iss_c + tbi + prior_cns',clinical_time='group + age_c + sex + iss_c + tbi + prior_cns + draw_hours_c',clinical_CIB='group + age_c + sex + iss_c + tbi + prior_cns + neut_CIB + tcell_CIB',clinical_disjoint='group + age_c + sex + iss_c + tbi + prior_cns + neut_disjoint + tcell_disjoint')
specs<-list(full48=c('unadjusted','age_sex'),clinical43=c('unadjusted','age_sex','clinical','clinical_CIB','clinical_disjoint'),time41=c('unadjusted','age_sex','clinical','clinical_time'))
for(n in names(sets)){
 d<-sets[[n]];stopifnot(!anyDuplicated(d$sample)); d[,group:=factor(group,levels=c('TC','SCI'))];d[,sex:=factor(sex,levels=c('F','M'))]
 for(v in c('tbi','prior_cns'))set(d,j=v,value=factor(d[[v]],levels=c('No','Yes','Unknown')))
 for(v in c('age','iss','draw_hours'))d[,(paste0(v,'_c')):=get(v)-mean(get(v),na.rm=TRUE)]
 for(cellname in marker_cells){z<-px[dataset=='GSE151371' & cell==cellname];d[,(if(cellname=='Neutrophils')'neut_disjoint' else 'tcell_disjoint'):=as.numeric(scale(z$value[match(sample,z$sample)]))]}
 d[,neut_CIB:=as.numeric(scale(neutrophil_fraction))];d[,tcell_CIB:=as.numeric(scale(tcell_fraction))];sets[[n]]<-d;w(d,paste0('fixed_input_',n))
}
hc3<-function(d){
 f<-lm(as.formula(paste('score~',rhs[d$spec[1]])),d,na.action=na.fail,x=TRUE,y=TRUE)
 X<-model.matrix(f);j<-match('groupSCI',colnames(X));B<-solve(crossprod(X));h<-hatvalues(f);e<-residuals(f);stopifnot(f$rank==ncol(X),all(h<1-1e-8))
 V<-B%*%crossprod(X,X*(e/(1-h))^2)%*%B;se<-sqrt(V[j,j]);est<-unname(coef(f)[j]);df<-df.residual(f)
 gr<-residuals(lm.fit(X[,-j,drop=FALSE],X[,j]));vif<-sum((X[,j]-mean(X[,j]))^2)/sum(gr^2)
 tab<-data.table(estimate=est,SE=se,df=df,lower_CL=est-qt(.975,df)*se,upper_CL=est+qt(.975,df)*se,p_value=2*pt(-abs(est/se),df),n_observations=nrow(d),n_participants=uniqueN(d$subject),SCI_n=sum(d$group=='SCI'),TC_n=sum(d$group=='TC'),rank=f$rank,columns=ncol(X),max_leverage=max(h),group_VIF=vif,shapiro_p=shapiro.test(e)$p.value,BP_p=lmtest::bptest(f)$p.value,model=paste(deparse(formula(f)),collapse=' '))
 inf<-data.table(sample=d$sample,gsm=d$gsm,group=d$group,fitted=fitted(f),residual=e,hat=h,cook=cooks.distance(f),dfbetas_group=dfbetas(f)[,j])
 list(tab=tab,inf=inf,fit=f)
}
res<-list(); infl<-list(); modeldata<-list(); loo<-list()
for(n in names(sets))for(spec in specs[[n]])for(v in variants)for(id in mods){
 d<-copy(sets[[n]]);z<-sc[dataset=='GSE151371' & variant==v & module_id==id];d[,score:=z$score[match(sample,z$sample)]];d[,spec:=spec]
 k<-paste(n,spec,v,id,sep='__');h<-hc3(d);h$tab[,`:=`(subset=n,spec=spec,variant=v,module_id=id)];res[[k]]<-h$tab;h$inf[,`:=`(subset=n,spec=spec,variant=v,module_id=id)];infl[[k]]<-h$inf;d[,`:=`(subset=n,variant=v,module_id=id)];modeldata[[k]]<-d
 if((n=='clinical43' & spec%in%c('clinical','clinical_CIB','clinical_disjoint'))|(n=='time41' & spec=='clinical_time'))for(samp in d[group=='TC',sample]){
  a<-hc3(d[sample!=samp])$tab;a[,`:=`(subset=n,spec=spec,variant=v,module_id=id,deleted_TC=samp,full_estimate=h$tab$estimate)];loo[[paste(k,samp)]]<-a
 }
}
cr<-rbindlist(res);stopifnot(nrow(cr)==88)
cr[,q_BH_within_spec:=p.adjust(p_value,'BH',n=4),by=.(subset,spec,variant)];cr[,q_BH_all88:=p.adjust(p_value,'BH',n=88)]
w(rbindlist(infl),'discovery_influence');w(rbindlist(modeldata,fill=TRUE),'discovery_model_inputs');ll<-rbindlist(loo);ll[,q_BH_within_deletion:=p.adjust(p_value,'BH',n=4),by=.(subset,spec,variant,deleted_TC)];w(ll,'TC_leave_one_out')
desc<-rbindlist(lapply(names(sets),function(n){d<-merge(sets[[n]][,.(sample,subject,group)],sc[dataset=='GSE151371',.(sample,variant,module_id,score)],by='sample');d[,.(n=.N,mean=mean(score),sd=sd(score),median=median(score),q25=quantile(score,.25),q75=quantile(score,.75)),by=.(variant,module_id,group)][,subset:=n]}));w(desc,'discovery_descriptives')
# Longitudinal: existing subject metadata, separate clinical sensitivities.
md<-meta[dataset=='GSE226238' & group=='SCI'];stopifnot(nrow(md)==37,uniqueN(md$subject)==10)
visitmap<-unique(sc[dataset=='GSE226238',.(sample,time)]);stopifnot(!anyDuplicated(visitmap$sample));md[,time:=visitmap$time[match(sample,visitmap$sample)]]
official<-rd(inputs['official']);md[,neut_official:=official[['Neutrophils']][match(sample,official$sample)]];md[,tcell_official:=official[['T cells']][match(sample,official$sample)]]
for(cellname in marker_cells){z<-px[dataset=='GSE226238' & cell==cellname];md[,(if(cellname=='Neutrophils')'neut_disjoint' else 'tcell_disjoint'):=z$value[match(sample,z$sample)]]}
md[,time:=factor(time,levels=c('acute','3m','6m','12m'))];md[,subject:=factor(subject)];md[,sex:=factor(sex)];md[,age_c:=age-mean(age)];md[,severity_AB:=factor(ifelse(baseline_ais%in%c('A','B'),'AB','CD'))]
complete<-md[,.(complete=all(c('acute','3m','6m','12m')%in%time)),by=subject][complete==TRUE,subject];stopifnot(length(complete)==7)
extras<-c(baseline='',official_both='+neut_official+tcell_official',disjoint_both='+neut_disjoint+tcell_disjoint',disjoint_neut='+neut_disjoint',disjoint_T='+tcell_disjoint',age='+age_c',sex='+sex',severity='+severity_AB')
lr<-list();diag<-list();longinput<-list();longresid<-list();fits<-list()
for(scope in c('all37','complete28'))for(spec in if(scope=='all37')names(extras)else c('baseline','disjoint_both'))for(v in variants)for(id in mods){
 d<-copy(if(scope=='all37')md else md[subject%in%complete]);d[,subject:=droplevels(subject)]
 for(z in c('neut_official','tcell_official','neut_disjoint','tcell_disjoint'))d[,(z):=as.numeric(scale(get(z)))]
 zz<-sc[dataset=='GSE226238' & variant==v & module_id==id];d[,score:=zz$score[match(sample,zz$sample)]]
 k<-paste(scope,spec,v,id,sep='__');f<-as.formula(paste('score~time',extras[spec],'+(1|subject)'))
 fit<-tryCatch(lmer(f,d,REML=TRUE,na.action=na.fail),error=function(e)e)
 if(inherits(fit,'error'))stop(paste(k,conditionMessage(fit)))
 t<-as.data.table(summary(contrast(emmeans(fit,~time,lmer.df='satterthwaite'),'trt.vs.ctrl',ref=1),infer=c(TRUE,TRUE),adjust='none'))
 setnames(t,c('lower.CL','upper.CL','p.value'),c('lower_CL','upper_CL','p_value'));t[,`:=`(scope=scope,spec=spec,variant=v,module_id=id,n_participants=uniqueN(d$subject),n_observations=nrow(d),model=paste(deparse(f),collapse=' '))];lr[[k]]<-t
 X<-model.matrix(fit);vc<-as.data.frame(VarCorr(fit)); dg<-data.table(scope=scope,spec=spec,variant=v,module_id=id,rank=qr(X)$rank,columns=ncol(X),condition_number=kappa(X),singular=isSingular(fit),random_variance=vc$vcov[1],residual_variance=sigma(fit)^2,shapiro_p=shapiro.test(residuals(fit))$p.value,convergence=paste(fit@optinfo$conv$lme4$messages,collapse=';'))
 diag[[k]]<-dg;d[,`:=`(scope=scope,spec=spec,variant=v,module_id=id)];longinput[[k]]<-d;longresid[[k]]<-data.table(scope=scope,spec=spec,variant=v,module_id=id,sample=d$sample,subject=d$subject,time=d$time,fitted=fitted(fit),residual=residuals(fit));fits[[k]]<-fit
}
lr<-rbindlist(lr);stopifnot(nrow(lr)==240)
lr[,q_BH_within_spec:=p.adjust(p_value,'BH',n=12),by=.(scope,spec,variant)];lr[,q_BH_all240:=p.adjust(p_value,'BH',n=240)]
cr[,q_BH_primary16:=NA_real_];lr[,q_BH_primary16:=NA_real_]
ci<-which(cr$subset=='clinical43' & cr$spec=='clinical' & cr$variant=='published');li<-which(lr$scope=='all37' & lr$spec=='disjoint_both' & lr$variant=='published');stopifnot(length(ci)==4,length(li)==12)
pq<-p.adjust(c(cr$p_value[ci],lr$p_value[li]),'BH',n=16);cr$q_BH_primary16[ci]<-pq[1:4];lr$q_BH_primary16[li]<-pq[5:16]
w(cr,'discovery_models');w(lr,'longitudinal_models');w(rbindlist(diag),'longitudinal_diagnostics');w(rbindlist(longinput,fill=TRUE),'longitudinal_model_inputs');w(rbindlist(longresid),'longitudinal_residuals')
saveRDS(fits,file.path(out,'logs/longitudinal_models.rds'))
ld<-sc[dataset=='GSE226238' & group=='SCI',.(n_observations=.N,n_participants=uniqueN(subject),mean=mean(score),sd=sd(score),median=median(score)),by=.(variant,module_id,time)];w(ld,'longitudinal_descriptives')
# Explicit stage changes, not a test of coefficient differences.
trans<-list(c('full48','unadjusted','clinical43','unadjusted'),c('clinical43','unadjusted','clinical43','age_sex'),c('clinical43','age_sex','clinical43','clinical'),c('clinical43','clinical','time41','clinical'),c('time41','clinical','time41','clinical_time'),c('clinical43','clinical','clinical43','clinical_CIB'),c('clinical43','clinical','clinical43','clinical_disjoint'))
changes<-rbindlist(lapply(trans,function(a)rbindlist(lapply(variants,function(v)rbindlist(lapply(mods,function(id){x<-cr[subset==a[1]&spec==a[2]&variant==v&module_id==id];y<-cr[subset==a[3]&spec==a[4]&variant==v&module_id==id];data.table(variant=v,module_id=id,from_set=a[1],from_spec=a[2],to_set=a[3],to_spec=a[4],before=x$estimate,after=y$estimate,change=y$estimate-x$estimate,CI_width_before=x$upper_CL-x$lower_CL,CI_width_after=y$upper_CL-y$lower_CL)}))))))
w(changes,'model_stage_changes')
w(data.table(role=names(inputs),path=unname(inputs),md5=unname(tools::md5sum(inputs))),'input_hashes_after')
capture.output(sessionInfo(),file=file.path(out,'logs/sessionInfo.txt'));cat('COMPLETE',nrow(cr),nrow(lr),'\n');sink()
