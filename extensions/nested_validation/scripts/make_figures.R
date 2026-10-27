args<-commandArgs(trailingOnly=TRUE);ROOT<-normalizePath(args[1],winslash='/',mustWork=TRUE)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(svglite);library(ragg)})
dir.create(file.path(ROOT,'figures'),recursive=TRUE,showWarnings=FALSE)
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7),strip.text=element_text(size=8,face='bold'),legend.position='bottom',plot.tag=element_text(size=10,face='bold'),plot.caption=element_text(size=7,hjust=0),axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.35)))
save_plot<-function(p,stem,width_mm=183,height_mm=155){
  base<-file.path(ROOT,'figures',stem);w<-width_mm/25.4;h<-height_mm/25.4
  svglite(paste0(base,'.svg'),width=w,height=h);print(p);dev.off()
  cairo_pdf(paste0(base,'.pdf'),width=w,height=h,family='Arial');print(p);dev.off()
  agg_tiff(paste0(base,'.tiff'),width=w,height=h,units='in',res=600,compression='lzw');print(p);dev.off()
  agg_png(paste0(base,'.png'),width=w,height=h,units='in',res=300);print(p);dev.off()
}
if(length(args)<2 || args[2]!='clinical_only') {
m<-fread(file.path(ROOT,'tables/nested_repeat_metrics.csv'));nul<-fread(file.path(ROOT,'tables/full_pipeline_permutation_null.csv'));res<-fread(file.path(ROOT,'tables/full_pipeline_permutation_results.csv'))
m[,comparison:=factor(contrast,levels=c('SCI_HC','SCI_TC'),labels=c('SCI vs healthy control','SCI vs trauma control'))]
m[,Source:=factor(source,levels=c('raw','processed'),labels=c('Raw counts','Archived processed'))]
meanm<-m[,.(AUC=mean(AUC)),by=.(comparison,Source)]
p1<-ggplot(m,aes(Source,AUC,color=Source))+geom_hline(yintercept=.5,linetype='dashed',color='grey65',linewidth=.4)+geom_point(position=position_jitter(width=.11,height=0,seed=42),alpha=.7,size=1.5)+geom_point(data=meanm,shape=95,size=8,color='black')+facet_wrap(~comparison,nrow=1)+scale_color_manual(values=c('#326889','#8a8a8a'))+coord_cartesian(ylim=c(0,1))+labs(x=NULL,y='Nested cross-validation AUC',title='Internal performance across 20 split repeats',subtitle='38 SCI and 10 controls per comparison; black ticks show means')+theme(legend.position='none')
null_plot<-function(cn,title){
  r<-res[contrast==cn]
  ggplot(nul[contrast==cn],aes(AUC))+geom_histogram(binwidth=.04,boundary=0,fill='#9cb9cb',color='white',linewidth=.2)+geom_vline(xintercept=r$observed_first_repeat_AUC,color='#1c3e54',linewidth=.8)+scale_x_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1))+labs(x='AUC after complete-pipeline label permutation',y='Number of permutations',title=title,subtitle=sprintf('Observed %.3f; P = %.3f; Holm2 = %.3f',r$observed_first_repeat_AUC,r$p_permutation,r$p_Holm2))
}
p2<-null_plot('SCI_HC','SCI vs healthy control');p3<-null_plot('SCI_TC','SCI vs trauma control')
fig10<-(p1/(p2+p3))+plot_layout(heights=c(1,1))+plot_annotation(tag_levels='a',caption='999 permutations per comparison; labels shuffled within sex. The observed statistic uses the same\nsingle nested-CV schedule as each permutation. Internal association does not establish clinical validity.')
save_plot(fig10,'Figure_S10_Nested_LASSO_and_full_pipeline_permutation')
}
c<-fread(file.path(ROOT,'tables/clinical_resampling_results.csv'));c[,Module:=factor(module_id,levels=rev(c('M37.1','M163','M7.0','M223')))]
eff<-c[method=='Freedman_Lane_within_sex']
q1<-ggplot(eff,aes(estimate,Module))+geom_vline(xintercept=0,linetype='dashed',color='grey55',linewidth=.4)+geom_errorbar(aes(xmin=lower_CL,xmax=upper_CL),orientation='y',width=.15,linewidth=.6,color='#326889')+geom_point(size=2,color='#326889')+labs(x='SCI minus trauma control estimate\n(pointwise 95% HC3 CI)',y=NULL,title='Frozen clinical effect estimates')
c[,Method:=factor(method,levels=c('Freedman_Lane_within_sex','null_wild_Rademacher_HC2'),labels=c('Residual permutation','Wild bootstrap'))]
c[,display_y:=as.numeric(Module)+ifelse(method=='Freedman_Lane_within_sex',-.12,.12)]
q2<-ggplot(c,aes(p_resampling,display_y,color=Method,shape=Method))+geom_point(size=2.2)+scale_y_continuous(breaks=1:4,labels=levels(c$Module),limits=c(.5,4.5))+scale_color_manual(values=c('#326889','#888888'))+scale_shape_manual(values=c(16,17))+scale_x_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1))+labs(x='Two-sided resampling P value',y=NULL,title='Alternative null resampling',subtitle=sprintf('Holm4-adjusted P range: %.3f to %.3f',min(c$p_Holm4),max(c$p_Holm4)))
fig11<-(q1+q2)+plot_layout(widths=c(1,1.15),guides='collect')+plot_annotation(tag_levels='a',caption='n = 43 independent participants (36 SCI / 7 trauma controls). Each method uses 9,999 resamples.\nThe model adjusts for age, sex, ISS, TBI and prior CNS pathology. Results do not establish equivalence.')
save_plot(fig11,'Figure_S11_Clinical_BTM_resampling',height_mm=125)
cat('FIGURES COMPLETE\n')
