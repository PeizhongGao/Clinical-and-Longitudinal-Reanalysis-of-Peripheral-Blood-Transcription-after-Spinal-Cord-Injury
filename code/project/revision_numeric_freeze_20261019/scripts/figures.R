# Public-package copy: see metadata/code_changes.csv for provenance.
args<-commandArgs(trailingOnly=TRUE)
root<-normalizePath(if(length(args))args[1] else getwd(),winslash='/',mustWork=TRUE)
out<-file.path(root,'revision_numeric_freeze_20261019')
stopifnot(!file.exists(file.path(out,'FREEZE.json')))
suppressPackageStartupMessages({library(ggplot2);library(data.table);library(svglite);library(ragg)})
theme_set(theme_classic(base_size=7,base_family='Arial')+
 theme(axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),
       legend.title=element_text(size=7),legend.text=element_text(size=7),
       plot.title=element_text(size=8,face='bold'),strip.text=element_text(size=7),
       plot.subtitle=element_text(size=7),plot.caption=element_text(size=7,hjust=0),
       plot.margin=margin(8,8,8,8),legend.position='bottom'))
save_pub<-function(p,name,width_mm=183,height_mm=120){
 stem<-file.path(out,'figures',name);w<-width_mm/25.4;h<-height_mm/25.4
 svglite::svglite(paste0(stem,'.svg'),width=w,height=h);print(p);dev.off()
 grDevices::cairo_pdf(paste0(stem,'.pdf'),width=w,height=h,family='Arial');print(p);dev.off()
 ragg::agg_tiff(paste0(stem,'.tiff'),width=w,height=h,units='in',res=600,compression='lzw');print(p);dev.off()
 ragg::agg_png(paste0(stem,'.png'),width=w,height=h,units='in',res=300);print(p);dev.off()
}
g<-fread(file.path(out,'tables','bridge43_fixed_eight_GO.csv'))
stopifnot(nrow(g)==24L,uniqueN(g$go_id)==8L,all(is.finite(g$FDR)))
terms<-unique(g$term)
g[,term_label:=factor(term,levels=rev(terms))]
g[,model:=factor(model,levels=c('M0','M1','M2'))]
fwrite(g,file.path(out,'figures','figure01_source.csv'))
wrapped<-function(x)vapply(x,function(z)paste(strwrap(z,35),collapse='\n'),character(1))
p1<-ggplot(g,aes(FDR,term_label,colour=model,shape=Direction,group=model))+
 geom_vline(xintercept=.05,linetype='dashed',linewidth=.35,colour='#777777')+
 geom_point(size=2.5,position=position_dodge(width=.55),stroke=.6)+
 scale_colour_manual(values=c(M0='#6C757D',M1='#4C87A1',M2='#A85F39'),
   labels=c('M0: group','M1: age + sex','M2: clinical'))+
 scale_shape_manual(values=c(Up=2,Down=6))+
 scale_y_discrete(labels=wrapped)+scale_x_continuous(limits=c(0,1),breaks=c(0,.05,.25,.5,.75,1))+
 labs(title='Fixed eight immune GO terms across the same 43 participants',
      subtitle='36 SCI / 7 trauma controls; complete 4,032-set BH family in each model',
      x='BH q value',y=NULL,colour='Model',shape='CAMERA direction',
      caption='All 24 displayed q values exceed 0.05. Direction is relative to the gene background.\nComplete GO results are retained; this panel was fixed before the bridge fits.')
save_pub(p1,'figure01_fixed_immune_GO',height_mm=125)

b<-fread(file.path(out,'tables','BTM_primary16.csv'))
stopifnot(nrow(b)==16L,all(is.finite(b$estimate)),all(is.finite(b$lower_CL)),all(is.finite(b$upper_CL)))
b[source_branch=='discovery',comparison:='SCI - TC\nClinical, n = 43']
b[source_branch=='longitudinal',comparison:=paste0(contrast,'\nComposition adjusted')]
b[,comparison:=factor(comparison,levels=c('SCI - TC\nClinical, n = 43','3m - acute\nComposition adjusted','6m - acute\nComposition adjusted','12m - acute\nComposition adjusted'))]
b[,module_id:=factor(module_id,levels=c('M223','M7.0','M163','M37.1'))]
b[,BH16:=factor(ifelse(q_BH_primary16<.05,'q < 0.05','q >= 0.05'),levels=c('q < 0.05','q >= 0.05'))]
fwrite(b,file.path(out,'figures','figure02_source.csv'))
p2<-ggplot(b,aes(estimate,module_id,colour=BH16))+
 geom_vline(xintercept=0,linetype='dashed',colour='#777777',linewidth=.35)+
 geom_segment(aes(x=lower_CL,xend=upper_CL,y=module_id,yend=module_id),linewidth=.5)+
 geom_point(size=2)+facet_wrap(~comparison,nrow=1)+
 scale_colour_manual(values=c('q < 0.05'='#356F8A','q >= 0.05'='#777777'))+
 scale_x_continuous(breaks=c(-.5,0,.5),limits=c(-.72,.93))+
 labs(title='External BTM: clinical uncertainty and composition-adjusted longitudinal changes',
      x='Coefficient and pointwise 95% CI (source-specific log-expression units)',y=NULL,colour='BH16',
      caption='Clinical: OLS / HC3, 43 participants. Longitudinal: random intercept, 10 people / 37 observations.\nLongitudinal models include both disjoint expression proxies; all 16 tests share BH correction.')
save_pub(p2,'figure02_BTM_primary16',height_mm=100)
cat('R FIGURES PASS: 24 fixed GO rows and 16 BTM rows exported without omission.\n')
