# Public-package copy: see metadata/code_changes.csv for provenance.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
o<-file.path(commandArgs(TRUE)[1],'revision_external_adjusted_20261016');d<-fread(file.path(o,'tables/discovery_models.csv'));l<-fread(file.path(o,'tables/longitudinal_models.csv'))
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(strip.background=element_blank(),strip.text=element_text(face='bold',size=9),axis.title=element_text(size=8),axis.text=element_text(size=7),legend.position='bottom',legend.title=element_blank(),plot.caption=element_text(hjust=0,size=7),plot.title=element_text(face='bold',size=10),panel.spacing=grid::unit(5,'mm')))
save_pub<-function(p,name,height_mm){
 path<-file.path(o,'figures',name);width_mm<-183
 svglite::svglite(paste0(path,'.svg'),width=width_mm/25.4,height=height_mm/25.4);print(p);dev.off()
 cairo_pdf(paste0(path,'.pdf'),width=width_mm/25.4,height=height_mm/25.4,family='Arial');print(p);dev.off()
 ragg::agg_tiff(paste0(path,'.tiff'),width=width_mm,height=height_mm,units='mm',res=600,compression='lzw');print(p);dev.off()
 ragg::agg_png(paste0(path,'.png'),width=width_mm,height=height_mm,units='mm',res=300);print(p);dev.off()
}
mods<-c('M37.1','M163','M7.0','M223');modlabs<-c('M37.1  Neutrophil-enriched','M163  Neutrophil-enriched','M7.0  T-cell-enriched','M223  T-cell-enriched')
order<-c('full48__unadjusted','full48__age_sex','clinical43__unadjusted','clinical43__age_sex','clinical43__clinical','clinical43__clinical_CIB','clinical43__clinical_disjoint','time41__unadjusted','time41__age_sex','time41__clinical','time41__clinical_time')
labels<-c('48  Unadjusted','48  Age + sex','43  Unadjusted','43  Age + sex','43  Clinical','43  Clinical + CIBERSORTx','43  Clinical + disjoint markers','41  Unadjusted','41  Age + sex','41  Clinical','41  Clinical + draw time')
x<-d[variant=='published'];x[,stage:=factor(paste(subset,spec,sep='__'),levels=rev(order),labels=rev(labels))];x[,module_id:=factor(module_id,levels=mods,labels=modlabs)]
x[,set:=factor(subset,levels=c('full48','clinical43','time41'),labels=c('48 participants','43 participants','41 participants'))]
fwrite(x,file.path(o,'tables/figure01_source.csv'))
p<-ggplot(x,aes(estimate,stage,colour=set))+geom_vline(xintercept=0,linetype='dashed',linewidth=.3,colour='gray50')+geom_errorbar(aes(xmin=lower_CL,xmax=upper_CL),orientation='y',width=.22,linewidth=.4)+geom_point(size=1.6)+facet_wrap(~module_id,ncol=2,scales='free_x')+scale_colour_manual(values=c('#566B82','#A27A42','#777777'))+labs(x='SCI minus trauma control (log-expression units)',y=NULL,title='External modules under fixed-participant clinical models',caption='Points: group coefficients; bars: pointwise 95% HC3 intervals.\n48 = 38 SCI / 10 TC; 43 = 36 / 7; 41 = 36 / 5. Published module definitions.')
save_pub(p,'figure01_external_fixed_sample',165)
y<-l[variant=='published' & scope=='all37' & spec%in%c('baseline','official_both','disjoint_both')];y[,module_id:=factor(module_id,levels=mods,labels=modlabs)];y[,visit:=factor(contrast,levels=c('3m - acute','6m - acute','12m - acute'),labels=c('3 months','6 months','12 months'))];y[,adjustment:=factor(spec,levels=c('baseline','official_both','disjoint_both'),labels=c('Time only','Official MCP-counter','Disjoint marker means'))]
fwrite(y,file.path(o,'tables/figure02_source.csv'))
pd<-position_dodge(width=.5)
p2<-ggplot(y,aes(visit,estimate,colour=adjustment,shape=adjustment))+geom_hline(yintercept=0,linetype='dashed',linewidth=.3,colour='gray50')+geom_errorbar(aes(ymin=lower_CL,ymax=upper_CL),width=.14,position=pd,linewidth=.4)+geom_point(position=pd,size=1.7)+facet_wrap(~module_id,ncol=2,scales='free_y')+scale_colour_manual(values=c('#777777','#A27A42','#566B82'))+scale_shape_manual(values=c(16,17,15))+labs(x=NULL,y='Later minus acute (log2 CPM units)',title='Longitudinal module changes depend on composition proxies',caption='10 SCI participants; 37 observations (10 / 10 / 10 / 7 by visit).\nRandom-intercept models; pointwise 95% intervals. Proxy models include neutrophil and T-cell terms.')
save_pub(p2,'figure02_external_composition',135)
capture.output(sessionInfo(),file=file.path(o,'logs/figure_sessionInfo.txt'))
