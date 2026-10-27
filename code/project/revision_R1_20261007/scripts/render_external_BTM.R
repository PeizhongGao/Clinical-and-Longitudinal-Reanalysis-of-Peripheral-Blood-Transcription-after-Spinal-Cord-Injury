# Public-package copy: see metadata/code_changes.csv for provenance.
# Figure contract: quantitative grid, all four locked BTM modules and contrasts.
# Claim: direction agreement varies in statistical support across discovery and follow-up.
# R-only; no observation exclusion. All exact statistics accompany the figure in CSV.
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(if(length(args))args[1] else '.',winslash='/');out<-file.path(root,'revision_R1_20261007')
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
d<-as.data.frame(fread(file.path(out,'tables/external_BTM_effects.csv')))
d<-d[d$variant=='published' & d$scope=='all_available',];stopifnot(nrow(d)==16,all(is.finite(d$estimate)))
d$module<-factor(d$module_id,levels=c('M223','M7.0','M163','M37.1'),labels=c('T cell II (M223)','T cell I (M7.0)','Neutrophil II (M163)','Neutrophil I (M37.1)'))
d$passed<-factor(ifelse(d$q_BH<.05,'BH q < 0.05','BH q >= 0.05'),levels=c('BH q < 0.05','BH q >= 0.05'))
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(legend.position='bottom',plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=8),axis.title=element_text(size=8),axis.text=element_text(size=7),plot.tag=element_text(size=10,face='bold')))
panels<-list()
for(con in c('SCI - TC','3m - acute','6m - acute','12m - acute')){
 dd<-d[d$contrast==con,]
 panels[[con]]<-ggplot(dd,aes(x=estimate,y=module,color=lineage,shape=passed))+
  geom_vline(xintercept=0,linetype='dashed',color='grey65',linewidth=.35)+
  geom_segment(aes(x=lower.CL,xend=upper.CL,yend=module),linewidth=.5,show.legend=FALSE)+geom_point(size=2.1,stroke=.7,show.legend=TRUE)+
  scale_color_manual(values=c(neutrophil='#B9564F',T_cell='#518AA0'),guide='none')+scale_shape_manual(values=c(16,1),name=NULL,drop=FALSE)+
  labs(title=con,subtitle=if(con=='SCI - TC')'38 SCI and 10 trauma controls' else '10 SCI participants; 37 observations',x=if(con=='SCI - TC')'Difference in mean log expression' else 'Change in mean log2 CPM',y=NULL)
}
p<-wrap_plots(panels,ncol=2,guides='collect')+plot_annotation(tag_levels='a',caption='External definitions; same study cohorts. Points: estimates; lines: pointwise 95% CI.\nBH correction over all 16 contrasts. No clinical covariate adjustment.')
p<-p&theme(plot.caption=element_text(size=7,hjust=0),legend.position='bottom')
stem<-file.path(out,'figures/Supp_external_BTM');width_mm<-183;height_mm<-143
grDevices::cairo_pdf(paste0(stem,'.pdf'),width=width_mm/25.4,height=height_mm/25.4,family='Arial');print(p);dev.off()
svglite::svglite(paste0(stem,'.svg'),width=width_mm/25.4,height=height_mm/25.4);print(p);dev.off()
ggsave(paste0(stem,'.png'),p,width=width_mm/25.4,height=height_mm/25.4,dpi=300,device=ragg::agg_png)
ggsave(paste0(stem,'.tiff'),p,width=width_mm/25.4,height=height_mm/25.4,dpi=600,device=ragg::agg_tiff,compression='lzw')
fwrite(d,file.path(out,'tables/external_BTM_figure_source.csv'))
cat('Exported 16 contrasts at 183 x 143 mm, PDF/SVG/600dpi TIFF/PNG.\n')
