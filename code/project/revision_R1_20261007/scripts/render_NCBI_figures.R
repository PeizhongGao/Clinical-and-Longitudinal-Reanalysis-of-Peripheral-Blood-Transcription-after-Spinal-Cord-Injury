# Public-package copy: see metadata/code_changes.csv for provenance.
# Figure contract: observed full-signature changes within SCI participants, not immune recovery.
# Quantitative grid. Structural adaptation of existing trajectory/forest design.
# All 37 SCI samples retained. Columns: sample/subject/time identify repeated observations;
# score is the control-referenced NCBI logCPM signature mean, estimate/CI are LMM contrasts.
# Main panels: trajectories and pointwise 95% CIs. Supplement: 12m source/coverage/adjustment effects.
# Export contract: 183 mm wide, Arial editable vector text; PDF/SVG/600dpi TIFF/300dpi PNG.
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(if(length(args))args[1] else '.',winslash='/');out<-file.path(root,'revision_R1_20261007')
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
state<-readRDS(file.path(out,'tables/NCBI_analysis_state.rds'));s<-state$sfull;mixed<-state$mixed
scores<-c('innate_score','adaptive_score','balance_score');labels<-c('Innate-associated','Lymphocyte/T-cell-associated','Balance (secondary)');visits<-c('acute','3m','6m','12m')
palette<-c('#B9564F','#518AA0','#756783')
theme_set(theme_classic(base_size=9,base_family='Arial')+theme(strip.text=element_text(face='bold'),plot.title=element_text(size=10,face='bold'),legend.position='bottom'))
save_plot<-function(p,name,w=183,h=160){
  stem<-file.path(out,'figures',name)
  grDevices::cairo_pdf(paste0(stem,'.pdf'),width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off()
  svglite::svglite(paste0(stem,'.svg'),width=w/25.4,height=h/25.4);print(p);dev.off()
  ggsave(paste0(stem,'.png'),p,width=w/25.4,height=h/25.4,dpi=300,device=ragg::agg_png)
  ggsave(paste0(stem,'.tiff'),p,width=w/25.4,height=h/25.4,dpi=600,device=ragg::agg_tiff,compression='lzw')
}
long<-melt(as.data.table(s[s$group=='SCI',]),id.vars=c('sample','subject','time'),measure.vars=scores,variable.name='score',value.name='value');long$time<-factor(long$time,levels=visits);long$score<-factor(long$score,levels=scores,labels=labels)
fwrite(long,file.path(out,'tables/NCBI_Figure5_trajectories.csv'))
p1<-ggplot(long,aes(time,value,group=subject))+geom_line(color='#A7ACB3',linewidth=.35,alpha=.8)+geom_point(aes(color=score),size=1.6)+facet_wrap(~score,scales='free_y',nrow=1)+scale_color_manual(values=palette)+labs(x='Visit: acute n=10; 3m n=10; 6m n=10; 12m n=7',y='Control-referenced signature score',title='Within-participant trajectories')+theme(legend.position='none')
forest<-mixed[mixed$source=='NCBI_full_15plus15' & mixed$scope=='all_available' & mixed$spec=='unadjusted',];forest$score<-factor(forest$score,levels=scores,labels=labels);forest$contrast<-factor(forest$contrast,levels=rev(c('3m - acute','6m - acute','12m - acute')))
fwrite(forest,file.path(out,'tables/NCBI_Figure5_contrasts.csv'))
p2<-ggplot(forest,aes(estimate,contrast,color=score))+geom_vline(xintercept=0,linetype=2,color='#888888')+geom_errorbar(aes(xmin=lower.CL,xmax=upper.CL),orientation='y',width=.1)+geom_point(size=2)+facet_wrap(~score,scales='free_x',nrow=1)+scale_color_manual(values=palette)+labs(x='Later minus acute estimate (pointwise 95% CI)',y=NULL,title='Random-intercept model; 10 participants / 37 samples')+theme(legend.position='none')
save_plot(p1/p2+plot_annotation(tag_levels='a',subtitle='NCBI-generated counts: full frozen coverage, 15 innate + 15 lymphocyte/T-cell genes'),'Figure5_revised')
z<-mixed[mixed$contrast=='12m - acute',];z$analysis<-ifelse(z$source=='NCBI_restricted_8plus13','NCBI: 8+13',ifelse(z$scope=='complete_cases','NCBI: 15+15, complete cases',ifelse(z$spec=='marker_adjusted','NCBI: 15+15, marker adjusted','NCBI: 15+15')))
z<-z[,c('score','estimate','lower.CL','upper.CL','analysis')]
old<-as.data.frame(fread(file.path(out,'tables/longitudinal_mixed_contrasts.csv')));old<-old[old$scope=='all_available' & old$contrast=='12m - acute',c('score','estimate','lower.CL','upper.CL')];old$analysis<-'Author matrix: 8+13'
z<-rbind(z,old);z$analysis<-factor(z$analysis,levels=rev(c('Author matrix: 8+13','NCBI: 8+13','NCBI: 15+15','NCBI: 15+15, complete cases','NCBI: 15+15, marker adjusted')));z$score<-factor(z$score,levels=scores,labels=labels)
fwrite(z,file.path(out,'tables/NCBI_source_sensitivity_figure.csv'))
p3<-ggplot(z,aes(estimate,analysis,color=score))+geom_vline(xintercept=0,linetype=2,color='#888888')+geom_errorbar(aes(xmin=lower.CL,xmax=upper.CL),orientation='y',width=.12)+geom_point(size=2)+facet_wrap(~score,scales='free_x',nrow=1)+scale_color_manual(values=palette)+labs(x='12 months minus acute estimate (pointwise 95% CI)',y=NULL,title='Source, gene coverage and sensitivity comparisons',subtitle='Same SCI participants; conditional on each source and control-reference transformation')+theme(legend.position='none',axis.text.y=element_text(size=7),strip.text=element_text(size=8))
save_plot(p3,'Supp_NCBI_source_sensitivity',183,100)
