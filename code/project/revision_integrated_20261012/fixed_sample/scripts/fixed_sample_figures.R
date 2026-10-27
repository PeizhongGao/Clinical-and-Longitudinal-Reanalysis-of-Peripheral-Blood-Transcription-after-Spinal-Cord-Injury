# Public-package copy: see metadata/code_changes.csv for provenance.
# All panels use the existing R backend and real source tables; no random sampling.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1],winslash='/',mustWork=TRUE) else normalizePath(getwd(),winslash='/',mustWork=TRUE)
out <- file.path(root,'revision_R2_20261008','fixed_sample_20261012')
tables <- file.path(out,'tables'); figures <- file.path(out,'figures')
read_tab <- function(n) fread(file.path(tables,paste0(n,'.csv')),na.strings=c('NA',''))
write_source <- function(x,n) fwrite(x,file.path(tables,paste0(n,'.csv')),na='NA')
theme_set(theme_classic(base_size=9,base_family='Arial')+
  theme(axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.3),
        strip.background=element_blank(),strip.text=element_text(size=9,face='bold'),
        plot.title=element_text(size=10,face='bold'),plot.subtitle=element_text(size=8),
        plot.caption=element_text(size=8,hjust=0),legend.text=element_text(size=8),legend.title=element_text(size=8),
        plot.margin=margin(6,8,5,5),panel.spacing=grid::unit(8,'pt')))
pal <- c(full48='#767676',clinical43='#367B89',time41='#B07B45')
groups <- c(TC='#367B89',SCI='#B07B45')
exports <- list()
save_pub <- function(plot,name,width_mm=183,height_mm=125,tiff=TRUE){
  path <- file.path(figures,name); w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(path,'.svg'),width=w,height=h,bg='white'); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(path,'.pdf'),width=w,height=h,family='Arial',bg='white'); print(plot); dev.off()
  ragg::agg_png(paste0(path,'.png'),width=w,height=h,units='in',res=300,background='white'); print(plot); dev.off()
  if(tiff){ragg::agg_tiff(paste0(path,'.tiff'),width=w,height=h,units='in',res=600,compression='lzw',background='white'); print(plot); dev.off()}
  exports[[name]] <<- data.table(figure=name,width_mm=width_mm,height_mm=height_mm,PNG_dpi=300,TIFF_dpi=if(tiff)600 else NA_integer_,
                                backend='R',vector_formats='PDF/SVG',font='Arial',minimum_text_pt=8)
}
models <- read_tab('fixed_sample_models')
order_keys <- c('full48__unadjusted','full48__age_sex','clinical43__unadjusted','clinical43__age_sex','clinical43__clinical',
                'time41__unadjusted','time41__age_sex','time41__clinical','time41__clinical_time')
labels <- c('48: unadjusted (38/10)','48: age + sex (38/10)','43: unadjusted (36/7)','43: age + sex (36/7)','43: clinical (36/7)',
            '41: unadjusted (36/5)','41: age + sex (36/5)','41: clinical (36/5)','41: clinical + time (36/5)')
models[,key:=paste(subset,spec,sep='__')]; models[,y:=match(key,rev(order_keys))]
models[,row_label:=labels[match(key,order_keys)]]
write_source(models,'figure01_source_all_scores')
panel_score <- function(outcome_name,title,show_y=TRUE){
  z <- models[outcome==outcome_name]
  stopifnot(nrow(z)==9L)
  p <- ggplot(z,aes(x=estimate,y=y,colour=subset,shape=subset))+
    geom_vline(xintercept=0,linewidth=.35,linetype='dashed',colour='grey60')+
    geom_hline(yintercept=c(4.5,7.5),colour='grey90',linewidth=.3)+
    geom_segment(aes(x=lower_CL,xend=upper_CL,yend=y),linewidth=.65)+geom_point(size=2.3)+
    scale_y_continuous(breaks=1:9,labels=rev(labels),limits=c(.5,9.5))+
    scale_colour_manual(values=pal,guide='none')+scale_shape_manual(values=c(full48=16,clinical43=17,time41=15),guide='none')+
    labs(x='SCI minus TC (frozen score units)',y=NULL,title=title)
  if(!show_y)p <- p+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  p
}
p1 <- panel_score('innate_score','Innate-associated score')+panel_score('adaptive_score','Lymphocyte-associated score',FALSE)+
  plot_layout(widths=c(1.5,1))+
  plot_annotation(title='Fixed-sample models: effect estimates and precision',
    subtitle='48-person anchor, 43-person clinical set and 41-person time set; counts = SCI/TC',
    caption='Points: OLS group coefficients. Bars: pointwise HC3 95% CIs (residual-df t reference).\nClinical = age, sex, ISS, TBI and prior CNS pathology. Post-selection revision diagnostics.',tag_levels='a')
save_pub(p1,'figure01_fixed_sample_effects',height_mm=145)
pb <- panel_score('balance_score','Secondary balance = innate minus lymphocyte')+
  labs(subtitle='Derived secondary outcome; same frozen score scale',caption='Points: OLS group coefficients; bars: pointwise HC3 95% CIs. Counts = SCI/TC.')
save_pub(pb,'figureS01_balance_effects',height_mm=125)

loo <- read_tab('TC_leave_one_out_models')
loo[,panel:=factor(paste(subset,spec,sep='__'),levels=c('clinical43__clinical','time41__clinical','time41__clinical_time'),
  labels=c('43: clinical (7 TC)','41: clinical (5 TC)','41: clinical + time (5 TC)'))]
loo[,score:=factor(outcome,levels=c('innate_score','adaptive_score','balance_score'),
                  labels=c('Innate-associated','Lymphocyte-associated','Secondary balance'))]
loo[,deleted_label:=paste0('Remove ',deleted_TC)]
write_source(loo,'figure02_source_TC_deletions')
refline <- unique(loo[,.(panel,score,full_estimate)])
p2 <- ggplot(loo[outcome!='balance_score'],aes(y=deleted_label,x=estimate))+
  geom_vline(xintercept=0,linetype='dashed',colour='grey65',linewidth=.3)+
  geom_vline(data=refline[score!='Secondary balance'],aes(xintercept=full_estimate),colour='#B07B45',linewidth=.55)+
  geom_segment(aes(x=lower_CL,xend=upper_CL,yend=deleted_label),colour='#367B89',linewidth=.55)+
  geom_point(colour='#367B89',size=1.9)+
  facet_grid(panel~score,scales='free_y',space='free_y')+
  labs(x='SCI minus TC after removing one TC',y=NULL,title='Influence of every trauma control',
       subtitle='Each deletion is a diagnostic refit; all participants remain in the main models',
       caption='Orange line: full-set estimate. Points and bars: deletion estimate and HC3 95% CI.\nPanel headings give the original TC count; each refit contains one fewer TC. No favorable deletion is selected.')+
  theme(strip.text.y=element_text(angle=0,size=8),axis.text.y=element_text(size=8))
save_pub(p2,'figure02_TC_leave_one_out',height_mm=190)

support <- rbindlist(lapply(c('clinical43','time41'),function(s){x<-read_tab(paste0('fixed_input_',s));x[,subset:=s];x}))
support <- melt(support,id.vars=c('sample','group','subset'),measure.vars=c('age','iss','draw_hours'),variable.name='variable',value.name='value')
support[,available:=is.finite(value)]
write_source(support,'figure03_source_covariate_support')
support <- support[available==TRUE]
support[,variable:=factor(variable,levels=c('age','iss','draw_hours'),labels=c('Age (years)','ISS','Blood draw (hours)'))]
support[,subset:=factor(subset,levels=c('clinical43','time41'),labels=c('43 participants (36 SCI / 7 TC)','41 participants (36 SCI / 5 TC)'))]
support[,group:=factor(group,levels=c('TC','SCI'))]
setorder(support,subset,variable,group,sample)
support[,offset:=if(.N==1)0 else seq(-.10,.10,length.out=.N),by=.(subset,variable,group)]
support[,y:=as.numeric(group)+offset]
overlap <- read_tab('covariate_range_overlap')[subset%in%c('clinical43','time41')]
overlap[,variable:=factor(variable,levels=c('age','iss','draw_hours'),labels=c('Age (years)','ISS','Blood draw (hours)'))]
overlap[,subset:=factor(subset,levels=c('clinical43','time41'),labels=c('43 participants (36 SCI / 7 TC)','41 participants (36 SCI / 5 TC)'))]
support_summary <- support[,.(median=as.numeric(median(value)),Q1=as.numeric(quantile(value,.25)),Q3=as.numeric(quantile(value,.75)),y=as.numeric(group[1]),n_available=.N),by=.(subset,variable,group)]
p3 <- ggplot(support,aes(x=value,y=y,colour=group,shape=group))+
  geom_rect(data=overlap,aes(xmin=overlap_low,xmax=overlap_high,ymin=.65,ymax=2.35),inherit.aes=FALSE,fill='grey94',colour=NA)+
  geom_segment(data=support_summary,aes(x=Q1,xend=Q3,y=y,yend=y),linewidth=1.3,alpha=.45)+
  geom_point(size=1.8,alpha=.9)+
  geom_point(data=support_summary,aes(x=median,y=y),shape=3,size=3,stroke=.8)+
  facet_grid(subset~variable,scales='free_x')+
  scale_y_continuous(breaks=1:2,labels=c('TC','SCI'),limits=c(.65,2.35))+
  scale_colour_manual(values=groups,guide='none')+scale_shape_manual(values=c(TC=16,SCI=17),guide='none')+
  labs(x='Observed covariate value',y=NULL,title='Available covariate support in the fixed sets',
    subtitle='Every available observation; no support-based exclusions',
    caption='Shading: intersection of observed group ranges, not multivariable common support. Cross: median; thick line: IQR.\nIn the 43-person set, blood-draw time is missing for 2 TC. Full availability and category counts accompany the figure.')+
  theme(strip.text.y=element_text(angle=0,size=8),axis.text.x=element_text(size=8))
save_pub(p3,'figure03_covariate_support',height_mm=125)

influence <- read_tab('all_participant_influence')[spec%in%c('clinical','clinical_time') & outcome!='balance_score']
influence[,panel:=factor(paste(subset,spec,sep='__'),levels=c('clinical43__clinical','time41__clinical','time41__clinical_time'),
 labels=c('43: clinical','41: clinical','41: clinical + time'))]
influence[,score:=factor(outcome,levels=c('innate_score','adaptive_score'),labels=c('Innate-associated','Lymphocyte-associated'))]
write_source(influence,'figureS02_source_residual_diagnostics')
pr <- ggplot(influence,aes(x=fitted,y=residual,colour=group,shape=group))+
  geom_hline(yintercept=0,colour='grey65',linetype='dashed',linewidth=.3)+geom_point(size=1.5)+
  facet_grid(panel~score)+scale_colour_manual(values=groups)+scale_shape_manual(values=c(TC=16,SCI=17))+
  labs(x='Fitted score',y='Residual',title='Residuals versus fitted values',colour=NULL,shape=NULL)+theme(legend.position='bottom',strip.text.y=element_text(angle=0,size=8))
pq <- ggplot(influence,aes(sample=residual))+
  stat_qq(size=1.3,colour='#367B89')+stat_qq_line(linewidth=.4,colour='grey45')+
  facet_grid(panel~score)+labs(x='Theoretical normal quantile',y='Residual quantile',title='Residual normal QQ diagnostics')+
  theme(strip.text.y=element_text(angle=0,size=8))
pd <- pr/pq+plot_annotation(title='Fuller-model assumption diagnostics',
 caption='All 43 or 41 participants in each panel. These plots diagnose model limitations;\nHC3 does not guarantee finite-sample calibration, linearity or adequate covariate support.')
save_pub(pd,'figureS02_residual_diagnostics',height_mm=250,tiff=FALSE)
write_source(rbindlist(exports),'figure_export_manifest')
capture.output(sessionInfo(),file=file.path(out,'logs','figure_R_sessionInfo.txt'))
write_json(list(status='RENDERED_PENDING_VISUAL_QA',backend='R',n_figures=length(exports),
                random_sampling=FALSE,full_score_source_rows=nrow(models),TC_deletion_source_rows=nrow(loo)),
           file.path(out,'logs','figure_render_QA.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Figures rendered:',length(exports),'\n')
