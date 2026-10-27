# Revised Figures 1-7. R-only rendering, 170 mm width, no biological model refits.
if(.Platform$OS.type=='windows') invisible(suppressWarnings(Sys.setlocale('LC_CTYPE','English_United States.utf8')))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
project<-normalizePath(if(length(args))args[1] else '.',winslash='/',mustWork=TRUE)
out<-project;src<-file.path(out,'source_data')
for(d in c('previews','qa','svg'))dir.create(file.path(out,d),showWarnings=FALSE,recursive=TRUE)
# Optional second argument renders only specified numbers (e.g. 1,5,7).
selected<-if(length(args)>1)as.integer(strsplit(args[2],',',fixed=TRUE)[[1]]) else 1:7
width_mm<-170
W<-width_mm
cohort_col<-c(HC='#64727C',TC='#367B89',SCI='#B9564F',CTL='#64727C')
score_keys<-c('innate_score','adaptive_score','balance_score')
score_labels<-c('Myeloid /\ninnate-associated','Lymphocyte /\nT-cell-associated','Balance\n(secondary)')
score_col<-setNames(c('#B9564F','#367B89','#756783'),score_labels)
group_shape<-c(HC=16,TC=17,SCI=15,CTL=16)
wrap<-function(x,n=105)paste(strwrap(x,n),collapse='\n')
base_theme<-theme_classic(base_size=8,base_family='Arial')+
 theme(axis.line=element_line(linewidth=.35,colour='#31373B'),
       axis.ticks=element_line(linewidth=.3,colour='#31373B'),
       axis.text=element_text(size=7.5,colour='#3B4145'),axis.title=element_text(size=8),
       plot.title=element_text(size=9,face='bold',colour='#202428'),
       plot.subtitle=element_text(size=7.5,colour='#515B63'),
       plot.caption=element_text(size=7,hjust=0,colour='#515B63'),
       legend.title=element_text(size=8),legend.text=element_text(size=7.5),
       strip.background=element_blank(),strip.text=element_text(size=8,face='bold'),
       plot.margin=margin(6,6,6,6),panel.spacing=grid::unit(8,'pt'),
       plot.background=element_rect(fill='white',colour=NA),
       legend.key.width=grid::unit(10,'pt'),legend.key.height=grid::unit(10,'pt'))
theme_set(base_theme)
tag_theme<-theme(plot.tag=element_text(size=10,face='bold'),plot.tag.position='topleft')
annotation_theme<-theme(plot.title=element_text(size=10,face='bold'),
 plot.subtitle=element_text(size=8),plot.caption=element_text(size=7,hjust=0),
 plot.margin=margin(8,8,8,8))
read_src<-function(n)fread(file.path(src,n))
exports<-list()
save_pub<-function(p,number,name,height_mm){
 if(!number %in% selected)return(invisible(NULL))
 w<-W/25.4;h<-height_mm/25.4;stem<-file.path(out,name)
 grDevices::cairo_pdf(paste0(stem,'.pdf'),width=w,height=h,family='Arial',bg='white');print(p);dev.off()
 svglite::svglite(file.path(out,'svg',paste0(name,'.svg')),width=w,height=h,bg='white');print(p);dev.off()
 ragg::agg_tiff(paste0(stem,'.tif'),width=w,height=h,units='in',res=600,compression='lzw',background='white');print(p);dev.off()
 ragg::agg_png(file.path(out,'previews',paste0(name,'.png')),width=w,height=h,units='in',res=300,background='white');print(p);dev.off()
 exports[[name]]<<-data.table(number=number,figure=name,width_mm=W,height_mm=height_mm,tiff_dpi=600,preview_dpi=300,
  font='Arial',minimum_designed_text_pt=7,backend='R',pdf='editable vector',tiff='RGB / LZW')
 cat('EXPORTED Figure',number,':',W,'x',height_mm,'mm\n')
}
sfactor<-function(x)factor(x,levels=score_keys,labels=score_labels)

# 1. Historical discovery: all participants, all genes, explicit thresholds.
pca<-read_src('Figure1_PCA.csv');pv<-read_src('Figure1_PCA_variance.csv')
pca[,group:=factor(group,levels=c('HC','TC','SCI'))]
p1a<-ggplot(pca,aes(PC1,PC2,colour=group,shape=group))+
 geom_point(size=1.9,alpha=.85)+scale_colour_manual(values=cohort_col)+scale_shape_manual(values=group_shape)+
 labs(title='Discovery transcriptomes',x=sprintf('PC1 (%.1f%%)',100*pv[component=='PC1',variance_fraction]),
      y=sprintf('PC2 (%.1f%%)',100*pv[component=='PC2',variance_fraction]),colour=NULL,shape=NULL)+
 theme(legend.position='bottom')
vol<-read_src('Figure1_SCI_vs_HC.csv');vol[,direction:=factor(direction,levels=c('Down','Other','Up'))]
p1b<-ggplot(vol,aes(logFC,-log10(pmax(adj.P.Val,1e-300)),colour=direction))+
 geom_point(size=.45,alpha=.65)+geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.3,colour='#68737B')+
 geom_vline(xintercept=c(-1,1),linetype=2,linewidth=.3,colour='#68737B')+
 scale_colour_manual(values=c(Down='#367B89',Other='#CDD1D4',Up='#B9564F'))+
 labs(title='SCI − HC differential expression',x='log2 fold change (SCI − HC)',y='−log10 BH q',colour=NULL)+theme(legend.position='bottom')
dc<-read_src('Figure1_DEG_counts.csv')
dc[,category:=factor(category,levels=c('Up','Down','Other'))]
sc<-dc[comparison=='SCI_vs_HC']
sc[,category:=factor(category,levels=c('Other','Down','Up'))]
p1c<-ggplot(sc,aes(n,category,fill=category))+geom_col(width=.6)+
 geom_text(aes(x=18000,label=sprintf('%s (%.1f%%)',format(n,big.mark=','),percent)),
           hjust=1,size=2.7,family='Arial',colour='#30363A')+
 scale_fill_manual(values=c(Up='#B9564F',Down='#367B89',Other='#D5D9DC'))+
 scale_x_continuous(limits=c(0,18500),breaks=c(0,5000,10000,15000),labels=c('0','5,000','10,000','15,000'),expand=expansion(mult=c(0,0)))+
 labs(title='SCI − HC threshold summary',x='Number of genes',y=NULL)+theme(legend.position='none')
d<-dc[category!='Other'];d[,comparison:=factor(comparison,levels=c('TC_vs_HC','SCI_vs_HC','SCI_vs_TC'),labels=c('TC − HC','SCI − HC','SCI − TC'))]
p1d<-ggplot(d,aes(comparison,n,fill=category))+geom_col(position=position_dodge(.72),width=.62)+
 geom_text(aes(label=n),position=position_dodge(.72),vjust=-.35,size=2.7,family='Arial')+
 scale_fill_manual(values=c(Up='#B9564F',Down='#367B89'))+
 scale_y_continuous(expand=expansion(mult=c(0,.13)),labels=scales::label_comma())+
 labs(title='Directed pairwise DEG counts',x=NULL,y='Genes passing both thresholds',fill=NULL)+theme(legend.position='bottom')
fig1<-((p1a|p1b)/(p1c|p1d))+plot_layout(heights=c(1.1,1))+
 plot_annotation(title='Discovery-cohort transcriptional comparisons',tag_levels='a',
 subtitle='HC n = 10  |  TC n = 10  |  SCI n = 38',
 caption=paste('Historical processed-expression analysis; BH q < 0.05 and |log2FC| > 1 define Up / Down.',
 'Other = not passing both thresholds. PCA and DEG counts are descriptive, not evidence beyond trauma.',sep='\n'),
 theme=annotation_theme) & tag_theme
save_pub(fig1,1,'Figure1_Discovery_transcriptomic_comparisons',173)

# 2. All eight fixed GO terms, all three models, same 43 participants.
go<-read_src('Figure2_fixed_eight_GO.csv')
go[,term_label:=factor(term,levels=rev(unique(term)))]
go[,model:=factor(model,levels=c('M0','M1','M2'))]
p2<-ggplot(go,aes(FDR,term_label,colour=model,shape=Direction,group=model))+
 geom_vline(xintercept=.05,linetype=2,linewidth=.35,colour='#68737B')+
 geom_point(size=2.4,position=position_dodge(width=.58),stroke=.65)+
 scale_colour_manual(values=c(M0='#767676',M1='#367B89',M2='#B9564F'),
 labels=c('M0: group','M1: age + sex','M2: clinical'))+
 scale_shape_manual(values=c(Up=2,Down=6))+
 scale_y_discrete(labels=function(x)vapply(x,wrap,character(1),n=34))+
 scale_x_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1))+
 labs(title='Fixed immune GO comparisons',subtitle='36 SCI / 7 TC; BH correction across 4,032 sets per model',
 x='CAMERA BH q value',y=NULL,colour='Model',shape='Competitive direction',
 caption=paste('All eight terms and all three models are displayed; dashed line: q = 0.05.',
 'Direction reflects competitive enrichment, not pathway activity. All displayed q values exceed 0.05.',sep='\n'))+
 theme(legend.position='bottom',plot.title=element_text(size=10))
save_pub(p2,2,'Figure2_Fixed_immune_GO_CAMERA',137)

# 3. All 30 frozen members; vector tiles and vector colorbar.
heat<-read_src('Figure3_heatmap.csv');membership<-read_src('Frozen_30_gene_membership.csv')
heat[,gene:=factor(gene,levels=rev(membership$gene))]
heat[,sample:=factor(sample,levels=unique(sample[order(sample_order)]))]
heat[,group:=factor(group,levels=c('HC','TC','SCI'),labels=c('HC (n = 10)','TC (n = 10)','SCI (n = 38)'))]
heat[,module:=factor(module,levels=unique(membership$module),labels=c('Myeloid / innate-associated','Lymphocyte / T-cell-associated'))]
p3<-ggplot(heat,aes(sample,gene,fill=z))+geom_tile(width=1,height=1)+
 facet_grid(module~group,scales='free',space='free')+
 scale_fill_gradient2(low='#367B89',mid='white',high='#B9564F',limits=c(-3,3),oob=scales::squish,
  breaks=c(-3,0,3),name='Row z',guide=guide_colourbar(display='rectangles',barheight=grid::unit(25,'mm')))+
 labs(title='Frozen 30-gene discovery signature',subtitle='Discovery-label-selected genes; descriptive display',
      x=NULL,y=NULL,caption=paste('All 30 frozen genes and all 58 participants are shown. Gene-wise z scores are display transformations.',
 'Color is saturated at ±3; original values are retained. This is not independent validation.',sep='\n'))+
 theme(axis.text.x=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),
       axis.text.y=element_text(size=8,face='italic'),panel.spacing=grid::unit(3,'pt'),
       legend.position='right',strip.text.y=element_text(size=7.5),
       plot.title=element_text(size=10))
save_pub(p3,3,'Figure3_Frozen_30_gene_heatmap',181)

# 4. Saved conditional diagnostics: 5 k values, 6 gap references, all 20 cluster summaries.
metrics<-read_src('Figure4_cluster_metrics.csv');gap<-read_src('Figure4_cluster_gap.csv')
boot<-read_src('Figure4_cluster_bootstrap_summary.csv');dec<-read_src('Figure4_cluster_decision.csv')
p4a<-ggplot(metrics,aes(k,silhouette))+geom_line(colour='#9BA3A8',linewidth=.45)+
 geom_point(size=2.2,colour='#367B89')+geom_point(data=metrics[k==dec$silhouette_k[1]],size=2.8,colour='#B9564F')+
 scale_x_continuous(breaks=2:6)+labs(title='Silhouette',subtitle='Largest value at k = 2',x='Number of gene clusters (k)',y='Mean silhouette')
p4b<-ggplot(gap,aes(k,gap))+geom_line(colour='#9BA3A8',linewidth=.45)+
 geom_errorbar(aes(ymin=gap-SE.sim,ymax=gap+SE.sim),width=.12,linewidth=.4)+geom_point(size=2.2,colour='#367B89')+
 geom_point(data=gap[k==dec$gap_k[1]],size=2.8,colour='#B9564F')+
 scale_x_continuous(breaks=1:6)+labs(title='Gap statistic',subtitle='One-SE choice: boundary k = 6',x='Number of gene clusters (k)',y='Gap statistic ± SE')
boot[,solution:=factor(k,levels=2:6,labels=paste0('k = ',2:6))]
p4c<-ggplot(boot,aes(cluster,mean_jaccard))+
 geom_hline(yintercept=.5,linetype=2,linewidth=.3,colour='#9BA3A8')+
 geom_errorbar(aes(ymin=q025,ymax=q975),width=.15,linewidth=.4,colour='#367B89')+
 geom_point(size=1.8,colour='#367B89')+facet_wrap(~solution,nrow=1,scales='free_x')+
 scale_x_continuous(breaks=1:6)+scale_y_continuous(limits=c(0,1.03),breaks=c(0,.5,1))+
 labs(title='Conditional bootstrap membership stability',x='Gene-cluster label within each fixed solution',y='Mean Jaccard and percentile range')+
 theme(panel.spacing=grid::unit(7,'pt'))
fig4<-((p4a|p4b)/p4c)+plot_layout(heights=c(1,1))+
 plot_annotation(title='Gene-clustering diagnostics',subtitle='2,509 selected gene profiles across 58 participants',tag_levels='a',
 caption=paste('Bootstrap points: mean best-matched Jaccard; bars: empirical 2.5th–97.5th percentiles over 200 resamples.',
 'Diagnostics condition on the saved gene set. The criteria do not consistently select k = 3.',
 'Clusters describe genes, not participant immune states.',sep='\n'),theme=annotation_theme) & tag_theme
save_pub(fig4,4,'Figure4_Gene_clustering_diagnostics',163)

# 5. Full frozen 15+15 longitudinal observations and all nine saved time-only contrasts.
long<-read_src('Figure5_longitudinal_scores.csv')
long[,time:=factor(time,levels=c('acute','3m','6m','12m'),labels=c('Acute','3m','6m','12m'))]
ll<-melt(long,id.vars=c('sample','subject','time'),measure.vars=score_keys,variable.name='score',value.name='value')
ll[,score:=sfactor(score)]
p5a<-ggplot(ll,aes(time,value,group=subject))+
 geom_line(colour='#A7AFB5',linewidth=.35)+geom_point(aes(colour=score),size=1.7)+
 facet_wrap(~score,scales='free_y',nrow=1)+scale_colour_manual(values=score_col)+
 labs(title='Participant-linked observations',x='Visit',y='Control-referenced score')+theme(legend.position='none')
cc<-read_src('Figure5_saved_mixed_model_contrasts.csv');cc[,score_label:=sfactor(score)]
cc[,contrast:=factor(contrast,levels=rev(c('3m - acute','6m - acute','12m - acute')),
 labels=rev(c('3m − acute','6m − acute','12m − acute')))]
p5b<-ggplot(cc,aes(estimate,contrast,colour=score_label))+
 geom_vline(xintercept=0,linetype=2,linewidth=.35,colour='#8D979F')+
 geom_segment(aes(x=lower.CL,xend=upper.CL,yend=contrast),linewidth=.55)+geom_point(size=2)+
 facet_wrap(~score_label,scales='free_x',nrow=1)+scale_colour_manual(values=score_col)+
 labs(title='Saved random-intercept model contrasts',x='Later minus acute estimate (pointwise 95% CI)',y=NULL)+
 theme(legend.position='none')
fig5<-(p5a/p5b)+plot_layout(heights=c(1.2,1))+
 plot_annotation(title='Longitudinal frozen-score patterns',subtitle='Full 15+15 coverage; 10 SCI participants, 37 observations',tag_levels='a',
 caption=paste('Visit counts: acute 10; 3 months 10; 6 months 10; 12 months 7. Lines connect the same participant.',
 'All nine time-only contrasts and pointwise CIs are shown; full model sensitivities remain in the supplement.',
 'Scores describe bulk-expression patterns, not immune or neurological recovery.',sep='\n'),theme=annotation_theme) & tag_theme
save_pub(fig5,5,'Figure5_Frozen_score_longitudinal',177)

# 6. Rebuild coordinates from frozen 15+15; no supervised continuum or transition arrows.
disc<-read_src('Figure6_frozen_discovery_scores.csv');disc[,group:=factor(group,levels=c('HC','TC','SCI'))]
p6a<-ggplot(disc,aes(innate_score,adaptive_score,colour=group,shape=group))+
 geom_hline(yintercept=0,linewidth=.3,colour='#D9DEE2')+geom_vline(xintercept=0,linewidth=.3,colour='#D9DEE2')+
 geom_point(size=2.1,alpha=.88)+scale_colour_manual(values=cohort_col)+scale_shape_manual(values=group_shape)+
 labs(title='Frozen component coordinates',x='Myeloid / innate-associated score',y='Lymphocyte / T-cell-associated score',colour=NULL,shape=NULL)+
 theme(legend.position='bottom')
dl<-melt(disc,id.vars=c('sample','group'),measure.vars=score_keys[1:2],variable.name='score',value.name='value')
dl[,score:=sfactor(score)]
p6b<-ggplot(dl,aes(group,value,colour=group))+
 geom_boxplot(width=.55,outlier.shape=NA,linewidth=.4,show.legend=FALSE)+
 geom_point(aes(shape=group),position=position_jitter(width=.1,height=0,seed=20261028),size=1.5,alpha=.8)+
 facet_wrap(~score,scales='free_y',nrow=1)+scale_colour_manual(values=cohort_col)+scale_shape_manual(values=group_shape)+
 labs(title='Component distributions',x=NULL,y='Control-referenced score')+theme(legend.position='none')
fig6<-(p6a/p6b)+plot_layout(heights=c(1.2,1))+
 plot_annotation(title='Frozen-score discovery distributions',subtitle='HC n = 10  |  TC n = 10  |  SCI n = 38',tag_levels='a',
 caption=paste('Every point is a participant. Boxes: median and IQR; whiskers: observations within 1.5 × IQR.',
 'These group labels contributed to signature selection. No postselection P values or disease-progression',
 'axis are added; this selected coordinate space is descriptive.',sep='\n'),theme=annotation_theme) & tag_theme
save_pub(fig6,6,'Figure6_Frozen_score_discovery_coordinates',170)

# 7. All 45 observations and all 12 final scope-specific effect estimates.
ext<-read_src('Figure7_additional_cohort_scores.csv');ext[,group:=factor(group,levels=c('CTL','SCI'))]
stopifnot(all(ext$overlap %in% c('new','reused')))
ext[,reuse:=factor(overlap,levels=c('new','reused'),labels=c('New SCI','Reused participant'))]
el<-melt(ext,id.vars=c('sample','subject_id','group','overlap','reuse'),measure.vars=score_keys,variable.name='score',value.name='value')
el[,score:=sfactor(score)]
p7a<-ggplot(el,aes(group,value,colour=group))+
 geom_boxplot(width=.5,outlier.shape=NA,linewidth=.4,show.legend=FALSE)+
 geom_point(aes(shape=reuse),position=position_jitter(width=.1,height=0,seed=20261028),size=1.65)+
 facet_wrap(~score,scales='free_y',nrow=1)+scale_colour_manual(values=cohort_col,guide='none')+
 scale_shape_manual(values=c('New SCI'=16,'Reused participant'=2))+
 scale_x_discrete(labels=c(CTL='CTL\nn = 9',SCI='SCI\nn = 36'))+
 labs(title='Observed additional-cohort scores',x=NULL,y='Supplied score',shape=NULL)+theme(legend.position='bottom')
ef<-read_src('Figure7_reuse_platform_effects.csv');ef[,score_label:=sfactor(score)]
scope_order<-c('all_available','exclude_reused_SCI','same_platform_all','same_platform_new_SCI')
stopifnot(setequal(unique(ef$scope),scope_order))
scope_text<-c('All samples\n36 / 9','Exclude reused SCI\n26 / 9','Same platform\n31 / 9','Same platform, new SCI\n21 / 9')
ef[,scope_label:=factor(scope,levels=rev(scope_order),labels=rev(scope_text))]
forest<-function(key,show_y=FALSE){
 z<-ef[score==key]
 p<-ggplot(z,aes(effect_g,scope_label))+geom_vline(xintercept=0,linetype=2,linewidth=.3,colour='#9BA4AB')+
  geom_segment(aes(x=CI_low,xend=CI_high,yend=scope_label),linewidth=.5,colour=score_col[[as.character(sfactor(key))]])+
  geom_point(size=1.9,colour=score_col[[as.character(sfactor(key))]])+
  labs(title=as.character(sfactor(key)),x='Hedges g',y=NULL)+theme(plot.title=element_text(size=8,face='bold'))
 if(!show_y)p<-p+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
 p
}
p7b<-forest('innate_score',TRUE)+forest('adaptive_score')+forest('balance_score')+plot_layout(widths=c(1.6,1,1))
fig7<-(p7a/p7b)+plot_layout(heights=c(1.2,1))+
 plot_annotation(title='Additional-cohort assessment with partial reuse',
 subtitle='GSE293559: 14+14 genes available; TLR5 and ITK absent',
 caption=paste('Triangles identify reused participants: 10 SCI and all 9 controls. Excluding reused SCI retains all controls.',
 'Lower panels: all four scopes (SCI / CTL n), saved Hedges g and participant-bootstrap 95% CIs.',
 'Scores retain the supplied transformation. Partial reuse and lack of trauma controls limit independence.',sep='\n'),
 tag_levels='a',theme=annotation_theme) & tag_theme
save_pub(fig7,7,'Figure7_Additional_cohort_reuse_sensitivity',187)

manifest<-rbindlist(exports)
mf<-file.path(out,'qa','export_manifest.csv')
if(file.exists(mf)&&length(selected)<7){
 old<-fread(mf);manifest<-rbind(old[!number %in% selected],manifest);setorder(manifest,number)
}
fwrite(manifest,mf)
capture.output(sessionInfo(),file=file.path(out,'qa','figure_R_sessionInfo.txt'))
write_json(list(status='RENDERED_PENDING_VISUAL_QA',backend='R',width_cap_mm=W,rendered=selected,
 biological_models_refit=FALSE,legacy_files_modified=FALSE),file.path(out,'qa','render_status.json'),pretty=TRUE,auto_unbox=TRUE)

