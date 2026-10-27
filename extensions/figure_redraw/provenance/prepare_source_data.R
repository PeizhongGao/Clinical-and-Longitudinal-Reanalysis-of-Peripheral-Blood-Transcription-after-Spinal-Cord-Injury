# Prepare compact, real source-data snapshots for the seven revised figures.
# All biological estimates are read from saved revision results. No model refits.
if(.Platform$OS.type=='windows') suppressWarnings(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(digest)})
args <- commandArgs(trailingOnly=TRUE)
project <- normalizePath(if(length(args)) args[1] else '.', winslash='/', mustWork=TRUE)
out <- file.path(project,'figures','Figures_new')
for (d in c('source_data','previews','qa')) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
provenance <- list()
record <- function(path,role,member='') {
  provenance[[length(provenance)+1L]] <<- data.table(source=sub(paste0(project,'/'),'',normalizePath(path,winslash='/'),fixed=TRUE),
    archive_member=member,role=role,sha256=digest(path,algo='sha256',file=TRUE))
}
put <- function(x,name) fwrite(as.data.table(x),file.path(out,'source_data',name),na='NA')
read_project <- function(path,role=path) {
  p<-file.path(project,path);record(p,role);fread(p)
}
read_archive <- function(prefix,member,role=member) {
  files<-list.files(file.path(project,'revision_signature_resolution_20261027','supplement'),recursive=TRUE,full.names=TRUE)
  archive<-files[basename(files)==prefix];stopifnot(length(archive)==1L)
  temp<-tempfile('figure_archive_');dir.create(temp)
  utils::unzip(archive,files=member,exdir=temp)
  d<-fread(file.path(temp,member));record(archive,role,member);d
}

legacy<-list.files(file.path(project,'figures'),pattern='^Figure[1-7].*\\.(pdf|tif)$',full.names=TRUE)
stopifnot(length(legacy)==14L)
legacy_hashes<-data.table(file=basename(legacy),sha256=vapply(legacy,function(p)digest(p,algo='sha256',file=TRUE),character(1)))
hash_file<-file.path(out,'qa','legacy_original_sha256.csv')
if(file.exists(hash_file)) stopifnot(identical(fread(hash_file)$sha256,legacy_hashes$sha256)) else fwrite(legacy_hashes,hash_file)

state_path<-file.path(project,'revision_R1_20261007','tables','analysis_state.rds')
record(state_path,'Saved processed-expression, frozen discovery scores and additional-cohort scores')
s<-readRDS(state_path); E<-s$discovery; m<-as.data.table(s$s151)
stopifnot(nrow(E)==17498L,ncol(E)==58L,identical(colnames(E),m$sample))
stopifnot(all(is.finite(E)),all(table(m$group)[c('HC','TC','SCI')]==c(10,10,38)))
variable <- apply(E,1,sd)>0
pca<-prcomp(t(E[variable,,drop=FALSE]),center=TRUE,scale.=TRUE)
v<-pca$sdev^2/sum(pca$sdev^2)
put(data.table(sample=m$sample,group=m$group,PC1=pca$x[,1],PC2=pca$x[,2]),'Figure1_PCA.csv')
put(data.table(component=paste0('PC',seq_along(v)),variance_fraction=v,n_genes=sum(variable),n_participants=58L),'Figure1_PCA_variance.csv')
counts<-list()
for(comp in c('TC_vs_HC','SCI_vs_HC','SCI_vs_TC')) {
  d<-read_archive('Supplementary_Data_2_Gene_GO_results.zip',paste0('Historical_',comp,'_processed_limma.csv'))
  stopifnot(nrow(d)==17498L,all(is.finite(d$logFC)),all(is.finite(d$adj.P.Val)),all(d$adj.P.Val>0))
  d[,direction := 'Other']
  d[adj.P.Val < .05 & logFC > 1, direction := 'Up']
  d[adj.P.Val < .05 & logFC < -1, direction := 'Down']
  put(d,paste0('Figure1_',comp,'.csv'))
  counts[[comp]]<-data.table(comparison=comp,category=c('Up','Down','Other'),n=as.integer(table(factor(d$direction,levels=c('Up','Down','Other')))),universe=nrow(d))
}
counts<-rbindlist(counts);counts[,percent:=n/universe*100]
stopifnot(counts[comparison=='SCI_vs_HC' & category=='Up',n]==1150L,
          counts[comparison=='TC_vs_HC' & category=='Up',n]==127L)
put(counts,'Figure1_DEG_counts.csv')

go<-read_archive('Supplementary_Data_2_Gene_GO_results.zip','Fixed_eight_GO.csv')
stopifnot(nrow(go)==24L,uniqueN(go$go_id)==8L,all(go$n_participants==43L),all(go$FDR>.05))
put(go,'Figure2_fixed_eight_GO.csv')

sigpath<-file.path(project,'revision_signature_resolution_20261027','evidence','locked_signature_v1.json')
record(sigpath,'Frozen 15+15 membership');sig<-fromJSON(sigpath)
genes<-data.table(gene=c(sig$innate_genes,sig$adaptive_genes),module=rep(c('Myeloid / innate-associated','Lymphocyte / T-cell-associated'),each=15),order=seq_len(30))
stopifnot(setequal(unlist(s$genes,use.names=FALSE),genes$gene),uniqueN(genes$gene)==30L)
put(genes,'Frozen_30_gene_membership.csv')
heat<-as.data.table(as.table(t(scale(t(E[genes$gene,,drop=FALSE])))))
setnames(heat,c('gene','sample','z'))
heat[,gene:=as.character(gene)];heat[,sample:=as.character(sample)]
heat[,group:=m$group[match(sample,m$sample)]]
heat[,module:=genes$module[match(gene,genes$gene)]]
heat[,sample_order:=match(sample,m$sample)]
stopifnot(nrow(heat)==1740L,all(is.finite(heat$z)))
put(heat,'Figure3_heatmap.csv')

for (name in c('Cluster_metrics.csv','Cluster_decision.csv','Cluster_bootstrap_summary.csv','cluster_gap.csv','cluster_bootstrap.csv')) {
  d<-read_archive('Supplementary_Data_7_Multiplicity_and_disposition.zip',name)
  put(d,paste0('Figure4_',tolower(name)))
}

ncbipath<-file.path(project,'revision_R1_20261007','tables','NCBI_analysis_state.rds')
record(ncbipath,'Full-count frozen longitudinal score observations')
n<-readRDS(ncbipath);long<-as.data.table(n$sfull)[group=='SCI',.(sample,subject,time,innate_score,adaptive_score,balance_score)]
stopifnot(nrow(long)==37L,uniqueN(long$subject)==10L,all(table(long$time)[c('acute','3m','6m','12m')]==c(10,10,10,7)))
put(long,'Figure5_longitudinal_scores.csv')
contrast<-read_archive('Supplementary_Data_4_Original_15plus15_scores.zip','Longitudinal_contrasts36.csv')
keep<-contrast[source=='NCBI_full_15plus15' & scope=='all_available' & spec=='unadjusted']
stopifnot(nrow(keep)==9L,all(keep$n_samples==37L),all(keep$n_participants==10L))
put(keep,'Figure5_saved_mixed_model_contrasts.csv')

put(m[,.(sample,group,innate_score,adaptive_score,balance_score)],'Figure6_frozen_discovery_scores.csv')
ref<-read_project('revision_R1_20261007/tables/standardization_parameters.csv','Archived control-reference parameters')
put(ref,'Frozen_score_reference_parameters.csv')

additional<-as.data.table(s$s293)[,.(sample,subject_id,group,overlap,innate_score,adaptive_score,balance_score)]
stopifnot(nrow(additional)==45L,uniqueN(additional$sample)==45L,
          additional[group=='SCI',.N]==36L,additional[group=='CTL',.N]==9L,
          additional[group=='SCI' & overlap=='reused',.N]==10L,
          additional[group=='CTL' & overlap=='reused',.N]==9L)
put(additional,'Figure7_additional_cohort_scores.csv')
effects<-read_project('revision_signature_resolution_20261027/supplement/Supplementary_Tables/Table_S5_GSE293559_reuse_sensitivities.csv','Final 12 additional-cohort sensitivity estimates')
stopifnot(nrow(effects)==12L,uniqueN(effects$scope)==4L,all(effects$n_control==9L))
put(effects,'Figure7_reuse_platform_effects.csv')

fwrite(rbindlist(provenance),file.path(out,'qa','source_provenance.csv'))
write_json(list(status='PASS',width_cap_mm=170,backend='R',no_biological_model_refit=TRUE,
  source_counts=list(discovery_participants=58,discovery_genes=17498,heatmap_genes=30,heatmap_tiles=1740,
  fixed_GO_rows=24,cluster_genes=2509,longitudinal_participants=10,longitudinal_observations=37,
  longitudinal_contrasts=9,additional_participants=45,additional_effects=12),
  display_rules=list(heatmap='row-wise z, color saturation only beyond [-3,3]; full values retained',
    PCA='all nonconstant genes, centered and unit-variance scaled; no participant exclusions',
    no_downsampling=TRUE)),file.path(out,'qa','source_data_checks.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'qa','prepare_R_sessionInfo.txt'))
cat('SOURCE DATA PASS: 7 figure contracts prepared; original 14 files preserved.\n')
