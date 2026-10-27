# Reconstruct and compare inputs before the sealed statistical script is executed.
# No model-selection or statistical-inference changes are made here.
library(data.table);library(jsonlite)
root<-normalizePath(commandArgs(TRUE)[1],winslash='/',mustWork=TRUE)
b<-file.path(root,'revision_scope_lock_20261018/bridge43')
cov<-fread(file.path(b,'inputs_local/clinical43_covariates.csv'))
genes<-fread(file.path(b,'inputs_local/gene_universe.csv'))$gene
raw<-fread(file.path(root,'explore/Report/GSE151371_normalized_expression_matrix_comBat_de-ID.csv'))
raw<-raw[!duplicated(raw[[1]])];ix<-match(genes,raw[[1]]);stopifnot(!anyNA(ix))
E<-as.matrix(raw[ix,cov$sample,with=FALSE]);storage.mode(E)<-'double';rownames(E)<-genes
saved<-readRDS(file.path(b,'inputs_local/expression43.rds'))
stopifnot(identical(dim(E),c(17498L,43L)),identical(colnames(E),colnames(saved)),max(abs(E-saved))<1e-12)
clinical<-fread(file.path(root,'revision_integrated_20261012/fixed_sample/tables/fixed_input_clinical43.csv'))
for(v in c('sample','group','sex','tbi','prior_cns'))stopifnot(identical(as.character(cov[[v]]),as.character(clinical[[v]])))
for(v in c('age','iss'))stopifnot(max(abs(cov[[v]]-clinical[[v]]))<1e-12)
cov[,group:=factor(group,levels=c('TC','SCI'))];cov[,sex:=factor(sex,levels=c('F','M'))]
for(v in c('tbi','prior_cns'))set(cov,j=v,value=factor(cov[[v]],levels=c('No','Yes','Unknown')))
cov[,age_c:=age-mean(age)];cov[,iss_c:=iss-mean(iss)]
spec<-fromJSON(file.path(b,'analysis_spec.json'),simplifyVector=FALSE)
designs<-readRDS(file.path(b,'designs_local/design_matrices.rds'))
for(m in spec$models){
 vars<-all.vars(as.formula(m$formula));fac<-intersect(vars,c('group','sex','tbi','prior_cns'));contr<-setNames(rep(list('contr.treatment'),length(fac)),fac)
 X<-model.matrix(as.formula(m$formula),data=cov,contrasts.arg=contr)
 stopifnot(identical(colnames(X),colnames(designs[[m$id]])),max(abs(unname(X)-unname(designs[[m$id]])))<1e-12)
}
membership<-fread(file.path(b,'inputs_local/GO_camera_membership.csv'))
idx<-lapply(strsplit(membership$genes,';',fixed=TRUE),match,table=genes);names(idx)<-membership$go_id
stopifnot(identical(idx,readRDS(file.path(b,'inputs_local/GO_indices.rds'))))
write_json(list(status='PASS',expression_source_match=TRUE,clinical_source_match=TRUE,design_reconstruction_match=TRUE,GO_index_reconstruction_match=TRUE),file.path(root,'bridge_input_reconstruction.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Source matrix, clinical membership, 3 designs and 4032 GO indices reconstructed and matched.\n')
