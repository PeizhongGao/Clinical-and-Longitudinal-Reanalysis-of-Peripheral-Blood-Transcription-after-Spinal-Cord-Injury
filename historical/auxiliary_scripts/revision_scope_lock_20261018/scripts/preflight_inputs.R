# Public-package copy: see metadata/code_changes.csv for provenance.
# Input selection and design checks only. No lm/lmFit/eBayes/CAMERA effect fit.
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(limma)})
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else getwd(), winslash='/', mustWork=TRUE)
out <- file.path(root,'revision_scope_lock_20261018')
stopifnot(!file.exists(file.path(out,'bridge43/protocol_lock.json')))
spec <- fromJSON(file.path(out,'bridge43/analysis_spec.json'),simplifyVector=FALSE)
src <- spec$sources
loc <- file.path(out,'bridge43/inputs_local')
des <- file.path(out,'bridge43/designs_local')
logs <- file.path(out,'bridge43/logs')
for(d in c(loc,des,logs)) dir.create(d,recursive=TRUE,showWarnings=FALSE)
checks <- list()
ck <- function(name,pass,detail='') {
  checks[[length(checks)+1L]] <<- list(check=name,passed=isTRUE(pass),detail=detail)
  if(!isTRUE(pass)) stop(paste(name,detail),call.=FALSE)
}
rd <- function(path) fread(path,check.names=FALSE,na.strings=c('','NA'))
d <- rd(src$clinical43_authority)
d2 <- rd(src$clinical43_reused_input)
d3 <- rd(src$clinical43_latest_covariates)
keys <- c('sample','gsm','subject','group')
cv <- c('age','sex','iss','tbi','prior_cns')
ck('exact_existing_43_membership',nrow(d)==43L && uniqueN(d$sample)==43L && uniqueN(d$gsm)==43L && uniqueN(d$subject)==43L)
ck('group_counts',sum(d$group=='SCI')==36L && sum(d$group=='TC')==7L)
ck('source_reuse_byte_content_fields',identical(d[,c(keys,cv),with=FALSE],d2[,c(keys,cv),with=FALSE]))
ck('latest_BTM_input_same_keys_and_covariates',identical(d[,c(keys,cv),with=FALSE],d3[,c(keys,cv),with=FALSE]))
cross <- rd(src$sample_crosswalk); j <- match(d$sample,cross$sample)
ck('crosswalk_43_found_unique',!anyNA(j) && !anyDuplicated(j))
ck('crosswalk_identity_group',all(vapply(keys,function(k)identical(as.character(d[[k]]),as.character(cross[[k]][j])),logical(1))))
ck('crosswalk_clinical_covariates',all(vapply(cv,function(k)isTRUE(all.equal(d[[k]],cross[[k]][j],check.attributes=FALSE,tolerance=1e-12)),logical(1))))
ck('complete_required_fields',all(complete.cases(d[,c(keys,cv),with=FALSE])))
ck('numeric_finite',all(is.finite(d$age)) && all(is.finite(d$iss)))

# Use raw clinical values; do not carry the old full48-centered columns or outcomes.
meta <- copy(d[,c(keys,cv),with=FALSE])
meta[,group:=factor(group,levels=c('TC','SCI'))]
meta[,sex:=factor(sex,levels=c('F','M'))]
for(v in c('tbi','prior_cns'))set(meta,j=v,value=factor(meta[[v]],levels=c('No','Yes','Unknown')))
ck('all_source_category_values_preserved',all(complete.cases(meta)))
centers <- data.table(variable=c('age','iss'),center=c(mean(meta$age),mean(meta$iss)),population='same_locked_clinical43')
meta[,age_c:=age-centers[variable=='age',center]]
meta[,iss_c:=iss-centers[variable=='iss',center]]
fwrite(meta,file.path(loc,'clinical43_covariates.csv'),na='NA')
fwrite(meta[,..keys],file.path(loc,'participant_order.csv'))
fwrite(centers,file.path(loc,'centering_constants.csv'))

g <- rd(src$tested_gene_universe)$gene
old_genes <- rd(src$original_gene_order)$gene
ck('17498_unique_genes_original_order',length(g)==17498L && !anyDuplicated(g) && identical(g,old_genes))
raw <- rd(src$expression_matrix)
gene_col <- names(raw)[1]
ck('raw_matrix_has_exact_58_sample_keys',ncol(raw)==59L && nrow(cross)==58L && setequal(names(raw)[-1],cross$sample))
duplicates <- sum(duplicated(raw[[gene_col]]))
raw <- raw[!duplicated(raw[[gene_col]])]
gi <- match(g,raw[[gene_col]])
ck('all_frozen_genes_present',!anyNA(gi))
ck('all_frozen_samples_present',all(meta$sample%in%names(raw)))
E <- as.matrix(raw[gi,meta$sample,with=FALSE]); storage.mode(E)<-'double';rownames(E)<-g
ck('matrix_17498_by43',identical(dim(E),c(17498L,43L)))
ck('matrix_finite',all(is.finite(E)))
ck('expression_and_clinical_order',identical(colnames(E),meta$sample))
saveRDS(E,file.path(loc,'expression43.rds'),version=3,compress='gzip')
ck('saved_matrix_exact_roundtrip',identical(E,readRDS(file.path(loc,'expression43.rds'))))
stopifnot(file.copy(src$tested_gene_universe,file.path(loc,'gene_universe.csv'),overwrite=TRUE))

go <- rd(src$go_membership)
members <- strsplit(go$genes,';',fixed=TRUE);names(members)<-go$go_id
ck('4032_unique_go_sets',nrow(go)==4032L && !anyDuplicated(go$go_id))
ck('go_set_sizes_and_uniqueness',all(lengths(members)==go$n_genes) && all(lengths(members)>=15L & lengths(members)<=500L) && all(vapply(members,function(z)!anyDuplicated(z),logical(1))))
ck('go_genes_in_frozen_universe',all(unlist(members,use.names=FALSE)%in%g))
idx <- lapply(members,match,table=g)
term <- rd(src$go_terms_and_old_camera)
ck('go_term_id_identity',identical(go$go_id,term$go_id))
ck('go_term_membership_counts',identical(as.integer(go$n_genes),as.integer(term$NGenes)))
focus <- rd(src$go_display_panel)
ck('existing_eight_term_display_panel',nrow(focus)==8L && !anyDuplicated(focus$go_id) && all(focus$go_id%in%go$go_id))
stopifnot(file.copy(src$go_membership,file.path(loc,'GO_camera_membership.csv'),overwrite=TRUE))
fwrite(term[,.(go_id,term)],file.path(loc,'GO_terms.csv'))
fwrite(focus[,.(go_id,term)],file.path(loc,'GO_display_panel.csv'))
stopifnot(file.copy(src$go_annotation_metadata,file.path(loc,'GO_annotation_metadata.csv'),overwrite=TRUE))
saveRDS(idx,file.path(loc,'GO_indices.rds'),version=3,compress='gzip')

design_summary <- list(); Xs <- list()
for(m in spec$models){
  vars <- all.vars(as.formula(m$formula))
  fac <- intersect(vars,c('group','sex','tbi','prior_cns'))
  cc <- setNames(rep(list('contr.treatment'),length(fac)),fac)
  X <- model.matrix(as.formula(m$formula),data=meta,contrasts.arg=cc)
  q <- qr(X); rank <- q$rank
  ck(paste0(m$id,'_full_rank_and_df'),rank==m$expected_rank && rank==ncol(X) && nrow(X)-rank==m$expected_residual_df)
  ck(paste0(m$id,'_groupSCI_estimable'),'groupSCI'%in%colnames(X))
  ck(paste0(m$id,'_identical_row_order'),nrow(X)==43L && all(is.finite(X)))
  Xs[[m$id]]<-X
  fwrite(data.table(sample=meta$sample,X,check.names=FALSE),file.path(des,paste0('design_',m$id,'.csv')))
  h <- rowSums(qr.Q(q)[,seq_len(rank),drop=FALSE]^2)
  design_summary[[m$id]] <- data.table(model_id=m$id,formula=paste('expression',m$formula),n=43L,SCI_n=36L,TC_n=7L,rank=rank,columns=ncol(X),residual_df=nrow(X)-rank,coef='groupSCI',max_design_leverage=max(h),condition_number_scaled=kappa(cbind(1,scale(X[,-1,drop=FALSE])),exact=TRUE))
}
saveRDS(Xs,file.path(des,'design_matrices.rds'),version=3,compress='gzip')
fwrite(rbindlist(design_summary),file.path(logs,'design_summary.csv'))

# Demonstrate the historical full48 vs current43 centering distinction at design level only.
old_centers <- rd(src$fixed_sample_centers)
old_X <- rd(src$existing_fixed_design)
oi <- match(meta$sample,old_X$sample)
expected <- copy(Xs$M2)
age_delta <- centers[variable=='age',center]-old_centers[variable=='age',center]
iss_delta <- centers[variable=='iss',center]-old_centers[variable=='iss',center]
expected[,'age_c']<-expected[,'age_c']+age_delta
expected[,'iss_c']<-expected[,'iss_c']+iss_delta
ck('old_score_design_only_differs_by_center_constants',!anyNA(oi) && max(abs(expected-as.matrix(old_X[oi,colnames(expected),with=FALSE])))<1e-10)

old_de <- rd(src$old_clinical_de_candidate)
ck('reuse_candidate_exact_gene_order',identical(old_de$gene,g) && nrow(old_de)==17498L)
ck('reuse_candidate_model_n_and_rank',all(old_de$n_samples==43L) && all(old_de$n_participants==43L) && all(old_de$case_n==36L) && all(old_de$reference_n==7L) && all(old_de$design_rank==9L) && all(old_de$residual_df==34L))
ck('reuse_candidate_formula',all(old_de$model=='expression ~ group + age_c + sex + iss_c + tbi + prior_cns'))
ck('reuse_candidate_full_finite_statistics',all(is.finite(as.matrix(old_de[,.(logFC,CI.L,CI.R,t,P.Value,adj.P.Val)]))))
ck('limma_version_locked',as.character(packageVersion('limma'))=='3.68.4')

capture.output(sessionInfo(),file=file.path(logs,'sessionInfo.txt'))
default_values <- list(eBayes=lapply(formals(limma::eBayes),deparse),CAMERA=lapply(formals(limma::camera.default),deparse))
write_json(default_values,file.path(logs,'limma_function_defaults.json'),auto_unbox=TRUE,pretty=TRUE)
report <- list(all_passed=TRUE,checks=checks,n_checks=length(checks),participants=list(n=43,SCI=36,TC=7),expression=list(n_genes=nrow(E),n_samples=ncol(E),source_duplicate_rows_retained_first=duplicates),GO=list(n_sets=nrow(go),minimum=min(go$n_genes),maximum=max(go$n_genes),fixed_display_terms=focus$go_id),centering=list(age=mean(meta$age),iss=mean(meta$iss),historical_center_difference_checked_at_design_only=TRUE),models=rbindlist(design_summary),old_M2_reuse_candidate=list(provenance_and_structure_checked=TRUE,numerical_reproduction_performed=FALSE),software=list(R=R.version.string,limma=as.character(packageVersion('limma')),data_table=as.character(packageVersion('data.table')),jsonlite=as.character(packageVersion('jsonlite'))),new_effect_models_fitted=0L,new_inferential_tests=0L,new_CAMERA_tests=0L,manuscript_or_reply_edited=FALSE,clinical_inputs_published=FALSE)
write_json(report,file.path(logs,'preflight_QA.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
cat('INPUT_AND_DESIGN_PREFLIGHT_PASS:',length(checks),'checks; 43 participants; 17498 genes; 4032 GO; no effect model fitted.\n')
