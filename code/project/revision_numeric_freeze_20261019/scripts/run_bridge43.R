# Public-package copy: see metadata/code_changes.csv for provenance.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1L) stop("Pass the absolute SCI project root as the first argument.")
root <- normalizePath(args[1], winslash="/", mustWork=TRUE)
if(file.exists(file.path(root,"revision_numeric_freeze_20261019/FREEZE.json"))) stop("Numeric results are frozen. Use the documented correction procedure.")
scope <- file.path(root,"revision_scope_lock_20261018")
bridge <- file.path(scope,"bridge43")
out <- file.path(root,"revision_followup_20261019","bridge43")
resdir <- file.path(out,"results")
logdir <- file.path(out,"logs")
dir.create(resdir, recursive=TRUE, showWarnings=FALSE)
dir.create(logdir, recursive=TRUE, showWarnings=FALSE)
logfile <- file.path(logdir,"bridge43_execution.log")
sink(logfile, split=TRUE)
on.exit(sink(), add=TRUE)
cat("Analysis: locked 43-participant fixed-set gene and GO bridge v1\n")
cat("Start Asia/Hong_Kong:",format(Sys.time(),tz="Asia/Hong_Kong",usetz=TRUE),"\n")
suppressPackageStartupMessages({library(data.table);library(limma);library(digest);library(jsonlite)})

lock <- fromJSON(file.path(bridge,"public_portability_lock.json"),simplifyVector=FALSE)
spec <- fromJSON(file.path(bridge,"analysis_spec.json"),simplifyVector=FALSE)
stopifnot(identical(lock$analysis_id,"bridge43_gene_GO_20261018_v1"),
          identical(spec$analysis_id,lock$analysis_id))
preflight <- fromJSON(file.path(bridge,"logs","preflight_QA.json"),simplifyVector=FALSE)
stopifnot(isTRUE(preflight$all_passed), preflight$participants$n==43,
          preflight$participants$SCI==36,preflight$participants$TC==7,
          preflight$expression$n_genes==17498,preflight$GO$n_sets==4032,
          identical(as.character(packageVersion("limma")),"3.68.4"))

## Preserve the sealed packet byte-for-byte before and after execution.
locked_checks <- rbindlist(lapply(lock$locked_files,function(x) {
  p <- file.path(scope,x$relative_path)
  if (!file.exists(p)) stop("Missing sealed packet input: ",x$relative_path)
  h <- digest(file=p,algo="sha256",serialize=FALSE)
  data.table(path=x$relative_path,expected_sha256=x$sha256,observed_sha256=h,
             matches=identical(h,x$sha256))
}))
if (!all(locked_checks$matches)) stop("Sealed packet SHA-256 mismatch; stopped before fitting.")

E <- readRDS(file.path(bridge,"inputs_local","expression43.rds"))
designs <- readRDS(file.path(bridge,"designs_local","design_matrices.rds"))
sets <- readRDS(file.path(bridge,"inputs_local","GO_indices.rds"))
genes <- fread(file.path(bridge,"inputs_local","gene_universe.csv"))$gene
cov <- fread(file.path(bridge,"inputs_local","clinical43_covariates.csv"))
terms <- fread(file.path(bridge,"inputs_local","GO_terms.csv"))
membership <- fread(file.path(bridge,"inputs_local","GO_camera_membership.csv"))
stopifnot(identical(dim(E),c(17498L,43L)),identical(rownames(E),genes),
          identical(colnames(E),as.character(cov$sample)),
          length(sets)==4032L,identical(names(sets),as.character(terms$go_id)),
          identical(membership$go_id,as.character(terms$go_id)),
          all(is.finite(E)),all(vapply(sets,function(i)length(i)>=15L&&length(i)<=500L,logical(1))))

coef_tabs <- list(); camera_tabs <- list(); fit_diags <- list(); repro <- NULL
old_candidate_path <- file.path(root, spec$sources$old_clinical_de_candidate)
old_candidate <- fread(old_candidate_path)
stopifnot(nrow(old_candidate)==17498L,identical(old_candidate$gene,genes),
          identical(unique(old_candidate$model),"expression ~ group + age_c + sex + iss_c + tbi + prior_cns"),
          all(old_candidate$n_participants==43),all(old_candidate$design_rank==9))

for (model_id in c("M0","M1","M2")) {
  design_csv <- fread(file.path(bridge,"designs_local",paste0("design_",model_id,".csv")),
                      check.names=FALSE)
  stopifnot("sample"%in%names(design_csv),
            identical(as.character(design_csv$sample),as.character(cov$sample)),
            identical(names(design_csv)[-1],colnames(designs[[model_id]])))
  D <- as.matrix(design_csv[,setdiff(names(design_csv),"sample"),with=FALSE])
  storage.mode(D) <- "double"
  rownames(D) <- as.character(design_csv$sample)
  storage.mode(designs[[model_id]]) <- "double"
  stopifnot(isTRUE(all.equal(unname(D),unname(designs[[model_id]]),
                             check.attributes=FALSE,tolerance=1e-12)),
            nrow(D)==43L,identical(rownames(D),colnames(E)),
            "groupSCI" %in% colnames(D),qr(D)$rank==ncol(D))
  expected <- spec$models[[match(model_id,vapply(spec$models,`[[`,"","id"))]]
  stopifnot(ncol(D)==expected$expected_rank,43L-ncol(D)==expected$expected_residual_df)
  fit <- lmFit(E, D, method="ls")
  fit <- eBayes(fit,proportion=spec$gene_test$ebayes$proportion,
                stdev.coef.lim=unlist(spec$gene_test$ebayes$stdev.coef.lim),
                trend=spec$gene_test$ebayes$trend,
                span=spec$gene_test$ebayes$span,
                robust=spec$gene_test$ebayes$robust,
                winsor.tail.p=unlist(spec$gene_test$ebayes$winsor.tail.p),
                legacy=spec$gene_test$ebayes$legacy)
  tt <- topTable(fit,coef="groupSCI",number=Inf,sort.by="none",confint=.95)
  tt$gene <- rownames(tt)
  stopifnot(nrow(tt)==17498L,identical(tt$gene,genes),
            all(is.finite(tt$logFC)),all(is.finite(tt$CI.L)),all(is.finite(tt$CI.R)),
            all(is.finite(tt$P.Value)),all(is.finite(tt$adj.P.Val)))
  tab <- as.data.table(tt)
  tab[,`:=`(model=model_id,formula=expected$formula,n_samples=43L,
            n_participants=43L,SCI_n=36L,TC_n=7L,
            design_rank=qr(D)$rank,residual_df=fit$df.residual[1],
            interval="moderated t pointwise 95% CI; not simultaneous")]
  tab[,q_family:=spec$gene_test$BH]
  fwrite(tab,file.path(resdir,paste0("gene_models_",model_id,".csv")),na="NA")
  coef_tabs[[model_id]] <- tab
  if (model_id=="M2") {
    cols <- c("logFC","CI.L","CI.R","t","P.Value","adj.P.Val")
    cmp <- merge(tab[,c("gene",cols),with=FALSE],
                 old_candidate[,c("gene",cols),with=FALSE],
                 by="gene",suffixes=c("_new","_old"),sort=FALSE)
    stopifnot(nrow(cmp)==17498L,identical(cmp$gene,genes))
    ctab <- rbindlist(lapply(cols,function(nm) {
      a <- cmp[[paste0(nm,"_new")]];b <- cmp[[paste0(nm,"_old")]]
      data.table(statistic=nm,max_abs_difference=max(abs(a-b)),
                 max_relative_difference=max(abs(a-b)/pmax(abs(b),.Machine$double.eps)),
                 tolerance=1e-8,passed=all(is.finite(a))&&all(is.finite(b))&&
                   all(abs(a-b)<=1e-8+1e-8*abs(b)))
    }))
    if (!all(ctab$passed)) {
      fwrite(ctab,file.path(resdir,"M2_reproduction_check_FAILED.csv"))
      stop("M2 reproduction differs from the frozen candidate outside tolerance.")
    }
    repro <- ctab
    fwrite(ctab,file.path(resdir,"M2_reproduction_check.csv"),na="NA")
  }
  ci <- camera(E,index=sets,design=D,contrast=match("groupSCI",colnames(D)),
               weights=NULL,use.ranks=FALSE,allow.neg.cor=FALSE,
               inter.gene.cor=.01,trend.var=FALSE,sort=FALSE,directional=TRUE)
  ci <- as.data.table(ci,keep.rownames="go_id")
  stopifnot(nrow(ci)==4032L,identical(ci$go_id,names(sets)),
            all(is.finite(ci$PValue)),all(is.finite(ci$FDR)),
            all(is.finite(ci$NGenes)),all(ci$NGenes>=15L),all(ci$NGenes<=500L))
  ci[,term:=terms$term[match(go_id,terms$go_id)]]
  ci[,`:=`(model=model_id,n_participants=43L,SCI_n=36L,TC_n=7L,
           design_rank=qr(D)$rank,residual_df=43L-qr(D)$rank,
           test="limma CAMERA; directional competitive gene-set test; inter.gene.cor=0.01",
           q_family=spec$pathway_test$BH)]
  setcolorder(ci,c("model","go_id","term","NGenes","Direction","PValue","FDR",
                   "n_participants","SCI_n","TC_n","design_rank","residual_df","test","q_family"))
  fwrite(ci,file.path(resdir,paste0("camera_models_",model_id,".csv")),na="NA")
  camera_tabs[[model_id]] <- ci
  fit_diags[[model_id]] <- data.table(model=model_id,formula=expected$formula,
       n_participants=43L,SCI_n=36L,TC_n=7L,design_rank=qr(D)$rank,
       residual_df=43L-qr(D)$rank,max_hat=max(diag(D%*%solve(crossprod(D),t(D)))),
       condition_number_scaled=preflight$models[[match(model_id,vapply(preflight$models,`[[`,"","model_id"))]]$condition_number_scaled,
       t2_prior_df=fit$df.prior,limma_version=as.character(packageVersion("limma")))
  cat("Completed",model_id,"43-participant genes:",nrow(tab),"GO sets:",nrow(ci),"\n")
}

fwrite(rbindlist(fit_diags),file.path(resdir,"model_diagnostics.csv"))
wide_effects <- lapply(c("M0","M1","M2"),function(m) {
  z <- copy(coef_tabs[[m]][,.(gene,logFC,P_value=P.Value,q_BH=adj.P.Val,
                              CI_lower=CI.L,CI_upper=CI.R)])
  setnames(z,setdiff(names(z),"gene"),paste0(setdiff(names(z),"gene"),"_",m))
  z
})
cmpm <- Reduce(function(x,y)merge(x,y,by="gene",sort=FALSE),wide_effects)
stopifnot(nrow(cmpm)==17498L,identical(cmpm$gene,genes))
cmpm[,`:=`(delta_M1_minus_M0=logFC_M1-logFC_M0,
           delta_M2_minus_M1=logFC_M2-logFC_M1,
           delta_M2_minus_M0=logFC_M2-logFC_M0,
           CI_zero_M0=CI_lower_M0<=0 & CI_upper_M0>=0,
           CI_zero_M1=CI_lower_M1<=0 & CI_upper_M1>=0,
           CI_zero_M2=CI_lower_M2<=0 & CI_upper_M2>=0)]
fwrite(cmpm,file.path(resdir,"gene_model_effects_descriptive_paired.csv"),na="NA")
gene_summary <- data.table(
  contrast=c("M1 minus M0","M2 minus M1","M2 minus M0"),
  median_gene_delta=c(median(cmpm$delta_M1_minus_M0),median(cmpm$delta_M2_minus_M1),median(cmpm$delta_M2_minus_M0)),
  median_absolute_gene_delta=c(median(abs(cmpm$delta_M1_minus_M0)),median(abs(cmpm$delta_M2_minus_M1)),median(abs(cmpm$delta_M2_minus_M0))),
  spearman_effect_rank=c(cor(cmpm$logFC_M0,cmpm$logFC_M1,method="spearman"),
                         cor(cmpm$logFC_M1,cmpm$logFC_M2,method="spearman"),
                         cor(cmpm$logFC_M0,cmpm$logFC_M2,method="spearman")),
  sign_reversals=c(sum((cmpm$logFC_M0>0&cmpm$logFC_M1<0)|(cmpm$logFC_M0<0&cmpm$logFC_M1>0)),
                   sum((cmpm$logFC_M1>0&cmpm$logFC_M2<0)|(cmpm$logFC_M1<0&cmpm$logFC_M2>0)),
                   sum((cmpm$logFC_M0>0&cmpm$logFC_M2<0)|(cmpm$logFC_M0<0&cmpm$logFC_M2>0))),
  CI_excludes_zero=c(sum(!cmpm$CI_zero_M0),sum(!cmpm$CI_zero_M1),sum(!cmpm$CI_zero_M2)),
  genes_BH_q_lt_0_05=c(sum(cmpm$q_BH_M0<.05),sum(cmpm$q_BH_M1<.05),sum(cmpm$q_BH_M2<.05)))
fwrite(gene_summary,file.path(resdir,"gene_model_change_summary.csv"))

panel_ids <- fread(file.path(bridge,"inputs_local","GO_display_panel.csv"))$go_id
panel <- rbindlist(lapply(names(camera_tabs),function(m)camera_tabs[[m]][go_id%chin%panel_ids]))
stopifnot(nrow(panel)==24L)
fwrite(panel,file.path(resdir,"fixed_eight_GO_panel_all_models.csv"))

## Fixed-sample clinical covariate and common-support summaries (no new fitting).
clinical <- copy(cov)
cat_summary <- rbindlist(lapply(c("sex","tbi","prior_cns"),function(v) {
  z <- clinical[, .N, by=c("group",v)]
  setnames(z,v,"category")
  z[,variable:=v]
  setcolorder(z,c("variable","group","category","N"))
  setnames(z,"N","n")
  z
}),fill=TRUE)
group_summary <- clinical[, .N, by="group"][,.(variable="group",group=as.character(group),category=as.character(group),n=N)]
cat_summary <- rbind(cat_summary,group_summary,fill=TRUE)
fwrite(cat_summary,file.path(resdir,"clinical43_category_counts.csv"),na="NA")

## Snapshot only; no later paths are changed or selected by p/q.
locked_checks_after <- rbindlist(lapply(lock$locked_files,function(x) {
  p <- file.path(scope,x$relative_path)
  h <- digest(file=p,algo="sha256",serialize=FALSE)
  data.table(path=x$relative_path,expected_sha256=x$sha256,observed_sha256=h,matches=identical(h,x$sha256))
}))
if (!all(locked_checks_after$matches)) stop("A sealed input changed during analysis.")
fwrite(locked_checks,file.path(logdir,"sealed_packet_hash_check.csv"))
fwrite(locked_checks_after,file.path(logdir,"sealed_packet_hash_check_after.csv"))
fwrite(data.table(source_role=c("old_M2_candidate","locked_protocol","sealed_input_packet"),
  source_path=c(old_candidate_path,file.path(bridge,"analysis_protocol_zh.md"),file.path(bridge,"public_portability_lock.json")),
  sha256=c(digest(file=old_candidate_path,algo="sha256",serialize=FALSE),
           digest(file=file.path(bridge,"analysis_protocol_zh.md"),algo="sha256",serialize=FALSE),
           digest(file=file.path(bridge,"public_portability_lock.json"),algo="sha256",serialize=FALSE))),
  file.path(logdir,"execution_references_sha256.csv"))

freeze <- list(analysis_id=lock$analysis_id,locked_protocol_id=lock$analysis_id,
  execution_date_hong_kong=format(Sys.time(),tz="Asia/Hong_Kong",usetz=TRUE),
  execution_stage="REVISION_STAGE_EXPLORATORY_EXECUTION_AFTER_PRIOR_RESULTS_KNOWN",
  prospectively_registered=FALSE,participants=list(n=43,SCI=36,TC=7),
  genes_tested=17498,GO_sets_per_model=4032,
  models=list(M0="expression ~ group",M1="expression ~ group + age_c + sex",
              M2="expression ~ group + age_c + sex + iss_c + tbi + prior_cns"),
  software=list(R=R.version.string,limma=as.character(packageVersion("limma")),
                data_table=as.character(packageVersion("data.table")),
                digest=as.character(packageVersion("digest"))),
  M2_reproduction=list(status="PASS",fields=repro$statistic,
       maximum_absolute_difference=max(repro$max_abs_difference),
       maximum_relative_difference=max(repro$max_relative_difference),tolerance=1e-8),
    corrections=list(gene="three separate BH families: 17,498 genes within each model",GO="three separate BH families: 4,032 sets within each model",CI="pointwise moderated 95% intervals; not simultaneous"),
  interpretation="Association and sensitivity to model specification only; not causal decomposition, independent validation, cell-intrinsic mechanism or prediction.",
  results=list(genes=list(M0=17498,M1=17498,M2=17498),GO=list(M0=4032,M1=4032,M2=4032)),
  sealed_packet_hashes_match=all(locked_checks_after$matches),new_LASSO_or_selection=FALSE,
  old_LASSO_and_old_nine_gene_dispositions="WITHDRAWN_HISTORICAL_NOT_CURRENT_EVIDENCE",
  invalid_supervised_label_permutation_disposition="WITHDRAWN_NO_VALID_P_INFERENCE")
write_json(freeze,file.path(out,"analysis_freeze.json"),auto_unbox=TRUE,pretty=TRUE,null="null")
capture.output(sessionInfo(),file=file.path(logdir,"sessionInfo.txt"))
cat("PASS: M2 prior-result reproduction (all six locked statistics, abs/rel <= 1e-8).\n")
cat("PASS: three gene families and three pathway families; all full result tables retained.\n")
cat("PASS: sealed pre-analysis protocol, input and design packet SHA-256 unchanged.\n")
cat("Finish Asia/Hong_Kong:",format(Sys.time(),tz="Asia/Hong_Kong",usetz=TRUE),"\n")





