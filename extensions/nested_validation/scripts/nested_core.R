suppressPackageStartupMessages({library(data.table);library(limma);library(glmnet)})
options(warn=1)
glmnet.control(fdev=0,devmax=1)
FRACTIONS <- exp(seq(log(1),log(.001),length.out=25))
INNATE_KEYS <- c('S100','FCGR','TLR','CLEC','NLR','IL1','CXCL','CCR','LILR','CD14','CD68','ITG','MMP','MPO','ELANE','AZU','CAMP','CEACAM','CSF','IRAK','MYD88','NFKB','ANXA','BST','BCL2A1','MIR223','SLC','PLIN','QPCT')
LYMPH_KEYS <- c('CD3','CD4','CD8','TRAC','TRBC','TRAV','TRBV','IL7R','CCR7','TCF7','LEF1','LCK','ZAP70','NKG7','GNLY','PRF1','GZMA','GZMB','KLR','FCGR3A')
INNATE_MANUAL <- c('S100A8','S100A9','S100A12','ANXA3','CLEC4D','FCGR1A','BST1','IRAK3','TLR5','NLRC4','BCL2A1','MIR223','QPCT','SLC37A3','SLC22A4','PLIN3')
LYMPH_MANUAL <- c('CD3D','CD3E','TRAC','TRBC1','TRBC2','IL7R','CCR7','TCF7','LEF1','ZAP70','NKG7','GNLY','PRF1','GZMA','GZMB','KLRD1','FCGR3A')

read_dataset <- function(root,source='raw',contrast='SCI_HC') {
  meta<-fread(file.path(root,'inputs/participants58.csv'))
  stopifnot(nrow(meta)==58,uniqueN(meta$subject)==58,all(meta$sex %in% c('F','M')))
  d<-fread(file.path(root,'inputs',if(source=='raw')'raw_counts.csv' else 'processed_expression.csv'))
  g<-d[[1]];E<-as.matrix(d[,-1,with=FALSE]);storage.mode(E)<-'double';rownames(E)<-g
  stopifnot(!anyDuplicated(g),all(is.finite(E)),setequal(colnames(E),meta$sample))
  ref<-if(contrast=='SCI_HC')'HC' else 'TC'
  meta<-meta[group %in% c('SCI',ref)];meta<-meta[order(sample)]
  E<-E[,meta$sample,drop=FALSE];CPM<-NULL
  if(source=='raw') {
    stopifnot(all(E>=0),all(E==floor(E)),all(colSums(E)>0))
    CPM<-sweep(E,2,colSums(E),'/')*1e6;E<-log2(CPM+.5)
  }
  y<-as.integer(meta$group=='SCI');stopifnot(sum(y)==38,sum(y==0)==10)
  list(E=E,CPM=CPM,meta=meta,y=y,source=source,contrast=contrast)
}

stratified_folds <- function(y,k,seed) {
  stopifnot(length(unique(y))==2,min(table(y))>=k)
  set.seed(seed);f<-integer(length(y))
  for(v in sort(unique(y))) {
    ix<-which(y==v);f[ix]<-sample(rep(seq_len(k),length.out=length(ix)))
  }
  for(i in seq_len(k))stopifnot(length(unique(y[f==i]))==2)
  f
}

training_state <- function(ds,idx,y) {
  E<-ds$E[,idx,drop=FALSE];n<-length(idx)
  v<-rowSums((E-rowMeans(E))^2)/(n-1)
  eligible<-is.finite(v)&v>1e-12
  if(ds$source=='raw')eligible<-eligible & rowSums(ds$CPM[,idx,drop=FALSE]>=1)>=ceiling(.2*n)
  E<-E[eligible,,drop=FALSE]
  stopifnot(nrow(E)>10,length(unique(y[idx]))==2)
  fit<-eBayes(lmFit(E,cbind('(Intercept)'=1,case=y[idx])),trend=FALSE,robust=FALSE)
  q<-p.adjust(fit$p.value[,2],method='BH');b<-fit$coefficients[,2];g<-rownames(E)
  innate<-g[q<.1 & b>0 & grepl(paste(INNATE_KEYS,collapse='|'),g,ignore.case=TRUE)]
  lymph<-g[q<.1 & b<0 & grepl(paste(LYMPH_KEYS,collapse='|'),g,ignore.case=TRUE)]
  candidate<-sort(unique(c(innate,lymph,intersect(c(INNATE_MANUAL,LYMPH_MANUAL),g))))
  x<-t(E[candidate,,drop=FALSE]);mu<-colMeans(x);sdv<-if(ncol(x))apply(x,2,sd) else numeric()
  stopifnot(all(is.finite(sdv)),all(sdv>0))
  x<-sweep(sweep(x,2,mu,'-'),2,sdv,'/')
  list(x=x,genes=candidate,mean=mu,sd=sdv,eligible_n=nrow(E),de_candidate_n=length(unique(c(innate,lymph))),prior=mean(y[idx]))
}

fit_path <- function(state,y) {
  x<-state$x
  if(!ncol(x))return(list(fit=NULL,genes=character(),lambda=rep(NA_real_,length(FRACTIONS)),prior=mean(y)))
  lm<-max(abs(drop(crossprod(x,y-mean(y)))))/length(y)
  if(lm<1e-12)return(list(fit=NULL,genes=state$genes,lambda=rep(NA_real_,length(FRACTIONS)),prior=mean(y)))
  if(ncol(x)==1)x<-cbind(x,.padding_zero=0)
  lam<-lm*FRACTIONS
  warns<-character()
  f<-withCallingHandlers(glmnet(x,y,family='binomial',alpha=1,lambda=lam,standardize=FALSE,intercept=TRUE,control=list(thresh=1e-7,maxit=100000)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
  if(any(!grepl('fewer than 8',warns)))stop('Unexpected glmnet warning: ',paste(warns,collapse=';'))
  if(f$jerr!=0)stop('glmnet nonzero jerr: ',f$jerr)
  list(fit=f,genes=state$genes,lambda=lam,prior=mean(y),small_class_warning_count=length(warns))
}

predict_path <- function(fit,state,ds,idx) {
  if(is.null(fit$fit))return(matrix(fit$prior,nrow=length(idx),ncol=length(FRACTIONS)))
  x<-t(ds$E[state$genes,idx,drop=FALSE]);x<-sweep(sweep(x,2,state$mean,'-'),2,state$sd,'/')
  if(ncol(x)==1)x<-cbind(x,.padding_zero=0)
  p<-as.matrix(predict(fit$fit,newx=x,s=fit$lambda,type='response'))
  # A numerically empty model must return the training prior exactly. Otherwise
  # ~1e-16 solver noise can create artificial ranks within an intercept-only fold.
  co<-as.matrix(coef(fit$fit,s=fit$lambda))
  empty<-colSums(abs(co[rownames(co)!='(Intercept)',,drop=FALSE])>1e-10)==0L
  if(any(empty))p[,empty]<-fit$prior
  stopifnot(all(is.finite(p)),all(p>=0 & p<=1),ncol(p)==length(FRACTIONS))
  p
}

write_predictions <- function(x,path) {
  out<-copy(x)
  for(nm in c('prediction','baseline'))out[,(nm):=sprintf('%.17g',get(nm))]
  fwrite(out,path,quote=FALSE)
}

logloss <- function(y,p) -mean(y*log(pmax(p,1e-15))+(1-y)*log(pmax(1-p,1e-15)))
auc_rank <- function(y,p) (sum(rank(p,ties.method='average')[y==1])-sum(y)*(sum(y)+1)/2)/(sum(y)*sum(y==0))
metrics <- function(y,p) {
  a<-data.table(y=y,p=p)[,.(pos=sum(y),n=.N),by=p][order(-p)]
  ap<-sum(a$pos/sum(y) * cumsum(a$pos)/cumsum(a$n))
  c(AUC=auc_rank(y,p),average_precision=ap,Brier=mean((p-y)^2),log_loss=logloss(y,p))
}

nested_once <- function(ds,y=ds$y,seed=20261027,detail=TRUE) {
  f<-stratified_folds(y,5,seed);p<-base<-rep(NA_real_,length(y));records<-list();splits<-list();inner_audit<-list()
  for(outer in 1:5) {
    train<-which(f!=outer);test<-which(f==outer);stopifnot(!length(intersect(train,test)))
    inn<-stratified_folds(y[train],3,seed+1000+outer)
    losses<-matrix(NA_real_,3,length(FRACTIONS));innern<-integer(3);innerwarnings<-0L
    for(j in 1:3) {
      tr<-train[inn!=j];va<-train[inn==j]
      stopifnot(!length(intersect(tr,va)),!length(intersect(tr,test)),!length(intersect(va,test)))
      st<-training_state(ds,tr,y);fit<-fit_path(st,y[tr]);pp<-predict_path(fit,st,ds,va)
      innerwarnings<-innerwarnings+sum(fit$small_class_warning_count)
      losses[j,]<-apply(pp,2,function(v)2*logloss(y[va],v));innern[j]<-length(st$genes)
      if(detail) {
        splits[[paste(outer,j)]]<-data.table(outer_fold=outer,inner_fold=j,sample=ds$meta$sample,role=ifelse(seq_along(y)%in%test,'outer_test',ifelse(seq_along(y)%in%va,'inner_validation','inner_train')))
        inner_audit[[paste(outer,j)]]<-data.table(outer_fold=outer,inner_fold=j,train_n=length(tr),validation_n=length(va),eligible_n=st$eligible_n,candidate_n=length(st$genes),candidates=paste(st$genes,collapse=';'))
      }
    }
    ml<-colMeans(losses);se<-apply(losses,2,sd)/sqrt(3);best<-which.min(ml);chosen<-which(ml<=ml[best]+se[best]+1e-12)[1]
    st<-training_state(ds,train,y);fit<-fit_path(st,y[train]);pp<-predict_path(fit,st,ds,test)
    p[test]<-pp[,chosen];base[test]<-mean(y[train]);selected<-character();coefvalues<-numeric()
    if(!is.null(fit$fit)) {
      co<-as.matrix(coef(fit$fit,s=fit$lambda[chosen]));nz<-rownames(co)[abs(co[,1])>1e-10];selected<-intersect(nz,st$genes);coefvalues<-co[selected,1]
    }
    records[[outer]]<-data.table(outer_fold=outer,train_n=length(train),test_n=length(test),eligible_n=st$eligible_n,candidate_n=length(st$genes),selected_n=length(selected),lambda_fraction=FRACTIONS[chosen],lambda=fit$lambda[chosen],inner_min_deviance=ml[best],inner_selected_deviance=ml[chosen],inner_candidate_counts=paste(innern,collapse=';'),inner_small_class_warnings=innerwarnings,outer_small_class_warnings=sum(fit$small_class_warning_count),selected=paste(selected,collapse=';'),coefficients=paste(signif(coefvalues,12),collapse=';'))
  }
  stopifnot(all(is.finite(p)),all(is.finite(base)))
  list(metrics=metrics(y,p),baseline=metrics(y,base),predictions=data.table(sample=ds$meta$sample,subject=ds$meta$subject,y=y,outer_fold=f,prediction=p,baseline=base),audit=rbindlist(records),splits=if(detail)rbindlist(splits) else NULL,inner_audit=if(detail)rbindlist(inner_audit) else NULL)
}

permute_sex <- function(y,sex,seed) {
  set.seed(seed);yp<-y
  for(s in sort(unique(sex))) {idx<-which(sex==s);yp[idx]<-sample(y[idx])}
  for(s in unique(sex))stopifnot(sum(yp[sex==s])==sum(y[sex==s]))
  yp
}
