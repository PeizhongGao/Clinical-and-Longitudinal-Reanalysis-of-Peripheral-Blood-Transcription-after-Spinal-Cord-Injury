root<-normalizePath(commandArgs(TRUE)[1],winslash='/',mustWork=TRUE)
required<-c('data.table','limma','readxl','jsonlite','ggplot2','patchwork','lmerTest','emmeans','AnnotationDbi','org.Hs.eg.db','edgeR','lmtest','digest','svglite','ragg','cluster')
missing<-required[!vapply(required,requireNamespace,quietly=TRUE,FUN.VALUE=logical(1))]
if(length(missing))stop('Missing packages: ',paste(missing,collapse=', '),'. See environment/README.md.')
if(as.character(packageVersion('limma'))!='3.68.4')stop('Exact locked gene/GO replication requires limma 3.68.4.')
if(as.character(packageVersion('org.Hs.eg.db'))!='3.23.1')stop('Locked gene mapping requires org.Hs.eg.db 3.23.1.')
write.csv(data.frame(package=required,version=vapply(required,function(p)as.character(packageVersion(p)),'')),file.path(root,'runtime_packages.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,'runtime_sessionInfo.txt'))
cat('Environment prerequisites passed. R:',R.version.string,'\n')
