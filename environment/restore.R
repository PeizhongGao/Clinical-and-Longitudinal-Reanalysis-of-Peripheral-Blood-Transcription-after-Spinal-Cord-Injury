# Install the recorded project environment into a separate, ignored library.
args<-commandArgs(TRUE)
root<-if(length(args))normalizePath(args[1],winslash='/',mustWork=TRUE)else normalizePath('.',winslash='/',mustWork=TRUE)
stopifnot(file.exists(file.path(root,'renv.lock')))
if(!requireNamespace('renv',quietly=TRUE))install.packages('renv',repos='https://cloud.r-project.org')
lib<-file.path(root,'environment','r_library');dir.create(lib,recursive=TRUE,showWarnings=FALSE)
renv::restore(project=root,lockfile=file.path(root,'renv.lock'),library=lib,prompt=FALSE)
cat('Restored library:',lib,'\nThe runner uses this directory through R_LIBS_USER.\n')
