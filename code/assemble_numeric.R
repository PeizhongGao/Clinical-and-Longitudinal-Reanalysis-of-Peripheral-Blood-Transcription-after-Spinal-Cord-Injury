# Assemble existing estimands for final reporting; no additional inference.
library(data.table)
root<-normalizePath(commandArgs(TRUE)[1],winslash='/',mustWork=TRUE)
out<-file.path(root,'revision_numeric_freeze_20261019/tables');dir.create(out,recursive=TRUE,showWarnings=FALSE)
base<-file.path(root,'revision_external_adjusted_20261016/tables')
d<-fread(file.path(base,'discovery_models.csv'));l<-fread(file.path(base,'longitudinal_models.csv'))
dp<-d[subset=='clinical43'&spec=='clinical'&variant=='published'];lp<-l[scope=='all37'&spec=='disjoint_both'&variant=='published']
dp[,source_branch:='discovery'];lp[,source_branch:='longitudinal']
b<-rbindlist(list(dp,lp),fill=TRUE);stopifnot(nrow(b)==16L,all(is.finite(b$q_BH_primary16)))
fwrite(b,file.path(out,'BTM_primary16.csv'));fwrite(d,file.path(out,'BTM_discovery_all88.csv'));fwrite(l,file.path(out,'BTM_longitudinal_all240.csv'))
g<-fread(file.path(root,'revision_followup_20261019/bridge43/results/fixed_eight_GO_panel_all_models.csv'))
fwrite(g,file.path(out,'bridge43_fixed_eight_GO.csv'))
