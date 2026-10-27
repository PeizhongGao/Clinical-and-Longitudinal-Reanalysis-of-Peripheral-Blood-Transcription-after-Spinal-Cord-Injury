# Public-package copy: see metadata/code_changes.csv for provenance.
"""Weighted preranked gene-set permutation test using sparse hit positions.

ES uses the same weighted hit/miss running-sum statistic as GSEA. Uniform
no-replacement sets are drawn by replacing duplicates until all indices are
unique. This is invariant under relabeling and therefore uniform over sets
of a fixed size. Independent null rows are used for each set. BH is our
reported correction; this does not implement the canonical pooled GSEA FDR.
"""
from pathlib import Path
import sys,json,time,importlib.util
import numpy as np,pandas as pd
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[2];O=ROOT/'revision_R2_20261008';R1=ROOT/'revision_R1_20261007'
seed=20261008;B=10000;rng=np.random.default_rng(seed)
def bh(p):
 ix=np.argsort(p);q=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=np.empty(len(p));out[ix]=np.minimum(q,1);return out
def statistic(h,weights,n):
 m=h.shape[-1];w=weights[h];den=w.sum(axis=-1,keepdims=True)
 if np.any(den==0):raise ValueError('zero-weight gene set')
 w=w/den;after=np.cumsum(w,axis=-1)-(h-np.arange(m))/(n-m)
 before=after-w;hi=np.max(after,axis=-1);lo=np.min(before,axis=-1)
 es=np.where(hi>=-lo,hi,lo);return es,after,before
def draw(n,m):
 x=rng.integers(0,n,size=(B,m),dtype=np.int32);x.sort(axis=1)
 while True:
  bad=np.zeros_like(x,dtype=bool);bad[:,1:]=x[:,1:]==x[:,:-1];k=bad.sum()
  if not k:return x
  x[bad]=rng.integers(0,n,size=k,dtype=np.int32);x.sort(axis=1)
sets=pd.read_csv(R1/'tables/GO_camera_membership.csv');terms=pd.read_csv(R1/'inputs/go_bp_terms.csv')
# Public package: import the installed cross-platform extension.
from gseapy import gse
audit=[]
for contrast in ['SCI_vs_TC','SCI_vs_HC']:
 start=time.time();tab=pd.read_csv(R1/f'tables/{contrast}_processed_limma.csv')
 ranks=tab[['gene','t']].sort_values(['t','gene'],ascending=[False,True]);assert ranks.gene.is_unique and np.isfinite(ranks.t).all()
 names=ranks.gene.to_numpy();vals=ranks.t.to_numpy();weight=np.abs(vals);n=len(names);lookup={g:i for i,g in enumerate(names)}
 ranks.to_csv(O/f'tables/{contrast}_GSEA_ranks.tsv',sep='\t',index=False)
 # Independently validate the observed ES against the installed reference Rust engine.
 check=sets.iloc[:8];ref=gse.prerank_rs(names.tolist(),vals.tolist(),{r.go_id:r.genes.split(';') for r in check.itertuples()},1.,15,500,1,1,seed)
 refes={z.term:z.es for z in ref.summaries};maxdiff=0.;rows=[]
 for j,r in enumerate(sets.itertuples()):
  h=np.array(sorted(lookup[g] for g in set(r.genes.split(';')) if g in lookup));m=len(h);assert 15<=m<=500
  es,a,b=statistic(h,weight,n);es=float(es)
  if r.go_id in refes:maxdiff=max(maxdiff,abs(es-refes[r.go_id]))
  x=draw(n,m);null,_,_=statistic(x,weight,n)
  same=null[null>=0] if es>=0 else null[null<0];ext=int(np.sum(same>=es if es>=0 else same<=es));p=(ext+1)/(len(same)+1)
  nes=es/np.mean(np.abs(same));peak=int(np.argmax(a) if es>=0 else np.argmin(b));lead=h[:peak+1] if es>=0 else h[peak:]
  rows.append(dict(go_id=r.go_id,ES=es,NES=nes,nominal_p=p,BH_q=np.nan,same_sign_permutations=len(same),extreme_permutations=ext,n_genes=m,leading_edge=';'.join(names[lead])))
  if (j+1)%500==0:print(contrast,j+1,'sets',round(time.time()-start,1),'seconds',flush=True)
 assert maxdiff<1e-10,maxdiff
 d=pd.DataFrame(rows);d.BH_q=bh(d.nominal_p.to_numpy());d=d.merge(terms,on='go_id',how='left');d['contrast']=contrast;d.to_csv(O/f'tables/{contrast}_GSEA.csv',index=False)
 audit.append(dict(contrast=contrast,n_genes=n,n_sets=len(d),BH_significant=int((d.BH_q<.05).sum()),permutations=B,seed=seed,weight=1,engine='numpy sparse weighted running sum; reference gseapy 1.3.1 observed ES check',reference_ES_max_absolute_difference=maxdiff,reference_checked_terms=8,p_rule='(same-sign exceedances+1)/(same-sign permutations+1)',correction='BH within all eligible sets; not canonical pooled GSEA FDR',tied_t=int(ranks.t.duplicated().sum()),elapsed_seconds=time.time()-start))
 (O/'logs/GSEA_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8');print('COMPLETE',audit[-1],flush=True)
