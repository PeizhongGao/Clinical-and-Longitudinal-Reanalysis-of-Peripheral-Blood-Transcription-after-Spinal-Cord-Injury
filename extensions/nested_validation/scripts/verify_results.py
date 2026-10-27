from pathlib import Path
import hashlib,json,csv
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
def read(name):return pd.read_csv(ROOT/'tables'/name,float_precision='round_trip')
checks={}
lock=json.loads((ROOT/'audit/protocol_input_lock.json').read_text())
checks['protocol_unchanged']=hashlib.sha256((ROOT/'PROTOCOL.md').read_bytes()).hexdigest()==lock['protocol_sha256']
checks['inputs_unchanged']=all(hashlib.sha256((ROOT/r['snapshot']).read_bytes()).hexdigest()==r['sha256'] for r in lock['inputs'])
m=read('nested_repeat_metrics.csv');p=read('nested_oof_predictions.csv');s=read('nested_split_assignments.csv');a=read('nested_outer_fit_audit.csv');inner=read('nested_inner_screen_audit.csv')
checks['all_expected_repeats_folds_predictions']=len(m)==80 and len(p)==3840 and len(a)==400 and len(s)==57600 and len(inner)==1200
checks['predictions_finite_in_range']=np.isfinite(p.prediction).all() and p.prediction.between(0,1).all()
checks['each_subject_once_per_outer_repeat']=not p.duplicated(['source','contrast','repeat_id','subject']).any()
checks['each_inner_role_partition_complete']=s.groupby(['source','contrast','repeat_id','outer_fold','inner_fold']).apply(lambda d: len(d)==48 and d['sample'].nunique()==48,include_groups=False).all().item()
fold_checks=[]
meta=pd.read_csv(ROOT/'inputs/participants58.csv').set_index('sample')
for key,d in s.groupby(['source','contrast','repeat_id','outer_fold','inner_fold']):
    pred=p[(p.source==key[0])&(p.contrast==key[1])&(p.repeat_id==key[2])].set_index('sample')
    fold_checks.append(set(d.loc[d.role=='outer_test','sample'])==set(pred.index[pred.outer_fold==key[3]]))
    for role in ['inner_train','inner_validation','outer_test']:
        fold_checks.append(pred.loc[d.loc[d.role==role,'sample'],'y'].nunique()==2)
checks['outer_test_never_inner_train_or_validation']=all(fold_checks)
def auc(y,v):
    y=np.asarray(y);v=np.asarray(v);n=y.sum();ranks=pd.Series(v).rank(method='average').to_numpy();return (ranks[y==1].sum()-n*(n+1)/2)/(n*(len(y)-n))
metric_ok=[]
for key,d in p.groupby(['source','contrast','repeat_id']):
    row=m[(m.source==key[0])&(m.contrast==key[1])&(m.repeat_id==key[2])].iloc[0]
    metric_ok.append(abs(auc(d.y,d.prediction)-row.AUC)<1e-12)
    metric_ok.append(abs(np.mean((d.prediction-d.y)**2)-row.Brier)<1e-12)
    grouped=d.groupby('prediction',sort=True).agg(pos=('y','sum'),n=('y','size')).sort_index(ascending=False)
    ap=np.sum((grouped.pos/d.y.sum())*(grouped.pos.cumsum()/grouped.n.cumsum()))
    loss=-np.mean(d.y*np.log(np.maximum(d.prediction,1e-15))+(1-d.y)*np.log(np.maximum(1-d.prediction,1e-15)))
    metric_ok.append(abs(ap-row.average_precision)<1e-12 and abs(loss-row.log_loss)<1e-12)
checks['independent_auc_brier_recalculation']=all(metric_ok)
checks['independent_AP_logloss_recalculation']=all(metric_ok)
empty=p.merge(a.loc[a.selected_n==0,['source','contrast','repeat_id','outer_fold']],on=['source','contrast','repeat_id','outer_fold'])
checks['empty_models_return_exact_training_prior']=bool((empty.prediction==empty.baseline).all())
null=read('full_pipeline_permutation_null.csv');r=read('full_pipeline_permutation_results.csv');labels=read('permutation_label_assignments.csv');pa=read('permutation_outer_fit_audit.csv')
checks['all_1998_permutations_and_9990_outer_fits']=len(null)==1998 and len(pa)==9990 and null.groupby('contrast').size().eq(999).all().item()
checks['every_null_statistic_finite']=np.isfinite(null[['AUC','average_precision','Brier','log_loss']].to_numpy()).all().item()
pc=[]
for row in r.itertuples():
    v=null.loc[null.contrast==row.contrast,'AUC'];obs=m.loc[(m.source=='raw')&(m.contrast==row.contrast)&(m.repeat_id==1),'AUC'].iloc[0]
    pc.append(abs(row.observed_first_repeat_AUC-obs)<1e-12)
    pc.append(abs(row.p_permutation-(1+(v>=obs-1e-12).sum())/1000)<1e-12)
checks['permutation_comparator_and_plus_one_correct']=all(pc)
def holm(values):
    p=np.asarray(values);order=np.argsort(p);out=np.empty_like(p,dtype=float)
    out[order]=np.minimum(np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))),1)
    return out
assert np.allclose(holm([.01,.04,.03]),[.03,.06,.06])
checks['classifier_holm2_recalculated']=np.allclose(r.p_Holm2,holm(r.p_permutation),atol=1e-12)
blocks=labels.merge(meta[['sex','group']],left_on='sample',right_index=True)
checks['every_permutation_preserves_sex_case_counts']=all(d.groupby('sex').permuted_label.sum().equals(d.assign(true=(d['group']=='SCI').astype(int)).groupby('sex').true.sum().rename('permuted_label')) for _,d in blocks.groupby(['contrast','permutation_id']))
cr=read('clinical_resampling_results.csv');ci=read('clinical_resampling_inputs.csv');X=read('clinical_design_matrix.csv').drop(columns='sample').to_numpy();j=1
iv=np.linalg.inv(X.T@X);A=iv@X.T;h=np.sum(X*A.T,axis=1);cr_ok=[]
for mid,d in ci.groupby('module_id',sort=False):
    y=d.score.to_numpy();beta=A@y;e=y-X@beta;se=np.sqrt(np.sum((e*A[j]/(1-h))**2));row=cr[cr.module_id==mid].iloc[0]
    cr_ok.append(abs(beta[j]-row.estimate)<1e-10 and abs(se-row.HC3_SE)<1e-10)
checks['independent_numpy_HC3_reproduction']=all(cr_ok)
cn=read('clinical_resampling_null.csv');checks['clinical_8_complete_null_distributions']=len(cn)==79992 and np.isfinite(cn.t).all().item()
checks['clinical_resampling_p_recalculated']=all(abs(row.p_resampling-(1+(cn.loc[(cn.module_id==row.module_id)&(cn.method==row.method),'t'].abs()>=abs(row.observed_t)-1e-12).sum())/10000)<1e-12 for row in cr.itertuples())
holm_ok=[]
for method,d in cr.groupby('method'):
    holm_ok.append(np.allclose(d.p_Holm4,holm(d.p_resampling),atol=1e-12))
checks['clinical_holm4_families_recalculated']=all(holm_ok)
report={'all_passed':all(checks.values()),'checks':{k:bool(v) for k,v in checks.items()},'unit':'participant; CV repeats/folds and resamples are not independent participants','scientific_boundary':'Internal evaluation and conditional resampling; no independent validation or restored historical nine-gene model'}
(ROOT/'audit/analysis_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));assert report['all_passed']
