"""Compare regenerated scientific outputs with the unchanged reference results."""
from pathlib import Path
import argparse,csv,json,math
BASE=Path(__file__).resolve().parents[1]
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def missing(x):return x is None or x in ('','NA','NaN','nan')
def compare(work,reference,keys,fields):
    a=read(work);b=read(reference)
    def index(rows):
        d={}
        for row in rows:
            key=tuple(row.get(k,'') for k in keys)
            if key in d:raise ValueError('Duplicate comparison key: '+str(key))
            d[key]=row
        return d
    aa=index(a);bb=index(b);errors=[];delta={k:0.0 for k in fields}
    if set(aa)!=set(bb):errors.append('Row key sets differ')
    compared=0
    for key in aa.keys()&bb.keys():
        for field in fields:
            if field not in aa[key] or field not in bb[key]:raise ValueError('Missing field '+field)
            x,y=aa[key][field],bb[key][field]
            if missing(x) and missing(y):continue
            if missing(x)!=missing(y):errors.append(f'{key}/{field}: missingness differs');continue
            try:
                fx,fy=float(x),float(y)
            except ValueError:
                if x!=y:errors.append(f'{key}/{field}: values differ')
                continue
            if not math.isfinite(fx) or not math.isfinite(fy):
                if fx!=fy:errors.append(f'{key}/{field}: nonfinite mismatch')
                continue
            gap=abs(fx-fy);delta[field]=max(delta[field],gap);compared+=1
            if gap>1e-8+1e-8*abs(fy):errors.append(f'{key}/{field}: gap={gap}')
    return dict(status='PASS' if not errors else 'FAIL',rows_rebuilt=len(a),rows_reference=len(b),numeric_values_compared=compared,maximum_absolute_differences=delta,error_count=len(errors),errors=errors[:12])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--work',type=Path,required=True);ap.add_argument('--report',type=Path,required=True);args=ap.parse_args()
    inp=BASE/'inputs/project';ref=BASE/'reference_results';tests=[]
    def add(name,rel,keys,fields,reference=None):tests.append((name,args.work/rel,reference or inp/rel,keys,fields))
    for module in ['innate','adaptive']:
        add('original_'+module+'_gene_members','code/01_discovery_signature/locked_'+module+'_gene_set_v1.csv',['gene'],['module'])
    add('discovery_15plus15_scores','revision_R1_20261007/tables/discovery_locked_scores_descriptive.csv',['sample'],['innate_score','adaptive_score','balance_score'])
    add('TMM_normalization','revision_R1_20261007/tables/NCBI_GSE226238_normalization.csv',['gsm'],['lib.size','norm.factors'])
    add('longitudinal_15plus15_scores','revision_R1_20261007/tables/NCBI_longitudinal_full_scores.csv',['sample'],['innate_score','adaptive_score','balance_score'])
    add('16_BTM_score_series','revision_R1_20261007/tables/external_BTM_scores.csv',['dataset','sample','variant','module_id'],['score'])
    add('BTM_gene_membership','revision_R1_20261007/tables/external_BTM_used_genes.csv',['dataset','variant','module_id','gene'],[])
    add('original_score_selection_sensitivity','revision_R1_20261007/tables/processed_matrix_topN_B_sensitivity.csv',['gene','module','top_n','iterations'],['selection_frequency','direction_consistency_500'])
    add('ODC_clinical_mapping','revision_R2_20261008/tables/ODC405_sample_crosswalk.csv',['sample'],['group','age','sex','iss','tbi','prior_cns'])
    add('fixed_sample_27_models','revision_R2_20261008/fixed_sample_20261012/tables/fixed_sample_models.csv',['subset','spec','outcome'],['estimate','lower_CL','upper_CL','p_value','robust_SE_HC3'])
    add('BTM_88_clinical_models','revision_external_adjusted_20261016/tables/discovery_models.csv',['subset','spec','variant','module_id'],['estimate','SE','lower_CL','upper_CL','p_value','q_BH_within_spec','q_BH_all88','q_BH_primary16'])
    add('BTM_240_longitudinal_contrasts','revision_external_adjusted_20261016/tables/longitudinal_models.csv',['scope','spec','variant','module_id','contrast'],['estimate','SE','lower_CL','upper_CL','p_value','q_BH_within_spec','q_BH_all240','q_BH_primary16'])
    add('original_score_longitudinal_sensitivities','revision_R2_20261008/tables/R2_longitudinal_sensitivities.csv',['spec','score','contrast'],['estimate','lower.CL','upper.CL','p.value','q_BH'])
    add('same_platform_reuse_sensitivity','revision_R2_20261008/tables/same_platform_sensitivity.csv',['scope','score'],['g','g_ci_low','g_ci_high','p.value','q_BH'])
    for model in ['M0','M1','M2']:
        for kind,key,fields in [('gene','gene',['logFC','CI.L','CI.R','t','P.Value','adj.P.Val']),('camera','go_id',['NGenes','Direction','PValue','FDR'])]:
            name=f'{kind}_models_{model}.csv'
            add(kind+'_'+model,'revision_followup_20261019/bridge43/results/'+name,[key],fields,ref/'bridge43'/name)
    add('main_Figure1_four_effects','revision_Day3_20261021/main_figures/Figure1_BTM_supported4_source_data.csv',['module_id','contrast'],['estimate','lower_CL','upper_CL','p_value','q_BH_primary16'],ref/'manuscript/main_figures/Figure1_BTM_supported4_source_data.csv')
    results=[]
    for name,work,reference,keys,fields in tests:
        try:r=compare(work,reference,keys,fields)
        except Exception as exc:r={'status':'FAIL','error':str(exc)}
        r.update(test=name,regenerated_file=work.relative_to(args.work).as_posix(),reference_file=reference.relative_to(BASE).as_posix());results.append(r);print(r['status'],name,flush=True)
    out={'all_passed':all(x['status']=='PASS' for x in results),'tests':results,'tolerance':'absolute 1e-8 + relative 1e-8 times reference magnitude','scope':'Rebuilt matrix-based standard pipeline; historical original selection, clustering and preranked GO are separate optional steps.'}
    args.report.write_text(json.dumps(out,indent=2),encoding='utf-8')
    raise SystemExit(0 if out['all_passed'] else 1)
if __name__=='__main__':main()
