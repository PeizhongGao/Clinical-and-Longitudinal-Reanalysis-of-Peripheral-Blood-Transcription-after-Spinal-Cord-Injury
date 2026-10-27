#!/usr/bin/env python3
"""Run the packaged SCI analyses in a new working directory, never in inputs/."""
from pathlib import Path
import argparse, csv, hashlib, json, os, shutil, subprocess, sys, time

BASE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def copyfiles(src,dst):
    dst.mkdir(parents=True,exist_ok=True)
    for p in src.glob('*'):
        if p.is_file():shutil.copy2(p,dst/p.name)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True,help='New, empty execution directory; reference inputs are never overwritten.')
    ap.add_argument('--rscript',default=shutil.which('Rscript'),help='Rscript executable, or put Rscript on PATH.')
    ap.add_argument('--profile',choices=['standard','extended'],default='standard',help='extended also recomputes historical clustering and preranked gene-set permutations')
    ap.add_argument('--original-selection',action='store_true',help='Also rerun the original count-based 15+15 selection block (exploratory, not nested validation).')
    ap.add_argument('--resume',action='store_true',help='Resume only steps from this same package whose completed scripts and outputs have unchanged hashes.')
    args=ap.parse_args()
    if not args.rscript:ap.error('Rscript not found. Supply --rscript or install R and add it to PATH.')
    output=args.output.expanduser().resolve()
    if output == BASE or (BASE in output.parents and not ((BASE/'outputs') == output or (BASE/'outputs') in output.parents)):
        ap.error('Inside the package, use outputs/ only; otherwise choose an external output directory.')
    if output in BASE.parents:
        ap.error('The output directory cannot be a parent of the release package.')
    if output.exists() and any(output.iterdir()) and not args.resume:ap.error('Output is not empty. Use a new directory or --resume for this run.')
    output.mkdir(parents=True,exist_ok=True)
    logs=output/'logs';logs.mkdir(exist_ok=True)
    work=output/'work';statuspath=output/'run_status.json'
    status=json.loads(statuspath.read_text()) if args.resume and statuspath.exists() else {'package_version':'1.0.0','profile':args.profile,'original_selection':args.original_selection,'steps':[]}
    def save():statuspath.write_text(json.dumps(status,indent=2),encoding='utf-8')
    if args.resume and status.get('profile')!=args.profile:ap.error('Resume profile must match the original run.')
    if args.resume and status.get('original_selection',False)!=args.original_selection:ap.error('Resume must use the same original-selection setting.')
    execution_files=[BASE/'reproduce.py']+[p for part in ['code','inputs'] for p in (BASE/part).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    fingerprint=hashlib.sha256('\n'.join(p.relative_to(BASE).as_posix()+':'+sha(p) for p in sorted(execution_files)).encode()).hexdigest()
    if args.resume and status.get('execution_fingerprint',fingerprint)!=fingerprint:ap.error('Package code or inputs changed; use a fresh output directory.')
    status['execution_fingerprint']=fingerprint;status['completed']=False;save()
    if not (output/'STAGED.json').exists():
        shutil.copytree(BASE/'inputs/project',work,dirs_exist_ok=True)
        shutil.copytree(BASE/'code/project',work,dirs_exist_ok=True)
        for branch in ['revision_R1_20261007','revision_R2_20261008','revision_integrated_20261012','revision_external_adjusted_20261016','revision_scope_lock_20261018','revision_numeric_freeze_20261019','revision_followup_20261019','revision_R2_20261008/fixed_sample_20261012']:
            for name in ['tables','logs','figures','supplement','inputs']:(work/branch/name).mkdir(parents=True,exist_ok=True)
        (output/'STAGED.json').write_text(json.dumps({'source':'package-relative inputs/project and code/project','timestamp':time.time()}),encoding='utf-8')
    reran_step=False
    def stage(name,script,expected,language='R',cwd=None):
        nonlocal reran_step
        scriptpath=BASE/script if script.startswith('code/') else work/script
        previous=next((s for s in status['steps'] if s['name']==name and s['status']=='PASS'),None)
        if previous and not reran_step:
            valid=previous.get('script_sha256')==sha(scriptpath) and all((work/x).is_file() and sha(work/x)==h for x,h in previous.get('output_sha256',{}).items())
            if valid:
                print('RESUME',name,flush=True);return
        reran_step=True
        # Expected output anchors must be regenerated, not inherited from checkpoints.
        for rel in expected:
            p=work/rel
            if p.exists():p.unlink()
        cmd=([args.rscript,'--vanilla',str(scriptpath),str(work)] if language=='R' else [sys.executable,'-X','utf8',str(scriptpath)])
        if name=='original_selection':cmd=[args.rscript,'--vanilla',str(scriptpath)]
        logfile=logs/(name+'.log');start=time.time();print('START',name,flush=True)
        env=os.environ.copy();env['PYTHONUTF8']='1'
        local_library=BASE/'environment/r_library'
        if local_library.is_dir():env['R_LIBS_USER']=str(local_library)
        with logfile.open('w',encoding='utf-8') as f:
            proc=subprocess.run(cmd,cwd=cwd or work,stdout=f,stderr=subprocess.STDOUT,env=env)
        missing=[x for x in expected if not (work/x).is_file()]
        item={'name':name,'script':script,'script_sha256':sha(scriptpath),'seconds':round(time.time()-start,2),'exit_code':proc.returncode,'missing_outputs':missing,'status':'PASS' if proc.returncode==0 and not missing else 'FAIL','output_sha256':{x:sha(work/x) for x in expected if (work/x).is_file()}}
        status['steps']=[x for x in status['steps'] if x['name']!=name]+[item];save()
        print(item['status'],name,item['seconds'],'seconds',flush=True)
        if item['status']!='PASS':
            print(logfile.read_text(encoding='utf-8',errors='replace')[-6000:]);raise SystemExit('Stopped; see '+str(logfile))

    stage('environment','code/check_environment.R',[])
    stage('metadata_R1','revision_R1_20261007/scripts/prepare_metadata.py',['revision_R1_20261007/tables/sample_manifest.csv'],'Python')
    if args.original_selection:
        stage('original_selection','code/original_signature_selection.R',['code/01_discovery_signature/locked_innate_gene_set_v1.csv','code/01_discovery_signature/locked_adaptive_gene_set_v1.csv'],'R',work/'code')
    stage('discovery_and_original_scores','revision_R1_20261007/scripts/run_R1.R',['revision_R1_20261007/tables/analysis_state.rds','revision_R1_20261007/tables/locked_signature_reused_S14.csv','revision_R1_20261007/tables/discovery_locked_scores_descriptive.csv','revision_R1_20261007/tables/processed_matrix_topN_B_sensitivity.csv','revision_R1_20261007/tables/SCI_vs_TC_GO_camera.csv'])
    stage('count_normalization_scores','revision_R1_20261007/scripts/NCBI_count_followup.R',['revision_R1_20261007/tables/NCBI_longitudinal_full_scores.csv','revision_R1_20261007/tables/NCBI_longitudinal_mixed_contrasts.csv'])
    stage('external_BTM_scores','revision_R1_20261007/scripts/external_BTM_analysis.R',['revision_R1_20261007/tables/external_BTM_scores.csv','revision_R1_20261007/tables/external_BTM_used_genes.csv'])
    # Explicitly wire newly regenerated upstream tables into downstream stages.
    for name in ['discovery_locked_scores_descriptive.csv','locked_signature_reused_S14.csv','standardization_parameters.csv']:
        shutil.copy2(work/'revision_R1_20261007/tables'/name,work/'revision_R2_20261008/supplement'/name)
        shutil.copy2(work/'revision_R1_20261007/tables'/name,work/'revision_integrated_20261012/supplement'/name)
    for p in (work/'revision_R1_20261007/tables').glob('external_BTM_*.csv'):
        shutil.copy2(p,work/'revision_integrated_20261012/supplement'/p.name)
    stage('metadata_R2','revision_R2_20261008/scripts/metadata_r2.py',['revision_R2_20261008/tables/sample_manifest_R2.csv'],'Python')
    for name in ['sample_manifest_R2.csv','GEO_metadata_full_audit.csv','covariate_availability_R2.csv','Table1_cohort_characteristics.csv']:
        p=work/'revision_R2_20261008/tables'/name
        if p.exists():shutil.copy2(p,work/'revision_integrated_20261012/tables'/name)
    stage('clinical_ODC405','revision_R2_20261008/scripts/odc405_adjusted_analysis.R',['revision_R2_20261008/tables/ODC405_sample_crosswalk.csv','revision_R2_20261008/tables/ODC405_locked_scores_with_clinical_and_cell_estimates.csv','revision_R2_20261008/supplement/ODC405_locked_score_models.csv'])
    for name in ['ODC405_sample_crosswalk.csv','ODC405_locked_scores_with_clinical_and_cell_estimates.csv']:
        src=work/'revision_R2_20261008/tables'/name
        if src.exists():shutil.copy2(src,work/'revision_integrated_20261012/tables'/name)
    for name in ['discovery_locked_scores_descriptive.csv','ODC405_locked_score_models.csv','standardization_parameters.csv']:
        src=work/'revision_R2_20261008/supplement'/name
        if src.exists():shutil.copy2(src,work/'revision_integrated_20261012/supplement'/name)
    stage('original_score_sensitivities','revision_R2_20261008/scripts/sensitivity_analysis.R',['revision_R2_20261008/tables/R2_longitudinal_sensitivities.csv','revision_R2_20261008/tables/GSE226238_official_MCPcounter.csv','revision_R2_20261008/tables/same_platform_sensitivity.csv'])
    stage('fixed_sample_models','revision_integrated_20261012/fixed_sample/scripts/fixed_sample_analysis.R',['revision_R2_20261008/fixed_sample_20261012/tables/fixed_sample_models.csv'])
    copyfiles(work/'revision_R2_20261008/fixed_sample_20261012/tables',work/'revision_integrated_20261012/fixed_sample/tables')
    copyfiles(work/'revision_R2_20261008/fixed_sample_20261012/inputs',work/'revision_integrated_20261012/fixed_sample/inputs')
    shutil.copy2(work/'revision_R2_20261008/tables/GSE226238_official_MCPcounter.csv',work/'revision_integrated_20261012/tables/GSE226238_official_MCPcounter.csv')
    stage('BTM_all_adjusted_models','revision_external_adjusted_20261016/scripts/analyze.R',['revision_external_adjusted_20261016/tables/discovery_models.csv','revision_external_adjusted_20261016/tables/longitudinal_models.csv'])
    stage('bridge_input_reconstruction','code/reconstruct_bridge_inputs.R',[])
    stage('bridge43_genes_GO','revision_numeric_freeze_20261019/scripts/run_bridge43.R',[f'revision_followup_20261019/bridge43/results/{kind}_models_{m}.csv' for kind in ['gene','camera'] for m in ['M0','M1','M2']])
    if args.profile=='extended':
        stage('clustering_diagnostics','revision_R2_20261008/scripts/cluster_validation.R',['revision_R2_20261008/tables/cluster_metrics.csv','revision_R2_20261008/tables/cluster_bootstrap_summary.csv'])
        stage('historical_preranked_GO','revision_R2_20261008/scripts/ranked_enrichment_sparse.py',['revision_R2_20261008/tables/SCI_vs_TC_GSEA.csv'],'Python')
    stage('assemble_final_numeric_tables','code/assemble_numeric.R',['revision_numeric_freeze_20261019/tables/BTM_primary16.csv','revision_numeric_freeze_20261019/tables/bridge43_fixed_eight_GO.csv'])
    stage('main_Figure1','revision_Day3_20261021/scripts/Figure1_BTM_supported4.R',['revision_Day3_20261021/main_figures/Figure1_BTM_supported4.png'])
    stage('supplement_GO_BTM_figures','revision_numeric_freeze_20261019/scripts/figures.R',['revision_numeric_freeze_20261019/figures/figure01_fixed_immune_GO.png','revision_numeric_freeze_20261019/figures/figure02_BTM_primary16.png'])
    stage('supplement_original_longitudinal','revision_R1_20261007/scripts/render_NCBI_figures.R',['revision_R1_20261007/figures/Figure5_revised.png'])
    stage('supplement_fixed_sample','revision_integrated_20261012/fixed_sample/scripts/fixed_sample_figures.R',['revision_R2_20261008/fixed_sample_20261012/figures/figure01_fixed_sample_effects.png','revision_R2_20261008/fixed_sample_20261012/figures/figure03_covariate_support.png'])
    stage('supplement_BTM_sensitivity','revision_external_adjusted_20261016/scripts/figures.R',['revision_external_adjusted_20261016/figures/figure02_external_composition.png'])
    stage('supplement_other_figures','revision_R2_20261008/scripts/render_R2.R',['revision_R2_20261008/figures/Figure4_revised.png','revision_R2_20261008/figures/Figure7_revised.png'])
    report=subprocess.run([sys.executable,'-X','utf8',str(BASE/'code/validate_results.py'),'--work',str(work),'--report',str(output/'validation.json')],cwd=work)
    status['validation_exit_code']=report.returncode;status['completed']=False
    status['not_recomputed']=[] if args.profile=='extended' else ['historical clustering resampling','historical preranked GSEA permutations']
    if not args.original_selection:status['not_recomputed'].append('original count-based signature selection; frozen gene lists and selection-frequency inputs reused')
    save()
    if report.returncode:raise SystemExit(report.returncode)
    exported=subprocess.run([sys.executable,'-X','utf8',str(BASE/'code/export_results.py'),'--work',str(work),'--output',str(output/'results')])
    status['export_exit_code']=exported.returncode;status['completed']=exported.returncode==0;save()
    if exported.returncode:raise SystemExit(exported.returncode)
    print('COMPLETE. Validation:',output/'validation.json',flush=True)

if __name__=='__main__':main()
