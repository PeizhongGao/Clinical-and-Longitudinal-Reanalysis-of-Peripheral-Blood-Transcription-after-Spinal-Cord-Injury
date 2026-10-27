"""Collect regenerated manuscript figures/tables without changing their numbers."""
from pathlib import Path
import argparse,csv,shutil
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--work',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    figures={
      'Figure1_BTM_supported4':'revision_Day3_20261021/main_figures/Figure1_BTM_supported4',
      'Figure_S1_Discovery_heatmap':'revision_R1_20261007/figures/Figure3_revised',
      'Figure_S2_Original_score_longitudinal':'revision_R1_20261007/figures/Figure5_revised',
      'Figure_S3_Original_score_fixed_models':'revision_R2_20261008/fixed_sample_20261012/figures/figure01_fixed_sample_effects',
      'Figure_S4_Clustering_diagnostics':'revision_R2_20261008/figures/Figure4_revised',
      'Figure_S5_Clinical_support':'revision_R2_20261008/fixed_sample_20261012/figures/figure03_covariate_support',
      'Figure_S6_BTM_composition_sensitivity':'revision_external_adjusted_20261016/figures/figure02_external_composition',
      'Figure_S7_Additional_cohort_reuse':'revision_R2_20261008/figures/Figure7_revised',
      'Figure_S8_Fixed_immune_GO':'revision_numeric_freeze_20261019/figures/figure01_fixed_immune_GO',
      'Figure_S9_BTM_primary16_all':'revision_numeric_freeze_20261019/figures/figure02_BTM_primary16',
    }
    rows=[];dest=args.output/'figures';dest.mkdir(parents=True,exist_ok=True)
    for name,stem in figures.items():
        found=False
        for ext in ['.png','.pdf','.svg','.tiff']:
            src=args.work/(stem+ext)
            if src.is_file():
                shutil.copy2(src,dest/(name+ext));rows.append({'item':name,'work_source':stem+ext,'exported':str(Path('figures')/(name+ext))});found=True
        if not found:raise FileNotFoundError(stem)
    for stage,relative in [('BTM','revision_external_adjusted_20261016/tables'),('bridge43','revision_followup_20261019/bridge43/results'),('original_fixed_sample','revision_R2_20261008/fixed_sample_20261012/tables'),('final_primary','revision_numeric_freeze_20261019/tables')]:
        d=args.output/'tables'/stage;d.mkdir(parents=True,exist_ok=True)
        for src in (args.work/relative).glob('*.csv'):
            shutil.copy2(src,d/src.name)
    with (args.output/'figure_index.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    (args.output/'README.txt').write_text('Generated files from the requested run. Figure S4 uses archived clustering checkpoints in standard mode; extended mode recomputes clustering. Source study CIBERSORTx values are reused, not re-estimated. Figure names follow the revision package; compare numeric results using validation.json.\n',encoding='utf-8')
if __name__=='__main__':main()
