# Public-package copy: see metadata/code_changes.csv for provenance.
from pathlib import Path
import pandas as pd,numpy as np,json,sys
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[2];O=ROOT/'revision_R2_20261008';R1=ROOT/'revision_R1_20261007'
m=pd.read_csv(R1/'tables/sample_manifest.csv');m['clinical_source']='GEO annotations; missing fields not imputed'
tb=pd.read_csv(R1/'inputs/followup/clinical_public_search/GSE226238_published_tb1.csv').set_index('ID')
check=pd.read_csv(R1/'inputs/followup/clinical_public_search/GSE226238_clinical_mapping_precheck.csv').set_index('gsm')
for i,row in m.iterrows():
    gsm=row.gsm if row.dataset=='GSE226238' else row.original_gsm if row.dataset=='GSE293559' and row.overlap=='reused_from_GSE226238' else None
    if pd.notna(gsm) and gsm in check.index:
        rr=check.loc[gsm]; a=tb.loc[int(rr.article_id)]
        assert a['Group (SCI or CTL)']==row.group
        m.loc[i,['age','sex']]=[a.Age,a.Gender]
        m.loc[i,'clinical_source']='Morrison 2023 published Table 1; verified GEO GSM/participant mapping'
        m.loc[i,'baseline_ais']=a['AIS grade (acute)'];m.loc[i,'baseline_nli']=a['NLI (acute)'];m.loc[i,'mechanism']=a.MOI
        if row.group=='SCI':m.loc[i,['ais','nli']]=[rr.ais,rr.nli]
    m.loc[i,'analysis_role']='metadata audit only; excluded from revised signature inference' if row.dataset=='GSE177034' else 'main analysis'
    # R1 treatment is a visit label, not treatment information.
    if row.dataset=='GSE226238' and pd.notna(row.treatment):
        m.loc[i,'original_treatment_annotation']=row.treatment;m.loc[i,'treatment']=np.nan
raw=[]
for ds in m.dataset.unique():
    f=R1/f'inputs/geo/{ds}_samples.json'
    if not f.exists():continue
    for rec in json.loads(f.read_text(encoding='utf-8')):
        raw.append({'dataset':ds,'gsm':rec['gsm'],**{k:' | '.join(v) if isinstance(v,list) else v for k,v in rec.items() if k.startswith('Sample_')}})
raw=pd.DataFrame(raw);raw.to_csv(O/'tables/GEO_metadata_full_audit.csv',index=False)
fields=['Sample_library_strategy','Sample_library_selection','Sample_library_source','Sample_extract_protocol_ch1','Sample_data_processing','Sample_submission_date']
for i,row in m.iterrows():
    rr=raw[raw.gsm.eq(row.gsm)]
    if rr.empty and pd.notna(row.original_gsm):rr=raw[raw.gsm.eq(row.original_gsm)]
    for f in fields:
        if f in rr and len(rr):m.loc[i,f]=rr.iloc[0][f]
    m.loc[i,'processing_date']=np.nan
m['sex']=m['sex'].replace({'Male':'M','Female':'F','male':'M','female':'F'})
m.to_csv(O/'tables/sample_manifest_R2.csv',index=False)
avail=[]
for ds,dd in m.groupby('dataset'):
    for col in ['age','sex','ais','nli','injury_time_hours','infection','treatment','batch','rin','platform','instrument','Sample_library_strategy','Sample_library_selection','Sample_library_source','processing_date']:
        x=dd[col] if col in dd else pd.Series(dtype=object)
        avail.append(dict(dataset=ds,variable=col,n_samples=len(dd),observed=int(x.notna().sum()),unique=int(x.nunique()),note='Submission date is not a processing date; repeated rows are not independent subjects'))
pd.DataFrame(avail).to_csv(O/'tables/covariate_availability_R2.csv',index=False)
rows=[]
for (ds,group),dd in m.groupby(['dataset','group']):
    participants=dd.drop_duplicates('subject');age=pd.to_numeric(participants.age,errors='coerce')
    rows.append(dict(dataset=ds,group=group,participants=participants.subject.nunique(),samples=len(dd),platform=';'.join(dd.platform.unique()),material=';'.join(dd.material.dropna().unique()),age_n=int(age.notna().sum()),age_mean=age.mean(),age_sd=age.std(),sex='; '.join(f'{k}:{v}' for k,v in participants.sex.value_counts().items()),AIS='; '.join(f'{k}:{v}' for k,v in participants.ais.value_counts().items()),NLI='; '.join(f'{k}:{v}' for k,v in participants.nli.value_counts().items()),visits=';'.join(dd.visit.dropna().astype(str).unique()),repeated_observations=len(dd)>len(participants),overlap='; '.join(f'{k}:{v}' for k,v in dd.overlap.value_counts().items()),infection='not sample-linked',treatment='not sample-linked',technical_QC='RIN/batch/processing date not sample-linked',role=dd.analysis_role.iloc[0]))
pd.DataFrame(rows).to_csv(O/'tables/Table1_cohort_characteristics.csv',index=False)
print(pd.DataFrame(rows)[['dataset','group','participants','samples','age_n','age_mean','sex']].to_string(index=False))
