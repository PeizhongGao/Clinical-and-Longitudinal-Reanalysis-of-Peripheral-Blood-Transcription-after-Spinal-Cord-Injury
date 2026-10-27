# Public-package copy: see metadata/code_changes.csv for provenance.
from pathlib import Path
import json,csv,re,time,hashlib,sys
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'revision_R1_20261007'
def load(a):return json.loads((OUT/'inputs/geo'/f'{a}_samples.json').read_text(encoding='utf-8'))
def val(z,k):return z.get(k,[''])[0]
def char(z):return {s.split(':',1)[0].strip().lower():s.split(':',1)[1].strip() for s in z.get('Sample_characteristics_ch1',[]) if ':' in s}
def save(rows,path):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
base=[]
for a in ['GSE151371','GSE226238','GSE177034']:
    for z in load(a):
        c=char(z);t=val(z,'Sample_title')
        if a=='GSE151371':
            group={'healthy control':'HC','trauma control':'TC','spinal cord injury':'SCI'}.get(c.get('disease state','').lower())
            if group is None:group='HC' if t.lower().startswith('healthy') else 'TC' if t.lower().startswith('trauma') else 'SCI'
            sample=group+f'{int(re.findall(r"\d+",t)[-1]):02d}';subject=sample;visit='acute' if group!='HC' else 'reference'
        elif a=='GSE226238':
            sample=t.split(',')[0].strip();group=c.get('group','');subject=re.match(r'ID\d+',sample).group();visit=c.get('treatment','')
        else:
            sample=t;subject=c.get('individual','');group=c.get('paingroup','');visit=c.get('timepoint','')
        base.append(dict(dataset=a,sample=sample,gsm=z['gsm'],subject=subject,group=group,visit=visit,age=c.get('age',''),sex=c.get('sex',''),smoking=c.get('smoking',''),ais=c.get('ais',''),nli=c.get('nli',''),injury_time_hours=c.get('time since injury',''),infection=c.get('infection',''),treatment=c.get('treatment',''),batch=c.get('batch',''),rin=c.get('rin',''),platform=val(z,'Sample_platform_id'),instrument=val(z,'Sample_instrument_model'),material=c.get('cell type',c.get('tissue',val(z,'Sample_source_name_ch1'))),include='yes',metadata_source='GEO SOFT retrieved 2026-10-07',mapping_rule='exact GEO title/characteristics; GSE151371 group prefix + title number',overlap='not_assessed_other_cohorts'))
new={val(z,'Sample_title').split(',')[0].strip():z for z in load('GSE293559')}
old={r['sample'].lower():r for r in base if r['dataset']=='GSE226238'}
rows293=[]
with (ROOT/'code/GSE293559_sample_metadata_curated.csv').open(encoding='utf-8-sig') as f:
    for r in csv.DictReader(f):
        sample=r['library_name'];z=new.get(sample);g=char(z) if z else {};orig=old.get(sample.lower());group=r['group']
        if not z and not orig:raise ValueError('No GEO provenance '+sample)
        rows293.append(dict(dataset='GSE293559',sample=sample,gsm=z['gsm'] if z else orig['gsm'],subject=r['subject_id'],group=group,visit='acute' if group=='SCI' else 'reference',age=g.get('age',''),sex=g.get('sex',''),ais=g.get('ais',''),nli=g.get('nli',''),platform=val(z,'Sample_platform_id') if z else orig['platform'],instrument=val(z,'Sample_instrument_model') if z else orig['instrument'],material='Whole blood',include='yes',metadata_source='GEO new-sample SOFT or exact GSE226238 reuse record',mapping_rule='exact library name and GEO series reanalysis relations',overlap='reused_from_GSE226238' if orig else 'new_GSE293559',original_gsm=orig['gsm'] if orig else '',infection='',batch='',rin='',injury_time_hours=''))
save(base+rows293,OUT/'tables/sample_manifest.csv')
print('manifest',len(base+rows293),'dataset counts', {a:sum(r['dataset']==a for r in base+rows293) for a in ['GSE151371','GSE226238','GSE177034','GSE293559']})
# Candidate screening is a revision-stage targeted search, not reconstructed original history.
j=json.loads((OUT/'tables/search_candidates.json').read_text(encoding='utf-8'));screen=[]
for a,z in j.items():
    text=(z.get('title','')+' '+z.get('summary','')).lower()
    if a in ['GSE151371','GSE226238','GSE293559']:decision='included';reason='Human clinical SCI bulk blood dataset; assigned distinct discovery/longitudinal/additional reproduction role'
    elif a=='GSE177034':decision='included_context';reason='Selected non-SCI longitudinal low-back-pain context; not identified as SCI cohort'
    elif a=='GSE309413':decision='excluded';reason='Circulating microRNA/plasma; not bulk blood mRNA transcriptome'
    elif any(x in text for x in ['cerebrospinal','csf','exosome']):decision='excluded';reason='CSF or exosomal material; not peripheral bulk blood mRNA'
    elif any(x in text for x in ['muscle','prostate','epidermis','colorectal','intervertebral disc degeneration','heterotopic ossification']):decision='excluded';reason='Non-blood tissue or a different clinical tissue question'
    else:decision='excluded';reason='GEO title/summary describes cultured, stem/progenitor, organoid, transplanted or nonclinical/non-blood material; record retained for audit'
    screen.append(dict(accession=a,title=z.get('title',''),geo_n_samples=z.get('n_samples',''),query_name=z['query_name'],decision=decision,reason=reason,screening_basis='GEO ESummary title and summary; included datasets also verified with full SOFT',url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+a))
save(screen,OUT/'tables/dataset_screening.csv')
save(json.loads((OUT/'tables/revision_search_log.json').read_text(encoding='utf-8')),OUT/'tables/revision_search_log.csv')

# Public package: pinned marker and GO term files are bundled; no network update.
for name in ['MCPcounter_genes.txt','go_bp_terms.csv']:
    assert (OUT/'inputs'/name).is_file(), name
