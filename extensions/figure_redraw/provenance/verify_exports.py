"""Non-visual export audit only. All graphics and previews are rendered in R."""
from pathlib import Path
import csv, hashlib, json, re, sys
import pymupdf as fitz
from PIL import Image
import pandas as pd

out=Path(__file__).resolve().parent
project=out.parent.parent
qa=out/'qa'
manifest=pd.read_csv(qa/'export_manifest.csv')
assert len(manifest)==7 and set(manifest.number)==set(range(1,8))
records=[]
for row in manifest.itertuples():
    p=out/(row.figure+'.pdf');t=out/(row.figure+'.tif')
    svg=out/'svg'/(row.figure+'.svg');png=out/'previews'/(row.figure+'.png')
    assert all(x.exists() and x.stat().st_size>1000 for x in [p,t,svg,png])
    with fitz.open(p) as doc:
        assert len(doc)==1
        page=doc[0]
        width=page.rect.width/72*25.4;height=page.rect.height/72*25.4
        # Cairo PDF pages round dimensions to integer points (~0.35 mm here).
        assert abs(width-170)<.5 and width<=170
        assert abs(height-row.height_mm)<.5
        images=page.get_images(full=True)
        assert len(images)==0,(p.name,images)
        text=page.get_text()
        assert len(text)>100
        fonts=[]
        for f in page.get_fonts(full=True):
            blob=doc.extract_font(f[0])[3]
            fonts.append({'name':f[3],'embedded':bool(blob)})
        assert fonts and all(f['embedded'] for f in fonts)
        spans=[s for b in page.get_text('dict')['blocks'] if b.get('type')==0 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside=[s['text'] for s in spans if s['bbox'][0]<-.5 or s['bbox'][1]<-.5 or s['bbox'][2]>page.rect.width+.5 or s['bbox'][3]>page.rect.height+.5]
        assert not outside,(p.name,outside)
        min_font=min(s['size'] for s in spans)
        assert min_font>=6.9
        (qa/(row.figure+'_pdf_text.txt')).write_text(text,encoding='utf-8')
    # Only read image metadata; no rasterization, resizing or image writing.
    with Image.open(t) as im:
        mode=im.mode;pixels=im.size;dpi=im.info.get('dpi');compression=im.info.get('compression')
        assert mode=='RGB' and dpi and abs(dpi[0]-600)<.01 and abs(dpi[1]-600)<.01
        assert compression=='tiff_lzw'
        assert abs(pixels[0]-170/25.4*600)<2
        assert abs(pixels[1]-row.height_mm/25.4*600)<2
    with Image.open(png) as im:
        assert abs(im.info['dpi'][0]-300)<.1
    svg_text=svg.read_text(encoding='utf-8')
    assert '<text' in svg_text and '<image' not in svg_text
    records.append({'number':row.number,'figure':row.figure,'width_mm':width,'height_mm':height,
      'pdf_pages':1,'pdf_raster_objects':0,'all_fonts_embedded':True,'minimum_font_pt':min_font,
      'text_outside_page':0,'tiff_width_px':pixels[0],'tiff_height_px':pixels[1],'tiff_dpi':float(dpi[0]),
      'tiff_mode':mode,'tiff_compression':compression,'pdf_bytes':p.stat().st_size,'tif_bytes':t.stat().st_size})

legacy=pd.read_csv(qa/'legacy_original_sha256.csv')
for row in legacy.itertuples():
    p=project/'figures'/row.file
    assert hashlib.sha256(p.read_bytes()).hexdigest()==row.sha256
assert len(legacy)==14

src=out/'source_data'
counts=pd.read_csv(src/'Figure1_DEG_counts.csv')
assert counts.set_index(['comparison','category']).loc[('SCI_vs_HC','Up'),'n']==1150
assert counts.set_index(['comparison','category']).loc[('TC_vs_HC','Up'),'n']==127
assert all(counts.groupby('comparison').n.sum()==17498)
go=pd.read_csv(src/'Figure2_fixed_eight_GO.csv')
assert len(go)==24 and go.go_id.nunique()==8 and set(go.model)=={'M0','M1','M2'}
hm=pd.read_csv(src/'Figure3_heatmap.csv');genes=pd.read_csv(src/'Frozen_30_gene_membership.csv')
assert len(hm)==1740 and hm.gene.nunique()==30 and hm['sample'].nunique()==58
assert set(hm.gene)==set(genes.gene)
boot=pd.read_csv(src/'Figure4_cluster_bootstrap_summary.csv');assert len(boot)==20
ls=pd.read_csv(src/'Figure5_longitudinal_scores.csv')
assert len(ls)==37 and ls.subject.nunique()==10
assert ls.time.value_counts().to_dict()=={'acute':10,'3m':10,'6m':10,'12m':7}
lc=pd.read_csv(src/'Figure5_saved_mixed_model_contrasts.csv')
assert len(lc)==9 and set(lc.n_participants)=={10} and set(lc.n_samples)=={37}
# BH9 independently reconstructed from the saved nine P values.
p=lc['p.value'];order=p.sort_values().index.tolist();q=[0.0]*len(p);running=1.0
for rank in range(len(order),0,-1):
    idx=order[rank-1];running=min(running,float(p[idx])*len(p)/rank);q[idx]=running
assert max(abs(lc.q_BH-pd.Series(q)))<1e-12
dc=pd.read_csv(src/'Figure6_frozen_discovery_scores.csv');assert len(dc)==58
ex=pd.read_csv(src/'Figure7_additional_cohort_scores.csv')
assert len(ex)==45 and ex['sample'].nunique()==45 and ex.subject_id.nunique()==45
assert sum((ex.group=='SCI')&(ex.overlap=='reused'))==10
assert sum((ex.group=='CTL')&(ex.overlap=='reused'))==9
ef=pd.read_csv(src/'Figure7_reuse_platform_effects.csv');assert len(ef)==12
assert set(ef.n_control)=={9}
assert all(ef.CI_low<ef.effect_g) and all(ef.CI_high>ef.effect_g)
preflight=json.loads((qa/'source_preflight.json').read_text(encoding='utf-8-sig'))
assert preflight['summary']['counts']['FAIL']==0
for name in ['README.md','Figure_legends.md','figure-contract.md','make_figures.R','prepare_source_data.R']:
    text=(out/name).read_text(encoding='utf-8')
    assert '\ufffd' not in text
pd.DataFrame(records).to_csv(qa/'export_verification.csv',index=False,encoding='utf-8-sig')
result={'status':'PASS','date':'2026-10-28','figures':7,'PDFs':7,'TIFFs':7,'SVGs':7,'R_PNG_previews':7,
 'width_cap_mm':170,'all_PDFs_vector_no_raster_objects':True,'all_fonts_embedded':True,
 'all_TIFFs_RGB_600dpi_LZW':True,'source_counts_and_BH9_checked':True,
 'all_14_legacy_files_unchanged':True,
 'source_preflight_failures':0,'source_preflight_warnings_reviewed':[
 'R delimiter warning resolved by successful full R execution',
 'Static width detector does not infer variable W; actual PDF widths are 169.686 mm, within the 170 mm design cap'],
 'visual_review':'All seven complete R-rendered previews inspected; no clipped labels or overlapping text found.',
 'remaining_boundaries':['No integration into manuscript numbering in this task','Saved biological models not refit',
 'Clustering diagnostics condition on saved 2509-gene universe',
 'Additional-cohort full score transform not independently reconstructed'],
 'outputs':records}
(qa/'final_qa.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
status=json.loads((qa/'render_status.json').read_text(encoding='utf-8'))
status['status']='PASS_RENDER_DATA_AND_VISUAL_QA'
(qa/'render_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# 绘图交付核验 — 2026-10-28','',
 '七幅图完成，全部以 R 原生绘制；旧图14个文件的SHA256均未改变。','',
 '| 图 | 宽 × 高 mm | PDF嵌入栅格数 | 最小文字 pt | TIFF |',
 '| --- | --- | --- | --- | --- |']
for r in records:
    lines.append(f"| Fig. {r['number']} | {r['width_mm']:.1f} × {r['height_mm']:.1f} | 0 | {r['minimum_font_pt']:.1f} | RGB / 600 dpi / LZW |")
lines+=['','逐图视觉检查已完成；字体嵌入、可提取文字、图面边界和全部源数据行数通过核验。静态脚本预检0项失败；两项提示分别以实际R运行及PDF页面尺寸核验解决。',
'','全体观察、已保存模型区间和完整指定比较均保留；热图只做颜色饱和，点图横向jitter只用于显示。没有重新拟合生物模型，也没有引入旧5+4分支。',
'','R绘图宽度设置为170 mm；Cairo PDF按整数point取整，实际宽169.686 mm。TIFF宽4015像素、600 dpi，约169.97 mm；两种格式均未超过170 mm。600 dpi为本包导出规格，不代表另行核实了期刊当日指南。',
'','正文和回复信仍需在采用这套图时同步编号；本次没有修改这两个文档。']
(qa/'QA_中文.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
files=[p for p in out.rglob('*') if p.is_file() and p.name!='SHA256SUMS.csv']
with (qa/'SHA256SUMS.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.writer(f);w.writerow(['file','sha256','bytes'])
    for p in sorted(files):
        w.writerow([p.relative_to(out).as_posix(),hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_size])
print(json.dumps({k:v for k,v in result.items() if k!='outputs'},ensure_ascii=False,indent=2))

