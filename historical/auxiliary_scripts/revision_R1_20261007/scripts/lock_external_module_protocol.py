# Public-package copy: see metadata/code_changes.csv for provenance.
"""Lock external definitions and revision-stage analysis before computing BTM outcomes."""
from pathlib import Path
import csv,json,hashlib,re,datetime
O=Path(__file__).resolve().parents[1];D=O/'inputs/external_modules'
lock=O/'tables/external_module_protocol_lock.json'
if lock.exists():raise SystemExit('Protocol already locked; not overwritten.')
gmt=D/'BTM_for_GSEA_20131008.gmt';commit=json.loads((D/'btm_tree.txt').read_text())['sha']
selected=[];allmodules=[]
for line in gmt.read_text(encoding='utf-8-sig').splitlines():
    fields=line.split('\t');title=fields[0];mid=re.search(r'\(([^()]+)\)$',title).group(1)
    selected_flag=bool(re.match(r'^enriched in (neutrophils|T cells) \(',title))
    allmodules.append(dict(module_id=mid,title=title,n_source_tokens=len(set(fields[2:])),selected=selected_flag,
                           reason='literal cell-enrichment title rule' if selected_flag else 'outside locked literal title rule'))
    if selected_flag:
        lineage='neutrophil' if 'neutrophils' in title else 'T_cell'
        selected.extend(dict(module_id=mid,title=title,lineage=lineage,source_token=g,source_url=fields[1]) for g in dict.fromkeys(fields[2:]))
assert set(r['module_id'] for r in selected)=={'M7.0','M223','M37.1','M163'}
for name,rows in [('external_BTM_source_members',selected),('external_BTM_module_selection',allmodules)]:
    with (O/'tables'/f'{name}.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
stamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
protocol=f'''# R1-11：外部血液转录模块补充分析方案（锁定版）

锁定UTC：{stamp}。这是本次修回新制定的方案，锁定发生在本次BTM评分与效应计算之前；此前已知原论文及既有修回结果，因此不是盲法、原研究事前方案或公开预注册。R1-10保持冻结。

## 目的与来源

检验原报告的中性粒细胞相关升高/T细胞相关降低及随访方向，是否同样出现在外部定义的模块中。保留原S14 15+15签名，不重新选取或优化。外部定义不等于独立人群验证，不消除细胞组成或临床混杂。

采用Li等发表于Nature Immunology 15,195–204 (2014)的BTM资源，DOI 10.1038/ni.2789，PMID24336226，PMC3946932（在线2013-12-15，正式刊期2014）。已读公开全文的模块构建和评分段落：整合超过30,000个人血液转录组，模块均值作为activity score。它是血液共表达/细胞富集模块，不是SCI功能恢复或细胞比例测量工具。论文的疫苗实验含PBMC，不能将其称为专门验证过的SCI全血分类器。

来源为作者仓库https://github.com/shuzhao-li/BTM，固定commit `{commit}` 下 `BTM/datasets/BTM_for_GSEA_20131008.gmt`（2013原始版，不使用2021 BTM_Plus）。该GMT有346条，包括M和S模块；论文正文报告334 BTMs，故明确具体文件而不混同数目。全文、README、许可、文件和hash本地保存。

按源文件标题固定选择所有以`enriched in neutrophils (`或`enriched in T cells (`开头的模块，得到M37.1、M163、M7.0、M223，全部等同报告。该规则针对原稿细胞相关表达主题，避免把激活/分化/细胞周期模块混为细胞富集指标。不按本项目效应/P值挑选模块，不扩大为全部先天或适应性免疫的功能判断。全346条纳排表保存。

## 基因处理、覆盖与两种分析

1. 原GMT每模块唯一token作为来源记录。含`///`的多基因token不拆分，作为歧义未映射；先匹配AnnotationDbi/org.Hs.eg.db的SYMBOL，其余仅在ALIAS可唯一映射到一个SYMBOL时转换。无法唯一映射保留缺失，不按结果人工修复。记录软件版本、所有映射与重复合并。
2. 对发现矩阵行名采用同一规则；多个行归一到一个symbol时在既有log表达尺度取均值，不重新做ComBat或差异筛选。使用既有S1的17,498基因分析范围。不得按显著性/方向过滤模块成员。
3. 对GSE226238保留既有全矩阵TMM归一化因子。使用官方GeneID注释，多个GeneID指向同一symbol时先合并原始counts，再用既有effective library size计算log2CPM，prior.count=0.5。不重新比对、不重新选择基因或重估归一化因子。
4. 主分析使用四个外部模块的全部可映射基因。可分析须至少5个有效唯一基因且覆盖≥80%；覆盖分母为源GMT唯一token数，歧义/未映射/映射折叠不增加分子。不按组差异或表达高低过滤，不插补缺失；模块不达标则原样报告失败，不换模块或放宽规则。
5. 另固定进行一次“排除原106候选基因”敏感性：以同一符号映射后的原候选池并集，从每个外部模块中排除重叠基因，再算分数；其覆盖分母扣除对应已明确映射的重叠source tokens，其他规则不变。此为BTM衍生集合，不能称为原文原样模块。即使集合太小也报告，不重新补基因。不只排除最终30基因，以检验依赖整个原候选池的程度。
6. 每模块评分为所用基因的既有log尺度表达算术均值，参考原文module activity方法；不使用对照组标准差重新缩放，不作模块间绝对丰度比较，不计算innate-minus-T balance。

## 样本、模型和报告

- GSE151371：既有处理后矩阵，SCI38与TC10，Welch双侧t检验，估计SCI−TC模块均值差及95%CI；独立单位为参与者。HC不进入这项主比较。保留全部SCI/TC，无新排除。
- GSE226238：既有同一队列，SCI10人37次观察（10/10/10/7），对每模块拟合score ~ factor(time) + (1|subject)，REML，Satterthwaite自由度；报告3/6/12月−急性期的估计、点态95%CI和双侧P。九名CTL不进入随访效应模型，表达导出仍保留全部46样本。不填补缺失随访。
- 主分析4个发现比较+12个随访比较为一族，BH按计划16个检验校正；未能估计者保留NA，校正n仍为16。排除106候选的16个比较作为另一个明确的敏感性族。完整7人28样本的12个随访比较作为第三个敏感性族，BH n=12。报告所有结果，不只报12月或显著项。
- 检查拟合收敛、随机效应奇异性、固定设计秩、残差及方差异常，保存诊断。奇异/失败模型不悄悄改公式；标为模型限制，预定完整病例结果作辅助，不据此挑选模型。
- 仅做上述表达模块对比，不加入临床协变量或细胞比例调整。新模块不用于“校正”原签名，不推进R1-10。没有可供本方案重算的独立GSE293559全表达矩阵时，不从S16的原签名分数伪造新模块分数。
- 判读：中性粒模块SCI−TC正、随访负；T细胞模块SCI−TC负、随访正为方向一致。分别报告方向、CI、BH结果；相反或不显著结果保留。不能靠多数显著模块投票证明机制，也不能视重叠模块为独立生物学重复。

## 预定输出与图合同

输出来源/选择/映射/coverage/原候选重叠/逐样本评分/全部效应与诊断CSV、R代码与sessionInfo、分析报告、R1-11补充回复和稿件段落。R制作183mm宽的补充森林图（主分析四模块×四对比），显示全部估计和95%CI，PDF/SVG可编辑文字、600dpi TIFF、PNG预览。图的任务是展示一致与不一致，不预设成功结论。原有主图和原签名不替换。
'''
p=O/'R1_11_external_module_protocol.md';p.write_text(protocol,encoding='utf-8')
tracked=[p,gmt,O/'tables/external_BTM_source_members.csv',O/'tables/external_BTM_module_selection.csv']
j={'locked_at_utc':stamp,'stage':'before new BTM score/effect calculation; previous study results known','commit':commit,
   'source_doi':'10.1038/ni.2789','modules':['M7.0','M223','M37.1','M163'],'files':[]}
for f in tracked:j['files'].append({'path':str(f.resolve()),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
lock.write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8')
print('Protocol locked before BTM outcomes:',stamp, j['modules'])
