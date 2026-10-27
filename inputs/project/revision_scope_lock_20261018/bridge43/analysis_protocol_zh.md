# 43人全基因与GO桥接分析：锁定协议 v1

日期：2026-10-18，香港时间。状态：仅协议、输入快照与设计预检；未计算新的桥接效应、P/q或43人CAMERA。原始结果、48人未调整结果及43人全临床结果已知；这是修回阶段探索性协议，不是事前预注册、盲法或独立验证。

## 研究问题与单位

在完全相同的36名SCI和7名TC中，比较未调整、年龄/性别调整、完整临床调整的全基因及GO关联。目的是把名单变化与调整规格的报告衔接开，不估计因果混杂贡献。43次横断面观察对应43名独立参与者。

## 锁定输入

1. 权威名单：[fixed_input_clinical43.csv](<PROJECT_ROOT/revision_integrated_20261012/fixed_sample/tables/fixed_input_clinical43.csv>)。与10月16日复用输入、最新组成模型输入、58人GSM交叉表逐一核对身份、组别及临床原值。精确复用已保存名单，不根据新结局重选完整病例。
2. 表达：[GSE151371_normalized_expression_matrix_comBat_de-ID.csv](<PROJECT_ROOT/explore/Report/GSE151371_normalized_expression_matrix_comBat_de-ID.csv>)。保持原处理后尺度，只按原规则保留首个重复symbol，再按锁定宇宙与43人名单取子矩阵；不再次log、归一化或ComBat，不用ODC长表替代。
3. 基因：[GO_tested_gene_universe.csv](<PROJECT_ROOT/revision_integrated_20261012/supplement/GO_tested_gene_universe.csv>)，精确17,498个，顺序与原limma表核对；不按新显著性、方差或覆盖筛选，不新增符号映射。
4. GO：[GO_camera_membership.csv](<PROJECT_ROOT/revision_integrated_20261012/supplement/GO_camera_membership.csv>)，精确4,032个BP集合，15–500成员；保留原去重成员、顺序、名称及注释来源快照，不联网更新GO。
5. 所有源文件、派生矩阵、临床字段、设计矩阵与协议保存SHA-256。临床关联和表达快照仅保存在inputs_local；本轮不发布，不默认附入公开仓库或投稿附件。

## 编码、缺失与三个模型

group按TC、SCI编码，系数groupSCI始终为SCI减TC；sex为F、M；TBI及既往CNS为No、Yes、Unknown，首类别为参照，使用contr.treatment。Unknown保持来源类别，并非插补，也不证明其属于已知No。数值字段不插补，锁定后不删人。

age_c与iss_c分别减去同43人的均值。旧原评分固定样本协议曾用48人中心；本次沿用旧全基因临床模型和最新BTM的集合内中心。在有截距、无交互的固定设计中，中心常数只改变截距、不改变groupSCI或拟合值；不回写历史文件，也不重标准化评分。

| 模型 | 公式 | 独立n / SCI / TC | 设计秩 / 残差df |
| --- | --- | --- | --- |
| M0 | expression ~ group | 43 / 36 / 7 | 2 / 41 |
| M1 | expression ~ group + age_c + sex | 43 / 36 / 7 | 4 / 39 |
| M2 | expression ~ group + age_c + sex + iss_c + tbi + prior_cns | 43 / 36 / 7 | 9 / 34 |

不加入采血时间、SCI特有AIS/NLI、种族、感染、治疗、batch/RIN或组成代理，不扩展41人基因/通路网格。已有其他模型保持其原角色。这一限定不代表缺失字段不重要。

## 全基因检验

使用limma 3.68.4的lmFit(method='ls')及eBayes；无权重，proportion=0.01、stdev.coef.lim=c(0.1,4)、trend=FALSE、span=NULL、robust=FALSE、winsor.tail.p=c(0.05,0.1)、legacy=NULL，与既有调用的默认设置一致。检验groupSCI，双侧moderated t；用topTable(confint=0.95, sort.by='none', number=Inf)保存点态moderated 95%区间。

逐模型对全部17,498基因进行BH，保存基因、logFC、CI、平均表达、t、原P、BH、B、公式、样本数、设计秩及df。此处使用limma推断，不能沿用原评分OLS的HC3标签。FDR<0.05计数只作描述；如保留历史|logFC|>1计数必须标明历史阈值，不引入新阈值。

## GO CAMERA

每个模型使用相同表达、同一GO索引和各自设计，contrast明确匹配groupSCI；inter.gene.cor=0.01、use.ranks=FALSE、allow.neg.cor=FALSE、trend.var=FALSE、directional=TRUE、weights=NULL、sort=FALSE。检验采用固定limma实现，不以新结果换富集算法。

每模型完整4,032集合BH，输出GO ID、原名称、NGenes、Direction、PValue、FDR和模型；不为CAMERA虚构效应CI。摘要展示只复用原面板八个GO ID，不因新显著性增删；完整结果无删减。旧48人未调整CAMERA不是本次43人结果。

## 多重比较与跨模型解释

共六个分别列示的校正族：M0/M1/M2各17,498基因，以及M0/M1/M2各4,032 GO。保留旧BTM16/88/240、原评分18及各规格内族，不重新计算旧q。失败项保留NA与计划族大小，不因失败缩小分母。点态95%CI非同时区间。

全部三模型同时报告；不以任一模型q<0.05作为跨三模型联合错误率受控的结论，不按最小q选“最佳模型”。调整后的估计可改变比较目标，ISS未必是无争议的前置混杂因素。DEG数量、方向及区间变化为路径依赖描述，不能解释为因果百分比；一个模型显著另一个不显著也不证明效应差异。

## 旧全临床结果复用及执行核验

已有候选为[ODC405_DE_SCI_vs_TC_clinical_complete_case.csv](<PROJECT_ROOT/revision_integrated_20261012/supplement/ODC405_DE_SCI_vs_TC_clinical_complete_case.csv>)。本轮只核对名单来源、字段、设计、宇宙和旧文件完整性，不宣称数值重现已经完成。执行时如M2需要生成拟合对象，可对同一M2作一次重现核验，不另增独立模型；按基因核对logFC、CI.L、CI.R、t、P.Value和adj.P.Val，绝对及相对容差均1e-8。匹配后复用既有M2表；有差异先记录并查明，不能从新旧结果择优。

## 诊断、停止与变更

本轮预检只检查输入与设计，不拟合效应：一人一行、组别/临床一致性、基因及GO、有限数值、完整秩、df和软件版本。执行时保存数值/方差警告及设计杠杆与条件数；不以诊断P或异常值规则删除样本换结果。有限共同支持沿用已有审计（ISS交集内仅2SCI/3TC等），HC3或limma都不能修复临床不可比性。

源/快照hash改变、人数/基因/GO不符、不可估计设计、M2无法重现或未解决数值失败均暂停受影响分析。阴性结果不算失败。不得悄悄换模型、减基因、改集合、重跑LASSO或扩展范围；如需变更，用另一个版本记录日期、原因、前后hash及是否已见新效应，不修改已封存协议。

## 执行交付与当前状态

执行后需交全部三份DE、三份CAMERA、模型/检验族说明、数值核验、描述性模型变化表及session信息；暂无这些新结果。当前inputs_local快照和design/preflight只证明拟合输入已锁定，不证明新的统计要求已完成。相同43人、阴性结果、受限数据与推断边界必须保持。
