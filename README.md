本仓库包含从 GEO/TCGA/Xena 数据准备 → 单基因/基因集单细胞分析 → 空间反卷积与通讯 → SMR 共定位 → Bulk 风险评分模型（Cox/LASSO/DeepSurv） 的完整、可复现实验脚本。

适配系统：Windows / Linux / macOS
主要语言：R、Python（可选）
典型场景：投稿复现、图件重绘、交叉验证与外部验证

目录结构
.
├─ bulk_pipeline.R               # Bulk 表达 + 生存建模/验证（单/多因素Cox、LASSO、ROC、临床分层、免疫/通路/检查点）
├─ geneset_pipeline.R            # 基因集驱动的单细胞分析（打分、分群、CellChat、拟时序、空间反卷积/MISTy）
├─ single_gene_pipeline.R        # 单基因（CASP9/8/10 等）驱动的单细胞分析与后续验证
├─ prep_data.R                   # GEO/TCGA/Xena 数据下载与预处理（表达矩阵、ID转换、TPM、分组等）
├─ smr_prep_and_plots.R          # SMR 输入整理、结果读取、候选基因筛选、区域关联图与效应图
├─ commot.py                     # （可选）空间通讯/共定位的 Python 实现或补充分析
├─ deepsurv.py                   # （可选）基于 PyTorch 的 DeepSurv 风险评分模型
└─ README.md

环境与依赖
1) R 包

建议使用 renv 或 conda + r-base 管理环境。核心依赖（按脚本聚合）：

通用/数据：data.table, tibble, dplyr, tidyr, stringr, readr, ggplot2, ggpubr, viridis, patchwork

单细胞/空间：Seurat, harmony, Nebulosa, slingshot, SingleCellExperiment, spacexr, SPOTlight, CellChat, ktplots, dbscan, Matrix, mistyR

GEO/注释/表达：GEOquery, limma, tinyarray, sva

生存/建模：survival, survminer, glmnet, timeROC, forestplot, caret

免疫与通路：IOBR, estimate

SMR/图件：vcfR, magick, ggplot2

安装示例：

install.packages(c(
  "data.table","tibble","dplyr","tidyr","stringr","readr",
  "Seurat","harmony","Nebulosa","slingshot","SingleCellExperiment",
  "spacexr","SPOTlight","CellChat","dbscan","Matrix","ggplot2","ggpubr",
  "viridis","patchwork","GEOquery","limma","tinyarray","sva",
  "survival","survminer","glmnet","timeROC","forestplot","caret",
  "IOBR","estimate","magick"
))
# ktplots & mistyR 可能需要 remotes：
remotes::install_github("zktuong/ktplots")
remotes::install_github("saezlab/mistyR")

2) Python 包（仅当使用 commot.py / deepsurv.py）

推荐 conda 新建环境：

conda create -n rcc-casp9 python=3.10 -y
conda activate rcc-casp9
pip install numpy pandas scikit-learn matplotlib seaborn
# DeepSurv（PyTorch / 或 pycox 二选一）
pip install torch torchaudio torchvision pycox lifelines
# 空间通讯/共定位（如 commot、scanpy、anndata）
pip install scanpy anndata commot


若脚本需要额外包，会在运行时报错提示，请按需补装。
