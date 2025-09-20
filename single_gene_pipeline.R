
options(stringsAsFactors = FALSE)
`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

# ----------------- 0) 依赖 -----------------
need <- c(
  # 基础
  "Seurat","data.table","dplyr","tibble","readr","stringr","Matrix","patchwork","ggplot2","viridis",
  # Single-cell 额外
  "SingleCellExperiment","Nebulosa","harmony","slingshot","RColorBrewer",
  # 通讯
  "CellChat","ktplots",
  # 空间 & 反卷积（可选）
  "png","cowplot","dbscan","spacexr","SPOTlight","mistyR"
)
inst <- rownames(installed.packages())
miss <- setdiff(need, inst)
if(length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(dplyr); library(tibble); library(readr); library(stringr)
  library(Matrix); library(patchwork); library(ggplot2); library(viridis)
  library(SingleCellExperiment); library(Nebulosa); library(harmony); library(slingshot); library(RColorBrewer)
  library(CellChat); library(ktplots)
  library(png); library(cowplot); library(dbscan); library(spacexr); library(SPOTlight); library(mistyR)
})

# ----------------- 1) 参数 -----------------
# 示例（只跑单细胞 + CellChat）：
# Rscript src/single_gene_pipeline.R \
#   gene=CASP9 outdir=./outputs sc_rds=./scRNA_anno.RDS
#
# 示例（从 10X 原始矩阵开始）：
# Rscript src/single_gene_pipeline.R \
#   gene=CASP9 outdir=./outputs \
#   set1_dir=./set1 set1_meta=./set1/GSE125449_Set1_samples.txt.gz \
#   set2_dir=./set2 set2_meta=./set2/GSE125449_Set2_samples.txt.gz
#
# 示例（附带空间/RCTD/mistyR）：
#   st_dir=./HRCC-1 lowres_img=./HRCC-1/spatial/tissue_lowres_image.png
#   （或已预处理好的） st_rds=./stRNA.RDS
args <- commandArgs(trailingOnly = TRUE)
kv   <- strsplit(args, "=")
opts <- setNames(vapply(kv, function(x) if(length(x)==2) x[2] else NA_character_, ""),
                 vapply(kv, function(x) x[1], ""))

GENE       <- opts[["gene"]]         %||% "CASP9"
OUTDIR     <- opts[["outdir"]]       %||% "./outputs"

# 单细胞输入二选一：已有 RDS，或者 set1/set2 原始矩阵
SC_RDS     <- opts[["sc_rds"]]       %||% ""
SET1_DIR   <- opts[["set1_dir"]]     %||% ""
SET1_META  <- opts[["set1_meta"]]    %||% ""
SET2_DIR   <- opts[["set2_dir"]]     %||% ""
SET2_META  <- opts[["set2_meta"]]    %||% ""

# 拟时序可选
DO_TRAJ    <- as.logical(opts[["do_traj"]] %||% "TRUE")   # 是否运行 Harmony+slingshot

# 空间 & RCTD（可选）
ST_DIR     <- opts[["st_dir"]]       %||% ""   # Visium 文件夹（含 filtered_feature_bc_matrix.h5）
LOWRES_IMG <- opts[["lowres_img"]]   %||% ""   # 低清图 PNG
ST_RDS     <- opts[["st_rds"]]       %||% ""   # 若已有 stRNA.RDS 可直接用
RUN_RCTD   <- as.logical(opts[["run_rctd"]] %||% "FALSE") # 是否运行 RCTD（需 ST & scRNA_<GENE>.RDS）
RCTD_OUT   <- opts[["rctd_out"]]     %||% file.path(OUTDIR,"weight_RCTD.csv")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(OUTDIR, sprintf("single_gene_pipeline_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
zz <- file(logf, open = "wt"); sink(zz, type = "output"); sink(zz, type = "message")
t0 <- Sys.time()
on.exit({
  cat("\n[INFO] Session info:\n"); print(sessionInfo())
  cat(sprintf("[INFO] Elapsed: %.2f secs\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))
  sink(type="message"); sink(type="output"); close(zz)
}, add = TRUE)
cat("[INFO]", as.character(Sys.time()), "Start single-gene pipeline for:", GENE, "\n")

# ----------------- 2) 读取单细胞 -----------------
if(nchar(SC_RDS) && file.exists(SC_RDS)){
  sc <- readRDS(SC_RDS)
  stopifnot(inherits(sc, "Seurat"))
  cat("[INFO] Loaded existing scRNA:", SC_RDS, "\n")
} else {
  stopifnot(nchar(SET1_DIR) && nchar(SET1_META) && nchar(SET2_DIR) && nchar(SET2_META))
  stopifnot(file.exists(SET1_META), file.exists(SET2_META))
  set1_mtx <- Read10X(SET1_DIR)
  set2_mtx <- Read10X(SET2_DIR)
  samp1 <- data.table::fread(SET1_META, data.table = FALSE); rownames(samp1) <- samp1[["Cell Barcode"]]
  samp2 <- data.table::fread(SET2_META, data.table = FALSE); rownames(samp2) <- samp2[["Cell Barcode"]]
  stopifnot(identical(rownames(samp1), colnames(set1_mtx)))
  stopifnot(identical(rownames(samp2), colnames(set2_mtx)))
  set1 <- CreateSeuratObject(set1_mtx, project = "set1", meta.data = samp1)
  set2 <- CreateSeuratObject(set2_mtx, project = "set2", meta.data = samp2)
  sc   <- merge(set1, set2)
  cat("[INFO] Built Seurat object from 10X matrices.\n")
}

# 清理潜在无用列
if("Cell.Barcode" %in% colnames(sc@meta.data)) sc$Cell.Barcode <- NULL

# 标准流程
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)
sc <- ScaleData(sc, features = VariableFeatures(sc))
sc <- RunPCA(sc, features = VariableFeatures(sc))
p1 <- DimPlot(sc, reduction = "pca", group.by = "orig.ident")
p2 <- ElbowPlot(sc, ndims = 20, reduction = "pca")
ggsave(file.path(OUTDIR,"sg_3.1_pca_elbow.pdf"), p1 + p2, width=14, height=8)

pc.num <- 1:15
sc <- FindNeighbors(sc, dims = pc.num)
sc <- FindClusters(sc)
sc <- RunTSNE(sc, dims = pc.num)
sc <- RunUMAP(sc, dims = pc.num)

# 细胞类型列统一（存在时）
if("Celltypes" %in% colnames(sc@meta.data)) sc$celltype <- sc$Celltypes
Idents(sc) <- if("celltype" %in% colnames(sc@meta.data)) sc$celltype else sc$seurat_clusters

saveRDS(sc, file = file.path(OUTDIR,"scRNA_anno.RDS"))
ggsave(file.path(OUTDIR,"sg_3.2_tsne_clusters.pdf"), DimPlot(sc, reduction="tsne", group.by="seurat_clusters", label=TRUE, label.size=5), width=10, height=8)
ggsave(file.path(OUTDIR,"sg_3.3_tsne_celltype.pdf"), DimPlot(sc, reduction="tsne", group.by="celltype", label=TRUE, label.size=5), width=10, height=8)
ggsave(file.path(OUTDIR,"sg_3.4_umap_celltype.pdf"), DimPlot(sc, reduction="umap", group.by="celltype", label=TRUE, label.size=4), width=10, height=8)

# ----------------- 3) 单基因表达可视化 -----------------
for(g in unique(c(GENE,"CASP8","CASP10"))){
  if(g %in% rownames(sc)){
    ggsave(file.path(OUTDIR, paste0("sg_3.5_dot_",g,".pdf")), DotPlot(sc, features=g), width=6, height=6)
    ggsave(file.path(OUTDIR, paste0("sg_3.6_umap_",g,".pdf")),
           FeaturePlot(sc, features=g, label=FALSE, order=TRUE, pt.size=0.5, cols=viridis::magma(10)),
           width=8, height=6)
    ggsave(file.path(OUTDIR, paste0("sg_3.7_density_",g,".pdf")),
           Nebulosa::plot_density(sc, g, reduction="umap"), width=8, height=6)
  } else {
    message("[WARN] gene not found in counts: ", g)
  }
}

# ----------------- 4) 单基因分组 & CellChat -----------------
# 癌细胞子集
if(!("celltype" %in% colnames(sc@meta.data)))
  stop("[ERROR] 缺少 celltype 列，请在 meta 中提供或修改脚本的标注列。")
sc_mali   <- subset(sc, celltype == "Cancer cells")
# 基于 raw counts 是否 >0 来判定阳性
if(!(GENE %in% rownames(sc_mali))) stop("[ERROR] 目标基因不存在：", GENE)
expr_vec <- as.numeric(sc_mali[["RNA"]]@counts[GENE, ])
sc_mali$gene_group <- ifelse(expr_vec > 0, paste0(GENE,"+Mali"), paste0(GENE,"-Mali"))

sc_other <- subset(sc, celltype != "Cancer cells")
sc_other$gene_group <- sc_other$celltype
sc_chat <- merge(sc_mali, sc_other)

saveRDS(sc_chat, file = file.path(OUTDIR, paste0("scRNA_",GENE,".RDS")))

# 控制 CellChat 细胞数（避免内存爆）
set.seed(123456)
if(ncol(sc_chat) > 4000) sc_chat <- sc_chat[, sample(colnames(sc_chat), 4000)]

meta <- sc_chat@meta.data
data_input <- as.matrix(sc_chat[["RNA"]]@data)
stopifnot(identical(colnames(data_input), rownames(meta)))

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "gene_group")
CellChatDB.use <- subsetDB(CellChatDB.human, search = c("Secreted Signaling","Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 0)
cellchat <- computeCommunProbPathway(cellchat, thresh = 1)
df.net   <- subsetCommunication(cellchat, thresh = 1)

# 转 cpdb 宽表
rt <- df.net; rt$cell_cell <- paste0(rt$source,"|",rt$target)
rt_means <- rt[, c("cell_cell","prob","interaction_name")]; colnames(rt_means)[3] <- "interacting_pair"
means_cellchat <- tidyr::pivot_wider(rt_means, names_from = "cell_cell", values_from = "prob", values_fill = 0)
means_cellchat$id_cp_interaction <- paste0("CPI-", seq_len(nrow(means_cellchat)))
rt_p <- rt[, c("cell_cell","pval","interaction_name")]; colnames(rt_p)[3] <- "interacting_pair"
pvalues_cellchat <- tidyr::pivot_wider(rt_p, names_from = "cell_cell", values_from = "pval", values_fill = 1)
pvalues_cellchat$id_cp_interaction <- paste0("CPI-", seq_len(nrow(pvalues_cellchat)))

# 热图
pheat <- plot_cpdb_heatmap(subset(pvalues_cellchat, select = -id_cp_interaction))
ggsave(file.path(OUTDIR, paste0("sg_3.8_cpdb_heatmap_",GENE,".pdf")), pheat, width=10, height=8)

# 典型配对图示例（可按需换 cell_type1/2）
sce <- as.SingleCellExperiment(sc_chat)
plot_one_cpdb <- function(ct1, ct2, fn){
  p <- plot_cpdb(
    scdata = sce, cell_type1 = ct1, cell_type2 = ct2,
    col_option = viridis::magma(10), max_highlight_size = 3,
    keep_significant_only = TRUE, max_size = 3, standard_scale = TRUE,
    celltype_key = "gene_group", means = means_cellchat, pvals = pvalues_cellchat
  ) + small_axis(7) + small_grid() + small_guide()
  ggsave(file.path(OUTDIR, fn), p, width=10, height=8)
}
suppressWarnings({
  if("Dendritic cells" %in% sc_chat$gene_group) plot_one_cpdb("Dendritic cells","Mali", paste0("sg_3.81_cpdb_",GENE,".pdf"))
  if("Macrophages"    %in% sc_chat$gene_group) plot_one_cpdb("Macrophages","Mali",    paste0("sg_3.82_cpdb_",GENE,".pdf"))
  if("Cycling cells"  %in% sc_chat$gene_group) plot_one_cpdb("Cycling cells","Mali",  paste0("sg_3.83_cpdb_",GENE,".pdf"))
})

# ----------------- 5) 拟时序（可选） -----------------
if(DO_TRAJ){
  sc_sub <- subset(sc_chat, celltype == "Cancer cells")
  sc_sub <- FindVariableFeatures(sc_sub, selection.method = "vst", nfeatures = 3000)
  sc_sub <- ScaleData(sc_sub, features = VariableFeatures(sc_sub))
  sc_sub <- RunPCA(sc_sub, features = VariableFeatures(sc_sub))
  ggsave(file.path(OUTDIR, paste0("sg_3.9_pca_elbow_",GENE,".pdf")),
         DimPlot(sc_sub, reduction = "pca", group.by = "orig.ident") + ElbowPlot(sc_sub, ndims=20, reduction="pca"),
         width=14, height=8)
  set.seed(1000)
  sc_sub <- RunHarmony(sc_sub, group.by.vars = "orig.ident")
  ggsave(file.path(OUTDIR, paste0("sg_3.91_harmony_",GENE,".pdf")),
         DimPlot(sc_sub, reduction="harmony", group.by="orig.ident"), width=7, height=6)
  ggsave(file.path(OUTDIR, paste0("sg_3.92_harmony_elbow_",GENE,".pdf")), ElbowPlot(sc_sub, reduction="harmony"), width=7, height=6)
  sc_sub <- FindNeighbors(sc_sub, reduction = "harmony", dims = 1:8)
  sc_sub <- FindClusters(sc_sub)
  sc_sub <- RunUMAP(sc_sub, reduction = "harmony", dims = 1:8)
  ggsave(file.path(OUTDIR, paste0("sg_3.93_umap_clusters_",GENE,".pdf")), DimPlot(sc_sub, label=TRUE), width=7, height=6)
  if(GENE %in% rownames(sc_sub)){
    ggsave(file.path(OUTDIR, paste0("sg_3.94_umap_",GENE,".pdf")),
           FeaturePlot(sc_sub, features=GENE, label=FALSE, order=TRUE, pt.size=0.5, cols=viridis::magma(10)), width=8, height=6)
    ggsave(file.path(OUTDIR, paste0("sg_3.95_density_",GENE,".pdf")), Nebulosa::plot_density(sc_sub, GENE), width=8, height=6)
  }
}

# ----------------- 6) 空间（Visium）& RCTD（可选） -----------------
as_matrix_sparse <- function(mat){
  # 将 dgCMatrix 转成普通 matrix，用于 RCTD & 自定义绘图
  as.matrix(mat)
}

# 构建或加载 stRNA
if(nchar(ST_RDS) && file.exists(ST_RDS)){
  st <- readRDS(ST_RDS)
  cat("[INFO] Loaded existing ST RDS:", ST_RDS, "\n")
} else if(nchar(ST_DIR)){
  stopifnot(file.exists(file.path(ST_DIR,"filtered_feature_bc_matrix.h5")))
  st <- Load10X_Spatial(data.dir = ST_DIR, filename = "filtered_feature_bc_matrix.h5", filter.matrix = FALSE)
  # on_tissue 过滤 + 常见 QC
  coldata <- GetTissueCoordinates(st, cols = c("row","col","tissue"), scale = NULL)
  st$tissue <- ifelse(coldata[rownames(st@meta.data),"tissue"] == 1, "on_tissue","not_on_tissue")
  st <- subset(st, subset = tissue == "on_tissue")
  st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "^MT-")
  # 去除线粒体/核糖体基因（只删基因，不删列）
  mt  <- grep("^MT-", rownames(st), value=TRUE)
  rps <- grep("^RPS", rownames(st), value=TRUE)
  mrp <- grep("^MRP", rownames(st), value=TRUE)
  rpl <- grep("^RPL", rownames(st), value=TRUE)
  st  <- st[!rownames(st) %in% c(mt,rps,mrp,rpl), ]
  st  <- st[rowSums(GetAssayData(st, assay="Spatial") > 0) > 10, ]
  st$nFeature_Spatial_filt <- colSums(GetAssayData(st, assay="Spatial") > 0)
  st$nCount_Spatial_filt   <- colSums(GetAssayData(st, assay="Spatial"))
  st <- subset(st, subset = nFeature_Spatial_filt > 300 & nCount_Spatial_filt > 500)
  st  <- SCTransform(st, assay = "Spatial", verbose = FALSE)
  st  <- NormalizeData(st, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  DefaultAssay(st) <- "SCT"
  st <- ScaleData(st, features = rownames(st), verbose = FALSE) %>% RunPCA() %>% RunUMAP(reduction="pca", dims=1:30, verbose=FALSE)
  saveRDS(st, file = file.path(OUTDIR,"stRNA.RDS"))
  ggsave(file.path(OUTDIR,"sg_spatial_umap.pdf"), SpatialDimPlot(st, label=TRUE, pt.size.factor=1), width=7, height=6)
} else {
  st <- NULL
}

# 单基因空间可视化（若提供 ST）
if(!is.null(st) && (GENE %in% rownames(st))){
  ggsave(file.path(OUTDIR, paste0("sg_spatial_",GENE,".pdf")), SpatialFeaturePlot(st, features=GENE, pt.size.factor=1.2), width=6, height=6)
}

# RCTD（需要 st & scRNA_<GENE>.RDS）
if(RUN_RCTD && !is.null(st)){
  sc_gene <- readRDS(file.path(OUTDIR, paste0("scRNA_",GENE,".RDS")))
  # 参考集：使用 gene_group
  sc_counts <- as_matrix_sparse(sc_gene[["RNA"]]@counts)
  sc_nUMI   <- colSums(sc_counts)
  cellType  <- data.frame(barcode=colnames(sc_gene), cell_type=sc_gene$gene_group, check.names = FALSE)
  ref <- Reference(sc_counts, as.factor(setNames(cellType$cell_type, cellType$barcode)), sc_nUMI)
  
  coords   <- GetTissueCoordinates(st, cols = c("row","col"), scale = NULL)
  colnames(coords) <- c("xcoord","ycoord")
  sp_counts <- as_matrix_sparse(st[["Spatial"]]@counts)
  sp_nUMI   <- colSums(sp_counts)
  puck <- SpatialRNA(coords, sp_counts, sp_nUMI)
  
  myRCTD <- create.RCTD(puck, ref)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
  save(myRCTD, file = file.path(OUTDIR, paste0("RCTD_",GENE,".Rdata")))
  w <- as.data.frame(myRCTD@results$weights)
  write.csv(w, file = RCTD_OUT, quote = FALSE)
  cat("[OK] RCTD finished:", RCTD_OUT, "\n")
}

# ----------------- 7) 同/异型空间网络与邻域富集（若有 RCTD 结果 & ST） -----------------
if(!is.null(st) && file.exists(RCTD_OUT) && nchar(LOWRES_IMG) && file.exists(LOWRES_IMG)){
  decon <- read.csv(RCTD_OUT, header = TRUE, row.names = 1, check.names = FALSE)
  # 对齐条形码
  ss <- intersect(colnames(st), rownames(decon))
  st <- st[, ss, drop=FALSE]; decon <- decon[ss, , drop=FALSE]
  
  slice <- names(st@images)[1]
  spatial_coord <- data.frame(st@images[[slice]]@coordinates) %>%
    tibble::rownames_to_column("barcodeID") %>%
    dplyr::mutate(imagerow_scaled = imagerow * st@images[[slice]]@scale.factors$lowres,
                  imagecol_scaled = imagecol * st@images[[slice]]@scale.factors$lowres)
  
  # knn 网络
  xys <- data.frame(x = spatial_coord$imagecol_scaled, y = spatial_coord$imagerow_scaled, row.names = spatial_coord$barcodeID)
  nNeighbours <- 6; maxdist <- 200; minK <- 4
  knn <- dbscan::kNN(as.matrix(xys), k = nNeighbours)
  net <- data.frame(from = rep(1:nrow(knn$id), nNeighbours),
                    to = as.vector(knn$id),
                    weight = 1/(1 + as.vector(knn$dist)),
                    distance = as.vector(knn$dist))
  spotnames <- rownames(xys); names(spotnames) <- seq_len(nrow(xys))
  net$from <- spotnames[net$from]; net$to <- spotnames[net$to]
  net <- net %>% group_by(from) %>% mutate(rnk = rank(distance)) %>% ungroup()
  spatnet <- subset(net, distance <= maxdist | rnk <= minK)
  spatnet <- cbind(spatnet, setNames(xys[spatnet$from, ], c("start_x","start_y")))
  spatnet <- cbind(spatnet, setNames(xys[spatnet$to, ],   c("end_x","end_y")))
  
  # ——同型网络：GENE+Mali——
  same_col <- paste0(GENE,"+Mali"); thr_same <- 0.1
  if(same_col %in% colnames(decon)){
    d_same <- decon[, same_col, drop=FALSE]; colnames(d_same) <- "prop"
    d_same <- subset(d_same, prop > thr_same)
    degree <- matrix(NA_real_, nrow=nrow(d_same), ncol=1, dimnames=list(rownames(d_same),"degree"))
    for(i in rownames(d_same)){
      su <- spatnet[spatnet$from == i, ]
      degree[i,1] <- sum(decon[su$to, same_col] > thr_same, na.rm = TRUE)
    }
    degree <- na.omit(as.data.frame(degree))
    img <- png::readPNG(LOWRES_IMG); img_grob <- grid::rasterGrob(img, interpolate=FALSE, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
    seg <- spatnet[spatnet$from %in% rownames(degree) & spatnet$to %in% rownames(degree), ]
    coord_deg <- spatial_coord[spatial_coord$barcodeID %in% rownames(degree), ]
    coord_deg$degree <- degree[coord_deg$barcodeID, "degree"]
    p_same <- ggplot() + annotation_custom(grob=img_grob, xmin=0,xmax=ncol(img),ymin=0,ymax=-nrow(img)) +
      geom_segment(data=seg, aes(x=start_x,y=start_y,xend=end_x,yend=end_y), color="white", linewidth=0.2) +
      geom_point(data=coord_deg, aes(x=imagecol_scaled,y=imagerow_scaled,color=degree), size=1.0) +
      ylim(nrow(img),0) + xlim(0,ncol(img)) + theme_void() + coord_fixed() + scale_y_reverse() +
      scale_color_gradient(low="#E5F4F8", high="#bc3c29") + ggtitle(paste0("Homotypic degree: ", same_col))
    ggsave(file.path(OUTDIR, paste0("sg_spatial_homo_",GENE,".pdf")), p_same, width=6, height=6)
  }
  
  # ——异型网络：GENE+Mali vs 单一免疫群（示例 Monocytes，可自行改）——
  cell2 <- if("Monocytes" %in% colnames(decon)) "Monocytes" else colnames(decon)[1]
  thr_het <- 0.2
  d1 <- decon[, paste0(GENE,"+Mali"), drop=FALSE]; colnames(d1) <- "prop1"; d1 <- subset(d1, prop1 > thr_het)
  d2 <- decon[, cell2, drop=FALSE]; colnames(d2) <- "prop2"; d2 <- subset(d2, prop2 > thr_het)
  seg_het <- spatnet[spatnet$from %in% rownames(d1) & spatnet$to %in% rownames(d2), ]
  if(nrow(seg_het) > 0){
    img <- png::readPNG(LOWRES_IMG); img_grob <- grid::rasterGrob(img, interpolate=FALSE, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
    c1 <- spatial_coord[spatial_coord$barcodeID %in% unique(seg_het$from), ]
    c2 <- spatial_coord[spatial_coord$barcodeID %in% unique(seg_het$to), ]
    p_het <- ggplot() + annotation_custom(grob=img_grob, xmin=0,xmax=ncol(img),ymin=0,ymax=-nrow(img)) +
      geom_segment(data=seg_het, aes(x=start_x,y=start_y,xend=end_x,yend=end_y), color="white", linewidth=0.2) +
      geom_point(data=c1, aes(x=imagecol_scaled,y=imagerow_scaled), size=1.0, color="yellow") +
      geom_point(data=c2, aes(x=imagecol_scaled,y=imagerow_scaled), size=1.0, color="purple") +
      ylim(nrow(img),0) + xlim(0,ncol(img)) + theme_void() + coord_fixed() + scale_y_reverse() +
      ggtitle(paste0("Heterotypic: ", GENE,"+Mali ↔ ", cell2))
    ggsave(file.path(OUTDIR, paste0("sg_spatial_hetero_",GENE,".pdf")), p_het, width=6, height=6)
  }
  
  # 邻域富集（以 cell2 为邻居）
  deg2 <- matrix(NA_real_, nrow=nrow(d1), ncol=1, dimnames=list(rownames(d1),"neighbor_cell2"))
  for(i in rownames(d1)){
    su <- spatnet[spatnet$from == i, ]
    deg2[i,1] <- sum(decon[c(su$to,su$from), cell2], na.rm = TRUE)
  }
  deg2 <- as.data.frame(deg2); deg2$enrich <- deg2$neighbor_cell2 / max(deg2$neighbor_cell2, na.rm = TRUE)
  coord_en <- spatial_coord[spatial_coord$barcodeID %in% rownames(deg2), ]
  coord_en$enrich <- deg2[coord_en$barcodeID, "enrich"]
  img <- png::readPNG(LOWRES_IMG); img_grob <- grid::rasterGrob(img, interpolate=FALSE, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
  p_en <- ggplot() + annotation_custom(grob=img_grob, xmin=0,xmax=ncol(img),ymin=0,ymax=-nrow(img)) +
    geom_point(data=coord_en, aes(x=imagecol_scaled,y=imagerow_scaled, color=enrich, alpha=enrich), size=1.0) +
    ylim(nrow(img),0) + xlim(0,ncol(img)) + theme_void() + coord_fixed() + scale_y_reverse() +
    scale_color_gradient2(low="blue", mid="white", high="yellow", midpoint=0.5) + guides(alpha="none") +
    ggtitle(paste0("Neighbor enrichment: ", cell2))
  ggsave(file.path(OUTDIR, paste0("sg_spatial_enrich_",GENE,".pdf")), p_en, width=6, height=6)
}

cat("[DONE] All finished. Outputs in: ", OUTDIR, "\n")
