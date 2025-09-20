

options(stringsAsFactors = FALSE)
`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

# ----------------- 0) 依赖 -----------------
need <- c(
  "Seurat","dplyr","readr","tibble","stringr","Matrix","SingleCellExperiment",
  "CellChat","viridis","ktplots","harmony","Nebulosa","slingshot","RColorBrewer",
  "Rcpp","spacexr","ggplot2","ggpubr","gridExtra","reshape2","png","patchwork",
  "dbscan","cowplot","mistyR","SPOTlight"
)
inst <- rownames(installed.packages())
miss <- setdiff(need, inst)
if(length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(readr); library(tibble); library(stringr)
  library(Matrix); library(SingleCellExperiment)
  library(CellChat); library(viridis); library(ktplots)
  library(harmony); library(Nebulosa)
  library(slingshot); library(RColorBrewer)
  library(Rcpp); library(spacexr)
  library(ggplot2); library(ggpubr); library(gridExtra); library(reshape2)
  library(png); library(patchwork); library(dbscan); library(cowplot)
  library(mistyR); library(SPOTlight)
})

# ----------------- 1) 命令行参数 -----------------
# 示例：
# Rscript src/geneset_pipeline.R \
#   sc_rds=./scRNA_anno.RDS \
#   geneset_file=./senescence.txt \
#   vector_script=./Vector.R \
#   st_rds=./stRNA.RDS \
#   lowres_img=./HRCC-1/spatial/tissue_lowres_image.png \
#   rctd_rdata=./RCTD_geneset.Rdata \
#   outdir=./outputs \
#   sample_n=3000 \
#   high_label=HighSet \
#   low_label=LowSet \
#   cell1_for_hetero=LowSet \
#   cell2_for_hetero='Dendritic cells'
args <- commandArgs(trailingOnly = TRUE)
kv   <- strsplit(args, "=")
opts <- setNames(vapply(kv, function(x) if(length(x)==2) x[2] else NA_character_, ""),
                 vapply(kv, function(x) x[1], ""))

SC_RDS       <- opts[["sc_rds"]]        %||% stop("缺少参数: sc_rds")
GENESET_FILE <- opts[["geneset_file"]]  %||% stop("缺少参数: geneset_file")
VECTOR_R     <- opts[["vector_script"]] %||% ""
ST_RDS       <- opts[["st_rds"]]        %||% ""
LOWRES_IMG   <- opts[["lowres_img"]]    %||% ""
RCTD_RDATA   <- opts[["rctd_rdata"]]    %||% ""
OUTDIR       <- opts[["outdir"]]        %||% "./outputs"
SAMPLE_N     <- as.integer(opts[["sample_n"]] %||% "3000")
HIGH_LAB     <- opts[["high_label"]]    %||% "SetHigh"
LOW_LAB      <- opts[["low_label"]]     %||% "SetLow"
CELL1_HET    <- opts[["cell1_for_hetero"]] %||% LOW_LAB
CELL2_HET    <- opts[["cell2_for_hetero"]] %||% "Dendritic cells"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(OUTDIR, sprintf("geneset_pipeline_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
zz <- file(logf, open = "wt"); sink(zz, type = "output"); sink(zz, type = "message")
t0 <- Sys.time()
on.exit({
  cat("\n[INFO] Session info:\n"); print(sessionInfo())
  cat(sprintf("[INFO] Elapsed: %.2f secs\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))
  sink(type="message"); sink(type="output"); close(zz)
}, add = TRUE)
cat("[INFO]", as.character(Sys.time()), "Start geneset pipeline\n")

# ----------------- 2) 读入单细胞 & 基因集 -----------------
stopifnot(file.exists(SC_RDS), file.exists(GENESET_FILE))
sc <- readRDS(SC_RDS)
stopifnot(inherits(sc, "Seurat"))

cat("[INFO] Celltype table:\n"); print(table(sc$celltype))
sc_mali <- subset(sc, celltype == "Cancer cells")

gs_vec <- read.table(GENESET_FILE, header = FALSE)[[1]]
stopifnot(length(gs_vec) > 0)
geneset <- list(geneset = unique(gs_vec))

# ----------------- 3) 基因集打分与分组 -----------------
sc_mali <- AddModuleScore(sc_mali, features = geneset, name = "geneset")
# AddModuleScore 会生成 "geneset1" 列
stopifnot("geneset1" %in% colnames(sc_mali@meta.data))
thr <- median(sc_mali$geneset1, na.rm = TRUE)
sc_mali$geneset_group <- ifelse(sc_mali$geneset1 > thr, paste0(HIGH_LAB,"Mali"), paste0(LOW_LAB,"Mali"))

# 其他细胞保持原类型名以便做 CellChat 的多群体
sc_other <- subset(sc, celltype != "Cancer cells")
sc_other$geneset_group <- sc_other$celltype

sc_chat <- merge(sc_mali, sc_other)
saveRDS(sc_chat, file = file.path(OUTDIR, "scRNA_geneset.RDS"))
cat("[OK] Saved scRNA_geneset.RDS\n")

# ----------------- 4) （可选）抽样控制细胞数 -----------------
set.seed(123456)
if (ncol(sc_chat) > SAMPLE_N) {
  idx <- sample(colnames(sc_chat), SAMPLE_N)
  sc_chat <- sc_chat[, idx]
  cat("[INFO] Downsample cells to:", SAMPLE_N, "\n")
}

# ----------------- 5) CellChat 管道 -----------------
# 使用规范输入：表达矩阵应与元数据条形码一致
meta <- sc_chat@meta.data
data_input <- as.matrix(sc_chat[["RNA"]]@data)
stopifnot(identical(colnames(data_input), rownames(meta)))

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "geneset_group")
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 0)
cellchat <- computeCommunProbPathway(cellchat, thresh = 1)
df.net <- subsetCommunication(cellchat, thresh = 1)

# 整理为 cpdb 风格矩阵以供 ktplots
rt <- df.net
rt$cell_cell <- paste0(rt$source, "|", rt$target)
rt_means <- rt[, c("cell_cell","prob","interaction_name")]
colnames(rt_means)[3] <- "interacting_pair"
means_cellchat <- tidyr::pivot_wider(rt_means, names_from = "cell_cell", values_from = "prob", values_fill = 0)
means_cellchat$id_cp_interaction <- paste0("CPI-", seq_len(nrow(means_cellchat)))
means_cellchat <- as.data.frame(means_cellchat)

rt_pvalues <- rt[, c("cell_cell","pval","interaction_name")]
colnames(rt_pvalues)[3] <- "interacting_pair"
pvalues_cellchat <- tidyr::pivot_wider(rt_pvalues, names_from = "cell_cell", values_from = "pval", values_fill = 1)
pvalues_cellchat$id_cp_interaction <- paste0("CPI-", seq_len(nrow(pvalues_cellchat)))
pvalues_cellchat <- as.data.frame(pvalues_cellchat)

# ktplots 绘图示例（注意 celltype_key 要对应元数据里的列）
sce <- as.SingleCellExperiment(sc_chat)
# 下面两个示例里，cell_type1/2 请使用你在 geneset_group 中的 **实际水平名称**
# 示例：巨噬 vs 癌细胞
try({
  p_cpdb1 <- plot_cpdb(
    scdata = sce,
    cell_type1 = "Macrophages",
    cell_type2 = paste0(HIGH_LAB,"Mali"),
    col_option = viridis::magma(10),
    max_highlight_size = 3,
    keep_significant_only = TRUE, max_size = 3,
    standard_scale = TRUE,
    celltype_key = "geneset_group",
    means = means_cellchat, pvals = pvalues_cellchat
  ) + small_axis(7) + small_grid() + small_guide()
  ggsave(file.path(OUTDIR,"cpdb_macrophage_vs_mali.pdf"), p_cpdb1, width=8, height=6)
}, silent = TRUE)

# ----------------- 6) 拟时序（Harmony + slingshot） -----------------
sc_sub <- subset(sc_chat, celltype == "Cancer cells")
sc_sub <- FindVariableFeatures(sc_sub, selection.method = "vst", nfeatures = 3000)
sc_sub <- ScaleData(sc_sub, features = VariableFeatures(sc_sub))
sc_sub <- RunPCA(sc_sub, features = VariableFeatures(sc_sub))
pca_dim <- ElbowPlot(sc_sub, ndims = 20)
ggsave(file.path(OUTDIR,"pca_elbow.pdf"), pca_dim, width=6, height=4)

set.seed(1000)
sc_sub <- RunHarmony(sc_sub, group.by.vars = "orig.ident")
sc_sub <- FindNeighbors(sc_sub, reduction = "harmony", dims = 1:10)
sc_sub <- FindClusters(sc_sub, resolution = 0.8)
sc_sub <- RunUMAP(sc_sub, reduction = "harmony", dims = 1:10)
sc_sub <- RunTSNE(sc_sub, reduction = "harmony", dims = 1:10)
ggsave(file.path(OUTDIR,"umap_clusters.pdf"), DimPlot(sc_sub, label=TRUE), width=6, height=5)

# 打分可视化
p_feat <- FeaturePlot(sc_sub, features = "geneset1", order = TRUE, pt.size = 0.4) +
  scale_color_viridis_c()
ggsave(file.path(OUTDIR,"umap_geneset1.pdf"), p_feat, width=6, height=5)
p_den <- Nebulosa::plot_density(sc_sub, "geneset1", reduction = "umap")
ggsave(file.path(OUTDIR,"umap_density_geneset1.pdf"), p_den, width=6, height=5)

# slingshot 使用 UMAP 坐标 + 聚类
emb <- sc_sub@reductions$umap@cell.embeddings
clu <- Idents(sc_sub)
sce_sl <- SingleCellExperiment(list(logcounts = sc_sub[["RNA"]]@data),
                               colData = data.frame(cluster = clu))
reducedDims(sce_sl)$UMAP <- emb
sce_sl <- slingshot(sce_sl, clusterLabels = 'cluster', reducedDim = 'UMAP')
saveRDS(sce_sl, file = file.path(OUTDIR,"slingshot_sce.rds"))

# ----------------- 7) 量子极化向量（可选） -----------------
if(nchar(VECTOR_R) && file.exists(VECTOR_R)){
  source(VECTOR_R)
  PCA <- sc_sub@reductions$harmony@cell.embeddings
  PCA <- vector.rankPCA(PCA)
  VEC <- sc_sub@reductions$umap@cell.embeddings
  rownames(VEC) <- colnames(sc_sub)
  
  OUT <- vector.buildGrid(VEC, N=30, SHOW=FALSE)
  OUT <- vector.buildNet(OUT, CUT=1, SHOW=FALSE)
  OUT <- vector.getValue(OUT, PCA, SHOW=FALSE)
  OUT <- vector.gridValue(OUT, SHOW=FALSE)
  OUT <- vector.autoCenter(OUT, UP=0.9, SHOW=FALSE)
  OUT <- vector.drawArrow(OUT, P=0.9, SHOW=FALSE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
  saveRDS(OUT, file = file.path(OUTDIR,"vector_OUT.rds"))
  cat("[OK] Saved vector_OUT.rds\n")
} else {
  cat("[INFO] 跳过 Vector.R（未提供或不存在）。\n")
}

# ----------------- 8) RCTD 反卷积（需要 ST_RDS） -----------------
if(nchar(ST_RDS) && file.exists(ST_RDS)){
  # --- as_matrix 工具（Rcpp） ---
  Rcpp::sourceCpp(code='
    #include <Rcpp.h>
    using namespace Rcpp;
    // [[Rcpp::export]]
    IntegerMatrix asMatrix(NumericVector rp, NumericVector cp, NumericVector z, int nrows, int ncols){
      int k = z.size();
      IntegerMatrix mat(nrows, ncols);
      for (int i = 0; i < k; i++){ mat(rp[i],cp[i]) = z[i]; }
      return mat;
    }')
  as_matrix <- function(mat){
    row_pos <- mat@i
    col_pos <- findInterval(seq_along(mat@x)-1, mat@p[-1])
    tmp <- asMatrix(rp=row_pos, cp=col_pos, z=mat@x, nrows=mat@Dim[1], ncols=mat@Dim[2])
    row.names(tmp) <- mat@Dimnames[[1]]; colnames(tmp) <- mat@Dimnames[[2]]
    tmp
  }
  
  sc_chat <- readRDS(file.path(OUTDIR,"scRNA_geneset.RDS"))
  # 可按需剔除未知分组
  sc_chat <- subset(sc_chat, geneset_group != "Unknown")
  sc_counts <- as_matrix(sc_chat[["RNA"]]@counts)
  sc_nUMI   <- colSums(sc_counts)
  cellType  <- data.frame(barcode=colnames(sc_chat), cell_type=sc_chat$geneset_group, check.names = FALSE)
  cell_types <- setNames(as.factor(cellType$cell_type), cellType$barcode)
  reference  <- Reference(sc_counts, cell_types, sc_nUMI)
  rm(sc_counts); gc()
  
  st <- readRDS(ST_RDS)
  coords <- GetTissueCoordinates(st, cols = c("row","col"), scale = NULL)
  colnames(coords) <- c("xcoord","ycoord")
  sp_counts <- as_matrix(st[["Spatial"]]@counts)
  sp_nUMI   <- colSums(sp_counts)
  puck <- SpatialRNA(coords, sp_counts, sp_nUMI)
  
  myRCTD <- create.RCTD(puck, reference)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  save(myRCTD, file = file.path(OUTDIR,"RCTD_geneset.Rdata"))
  
  weight <- as.data.frame(myRCTD@results$weights)
  write.csv(weight, file = file.path(OUTDIR,"weight_RCTD_geneset.csv"), quote = FALSE)
  cat("[OK] RCTD finished & saved.\n")
  
  # 可视化（若需要把权重拼回 ST Seurat 对象）
  st <- st[, rownames(weight)]
  st@meta.data <- cbind(st@meta.data, weight)
  saveRDS(st, file = file.path(OUTDIR,"stRNA_with_weights.rds"))
} else {
  cat("[INFO] 跳过 RCTD（未提供 st_rds）。\n")
}

# ----------------- 9) 空间网络（同型/异型） -----------------
if(nchar(ST_RDS) && file.exists(ST_RDS)){
  my_sp <- readRDS(ST_RDS)
  decon_mtrx <- read.csv(file.path(OUTDIR,"weight_RCTD_geneset.csv"), header = TRUE, row.names = 1, check.names = FALSE)
  if("celltype" %in% colnames(decon_mtrx)) decon_mtrx$celltype <- NULL
  
  # 对齐条形码
  ss <- intersect(colnames(my_sp), rownames(decon_mtrx))
  my_sp <- my_sp[, ss, drop=FALSE]
  decon_mtrx <- decon_mtrx[ss, , drop=FALSE]
  stopifnot(identical(colnames(my_sp), rownames(decon_mtrx)))
  
  # 取第一个切片
  slice <- names(my_sp@images)[1]
  spatial_coord <- data.frame(my_sp@images[[slice]]@coordinates) %>%
    tibble::rownames_to_column("barcodeID") %>%
    dplyr::mutate(imagerow_scaled = imagerow * my_sp@images[[slice]]@scale.factors$lowres,
                  imagecol_scaled = imagecol * my_sp@images[[slice]]@scale.factors$lowres)
  
  # knn 空间网络
  xys <- data.frame(x = spatial_coord$imagecol_scaled, y = spatial_coord$imagerow_scaled, row.names = spatial_coord$barcodeID)
  spotnames <- rownames(xys); names(spotnames) <- seq_len(nrow(xys))
  nNeighbours <- 6; maxdist <- 200; minK <- 4
  knn <- dbscan::kNN(as.matrix(xys[, c("x","y")]), k = nNeighbours)
  net <- data.frame(from = rep(1:nrow(knn$id), nNeighbours),
                    to = as.vector(knn$id),
                    weight = 1/(1 + as.vector(knn$dist)),
                    distance = as.vector(knn$dist))
  net$from <- spotnames[net$from]; net$to <- spotnames[net$to]
  net <- net %>% group_by(from) %>% mutate(rnk = rank(distance)) %>% ungroup()
  spatnet <- subset(net, distance <= maxdist | rnk <= minK)
  spatnet <- cbind(spatnet, setNames(xys[spatnet$from, 1:2], c("start_x","start_y")))
  spatnet <- cbind(spatnet, setNames(xys[spatnet$to,   1:2], c("end_x","end_y")))
  
  # ——同型网络（基于某一细胞占比阈值）——
  thr_same <- 0.1
  same_name <- paste0(LOW_LAB,"Mali")
  if(!(same_name %in% colnames(decon_mtrx))){
    message("[WARN] 找不到列：", same_name, "，将尝试使用列名包含 'Mali' 的第一列。")
    mali_cols <- grep("Mali", colnames(decon_mtrx), value = TRUE)
    stopifnot(length(mali_cols) >= 1); same_name <- mali_cols[1]
  }
  data_same <- decon_mtrx[, same_name, drop = FALSE]
  colnames(data_same) <- "prop"
  data_same <- subset(data_same, prop > thr_same)
  degree <- matrix(NA_real_, nrow = nrow(data_same), ncol = 1, dimnames = list(rownames(data_same), "degree"))
  for(i in rownames(data_same)){
    su <- spatnet[spatnet$from == i, ]
    degree[i,1] <- sum(decon_mtrx[su$to, same_name] > thr_same, na.rm = TRUE)
  }
  degree <- na.omit(as.data.frame(degree))
  saveRDS(degree, file = file.path(OUTDIR,"homotypic_degree.rds"))
  
  if(nchar(LOWRES_IMG) && file.exists(LOWRES_IMG)){
    img <- png::readPNG(LOWRES_IMG)
    img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1,"npc"), height = grid::unit(1,"npc"))
    seg <- spatnet[spatnet$from %in% rownames(degree) & spatnet$to %in% rownames(degree), ]
    coord_deg <- spatial_coord[spatial_coord$barcodeID %in% rownames(degree), ]
    coord_deg$degree <- degree[coord_deg$barcodeID, "degree"]
    p_same <- ggplot() +
      annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
      geom_segment(data = seg, aes(x = start_x, y = start_y, xend = end_x, yend = end_y), color = "white", linewidth = 0.2) +
      geom_point(data = coord_deg, aes(x = imagecol_scaled, y = imagerow_scaled, color = degree), size = 0.9) +
      ylim(nrow(img), 0) + xlim(0, ncol(img)) + theme_void() + coord_fixed() +
      scale_y_reverse() + ggtitle("Homotypic degree") +
      scale_color_gradient(low="#E5F4F8", high="#bc3c29")
    ggsave(file.path(OUTDIR,"homotypic_degree.pdf"), p_same, width=6, height=6)
  }
  
  # ——异型网络（CELL1_HET vs CELL2_HET）——
  stopifnot(CELL1_HET %in% colnames(decon_mtrx))
  stopifnot(CELL2_HET %in% colnames(decon_mtrx))
  thr_het <- 0.2
  data1 <- decon_mtrx[, CELL1_HET, drop=FALSE] %>% setNames("prop1") %>% subset(prop1 > thr_het)
  data2 <- decon_mtrx[, CELL2_HET, drop=FALSE] %>% setNames("prop2") %>% subset(prop2 > thr_het)
  
  seg_het <- spatnet[spatnet$from %in% rownames(data1) & spatnet$to %in% rownames(data2), ]
  if(nrow(seg_het) > 0 && nchar(LOWRES_IMG) && file.exists(LOWRES_IMG)){
    img <- png::readPNG(LOWRES_IMG)
    img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1,"npc"), height = grid::unit(1,"npc"))
    c1 <- spatial_coord[spatial_coord$barcodeID %in% unique(seg_het$from), ]
    c2 <- spatial_coord[spatial_coord$barcodeID %in% unique(seg_het$to), ]
    p_het <- ggplot() +
      annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
      geom_segment(data = seg_het, aes(x = start_x, y = start_y, xend = end_x, yend = end_y), color = "white", linewidth = 0.2) +
      geom_point(data = c1, aes(x = imagecol_scaled, y = imagerow_scaled), size = 1.0, color = "yellow") +
      geom_point(data = c2, aes(x = imagecol_scaled, y = imagerow_scaled), size = 1.0, color = "purple") +
      ylim(nrow(img), 0) + xlim(0, ncol(img)) + theme_void() + coord_fixed() + scale_y_reverse() +
      ggtitle(paste0("Heterotypic: ", CELL1_HET, " ↔ ", CELL2_HET))
    ggsave(file.path(OUTDIR,"heterotypic_links.pdf"), p_het, width=6, height=6)
  }
  
  # 邻居富集图（以 CELL2_HET 为邻居）
  deg2 <- matrix(NA_real_, nrow = nrow(data1), ncol = 1, dimnames = list(rownames(data1), "neighbor_cell2"))
  for (i in rownames(data1)){
    su <- spatnet[spatnet$from == i, ]
    deg2[i,1] <- sum(decon_mtrx[c(su$to, su$from), CELL2_HET], na.rm = TRUE)
  }
  deg2 <- as.data.frame(deg2)
  deg2$enrich <- deg2$neighbor_cell2 / max(deg2$neighbor_cell2, na.rm = TRUE)
  coord_en <- spatial_coord[spatial_coord$barcodeID %in% rownames(deg2), ]
  coord_en$enrich <- deg2[coord_en$barcodeID, "enrich"]
  if(nchar(LOWRES_IMG) && file.exists(LOWRES_IMG)){
    img <- png::readPNG(LOWRES_IMG)
    img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1,"npc"), height = grid::unit(1,"npc"))
    p_en <- ggplot() +
      annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
      geom_point(data = coord_en, aes(x = imagecol_scaled, y = imagerow_scaled, color = enrich, alpha = enrich), size = 1.0) +
      ylim(nrow(img), 0) + xlim(0, ncol(img)) + theme_void() + coord_fixed() + scale_y_reverse() +
      scale_color_gradient2(low = "blue", mid = "white", high = "yellow", midpoint = 0.5) +
      guides(alpha = "none") + ggtitle(paste0("Neighbor enrichment: ", CELL2_HET))
    ggsave(file.path(OUTDIR,"neighbor_enrichment.pdf"), p_en, width=6, height=6)
  }
} else {
  cat("[INFO] 跳过空间网络（未提供 st_rds）。\n")
}

# ----------------- 10) mistyR -----------------
if(nchar(RCTD_RDATA) && file.exists(RCTD_RDATA) && nchar(ST_RDS) && file.exists(ST_RDS)){
  load(RCTD_RDATA)  # myRCTD
  st <- readRDS(ST_RDS)
  
  prop <- as.matrix(myRCTD@results$weights)
  colnames(prop) <- stringr::str_replace(colnames(prop), "\\+", "_pos_")
  st[["RCTD_props"]] <- CreateAssayObject(t(prop))
  DefaultAssay(st) <- "RCTD_props"
  useful_features <- rownames(st)
  
  # 并行
  future::plan(future::multisession)
  
  run_colocalization <- function(slide, assay, useful_features, out_label, misty_out_alias="./misty_geneset/cm_"){
    view_assays   <- list(main=assay, juxta=assay, para=assay)
    view_features <- list(main=useful_features, juxta=useful_features, para=useful_features)
    view_types    <- list(main="intra", juxta="juxta", para="para")
    view_params   <- list(main=NULL, juxta=5, para=15)
    misty_out <- paste0(misty_out_alias, out_label, "_", assay)
    dir.create(dirname(misty_out), showWarnings = FALSE, recursive = TRUE)
    mistyR::run_misty_seurat(
      visium.slide = slide,
      view.assays  = view_assays,
      view.features= view_features,
      view.types   = view_types,
      view.params  = view_params,
      spot.ids     = NULL,
      out.alias    = misty_out
    )
    misty_out
  }
  
  mout <- run_colocalization(slide = st, assay = "RCTD_props", useful_features = useful_features,
                             out_label = "stRNA", misty_out_alias = file.path(OUTDIR,"misty/"))
  res <- mistyR::collect_results(mout)
  ggsave(file.path(OUTDIR,"misty_improvement.pdf"), mistyR::plot_improvement_stats(res), width=7, height=5)
  ggsave(file.path(OUTDIR,"misty_view_contrib.pdf"), mistyR::plot_view_contributions(res), width=7, height=5)
  ggsave(file.path(OUTDIR,"misty_intra_heatmap.pdf"), mistyR::plot_interaction_heatmap(res, "intra", cutoff=0), width=7, height=6)
  ggsave(file.path(OUTDIR,"misty_intra_communities.pdf"), mistyR::plot_interaction_communities(res, "intra", cutoff=0.5), width=7, height=6)
  ggsave(file.path(OUTDIR,"misty_juxta_heatmap.pdf"), mistyR::plot_interaction_heatmap(res, "juxta_5", cutoff=0), width=7, height=6)
  ggsave(file.path(OUTDIR,"misty_para_heatmap.pdf"),  mistyR::plot_interaction_heatmap(res, "para_15", cutoff=0), width=7, height=6)
} else {
  cat("[INFO] 跳过 mistyR（缺少 RCTD 或 ST 数据）。\n")
}

cat("[DONE] All finished. Outputs in: ", OUTDIR, "\n")
