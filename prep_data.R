
# -------------------- 0) 依赖与工具 ----------------------
required_cran <- c("data.table","tibble","stringr","tidyverse","limma")
required_bioc <- c("GEOquery","BiocManager")
opt_pkg       <- c("tinyarray")  # 可选：用于 Ensembl → Symbol 转换

install_if_missing <- function(pkgs, bioc = FALSE){
  inst <- rownames(installed.packages())
  need <- setdiff(pkgs, inst)
  if(length(need)){
    if(bioc) {
      if(!"BiocManager" %in% inst) install.packages("BiocManager", repos="https://cloud.r-project.org")
      BiocManager::install(need, update = FALSE, ask = FALSE)
    } else {
      install.packages(need, repos="https://cloud.r-project.org")
    }
  }
}
install_if_missing(required_cran, bioc = FALSE)
install_if_missing(required_bioc, bioc = FALSE)
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(data.table); library(tibble)
  library(stringr); library(tidyverse)
})
# 可选包
if(!"tinyarray" %in% rownames(installed.packages())){
  message("[WARN] 'tinyarray' 未安装，将跳过 TCGA ID 转换。建议安装后重跑。")
} else {
  suppressPackageStartupMessages(library(tinyarray))
}

`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

# -------------------- 1) CLI 参数 ------------------------
# 示例：
# Rscript src/prep_data.R \
#   gse_dir=./data/geo  tcga_file=./data/TCGA-LIHC.htseq_fpkm.tsv.gz \
#   outdir=./outputs    gpl_filter=GPL571
args <- commandArgs(trailingOnly = TRUE)
kv   <- strsplit(args, "=")
opts <- setNames(vapply(kv, function(x) if(length(x)==2) x[2] else NA_character_, ""),
                 vapply(kv, function(x) x[1], ""))
OUTDIR     <- opts[["outdir"]]     %||% "./outputs"
GSE_DIR    <- opts[["gse_dir"]]    %||% "./data/geo"
TCGA_FILE  <- opts[["tcga_file"]]  %||% "./data/TCGA-LIHC.htseq_fpkm.tsv.gz"  # TODO: 确认你的真实文件名
GPL_FILTER <- opts[["gpl_filter"]] %||% "GPL571"  # 仅保留该平台；设为空字符串可保留全部平台
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# 日志
logf <- file.path(OUTDIR, sprintf("prep_data_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
zz <- file(logf, open = "wt"); sink(zz, type = "output"); sink(zz, type = "message")
time_start <- Sys.time()
on.exit({
  cat("\n[INFO] Session info:\n"); print(sessionInfo())
  cat(sprintf("[INFO] Elapsed: %.2f secs\n", as.numeric(difftime(Sys.time(), time_start, units="secs"))))
  sink(type="message"); sink(type="output"); close(zz)
}, add = TRUE)
cat("[INFO]", as.character(Sys.time()), "Start prep_data\n")

# -------------------- 2) GEO: GSE14323 -------------------
# 备注：GSE14323 同时含多个平台（如 GPL571、GPL96）。为避免探针注释混用，
# 默认仅保留 GPL571（可通过 gpl_filter 参数修改）。
Sys.setenv('VROOM_CONNECTION_SIZE' = 131072 * 2)

gse_id <- "GSE14323"
dir.create(GSE_DIR, showWarnings = FALSE, recursive = TRUE)
cat("[INFO] Fetch GEO series:", gse_id, "into:", GSE_DIR, "\n")
gset <- getGEO(gse_id, destdir = GSE_DIR, getGPL = FALSE)
stopifnot(length(gset) >= 1)

eset   <- gset[[1]]
expr0  <- exprs(eset)                 # probe-level
pdata  <- pData(eset)
fvars  <- fData(eset)                 # 可能为空（getGPL=FALSE 时）

# 可选：按平台筛选
if(nchar(GPL_FILTER)){
  keep <- which(pdata$platform_id == GPL_FILTER | pdata$"platform_id" == GPL_FILTER)
  if(length(keep) == 0) stop("[ERROR] 在 GSE14323 中找不到平台：", GPL_FILTER)
  pdata <- pdata[keep, , drop=FALSE]
  expr0 <- expr0[, rownames(pdata), drop=FALSE]
  cat("[INFO] Keep platform:", GPL_FILTER, " Samples:", nrow(pdata), "\n")
} else {
  cat("[INFO] Keep all platforms. Samples:", nrow(pdata), "\n")
}

# 自动分组（不再手写 19/38）
# 依据 source_name_ch1 的标签归类 normal / tumor
norm_flag  <- grepl("normal", pdata$source_name_ch1, ignore.case = TRUE)
tumor_flag <- grepl("hcc|tumou?r|cancer|carcinoma", pdata$source_name_ch1, ignore.case = TRUE)
grp <- ifelse(norm_flag, "normal",
              ifelse(tumor_flag, "tumor", "other"))
# 仅保留 normal/tumor，去除 other
sel <- grp %in% c("normal","tumor")
pdata <- pdata[sel, , drop=FALSE]
expr0 <- expr0[, rownames(pdata), drop=FALSE]
group_geo <- factor(grp[sel], levels = c("normal","tumor"))
stopifnot(ncol(expr0) == length(group_geo))
cat("[INFO] GEO groups:", table(group_geo), "\n")

# 标准化
# 注意：series matrix 通常已 log 处理；此处示范用 limma::normalizeBetweenArrays 做分位数归一化
expr_geo <- normalizeBetweenArrays(expr0)
expr_geo <- as.data.frame(expr_geo)

# 注释（平台注释文件 *.txt）
# TODO: 请把 GPL571 注释文件放到 GSE_DIR，并确认列名（ID 与 Gene Symbol）
anno_file <- file.path(GSE_DIR, "GPL571-17391.txt")
if(!file.exists(anno_file)){
  stop("[ERROR] 缺少平台注释文件：", anno_file,
       "\n请到 GEO 平台页面下载 full table（TXT）并放置到：", GSE_DIR)
}
anno <- fread(anno_file, data.table = FALSE, sep = "\t")
# 兼容列名（不同版本可能大小写或空格不同）
id_col    <- c("ID","Id","ProbeID")
sym_col   <- c("Gene Symbol","Gene Symbol(s)","GeneSymbol")
id_col    <- id_col[id_col %in% colnames(anno)][1]
sym_col   <- sym_col[sym_col %in% colnames(anno)][1]
if(is.na(id_col) || is.na(sym_col))
  stop("[ERROR] 在注释文件中未找到 ID 与 Gene Symbol 列，请检查列名。")
anno <- anno[, c(id_col, sym_col)]
colnames(anno) <- c("ID","Symbol")
anno <- anno[anno$Symbol != "" & !is.na(anno$Symbol), ]
rownames(anno) <- anno$ID

# 对齐并聚合到 gene symbol
common <- intersect(rownames(anno), rownames(expr_geo))
expr_geo <- expr_geo[common, , drop=FALSE]
anno     <- anno[common, , drop=FALSE]
expr_geo$Symbol <- anno$Symbol
expr_geo <- aggregate(x = expr_geo[, setdiff(colnames(expr_geo),"Symbol")],
                      by = list(Symbol = expr_geo$Symbol),
                      FUN = mean)
expr_geo <- column_to_rownames(expr_geo, var = "Symbol")

# 保存
geo_rdata <- file.path(OUTDIR, "GSE14323_GPL571_expr_group.Rdata")
exprSet2  <- expr_geo
group_list2 <- group_geo
save(exprSet2, group_list2, file = geo_rdata)
cat("[OK] Saved GEO:", geo_rdata, "\n")

# -------------------- 3) TCGA: LIHC ---------------------
# 期望输入：Xena 导出的 HTSeq-FPKM 表（行=基因，列=样本；第一列为 Ensembl_ID）
# 注意：你的原代码文件名是 “TCGA-KIRC.htseq_fpkm.tsv.gz”，但对象命名为 LIHC。
# 这里统一改为 LIHC；请确认并把文件名改成真实路径。
if(!file.exists(TCGA_FILE)){
  stop("[ERROR] 未找到 TCGA FPKM 文件：", TCGA_FILE,
       "\n请从 Xena 下载 LIHC 的 HTSeq-FPKM 数据并修正 tcga_file 参数。")
}
tcga_dt <- fread(file = TCGA_FILE, header = TRUE, sep = "\t")
if(!"Ensembl_ID" %in% colnames(tcga_dt))
  stop("[ERROR] TCGA 文件缺少列：Ensembl_ID")
expr_tcga <- column_to_rownames(tcga_dt, var = "Ensembl_ID")
expr_tcga <- as.matrix(expr_tcga)

# FPKM → TPM
fpkmToTpm <- function(fpkm) exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
expr_tcga <- 2^expr_tcga - 1  # 若原表为 log2(FPKM+1)
expr_tpm  <- apply(expr_tcga, 2, fpkmToTpm)
expr_tpm  <- as.data.frame(expr_tpm)

# （可选）ID 转换：Ensembl → Gene Symbol
if("tinyarray" %in% rownames(installed.packages())){
  # 仅保留 mRNA，自动转符号
  expr_tpm_sym <- tryCatch({
    trans_exp(expr_tpm, mrna_only = TRUE)
  }, error = function(e){
    message("[WARN] tinyarray::trans_exp 失败，保留 Ensembl 行名。Error: ", e$message)
    expr_tpm
  })
} else {
  expr_tpm_sym <- expr_tpm
}

# 分组：TCGA barcode 第 14–15 位 <10 为肿瘤，>=10 为正常
k1 <- str_starts(colnames(expr_tpm_sym), "TCGA")
k2 <- suppressWarnings(as.numeric(str_sub(colnames(expr_tpm_sym), 14, 15))) < 10
group_tcga <- ifelse(k1 & k2, "tumor", "normal")
group_tcga <- factor(group_tcga, levels = c("normal","tumor"))
cat("[INFO] TCGA groups:", table(group_tcga), "\n")

# 保存
tcga_rdata <- file.path(OUTDIR, "LIHC_TPM_expr_group.Rdata")
exprSet <- expr_tpm_sym
group_list <- group_tcga
save(exprSet, group_list, file = tcga_rdata)
cat("[OK] Saved TCGA:", tcga_rdata, "\n")

cat("[DONE] All finished. Outputs in:", OUTDIR, "\n")
