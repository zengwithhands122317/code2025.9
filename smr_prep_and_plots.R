
options(stringsAsFactors = FALSE)

# ------------ 0) 依赖安装与加载 --------------
required_cran <- c("vcfR","tidyverse","data.table","tibble","stringr","ggplot2","magick")
install_if_missing <- function(pkgs){
  inst <- rownames(installed.packages())
  need <- setdiff(pkgs, inst)
  if(length(need)) install.packages(need, repos = "https://cloud.r-project.org")
}
install_if_missing(required_cran)
suppressPackageStartupMessages({
  library(vcfR); library(tidyverse); library(data.table); library(tibble)
  library(stringr); library(ggplot2); library(magick)
})

`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

# ------------ 1) 命令行参数 -------------------
# 示例：
# Rscript src/smr_prep_and_plots.R \
#   vcf=./input/ukb-b-1316.vcf \
#   smr_out=./smr-1.3.1-win-x86_64/output/output.smr \
#   gtf_rdata=./gtf_df.Rdata \
#   apoptosis_gmt=./APOPTOSIS.txt \
#   plot_helper=./plot/plot_SMR.r \
#   locus_file=./smr-1.3.1-win-x86_64/output/plot/myplot.ENSG00000132906.txt \
#   outdir=./outputs \
#   n_const=1673211
args <- commandArgs(trailingOnly = TRUE)
kv   <- strsplit(args, "=")
opts <- setNames(vapply(kv, function(x) if(length(x)==2) x[2] else NA_character_, ""),
                 vapply(kv, function(x) x[1], ""))
VCF_FILE    <- opts[["vcf"]]          %||% stop("参数缺失：vcf")
SMR_RESULT  <- opts[["smr_out"]]      %||% stop("参数缺失：smr_out")
GTF_RDATA   <- opts[["gtf_rdata"]]    %||% stop("参数缺失：gtf_rdata")
APOP_FILE   <- opts[["apoptosis_gmt"]]%||% ""
PLOT_HELPER <- opts[["plot_helper"]]  %||% ""
LOCUS_FILE  <- opts[["locus_file"]]   %||% ""
OUTDIR      <- opts[["outdir"]]       %||% "./outputs"
N_CONST     <- opts[["n_const"]]      %||% NA_character_

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# 日志
logf <- file.path(OUTDIR, sprintf("smr_prep_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
zz <- file(logf, open = "wt"); sink(zz, type = "output"); sink(zz, type = "message")
t0 <- Sys.time()
on.exit({
  cat("\n[INFO] Session info:\n"); print(sessionInfo())
  cat(sprintf("[INFO] Elapsed: %.2f secs\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))
  sink(type="message"); sink(type="output"); close(zz)
}, add = TRUE)
cat("[INFO]", as.character(Sys.time()), "Start SMR prep+plots\n")

# ------------ 2) 读取 VCF 并整理为 SMR .ma ----------
stopifnot(file.exists(VCF_FILE))
vcf <- vcfR::read.vcfR(VCF_FILE, verbose = FALSE)
cat("[INFO] VCF loaded. Variants:", nrow(vcf@fix), " Samples:", vcf@n samples, "\n")

# FIX 区（SNP/ALT/REF）
fix <- as.data.frame(vcf@fix)
stopifnot(all(c("ID","ALT","REF") %in% colnames(fix)))
fix <- fix %>% dplyr::select(ID, ALT, REF)  # A1=ALT, A2=REF

# GT 区：自动找到唯一样本列；VCF GT 第一列是 "FORMAT"
gt <- as.data.frame(vcf@gt, check.names = FALSE)
if(!"FORMAT" %in% colnames(gt)) stop("[ERROR] VCF 缺少 FORMAT 列")
sample_cols <- setdiff(colnames(gt), "FORMAT")
if(length(sample_cols) != 1){
  stop("[ERROR] 预期单样本 VCF，但检测到样本列：", paste(sample_cols, collapse=", "))
}
sample_col <- sample_cols[1]
cat("[INFO] Using sample column:", sample_col, "\n")

# 分离 sample 列；FORMAT 如 "ES:SE:LP:AF:ID(:SS?)"
fmt_fields <- unique(gt$FORMAT)
fmt_ref <- strsplit(fmt_fields[1], ":", fixed = TRUE)[[1]]
# 兼容字段顺序变动：我们把该顺序作为参考
gt_sep <- tidyr::separate(gt, col = all_of(sample_col), into = paste0("F", seq_along(fmt_ref)), sep = ":", fill = "right", remove = FALSE)

# 建立字段名映射（在常见的 GWAS-SMR导出中，常见标签为 ES,SE,LP,AF,ID,SS）
canonical <- c("ES","SE","LP","AF","ID","SS","N")
names(canonical) <- canonical
# 用 FORMAT 的位置推断字段名（若包含这些关键字）
name_guess <- tolower(fmt_ref)
map_field <- function(target){
  # 在 fmt_ref 中寻找包含目标名的字段位置
  ix <- which(grepl(tolower(target), name_guess))
  if(length(ix) >= 1) paste0("F", ix[1]) else NA_character_
}
F_ES <- map_field("ES"); F_SE <- map_field("SE"); F_LP <- map_field("LP")
F_AF <- map_field("AF"); F_ID <- map_field("ID"); F_SS <- map_field("SS")
F_N  <- map_field("N")   # 偶尔用 N 表示样本量

# 兜底：如果没检测到 ES/SE/LP/AF/ID，则按常见 1:2:3:4:5 顺序
F_ES <- F_ES %||% "F1"; F_SE <- F_SE %||% "F2"; F_LP <- F_LP %||% "F3"; F_AF <- F_AF %||% "F4"; F_ID <- F_ID %||% "F5"

# 组装 MR 表
MR_rt <- data.frame(
  ID   = gt_sep[[F_ID]],
  freq = gt_sep[[F_AF]],
  beta = gt_sep[[F_ES]],
  se   = gt_sep[[F_SE]],
  p    = gt_sep[[F_LP]],
  stringsAsFactors = FALSE
)

# n 的来源：优先 SS/N 字段，否则用 n_const 参数（若仍缺失则报错）
if(!is.na(F_SS) && F_SS %in% colnames(gt_sep)){
  MR_rt$n <- suppressWarnings(as.numeric(gt_sep[[F_SS]]))
} else if(!is.na(F_N) && F_N %in% colnames(gt_sep)){
  MR_rt$n <- suppressWarnings(as.numeric(gt_sep[[F_N]]))
} else if(!is.na(N_CONST)){
  MR_rt$n <- as.numeric(N_CONST)
} else {
  stop("[ERROR] 未找到 SS/N 字段且未提供 n_const 参数。请通过 n_const=... 传入样本量。")
}

# 类型转换 & LP→p
num_cols <- c("freq","beta","se","p","n")
for(cl in num_cols) MR_rt[[cl]] <- suppressWarnings(as.numeric(MR_rt[[cl]]))
# p 可能是 -log10(p)（LP），若范围>1 则判定为 LP 并转换
if(max(MR_rt$p, na.rm = TRUE) > 1){
  MR_rt$p <- 10^(-MR_rt$p)
  cat("[INFO] Converted LP to p using p=10^(-LP)\n")
}
stopifnot(all(MR_rt$p >= 0 & MR_rt$p <= 1, na.rm = TRUE))

# 合并 FIX 与效应列；A1=ALT, A2=REF
MR <- cbind(fix, MR_rt) %>%
  dplyr::select(ID, ALT, REF, freq, beta, se, p, n) %>%
  stats::na.omit()
colnames(MR) <- c("SNP","A1","A2","freq","b","se","p","n")

# 导出 .ma
ma_file <- file.path(OUTDIR, "mygwas.ma")
write.table(MR, ma_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
cat("[OK] SMR .ma written:", ma_file, "\n")

# ------------ 3) 读取 SMR 结果并做注释 --------------
stopifnot(file.exists(GTF_RDATA))
load(GTF_RDATA)  # 需要提供对象 gtf_df: 含 gene_id, gene_biotype 等列
if(!exists("gtf_df")) stop("[ERROR] gtf_df 不存在于 ", GTF_RDATA)

stopifnot(file.exists(SMR_RESULT))
gene <- read.table(SMR_RESULT, header = TRUE, sep = "", check.names = FALSE)
if(!all(c("probeID","p_SMR") %in% colnames(gene)))
  stop("[ERROR] SMR 结果缺少必要列：probeID / p_SMR")

# 去重并按 gtf_df 对齐
gtf_df <- gtf_df[!duplicated(gtf_df$gene_id), ]
gene   <- gene[!duplicated(gene$probeID), ]
rownames(gene) <- gene$probeID
common_ids <- intersect(gtf_df$gene_id, rownames(gene))
gene <- gene[common_ids, , drop = FALSE]
gtf_df_sub <- gtf_df[match(rownames(gene), gtf_df$gene_id), , drop=FALSE]

# 添加 biotype 并筛选
gene$type <- gtf_df_sub$gene_biotype
gene_sig  <- subset(gene, p_SMR < 0.05)
gene_sig  <- subset(gene_sig, type == "protein_coding")

csv_out <- file.path(OUTDIR, "causal_gene_protein.csv")
write.csv(gene_sig, file = csv_out, row.names = FALSE, quote = FALSE)
cat("[OK] Significant protein_coding genes saved:", csv_out, "\n")

# ------------ 4) 简单富集：与凋亡/凋亡相关基因交集 ------------
if(nchar(APOP_FILE) && file.exists(APOP_FILE)){
  senes <- read.table(APOP_FILE, header = FALSE)
  senes <- senes[[1]]
  # 兼容列名：SMR 输出列中常见为 Gene 或 Symbol；若没有则尝试 probeID→gtf_df 衔接
  gene_col <- if("Gene" %in% colnames(gene_sig)) "Gene" else if("Symbol" %in% colnames(gene_sig)) "Symbol" else NA
  if(is.na(gene_col)){
    message("[WARN] SMR 结果无 Gene/Symbol 列，将尝试从 gtf_df 中映射 symbol（若存在）。")
    if("gene_name" %in% colnames(gtf_df_sub)){
      gene_names <- gtf_df_sub$gene_name
    } else if("symbol" %in% colnames(gtf_df_sub)){
      gene_names <- gtf_df_sub$symbol
    } else {
      gene_names <- rownames(gene_sig)
    }
  } else {
    gene_names <- gene_sig[[gene_col]]
  }
  inter_genes <- intersect(senes, gene_names)
  cat("[INFO] Intersection with APOPTOSIS set (n):", length(inter_genes), "\n")
  if(length(inter_genes)){
    write.table(inter_genes, file = file.path(OUTDIR, "intersection_APOPTOSIS.txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("[OK] Saved intersection list.\n")
  }
} else {
  cat("[INFO] Skipped APOPTOSIS intersection (file missing or not provided).\n")
}

# ------------ 5) Locus/Effect 绘图 ---------------------
if(nchar(PLOT_HELPER) && file.exists(PLOT_HELPER)){
  source(PLOT_HELPER)
  if(nchar(LOCUS_FILE) && file.exists(LOCUS_FILE)){
    smrdata <- ReadSMRData(LOCUS_FILE)
    # Locus plot
    p1 <- SMRLocusPlot(data = smrdata,
                       smr_thresh = 0.05,
                       heidi_thresh = 0.05,
                       plotWindow = 1000,
                       anno_selfdef = FALSE)
    ggsave(filename = file.path(OUTDIR, "SMR_locus.pdf"), plot = p1, width = 12, height = 12, dpi = 300)
    cat("[OK] Saved locus plot: SMR_locus.pdf\n")
    
    # Effect plot（若脚本支持）
    p2 <- tryCatch({
      SMREffectPlot(data = smrdata, trait_name = "HCC")
    }, error = function(e){
      message("[WARN] SMREffectPlot 报错：", e$message); NULL
    })
    if(!is.null(p2)){
      ggsave(filename = file.path(OUTDIR, "SMR_effect.pdf"), plot = p2, width = 12, height = 12, dpi = 300)
      cat("[OK] Saved effect plot: SMR_effect.pdf\n")
    }
  } else {
    cat("[INFO] 跳过绘图：未提供或找不到 locus_file。\n")
  }
} else {
  cat("[INFO] 跳过绘图：未提供或找不到 plot_helper。\n")
}

cat("[DONE] All finished. Outputs in:", OUTDIR, "\n")
