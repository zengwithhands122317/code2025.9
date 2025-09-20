

options(stringsAsFactors = FALSE)
`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

needed <- c(
  "data.table","tibble","dplyr","tidyr","stringr","readr",
  "Seurat","tinyarray","sva","survival","survminer","caret",
  "glmnet","ggplot2","ggpubr","viridis","scales","forestplot",
  "timeROC","ggridges","IOBR","estimate"
)
inst <- rownames(installed.packages())
miss <- setdiff(needed, inst)
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(data.table); library(tibble); library(dplyr); library(tidyr); library(stringr); library(readr)
  library(Seurat); library(tinyarray); library(sva); library(survival); library(survminer); library(caret)
  library(glmnet); library(ggplot2); library(ggpubr); library(viridis); library(scales); library(forestplot)
  library(timeROC); library(ggridges); library(IOBR); library(estimate)
})

# -------------------- 参数 --------------------
# 关键输入（训练集+生存）：
#   train_expr_rdata: 包含 exprSet, group_list 的 Rdata（如 LIHC_tpm.Rdata / KIRC_tpm.Rdata）
#   train_survival_tsv: 生存表，含 OS.time / OS 或 futime / fustat
# 验证集（表达+生存）：
#   vali_expr_file: 行为 gene，列为样本（或 geneNames + 样本），脚本会自动汇总平均重复基因
#   vali_survival_file: 生存信息，含列名 Sample、OS.time、OS（脚本会自动对齐校正）
#
# 上游“差异基因候选”来源二选一：
#   marker_rds: 直接提供 .RData（包含对象 df，行名为基因）
#   or 从 Seurat RDS + 分组方式自动计算：
#     sc_mode = "gene" / "geneset"
#     sc_rds_gene / sc_rds_geneset
#     group_col_gene = "gene_group"（如 CASP9+Mali vs CASP9-Mali）
#     group_col_geneset = "geneset_group"（如 ApoptosishighMali vs ApoptosislowMali）
#     ident1_name = 目标分组名（如 "CASP9-Mali" 或 "ApoptosislowMali"）
#
# 其它：
#   outdir: 输出目录
#   tumor_only: TRUE 仅取肿瘤样本（从 train expr 的 group_list 过滤）
#   log2_add1: TRUE 对 FPKM/TPM 取 log2(x+1)
#   lfc_thr: FindMarkers 的 logFC 阈值
#   remove_ribo: TRUE 去除 RPS/RPL/MRP 核糖体基因
args <- commandArgs(trailingOnly = TRUE)
kv   <- strsplit(args, "=")
opt  <- setNames(vapply(kv, function(x) if(length(x)==2) x[2] else NA_character_, ""),
                 vapply(kv, function(x) x[1], ""))

OUTDIR             <- opt[["outdir"]]              %||% "./bulk_outputs"
TRAIN_EXPR_RDATA   <- opt[["train_expr_rdata"]]    %||% "E:/article4/smrsp1/LIHC_tpm.Rdata"
TRAIN_SURVIVAL_TSV <- opt[["train_survival_tsv"]]  %||% "TCGA-KIRC.survival.tsv"  # 示例
VALI_EXPR_FILE     <- opt[["vali_expr_file"]]      %||% "E:/article4/smrsp2/ccRCC数据集/E-MTAB-1980/EMT-symbol.txt"
VALI_SURV_FILE     <- opt[["vali_survival_file"]]  %||% "E:/article4/smrsp2/ccRCC数据集/E-MTAB-1980/EMT1980-sur.txt"

# 差异基因候选来源
MARKER_RDATA       <- opt[["marker_rdata"]]        %||% ""   # 若提供则优先；否则用 sc_* 自动算
SC_MODE            <- opt[["sc_mode"]]             %||% "geneset"  # "gene" 或 "geneset"
SC_RDS_GENE        <- opt[["sc_rds_gene"]]         %||% "E:/article4/smrsp1/scRNA_CASP9.RDS"
SC_RDS_GENESET     <- opt[["sc_rds_geneset"]]      %||% "E:/article4/smrsp1/scRNA_geneset.RDS"
GROUP_COL_GENE     <- opt[["group_col_gene"]]      %||% "gene_group"
GROUP_COL_GENESET  <- opt[["group_col_geneset"]]   %||% "geneset_group"
IDENT1_NAME        <- opt[["ident1_name"]]         %||% "ApoptosislowMali"  # 或 "CASP9-Mali"

# 其它开关
TUMOR_ONLY         <- as.logical(opt[["tumor_only"]]  %||% "TRUE")
LOG2_ADD1          <- as.logical(opt[["log2_add1"]]   %||% "TRUE")
LFC_THR            <- as.numeric(opt[["lfc_thr"]]     %||% "0.5")
REMOVE_RIBO        <- as.logical(opt[["remove_ribo"]] %||% "TRUE")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(OUTDIR, sprintf("bulk_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
zz <- file(logf, open="wt"); sink(zz, type="output"); sink(zz, type="message")
t0 <- Sys.time()
on.exit({
  cat("\n[INFO] Session:\n"); print(sessionInfo())
  cat(sprintf("[INFO] Elapsed: %.1f sec\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))
  sink(type="message"); sink(type="output"); close(zz)
}, add=TRUE)
cat("[INFO] Bulk pipeline started.\n")

# -------------------- 工具函数 --------------------
stopifnot_msg <- function(ok, msg){ if(!ok) stop(msg, call.=FALSE) }

read_train_expr <- function(rdata_file, tumor_only = TRUE, log2_add1 = TRUE){
  stopifnot_msg(file.exists(rdata_file), paste("Rdata not found:", rdata_file))
  env <- new.env(); load(rdata_file, envir = env)
  stopifnot_msg(exists("exprSet", envir = env), "Rdata must contain exprSet")
  expr <- get("exprSet", envir = env)
  stopifnot_msg(is.matrix(expr) || is.data.frame(expr), "exprSet must be matrix/data.frame")
  expr <- as.data.frame(expr)
  if(exists("group_list", envir = env) && tumor_only){
    gl <- get("group_list", envir = env)
    md <- data.frame(id = colnames(expr), group_list = gl)
    tumor_ids <- md$id[md$group_list == "tumor"]
    expr <- expr[, intersect(colnames(expr), tumor_ids), drop = FALSE]
  }
  if(log2_add1){
    # 若 expr 已经是 TPM/FPKM 的 log2 值，重复 log2 可能会压缩；这里按你的原流程保留 log2(x+1)
    expr <- log2(expr + 1)
  }
  expr
}

read_train_survival <- function(surv_file){
  stopifnot_msg(file.exists(surv_file), paste("Survival file not found:", surv_file))
  surv <- read.table(surv_file, header = TRUE, sep = ifelse(grepl("\\.tsv$", surv_file),"\\t","\t"), check.names = FALSE, row.names = 1)
  # 兼容 OS/OS.time 或 futime/fustat
  if(all(c("OS.time","OS") %in% colnames(surv))){
    out <- dplyr::select(surv, "OS.time","OS")
    colnames(out) <- c("futime","fustat")
  } else if(all(c("futime","fustat") %in% colnames(surv))){
    out <- surv[, c("futime","fustat")]
  } else stop("Survival table must contain OS.time/OS or futime/fustat")
  out
}

calc_markers_from_sc <- function(mode, sc_rds, group_col, ident1, lfc_thr, remove_ribo){
  stopifnot_msg(file.exists(sc_rds), paste("scRNA RDS not found:", sc_rds))
  sc <- readRDS(sc_rds)
  stopifnot_msg("celltype" %in% colnames(sc@meta.data), "Seurat meta must contain 'celltype'")
  sc_sub <- subset(sc, celltype == "Cancer cells")
  stopifnot_msg(group_col %in% colnames(sc_sub@meta.data), sprintf("Meta must contain '%s'", group_col))
  Idents(sc_sub) <- sc_sub@meta.data[[group_col]]
  df <- FindMarkers(sc_sub, ident.1 = ident1, only.pos = TRUE, logfc.threshold = lfc_thr)
  if(remove_ribo){
    ribo_idx <- grep("^RP[SL]", rownames(df))
    if(length(ribo_idx) > 0) df <- df[-ribo_idx, , drop = FALSE]
  }
  save(df, file = file.path(OUTDIR, "marker_gene.Rdata"))
  rownames(df)
}

read_vali_expr_surv <- function(expr_file, surv_file){
  stopifnot_msg(file.exists(expr_file), paste("Validation expr not found:", expr_file))
  stopifnot_msg(file.exists(surv_file),  paste("Validation surv not found:", surv_file))
  expr <- read.table(expr_file, header = TRUE, sep = "\t", check.names = FALSE)
  if("geneNames" %in% colnames(expr)){
    expr <- as.data.frame(expr)
    expr <- aggregate(x = expr[, 2:ncol(expr), drop=FALSE], by = list(expr$geneNames), FUN = mean)
    expr <- column_to_rownames(expr, var = "Group.1")
  } else {
    # 假设第一列为基因名
    first <- colnames(expr)[1]
    expr <- column_to_rownames(expr, var = first)
  }
  surv <- read.table(surv_file, header = TRUE, sep = "\t", check.names = FALSE)
  if("Sample" %in% colnames(surv)){
    # 你的原始替换规则：ccRCC-### -> ccRCC.###
    surv$Sample <- gsub("ccRCC-(\\d+)", "ccRCC.\\1", surv$Sample)
    rownames(surv) <- surv$Sample
  } else if(rownames(surv) %in% colnames(expr)) {
    # ok
  } else stop("Validation survival must contain 'Sample' or rownames match expression")
  # 兼容 OS/OS.time → futime/fustat
  if(all(c("OS.time","OS") %in% colnames(surv))){
    surv2 <- dplyr::select(surv, "OS.time","OS"); colnames(surv2) <- c("futime","fustat")
    rownames(surv2) <- rownames(surv)
  } else if(all(c("futime","fustat") %in% colnames(surv))){
    surv2 <- surv[, c("futime","fustat")]
  } else stop("Validation survival must contain OS.time/OS or futime/fustat")
  list(expr = expr, surv = surv2)
}

draw_forest_from_table <- function(tab, outfile){
  dat <- tab
  colnames(dat)[1] <- "gene"
  d1 <- dat %>% dplyr::filter(HR > 1)
  d2 <- dat %>% dplyr::filter(HR < 1)
  hrtable <- rbind(
    c("Increase in Hazard", NA, NA, NA, NA, NA),
    d1,
    c("Reduce in Hazard",  NA, NA, NA, NA, NA),
    d2
  )
  tabletext <- cbind(
    c("Gene","Increase in Hazard", d1$gene, "Reduce in Hazard", d2$gene),
    c("P-value", NA, d1$pvalue, NA, d2$pvalue),
    c("Hazard Ratio", NA, d1$HR, NA, d2$HR)
  )
  pdf(outfile, width=8, height=6)
  forestplot(
    labeltext = tabletext,
    mean  = c(NA, as.numeric(hrtable$HR)),
    lower = c(NA, as.numeric(hrtable$HR.95L)),
    upper = c(NA, as.numeric(hrtable$HR.95H)),
    graph.pos = 4,
    graphwidth = unit(.3, "npc"),
    fn.ci_norm = "fpDrawDiamondCI",
    col = fpColors(box = "#CD534CFF", lines = "#0073C2FF", zero = "#0073C2FF"),
    boxsize = 0.2, lwd.ci = 3, ci.vertices.height = 0.125, ci.vertices = TRUE,
    zero = 1, lwd.zero = 2, xticks = seq(0,4,1),
    lwd.xaxis = 2, xlab = "Hazard Ratio",
    lineheight = unit(0.8,"cm"),
    txt_gp = fpTxtGp(label=gpar(cex=1.2), ticks=gpar(cex=1), xlab=gpar(cex=1.3))
  )
  dev.off()
}

# -------------------- 读取训练集 --------------------
expr_tr <- read_train_expr(TRAIN_EXPR_RDATA, tumor_only = TUMOR_ONLY, log2_add1 = LOG2_ADD1)
surv_tr <- read_train_survival(TRAIN_SURVIVAL_TSV)

# 对齐样本
common <- intersect(colnames(expr_tr), rownames(surv_tr))
expr_tr <- expr_tr[, common, drop=FALSE]
surv_tr <- surv_tr[common, , drop=FALSE]

# -------------------- 读取/计算 marker 列表 --------------------
if(nchar(MARKER_RDATA) && file.exists(MARKER_RDATA)){
  load(MARKER_RDATA)  # 需包含对象 df
  stopifnot_msg(exists("df"), "marker_rdata 必须包含对象 df")
  marker_genes <- rownames(df)
} else {
  if(tolower(SC_MODE) == "gene"){
    marker_genes <- calc_markers_from_sc("gene", SC_RDS_GENE, GROUP_COL_GENE, IDENT1_NAME, LFC_THR, REMOVE_RIBO)
  } else {
    marker_genes <- calc_markers_from_sc("geneset", SC_RDS_GENESET, GROUP_COL_GENESET, IDENT1_NAME, LFC_THR, REMOVE_RIBO)
  }
}
stopifnot_msg(length(marker_genes) > 0, "候选差异基因为空")

# -------------------- 拼接验证集 + 去批次（train vs validation） --------------------
vali <- read_vali_expr_surv(VALI_EXPR_FILE, VALI_SURV_FILE)
expr_va <- vali$expr; surv_va <- vali$surv

# 只保留共同基因
ss <- intersect(rownames(expr_tr), rownames(expr_va))
expr_tr2 <- expr_tr[ss, , drop=FALSE]
expr_va2 <- expr_va[ss, , drop=FALSE]
expr_all <- cbind(expr_tr2, expr_va2)
batch <- factor(c(rep("train", ncol(expr_tr2)), rep("validation", ncol(expr_va2))))

# 批次前PCA
pdf(file.path(OUTDIR, "pca_before_combat.pdf"), width=7, height=6)
tinyarray::draw_pca(exp = expr_all, group_list = batch)
dev.off()

expr_combat <- ComBat(dat = as.matrix(expr_all), batch = batch, par.prior = TRUE, prior.plots = FALSE)
expr_combat <- as.data.frame(expr_combat)

# 批次后PCA
pdf(file.path(OUTDIR, "pca_after_combat.pdf"), width=7, height=6)
tinyarray::draw_pca(exp = expr_combat, group_list = batch)
dev.off()

train_mat <- expr_combat[, colnames(expr_tr2), drop=FALSE]
vali_mat  <- expr_combat[, colnames(expr_va2), drop=FALSE]

save(expr_va2, surv_va, file = file.path(OUTDIR, "validation_inputs.Rdata"))

# -------------------- 单因素 Cox --------------------
train_use <- t(train_mat[rownames(train_mat) %in% marker_genes, , drop=FALSE])
rt <- cbind(surv_tr, as.data.frame(train_use))
rt$futime <- rt$futime / 365  # 单位转为年

pFilter <- 0.05
outTab <- data.frame()
sigGenes <- c("futime","fustat")

for(g in colnames(rt)[-(1:2)]){
  cox <- coxph(Surv(futime, fustat) ~ rt[, g], data = rt)
  sm  <- summary(cox)
  p   <- sm$coefficients[,"Pr(>|z|)"]
  if(!is.na(p) && p < pFilter){
    sigGenes <- c(sigGenes, g)
    outTab <- rbind(outTab, cbind(
      id = g,
      HR = sm$conf.int[,"exp(coef)"],
      `HR.95L` = sm$conf.int[,"lower .95"],
      `HR.95H` = sm$conf.int[,"upper .95"],
      pvalue = p
    ))
  }
}
write.table(outTab, file = file.path(OUTDIR,"uniCox_train.txt"), sep="\t", row.names=FALSE, quote=FALSE)
uniSigExp <- rt[, sigGenes]; uniSigExp <- cbind(id = rownames(uniSigExp), uniSigExp)
write.table(uniSigExp, file = file.path(OUTDIR,"uniSigExp_gene.txt"), sep="\t", row.names=FALSE, quote=FALSE)

# 森林图
uni_tab <- read.table(file.path(OUTDIR,"uniCox_train.txt"), header=TRUE, sep="\t", check.names=FALSE)
if(nrow(uni_tab) > 0){
  draw_forest_from_table(uni_tab, outfile = file.path(OUTDIR,"forest_uniCox.pdf"))
}

# -------------------- LASSO Cox 构建风险评分 --------------------
rt_uni <- read.table(file.path(OUTDIR,"uniSigExp_gene.txt"), header=TRUE, sep="\t", check.names=FALSE, row.names=1)
set.seed(2)
x <- as.matrix(rt_uni[, -(1:2)])
y <- data.matrix(Surv(rt_uni$futime, rt_uni$fustat))
fit   <- glmnet(x, y, family="cox", alpha=1)
cvfit <- cv.glmnet(x, y, family="cox", nfolds=10, alpha=1)

pdf(file.path(OUTDIR,"lasso_path_and_cv.pdf"), width=10, height=4)
par(mfrow=c(1,2)); plot(fit, xvar="lambda", label=FALSE); plot(cvfit)
dev.off()

coef  <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
act   <- coef[index]
lasso_genes <- row.names(coef)[index]
geneCoef <- data.frame(Gene = lasso_genes, Coef = as.numeric(act))
write.table(geneCoef, file = file.path(OUTDIR,"lasso_output_final.txt"), sep="\t", quote=FALSE, row.names=FALSE)

# 准备导出给 Python（可选）
rt_train_for_py <- cbind(as.data.frame(rt_uni[, lasso_genes, drop=FALSE]), Event = rt_uni$fustat, Time = rt_uni$futime)
write.csv(rt_train_for_py, file = file.path(OUTDIR,"train_data.csv"), row.names=FALSE, quote=FALSE)

vali_use <- t(vali_mat[lasso_genes, , drop=FALSE])
rt_vali_for_py <- cbind(as.data.frame(vali_use), Event = surv_va$fustat[rownames(vali_use)], Time = surv_va$futime[rownames(vali_use)])
write.csv(rt_vali_for_py, file = file.path(OUTDIR,"vali_data.csv"), row.names=FALSE, quote=FALSE)

cat("[NOTE] 你也可以用 R 内部直接计算风险分：risk = X %*% Coef。\n")

# -------------------- 假定已从 Python 返回 risk（或直接在 R 里算） --------------------
# 为方便，这里直接在 R 里算训练集/验证集的 risk：
risk_train <- as.numeric(as.matrix(rt_uni[, lasso_genes, drop=FALSE]) %*% geneCoef$Coef)
write.csv(data.frame(risk=risk_train), file = file.path(OUTDIR,"risk_train.csv"), row.names=FALSE)

risk_vali  <- as.numeric(as.matrix(vali_use[, lasso_genes, drop=FALSE]) %*% geneCoef$Coef)
write.csv(data.frame(risk=risk_vali),  file = file.path(OUTDIR,"risk_vali.csv"),  row.names=FALSE)

# -------------------- Kaplan-Meier & timeROC（训练） --------------------
rt1 <- rt_uni
rt1$risk <- risk_train
rt1$group <- factor(ifelse(rt1$risk > median(rt1$risk), "High", "Low"), levels = c("Low","High"))

fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)
pdf(file.path(OUTDIR,"KM_train.pdf"), width=6, height=5)
print(ggsurvplot(fit, data=rt1, pval=TRUE, conf.int=FALSE, pval.size=4, risk.table=FALSE,
                 legend=c(0.88,0.88), legend.labs=c("Low","High"), legend.title="risk",
                 xlab="Years", surv.median.line='hv', palette=c('#0072B5','#BC3C29')))
dev.off()

roc_tr <- with(rt1, timeROC(T = futime, delta = fustat, marker = risk, cause = 1, times = c(1,3,5), ROC=TRUE, iid=TRUE))
df_tr <- data.frame(FPR = as.numeric(roc_tr$FP), TPR = as.numeric(roc_tr$TP), time = rep(as.factor(c(1,3,5)), each = nrow(roc_tr$TP)))
pdf(file.path(OUTDIR,"timeROC_train.pdf"), width=5.5, height=5)
print(ggplot() +
        geom_line(data=df_tr, aes(x=FPR, y=TPR, color=time), linewidth=1) +
        geom_abline(slope=1, intercept=0, color="grey") +
        scale_color_manual(NULL, values=c("#8AD879","#FA9F42","#F3533A"), labels=c("1 year","3 years","5 years")) +
        theme_bw() + labs(x="1 - Specificity (FPR)", y="Sensitivity (TPR)"))
dev.off()

# -------------------- Kaplan-Meier & timeROC（验证，限制时长可选） --------------------
rt2 <- data.frame(surv_va); rt2$risk <- risk_vali
rt2$group <- factor(ifelse(rt2$risk > median(rt2$risk), "High","Low"), levels=c("Low","High"))
rt2$futime <- rt2$futime / 365
# 你原代码限制到 6 年；可按需切片：
# rt2 <- rt2[rt2$futime <= 6, ]

fit2 <- survfit(Surv(futime, fustat) ~ group, data = rt2)
pdf(file.path(OUTDIR,"KM_vali.pdf"), width=6, height=5)
print(ggsurvplot(fit2, data=rt2, pval=TRUE, conf.int=FALSE, pval.size=4, risk.table=FALSE,
                 legend=c(0.88,0.88), legend.labs=c("Low","High"), legend.title="risk",
                 xlab="Years", surv.median.line='hv', palette=c('#0072B5','#BC3C29')))
dev.off()

roc_va <- with(rt2, timeROC(T = futime, delta = fustat, marker = risk, cause = 1, times = c(1,3,5), ROC=TRUE, iid=TRUE))
df_va <- data.frame(FPR = as.numeric(roc_va$FP), TPR = as.numeric(roc_va$TP), time = rep(as.factor(c(1,3,5)), each = nrow(roc_va$TP)))
pdf(file.path(OUTDIR,"timeROC_vali.pdf"), width=5.5, height=5)
print(ggplot() +
        geom_line(data=df_va, aes(x=FPR, y=TPR, color=time), linewidth=1) +
        geom_abline(slope=1, intercept=0, color="grey") +
        scale_color_manual(NULL, values=c("#8AD879","#FA9F42","#F3533A"), labels=c("1 year","3 years","5 years")) +
        theme_bw() + labs(x="1 - Specificity (FPR)", y="Sensitivity (TPR)"))
dev.off()

# -------------------- 临床协变量（需 GDC phenotype）单/多因素 + 森林图 --------------------
# 可选：只在提供临床表时运行
CLIN_PHENO_GZ <- opt[["gdc_pheno_tsv_gz"]] %||% "TCGA-KIRC.GDC_phenotype.tsv.gz"
if(file.exists(CLIN_PHENO_GZ)){
  cl  <- fread(CLIN_PHENO_GZ, data.table = FALSE)
  phe <- cl[, c("submitter_id.samples","tumor_stage.diagnoses","pathologic_M","pathologic_N","pathologic_T","age_at_initial_pathologic_diagnosis","gender.demographic")]
  colnames(phe) <- c("ID","stage","M","N","T","age","gender")
  # 清洗
  phe$age <- suppressWarnings(as.numeric(phe$age))
  # stage
  phe <- phe[phe$stage != "not reported",]
  phe$stage <- toupper(gsub("[a-c]|stage|\\s","", phe$stage))
  phe$stage <- as.numeric(gsub("IV","4", gsub("III","3", gsub("II","2", gsub("I","1", phe$stage)))))
  # T/N/M
  phe$T <- as.numeric(gsub("T","", gsub("[a-c]","", phe$T)));  phe <- phe[!is.na(phe$T),]
  phe$N <- gsub("NX","", phe$N); phe$N <- as.numeric(gsub("N","", phe$N)); phe <- phe[!is.na(phe$N),]
  phe$M <- gsub("MX","", phe$M); phe$M <- as.numeric(gsub("M","", gsub("[a-b]","", phe$M))); phe <- phe[!is.na(phe$M),]
  # 性别
  phe$gender <- ifelse(phe$gender == "male", 1, 0)
  phe <- phe[!duplicated(phe$ID),]; rownames(phe) <- phe$ID; phe$ID <- NULL
  save(phe, file = file.path(OUTDIR,"clinical_clean.Rdata"))
  
  # 对齐训练样本
  uniTab <- read.table(file.path(OUTDIR,"uniSigExp_gene.txt"), header=TRUE, sep="\t", check.names=FALSE, row.names=1)
  risk  <- risk_train; names(risk) <- rownames(uniTab)
  rt1c  <- cbind(uniTab[, c("futime","fustat")], risk = risk)
  ss    <- intersect(rownames(phe), rownames(rt1c))
  rt1c  <- rt1c[ss,]; phe <- phe[ss,]
  rt1c  <- cbind(rt1c, phe)
  
  # 单因素
  uni_vars <- c("age","T","N","M","stage","gender","risk")
  uni_res <- lapply(uni_vars, function(v){
    sm <- summary(coxph(Surv(futime, fustat) ~ rt1c[, v], data = rt1c))
    data.frame(Names=v,
               HR = round(sm$conf.int[,"exp(coef)"],3),
               lower.95 = round(sm$conf.int[,"lower .95"],3),
               upper.95 = round(sm$conf.int[,"upper .95"],3),
               p.value = round(sm$coefficients[,"Pr(>|z|)"],3))
  })
  data.sig <- do.call(rbind, uni_res)
  rownames(data.sig) <- c("Age","T_stage","N_stage","M_stage","Stage","Gender","risk_score")
  data.sig$Names <- rownames(data.sig)
  data.sig$p.value[data.sig$p.value == 0] <- "<0.001"
  write.table(data.sig, file = file.path(OUTDIR,"clinical.sig_unicox.txt"), sep="\t", quote=FALSE)
  
  # 多因素（根据你的原式样）
  muti <- summary(coxph(Surv(futime, fustat) ~ age + T + N + M + stage + risk, data = rt1c))
  data.muti <- data.frame(Names = rownames(muti$conf.int),
                          HR = round(muti$conf.int[,"exp(coef)"],3),
                          lower.95 = round(muti$conf.int[,"lower .95"],3),
                          upper.95 = round(muti$conf.int[,"upper .95"],3),
                          p.value = round(muti$coefficients[,"Pr(>|z|)"],3))
  rownames(data.muti) <- c("Age","T_stage","N_stage","M_stage","Stage","risk_score")
  data.muti$Names <- rownames(data.muti)
  data.muti$p.value[data.muti$p.value == 0] <- "<0.001"
  write.table(data.muti, file = file.path(OUTDIR,"clinical.muti.txt"), sep="\t", quote=FALSE)
  
  # 森林图（单/多）
  draw_forest_from_table(transform(data.sig, gene=Names), file.path(OUTDIR,"forest_clinical_uni.pdf"))
  draw_forest_from_table(transform(data.muti, gene=Names), file.path(OUTDIR,"forest_clinical_multi.pdf"))
}

# -------------------- 免疫浸润（MCPcounter/ESTIMATE） --------------------
# 复用训练表达（log2+1）作为 IOBR 的表达输入（行=gene，列=sample）
immune_eset <- expr_tr
# MCP
immune_mcp <- IOBR::deconvo_mcpcounter(eset = immune_eset)
rownames(immune_mcp) <- immune_mcp$ID; immune_mcp$ID <- NULL
immune_mcp <- immune_mcp[rownames(rt1), , drop=FALSE]
immune_mcp$group <- rt1$group

# 箱线图
rt_long <- tidyr::pivot_longer(immune_mcp, cols = -group, names_to = "cell", values_to = "score")
pdf(file.path(OUTDIR,"MCPcounter_box.pdf"), width=9, height=5)
print(ggboxplot(rt_long, x="cell", y="score", color="black", notch=TRUE, fill="group",
                xlab="risk_group", ylab="score", palette=c("skyblue","orange")) +
        stat_compare_means(aes(group=group), label="p.signif", method="wilcox.test", hide.ns=TRUE, size=4) +
        theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1)))
dev.off()

# ESTIMATE
immune_est <- IOBR::deconvo_estimate(eset = immune_eset)
rownames(immune_est) <- immune_est$ID; immune_est$ID <- NULL
immune_est <- immune_est[rownames(rt1), , drop=FALSE]
immune_est$group <- rt1$group
rt_long2 <- tidyr::pivot_longer(immune_est, cols = -group, names_to = "cell", values_to = "score")
pdf(file.path(OUTDIR,"ESTIMATE_box.pdf"), width=8, height=5)
print(ggboxplot(rt_long2, x="cell", y="score", color="black", notch=TRUE, fill="group",
                xlab="risk_group", ylab="score", palette=c("skyblue","orange")) +
        stat_compare_means(aes(group=group), label="p.signif", method="wilcox.test", hide.ns=TRUE, size=4) +
        theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1)))
dev.off()

# -------------------- ssGSEA pathway --------------------
pathway <- IOBR::calculate_sig_score(eset = immune_eset, method = "ssgsea")
pathway <- as.data.frame(pathway)
rownames(pathway) <- pathway$ID; pathway$ID <- NULL
pathway <- pathway[rownames(rt1), , drop=FALSE]
pathway$group <- rt1$group

# 仅 TIP_* 通路（如有）
tip_cols <- grep("^TIP_", colnames(pathway), value = TRUE)
if(length(tip_cols)){
  tip <- pathway[, c(tip_cols, "group"), drop=FALSE]
  rt_tip <- tidyr::pivot_longer(tip, cols = -group, names_to = "Pathway", values_to = "score")
  pdf(file.path(OUTDIR,"TIP_box.pdf"), width=10, height=5)
  print(ggboxplot(rt_tip, x="Pathway", y="score", color="black", notch=TRUE, fill="group",
                  xlab="group", ylab="score", palette=c("#EFC000","#0073C2")) +
          stat_compare_means(aes(group=group), label="p.signif", method="wilcox.test", hide.ns=TRUE, size=4) +
          theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1)))
  dev.off()
}

# -------------------- checkpoint 基因表达差异 --------------------
CHECK_FILE <- opt[["checkpoint_file"]] %||% ""
if(nchar(CHECK_FILE) && file.exists(CHECK_FILE)){
  ck <- read.table(CHECK_FILE, header = FALSE)
  genes_ck <- ck[[1]]
  expr_ck <- immune_eset[rownames(immune_eset) %in% genes_ck, , drop=FALSE]
  expr_ck <- as.data.frame(t(expr_ck))
  expr_ck <- expr_ck[rownames(rt1), , drop=FALSE]
  expr_ck$group <- rt1$group
  rt_ck <- tidyr::pivot_longer(expr_ck, cols = -group, names_to = "gene", values_to = "expression")
  pdf(file.path(OUTDIR,"checkpoint_box.pdf"), width=10, height=5)
  print(ggboxplot(rt_ck, x="gene", y="expression", color="black", notch=TRUE, fill="group",
                  xlab="group", ylab="expression", palette=c("#EFC000","#0073C2")) +
          stat_compare_means(aes(group=group), label="p.signif", method="wilcox.test", hide.ns=TRUE, size=4) +
          theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1)))
  dev.off()
}

cat("[DONE] All outputs at:", OUTDIR, "\n")
