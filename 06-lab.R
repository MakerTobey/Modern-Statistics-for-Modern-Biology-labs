pkgs_needed = c("dplyr", "ggplot2", "DESeq2", "pasilla", "genefilter",
                "pheatmap", "readr", "tibble", "apeglm",
                "TENxPBMCData", "MASS")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(setdiff(pkgs_needed, installed.packages()))

library("tidyverse")
library("ggplot2")
library("DESeq2")
library("pasilla")
library("genefilter")
library("pheatmap")
library("MASS")
library("TENxPBMCData")

fn = system.file("extdata", "pasilla_gene_counts.tsv",
                 package = "pasilla", mustWork = TRUE)
counts = as.matrix(read.csv(fn, sep = "\t", row.names = "gene_id"))

dim(counts)
counts[ 2000+(0:3), ]

annotationFile = system.file("extdata", "pasilla_sample_annotation.csv",
                             package = "pasilla", mustWork = TRUE)
pasillaSampleAnno = readr::read_csv(annotationFile)
pasillaSampleAnno

pasillaSampleAnno = mutate(
  pasillaSampleAnno,
  condition = factor(condition, levels = c("untreated", "treated")),
  type      = factor(type, levels = c("single-read", "paired-end")))

mt = match(colnames(counts), sub("fb$", "", pasillaSampleAnno$file))
pasilla = DESeqDataSetFromMatrix(
  countData = counts,
  colData   = pasillaSampleAnno[mt, ],
  design    = ~ condition)
class(pasilla)

is(pasilla, "SummarizedExperiment")

library("ggplot2")
library("matrixStats")
sf = estimateSizeFactorsForMatrix(counts)
ncounts  = counts / matrix(sf, 
                           byrow = TRUE, ncol = ncol(counts), nrow = nrow(counts))
uncounts = ncounts[, grep("^untreated", colnames(ncounts)), drop = FALSE]

ggplot(data.frame( mean = rowMeans(uncounts), var  = rowVars( uncounts)),
       aes(x = log(mean), y = log(var))) +
  geom_hex() + 
  geom_abline(slope = c(1, 2), intercept = c(0, -4),
              color = c("forestgreen", "red")) +
  coord_fixed() + theme(legend.position = "none") 

ggplot(tibble(
  `size factor` = estimateSizeFactorsForMatrix(counts),
  `sum` = colSums(counts)), aes(x = `size factor`, y = `sum`)) +
  geom_point()

pasilla = DESeq(pasilla)

res = results(pasilla)
res[order(res$padj), ] %>% head

ggplot(as(res, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)

pv_hist = hist(res$pvalue, breaks = seq(0, 1, length = 100), plot = FALSE)
pvals = na.omit(res$pvalue)

set.seed(0xdada2)
y = cbind(rnorm(10000, 0, 1),
          rnorm(10000, 0, 1),
          rnorm(10000, 0, 1),
          rnorm(10000, 0, 1))
library(genefilter)
pvalue = rowttests(y, factor(c("C","C","T","T")))$p.value
ggplot(tibble(pvalue), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)

plotMA(pasilla, ylim = c( -2, 2))

lfc_shrink_res = lfcShrink(pasilla, coef="condition_treated_vs_untreated", type="apeglm")
plotMA(lfc_shrink_res, ylim = c( -2, 2))

pas_rlog = rlogTransformation(pasilla)
plotPCA(pas_rlog, intgroup=c("condition", "type")) + coord_fixed()

library("pheatmap")
select = order(rowMeans(assay(pas_rlog)), decreasing = TRUE)[1:30]
pheatmap(
  assay(pas_rlog)[select, ],
  scale = "row",
  annotation_col = as.data.frame(colData(pas_rlog)[, c("condition", "type")]))

#BiocManager::install("ComplexHeatmap")
#library("ComplexHeatmap")

pasillaTwoFactor = pasilla
design(pasillaTwoFactor) = formula(~ type + condition)
pasillaTwoFactor = DESeq(pasillaTwoFactor)

res2 = results(pasillaTwoFactor)
head(res2, n = 3)

resType = results(pasillaTwoFactor, 
                  contrast = c("type", "single-read", "paired-end"))
head(resType, n = 3)

trsf = function(x) ifelse(is.na(x), 0, (-log10(x)) ^ (1/6))
ggplot(tibble(pOne = res$pvalue,
              pTwo = res2$pvalue),
       aes(x = trsf(pOne), y = trsf(pTwo))) +
  geom_hex(bins = 75) + coord_fixed() +
  xlab("Single factor analysis (condition)") +
  ylab("Two factor analysis (type + condition)") +
  geom_abline(col = "orange")

compareRes = table(
  `simple analysis` = res$padj < 0.1,
  `two factor` = res2$padj < 0.1 )
addmargins( compareRes )

# I don't have this file!
sce = TENxPBMCData::TENxPBMCData("frozen_pbmc_donor_a")
sce

sce_counts = as.matrix(counts(sce))
sce_vars  = rowVars(sce_counts)
sce_means = rowMeans(sce_counts)
sce_df = tibble(mean = sce_means, var=sce_vars)

ggplot(sce_df, aes(x = mean, y = var)) + 
  geom_point(alpha = 0.6, size = 0.3) +
  xlab("Mean") + ylab("Variance")

rlm_fit = rlm(var-mean ~ 0  + I(mean^2), data=sce_df_filt , maxit=100)

alpha_dispersion = coef(rlm_fit)
alpha_dispersion

sce_df_filt  = dplyr::mutate(sce_df_filt , fitted_var = mean + alpha_dispersion*mean^2)
ggplot(sce_df_filt, aes(x = mean, y = var)) + 
  geom_point(alpha = 0.6, size = 0.3) + 
  geom_line(aes(y = fitted_var), col = "red")

zero_prop = rowMeans( sce_counts_filt == 0 )

ggplot(tibble(zero_prop = zero_prop), aes(x = zero_prop)) +
  geom_histogram(binwidth = 0.05, fill = "Royalblue", boundary = 0) + 
  xlab("Proportion of zeros")

sce_df_filt  = dplyr::mutate(sce_df_filt, 
                             zero_prop = zero_prop,
                             theoretical_zero_prop = dnbinom(0, mu = mean, size = 1/alpha_dispersion))

ggplot(sce_df_filt, aes(x = mean, y = zero_prop))  +
  geom_point(alpha = 0.6, size = 0.3) +
  geom_line(aes(y=theoretical_zero_prop), col="red") +
  scale_x_log10() + ylab("proportion of zeros")
