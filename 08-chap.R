## ----initialize, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
source("../chapter-setup.R"); chaptersetup("/Users/Susan/Courses/CUBook-html/CUBook/Chap7-RNAseq/rnaseqanalysis.Rnw", "7")
knitr::opts_chunk$set(dev = 'png', dpi = 100, fig.margin = TRUE, fig.show = 'hold', fig.keep = 'none')

## ---- xkcd-1725-linear-regression-2x, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/xkcd-1725-linear_regression_2x.png'))


## ----loadpas, results="hide", error=FALSE--------------------------------
fn = system.file("extdata", "pasilla_gene_counts.tsv",
                  package = "pasilla", mustWork = TRUE)
counts = as.matrix(read.csv(fn, sep = "\t", row.names = "gene_id"))


## ----counts--------------------------------------------------------------
dim(counts)
counts[ 2000+(0:3), ]


## ----checkClaimMadeAbove, echo = FALSE-----------------------------------
conditionNames = (sub("[[:digit:]]$", "", colnames(counts))) #$
stopifnot(length(unique(conditionNames)) == 2,
  sum(conditionNames=="untreated") == 4,
  sum(conditionNames=="treated")   == 3)


## ----rnaseq-normalization, fig.keep = 'high', fig.cap = "Size factor estimation. The points correspond to hypothetical genes whose counts in two samples are indicated by their $x$- and $y$-coordinates. The lines indicate two different ways of size factor estimation explained in the text.", fig.width = 2.5, fig.height = 2.5, echo = FALSE----
szfcDemo = data.frame(
  x = c(2, 4, 6, 6,  8) * 10,
  y = c(3, 6, 2, 9, 12) * 10,
  name = LETTERS[1:5],
  check.names = FALSE)
slopes =  c(
  blue = with(szfcDemo, sum(y) / sum(x)),
  red = szfcDemo[, c("x", "y")] %>% as.matrix %>%
    (DESeq2::estimateSizeFactorsForMatrix) %>% (function(x) x[2]/x[1]) %>% as.vector)
ggplot(szfcDemo, aes(x = x, y = y, label = name)) + geom_point() +
  coord_fixed() + xlim(c(0, 128)) + ylim(c(0, 128)) + xlab("sample 1") + ylab("sample 2") +
  geom_text(hjust= 0.5, vjust = -0.6) +
  geom_abline(slope = slopes[1], col = names(slopes)[1]) +
  geom_abline(slope = slopes[2], col = names(slopes)[2])

## ----rnaseq-sfvssum, fig.keep = 'high', fig.cap = "Size factors versus sums for the pasilla data.", fig.width = 2.5, fig.height = 2.5----
library("tibble")
ggplot(tibble(
  `size factor` = estimateSizeFactorsForMatrix(counts),
  `sum` = colSums(counts)), aes(x = `size factor`, y = `sum`)) +
  geom_point()

## ----rnaseq-varmean, fig.keep = 'high', fig.cap = "Variance versus mean for the (size factor adjusted) `counts` data. The axes are logarithmic. Also shown are lines through the origin with slopes 1 (green) and 2 (red). ", fig.width = 2.5, fig.height = 3, warning = FALSE----
library("ggplot2")
library("matrixStats")
sf = estimateSizeFactorsForMatrix(counts)
ncounts  = counts / matrix(sf,
   byrow = TRUE, ncol = ncol(counts), nrow = nrow(counts))
uncounts = ncounts[, grep("^untreated", colnames(ncounts)),
                     drop = FALSE]
ggplot(tibble(
        mean = rowMeans(uncounts),
        var  = rowVars( uncounts)),
     aes(x = log(mean), y = log(var))) +
  geom_hex() + coord_fixed() + theme(legend.position = "none") +
  geom_abline(slope = 1:2, color = c("forestgreen", "red"))


## ----annotationFile------------------------------------------------------
annotationFile = system.file("extdata",
  "pasilla_sample_annotation.csv",
  package = "pasilla", mustWork = TRUE)
pasillaSampleAnno = readr::read_csv(annotationFile)
pasillaSampleAnno

## ----factors-------------------------------------------------------------
library("dplyr")
pasillaSampleAnno = mutate(pasillaSampleAnno,
condition = factor(condition, levels = c("untreated", "treated")),
type = factor(sub("-.*", "", type), levels = c("single", "paired")))

## ----checkfacs,  echo=FALSE----------------------------------------------
stopifnot(
  !any(is.na(pasillaSampleAnno$condition)),
  !any(is.na(pasillaSampleAnno$type)),
  sum(pasillaSampleAnno$type == "single") == 3,
  sum(pasillaSampleAnno$type == "paired") == 4)

## ----condvstype----------------------------------------------------------
with(pasillaSampleAnno,
       table(condition, type))

## ----DESeq2, message = FALSE, warning = FALSE----------------------------
mt = match(colnames(counts), sub("fb$", "", pasillaSampleAnno$file))
stopifnot(!any(is.na(mt)))

library("DESeq2")
pasilla = DESeqDataSetFromMatrix(
  countData = counts,
  colData   = pasillaSampleAnno[mt, ],
  design    = ~ condition)
class(pasilla)
is(pasilla, "SummarizedExperiment")

## ----ispasillaSummarizedExperiment, echo = FALSE-------------------------
stopifnot(is(pasilla, "SummarizedExperiment"))

## ----deseq---------------------------------------------------------------
pasilla = DESeq(pasilla)

## ----theresults----------------------------------------------------------
res = results(pasilla)
res[order(res$padj), ] %>% head

## ----rnaseq-hist1, fig.keep = 'high', fig.cap = "Histogram of p-values of a differential expression analysis.", fig.width = 4.5, fig.height = 4.5, warning = FALSE----
ggplot(as(res, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)

## ----hist2, echo=FALSE---------------------------------------------------
thehist = hist(res$pvalue, breaks = 100, plot=FALSE)
thehist$bgl = median(thehist$counts)

## ----rnaseq-MA, fig.keep = 'high', fig.cap = "\\href{https://en.wikipedia.org/wiki/MA_plot}{MA plot}: fold change versus mean of size-factor normalized counts. Logarithmic scaling is used for both axes. By default, points are colored red if the adjusted p-value is less than 0.1. Points which fall out of the $y$-axis range are plotted as triangles.", fig.width = 4.5, fig.height = 4.5----
plotMA(pasilla, ylim = c( -2, 2))

## ----rnaseq-PCA, fig.keep = 'high', fig.cap = "PCA plot. The (ref:rnaseq-PCA-1) samples are shown in the 2D plane spanned by their first two principal components.", fig.width = 4, fig.height = 4----
pas_rlog = rlogTransformation(pasilla)
plotPCA(pas_rlog, intgroup=c("condition", "type")) + coord_fixed()

## ----chap7-r-figHeatmap-1, fig.keep = 'high', fig.cap = "Heatmap of regularized log transformed data of the top (ref:chap7-r-figHeatmap-1-1) genes.", fig.width = 4, fig.height = 6----
library("pheatmap")
select = order(rowMeans(assay(pas_rlog)), decreasing = TRUE)[1:30]
pheatmap( assay(pas_rlog)[select, ],
     scale = "row",
     annotation_col = as.data.frame(
        colData(pas_rlog)[, c("condition", "type")] ))

## ----writecsv------------------------------------------------------------
write.csv(as.data.frame(res), file = "treated_vs_untreated.csv")






## ----rnaseq-mestimator, fig.keep = 'high', fig.cap = "Graph of $\\rho_s(\\varepsilon)$, for a choice of $s=2$.", fig.width = 3, fig.height = 2.5----
rho = function(x, s)
  ifelse(abs(x) < s, x^2 / 2,  s * abs(x) - s^2 / 2)

df = tibble(
  x        = seq(-7, 7, length.out = 100),
  parabola = x ^ 2 / 2,
  Huber    = rho(x, s = 2))

ggplot(reshape2::melt(df, id.vars = "x"),
  aes(x = x, y = value, col = variable)) + geom_line()


## ----replaceDesign-------------------------------------------------------
pasillaTwoFactor = pasilla
design(pasillaTwoFactor) = formula(~ type + condition)
pasillaTwoFactor = DESeq(pasillaTwoFactor)

## ----multiResults--------------------------------------------------------
res2 = results(pasillaTwoFactor)
head(res2, n = 3)

## ----multiTypeResults----------------------------------------------------
resType = results(pasillaTwoFactor,
  contrast = c("type", "single", "paired"))
head(resType, n = 3)

## ----rnaseq-scpres1res2, fig.keep = 'high', fig.cap = "Comparison of p-values from the models with a single factor (condition) and with two factors (type + condition). The axes correspond to $(-\\log_{10}p)^{\\frac{1}{6}}$, an arbitrarily chosen monotonically decreasing transformation that compresses the dynamic range of the p-values for the purpose of visualization. We can see a trend for the joint distribution to lie above the bisector, indicating that the small p-values in the two-factor analysis are generally smaller than those in the one-factor analysis.", fig.width = 3.8, fig.height = 3.8, warning = FALSE----
trsf = function(x) ifelse(is.na(x), 0, (-log10(x)) ^ (1/6))
ggplot(tibble(pOne = res$pvalue,
              pTwo = res2$pvalue),
    aes(x = trsf(pOne), y = trsf(pTwo))) +
    geom_hex(bins = 75) + coord_fixed() +
    xlab("Single factor analysis (condition)") +
    ylab("Two factor analysis (type + condition)") +
    geom_abline(col = "orange")

## ----compareRes----------------------------------------------------------
compareRes = table(
   `simple analysis` = res$padj < 0.1,
   `two factor` = res2$padj < 0.1 )
addmargins( compareRes )

## ----checkClaim, echo = FALSE, results = "hide", error = TRUE------------
stopifnot(compareRes[1, 2] > compareRes[2, 1])

## ----shrink0, echo = FALSE-----------------------------------------------
.suitableDESeq2 = (packageVersion("DESeq2") >= "1.17.8") #$

## ----selgenes, echo = FALSE, eval = FALSE--------------------------------
## ## This chunk reproduces Mike's way of selecting the two genes: they should have similar
## ## 'intercepts', large unshrunken FC and very different Wald statistic (i.e. small/large dispersion)
## with(res1,
##   plot(baseMean, log2FoldChange, log = "x", ylim = c(0, 3), xlim = c(10, 1e5),
##        col = ifelse(padj < 0.1, "red", "black"), cex = log(abs(stat))))
## rownames(res1)[with(res1, identify(baseMean, log2FoldChange))]

## ----defvsp--------------------------------------------------------------
vsp = varianceStabilizingTransformation(pasilla)

## ----rnaseq-plotvst, fig.keep = 'high', fig.cap = "Graph of variance-stabilizing transformation for the data of one of the samples, and for comparison also of the $\\log_2$ transformation. The variance-stabilizing transformation has finite values and finite slope even for counts close to zero, whereas the slope of $\\log_2$ becomes very steep for small counts and is undefined for counts of zero. For large counts, the two transformation are essentially the same.", fig.width = 3.5, fig.height = 3, warning = FALSE----
j = 1
ggplot(tibble(
         x    = assay(pasilla)[, j],
         VST  = assay(vsp)[, j],
         log2 = log2(assay(pasilla)[, j])) %>%
             reshape2::melt(id.vars = "x"),
       aes(x = x, y = value, col = variable)) +
  geom_line() + xlim(c(0, 600)) + ylim(c(0, 9)) +
  xlab("counts") + ylab("transformed")

## ----rnaseq-meansd, fig.keep = 'high', fig.cap = "Per-gene standard deviation (sd, taken across samples) against the rank of the mean, for the shifted logarithm $\\log_2(n+1)$, the variance-stabilizing transformation (vst) and the rlog. Note that for the leftmost $\\approx$ 2,500 genes, the counts are all zero, and hence their standard deviation is zero. The mean-sd dependence becomes more interesting for genes with non-zero counts. Note also the high value of the standard deviation for genes that are weakly detected (but not with all zero counts) when the shifted logarithm is used, and compare to the relatively flat shape of the mean-sd relationship for the variance-stabilizing transformation. ", fig.margin = FALSE, fig.width = 10, fig.height = 4, warning = FALSE----
library("vsn")
rlp = rlogTransformation(pasilla)

msd = function(x)
  meanSdPlot(x, plot = FALSE)$gg + ylim(c(0, 1)) +
     theme(legend.position = "none")

gridExtra::grid.arrange(
  msd(log2(counts(pasilla, normalized = TRUE) + 1)) +
    ylab("sd(log2)"),
  msd(assay(vsp)) + ylab("sd(vst)"),
  msd(assay(rlp)) + ylab("sd(rlog)"),
  ncol = 3
)

## ----rnaseq-lfcThresh, fig.keep = 'high', fig.cap = "MA-plots of tests of $\\log_2$ fold change with respect to a threshold value. From top to bottom, the tests are for \\texttt{altHypothesis = \"greaterAbs\"}, \\texttt{\"lessAbs\"}, \\texttt{\"greater\"}, and \\texttt{\"less\"}.", fig.width = 2.7, fig.height = 9----
par(mfrow = c(4, 1), mar = c(2, 2, 1, 1))
myMA = function(h, v, theta = 0.5) {
  plotMA(pasilla, lfcThreshold = theta, altHypothesis = h,
         ylim = c(-2.5, 2.5))
  abline(h = v * theta, col = "dodgerblue", lwd = 2)
}
myMA("greaterAbs", c(-1, 1))
myMA("lessAbs",    c(-1, 1))
myMA("greater",          1)
myMA("less",         -1   )

## ----ui.R, eval = FALSE--------------------------------------------------
## library("shiny")
## shinyUI(fluidPage(
##   titlePanel("Breakdown"),
##   sidebarLayout(
##     sidebarPanel(     # select oulier shift
##       sliderInput("shift", "Outlier:", min = 0, max = 100, value = 0),
##       radioButtons("method", "Method:",
##                    c("Non-robust least squares" = "lm",
##                      "M-estimation" = "rlm"))
##     ),
##     mainPanel(       # show fit
##       plotOutput("regPlot")
##     )
##   )
## ))

## ----server.R, eval = FALSE----------------------------------------------
## library("shiny")
## library("ggplot2")
## library("MASS")
## shinyServer(function(input, output) {
##   output$regPlot = renderPlot({
##     whpt = 15
##     mtcars_new = mtcars
##     mtcars_new$mpg[whpt] = mtcars_new$mpg[whpt] + input$shift
##     reg = switch(input$method,
##       lm = lm(mpg ~ disp, data = mtcars_new),
##       rlm = rlm(mpg ~ disp, data = mtcars_new),
##       stop("Unimplemented method:", input$method)
##     )
##     ggplot(mtcars_new, aes(x = disp, y = mpg)) + geom_point() +
##       geom_abline(intercept = reg$coefficients["(Intercept)"],
##                   slope = reg$coefficients["disp"], col = "blue")
##   })
## })

