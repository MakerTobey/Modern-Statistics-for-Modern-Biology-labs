## ----initialize, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
source("../chapter-setup.R"); chaptersetup("/Users/Susan/Courses/CUBook-html/CUBook/Chap10-AboutDesign/ExperimentalDesign.Rnw", "10")
knitr::opts_chunk$set(dev = 'png', dpi = 100, fig.margin = TRUE, fig.show = 'hold', fig.keep = 'none')

## ---- RAFisherSmoking, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/RAFisherSmoking.png'))

## ---- dailies-icon, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/dailies_icon.png'))




## ---- cointosser3, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "A carefully contructed coin tossing machine can be made to provide deterministic coin flips."----
knitr::include_graphics(c('images/cointosser3.png'))


## ----generatedata, echo=FALSE--------------------------------------------
library("dplyr")
bat6 = tibble(
  state  = factor(c("healthy", "disease")[rep(1:2, 6)]),
  time  = factor(rep(1:3, each = 4)),
  exprst = 0.5 * as.integer(time) + rep(c(0.5, 1), 6) + rnorm(12, 0, 0.3),
  exprs0 = rep(c(1.5, 2), 6) + rnorm(12,0,0.1),
  batch  = factor(c("Batch 1", "Batch 2")[rep(c(1, 2), 6)]))

ms0 = group_by(bat6, state) %>% summarize(y = median(exprs0))

bat60 = tibble(
  state = factor(c("healthy", "disease")[rep(c(1, 2), 60)], levels = c("healthy", "disease")),
  exprs = rep(c(1.5, 2), 60) + rnorm(120, 0, 0.3))
## save(bat60, bat6, ms0, file = "../data/designI.rda")

## ----chap10-r-confounding-1, fig.keep = 'high', fig.cap = "Comparison of a (hypothetical) biomarker between samples from disease and healthy states. If we are only given the information shown in the left panel, we might conclude that this biomarker performs well in detecting the disease. If, in addition, we are told that the data were acquired in two separate batches (e.g., different labs, different machines, different time points) as indicated in the panel on the right hand side, the conclusion will be different.", fig.width = 4.25, fig.height = 3, echo = FALSE----
library("ggplot2")
library("gridExtra")
library("ggbeeswarm")
## load("../data/designI.rda")
p0 = ggplot(bat6, aes(x = state, y = exprs0)) +
       geom_boxplot(alpha = 0.5, col="blue") + geom_beeswarm(size = 2, cex = 6) + # geom_point(size = 2) +
       ylab("biomarker level")
grid.arrange(p0, p0 + geom_beeswarm(aes(col = batch), size = 2, cex = 6),  # geom_point(aes(col = batch), size = 2),
  ncol = 2, widths = c(1.3, 2))

## ---- Avicenna, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Confounding is the reason that one of the seven rules of experimental design listed by the Persian physician-scientist \\href{https://en.wikipedia.org/wiki/Avicenna}{Abu \'Ali al-Husayn ibn Sina (Avicenna)} around AD 1020 was \"to study one possible cause of a disease at a time\"  [@Stigler:sevenpillars], ."----
knitr::include_graphics(c('images/Avicenna.png'))

## ----Design-effectsize, fig.keep = 'high', fig.cap = "The red arrow shows the effect size, as measured by the difference between the centers of the two groups. Here we locate the centers by the medians; sometimes the mean is used.", echo = FALSE, fig.width = 2, fig.height = 3----
p0 + geom_segment(data = ms0, aes(y = ms0$y[1], yend = ms0$y[2]),
    x = 1.5, xend = 1.5, col = "red", arrow = arrow(length = unit(0.5, "cm"),
    ends = "both", type = "closed"))

## ----Design-comparesamplesize, fig.keep = 'high', fig.cap = "On the left, the boxplot was created with samples of size 6. On the right the sample sizes are 60. The measurements have the same underlying error distribution in both cases.", echo = FALSE, fig.width = 3.7, fig.height = 3----
p = ggplot(bat6, aes(x = state, y = exprst)) + geom_boxplot(alpha = 0.5, col = "blue") +
    ylim(c(0.5, 3)) + ylab("biomarker level")
p1 = p + geom_beeswarm(size = 2, cex = 6) # geom_point(size = 2)
p2 = p + geom_beeswarm(aes(col = time), size = 2, cex = 6) # geom_point(aes(col = time), size = 2)

mN = summarise(group_by(bat60, state), med = median(exprs))
pN = ggplot(bat60, aes(x = state, y = exprs)) + geom_boxplot(alpha = 0.5, col="blue") +
  ylab("biomarker level") + ylim(c(0.5,3)) + geom_beeswarm(size = 2, cex = 2)  +
  geom_segment(data = mN, aes(y = med[1],yend=med[2]), x = 1.5, xend = 1.5,
               col = "red", arrow = arrow(length = unit(0.5, "cm"), ends = "both", type = "closed"))
grid.arrange(p1, pN, ncol = 2, widths = c(1.6, 2.5))

## ---- Design-balance, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The example in this section uses the pharmacist\'s balance weighing analogy introduced by Yates and developed by  @Hotelling1944 and  @Mood1946."----
knitr::include_graphics(c('images/balancechem.png'))

## ----HotellingsExpt------------------------------------------------------
theta = round((2 * sample(8, 8) + rnorm(8)), 1)
theta

## ----SimpleWeighing------------------------------------------------------
X = theta + rnorm(length(theta), 0, 0.1)
X
errors1 = X - theta
errors1
sum(errors1^2)

## ----HotellingsMethod----------------------------------------------------
library("survey")
h8 = hadamard(6)
coef8 = 2*h8 - 1
coef8

## ------------------------------------------------------------------------
Y = theta  %*% coef8 + rnorm(length(theta), 0, 0.1)

## ----coef8---------------------------------------------------------------
coef8 %*% t(coef8)
theta %*% coef8 %*% t(coef8) / ncol(coef8)

## ----thetahat------------------------------------------------------------
thetahat = Y %*% t(coef8) / ncol(coef8)

## ----Hoterrors-----------------------------------------------------------
errors2 = as.vector(thetahat) - theta
errors2
sum(errors2^2)

## ----bootstrapHotelling--------------------------------------------------
B  = 10000
tc = t(coef8) / ncol(coef8)
sse = replicate(B, {
  theta = round((2 * sample(8, 8)) + rnorm(8), 1)
  X = theta + rnorm(length(theta), 0, 0.1)
  err1 = sum((X - theta)^2)
  Y = coef8 %*% theta + rnorm(length(theta), 0, 0.1)
  thetahat = tc %*% Y
  err2 = sum((thetahat - theta)^2)
  c(err1, err2)
})
rowMeans(sse)

## ----Design-logsseratios, fig.keep = 'high', fig.cap = "Logarithm (base 2) of the ratios of sum of squared error for the two methods. The vertical orange line corresponds to 8.", fig.width = 3.2, fig.height = 2----
ggplot(tibble(lr = log2(sse[1, ] / sse[2, ])), aes(x = lr)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = log2(8), col = "orange") +
  xlab("log2 ratio of SSE, Method 1 vs 2")

## ----Design-blockbox, fig.keep = 'high', fig.cap = "On the left, two samples each of size 6 are being compared. On the right, the same data are shown, but colored by the time of data collection. We note a tendency of the data to fall into blocks according to these times. Because of this, comparison between the groups is diluted. This effect can be mitigated by comparing within times, i.\\,e., by blocking into three groups. Paired analysis, such as demonstrated in Questions \\@ref(ques:Design-ques-paired)-\\@ref(ques:Design-ques-powerPairedUnpaired), is a special case of blocking.", echo = FALSE, fig.width = 4.25, fig.height = 3----
grid.arrange(p1, p2, ncol = 2, widths = c(1.3, 2))

## ---- design-maizedarwin, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "A paired experiment is the simplest case of blocking."----
knitr::include_graphics(c('images/maizeDarwin.png'))

## ----zeamays-------------------------------------------------------------
n = 15
effect = 0.2
pots   = rnorm(n, 0, 1)
noiseh = rnorm(n, 0, 0.25)
noisea = rnorm(n, 0, 0.25)
hybrid = pots + effect + noiseh
autoz  = pots + noisea

## ----ttestpairedornot----------------------------------------------------
t.test(hybrid, autoz, paired = FALSE)
t.test(hybrid, autoz, paired = TRUE)

## ----bootstrapPower------------------------------------------------------
B     = 1000
alpha = 0.05
what  = c(FALSE, TRUE)
pvs = replicate(B, {
  pots   = rnorm(n, 0, 1)
  noiseh = rnorm(n, 0, 0.25)
  noisea = rnorm(n, 0, 0.25)
  hybrid = pots + effect + noiseh
  autoz  = pots + noisea
  vapply(what,
    function(paired)
      t.test(hybrid, autoz, paired = paired)$p.value,
    double(1)) %>% setNames(paste(what))
})
rowMeans(pvs <= alpha)

## ----chap10-r-pvaluescompare-1, fig.keep = 'high', fig.cap = "Results from the power calculation, comparing the p-value distributions from the ordinary unpaired and the paired $t$-test.", fig.height = 2, fig.width = 3----
library("reshape2")
ggplot(melt(pvs), aes(x = value, fill = Var1)) +
  geom_histogram(binwidth = 0.01, boundary = 0, alpha = 1/3) +
  scale_fill_discrete(name = "Paired")

## ----powerPairedUnpaired-------------------------------------------------
powercomparison = function(effect = 0.2, n = 15, alpha = 0.05,
                sdnoise, sdpots, B = 1000) {
  what = c(FALSE, TRUE)
  pvs = replicate(B, {
    pots   = rnorm(n, 0, sdpots)
    noiseh = rnorm(n, 0, sdnoise)
    noisea = rnorm(n, 0, sdnoise)
    hybrid = pots + effect + noiseh
    autoz  = pots + noisea
    vapply(what,
      function(paired)
        t.test(hybrid, autoz, paired = paired)$p.value,
      double(1)) %>% setNames(paste(what))
  })
  rowMeans(pvs <= alpha)
}

## ------------------------------------------------------------------------
powercomparison(sdpots = 0.5,  sdnoise = 0.25)
powercomparison(sdpots = 0.25, sdnoise = 0.25)
powercomparison(sdpots = 0.1,  sdnoise = 0.25)

## ------------------------------------------------------------------------
powercomparison(sdpots = 0.5, sdnoise = 0.5, n = 100)


## ---- Design-elephant, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The elephant in the room with power calculations is the effect size. Especially in \'omics studies, when we are screening thousands of genes (or other features) for differences, we rarely have a precise idea of what effect size to expect. However, even so, power calculations are useful for order-of-magnitude calculations, or for qualitative comparisons such as shown in this section for paired versus unpaired tests. Source: [Wikimedia CH](https://en.wikipedia.org/wiki/Elephant)."----
knitr::include_graphics(c('images/African_Bush_Elephant.jpg'))

## ----pwr.t.test----------------------------------------------------------
library("pwr")
str(pwr.t.test)

## ------------------------------------------------------------------------
pwr.t.test(n = 15, d = 0.4, sig.level = 0.05, type = "two.sample")
pwr.t.test(n = 15, d = 0.4, sig.level = 0.05, type = "paired")

## ------------------------------------------------------------------------
pwr.t.test(d = 0.4, sig.level = 0.05, type = "two.sample", power=0.8)
pwr.t.test(d = 0.4, sig.level = 0.05, type = "paired", power=0.8)

## ----Design-r-effective-sample-size-sim-1, fig.keep = 'high', fig.cap = "Density estimates for the polling result using the two sampling methods. The correlated method has higher spread. The truth is indicated by the vertical line.", fig.width = 4, fig.height = 2.5----
doPoll = function(n = 100, numPeoplePolled = 12) {
  opinion = sort(rnorm(n))
  i1 = sample(n, numPeoplePolled)
  i2 = sample(seq(3, n, by = 3), numPeoplePolled / 3)
  i2 = c(i2, i2 - 1, i2 - 2)
  c(independent = mean(opinion[i1]), correlated = mean(opinion[i2]))
}
responses = replicate(5000, doPoll())
ggplot(melt(responses), aes(x = value, col = Var1)) + geom_density() +
  geom_vline(xintercept = 0) + xlab("Opinion poll result")

## ---- 1896-Ford-Quadricycle, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Henry Ford\'s (possibly apocryphal) quote: \"If I had asked people what they wanted, they would have said faster horses.\" expresses the view of quality as **fitness for purpose**, versus adherence to specifications. ([Source: Ford](https://corporate.ford.com/history.html))"----
knitr::include_graphics(c('images/1896_Ford_Quadricycle.jpeg'))


## ----dt2-----------------------------------------------------------------
library("Hiiragi2013")
library("magrittr")
data("x")
xwdf = tibble(
  probe  = c("1420085_at", "1418863_at", "1425463_at", "1416967_at"),
  symbol = c(      "Fgf4",      "Gata4",      "Gata6",       "Sox2"))
xwdf %<>% bind_cols(as_tibble(Biobase::exprs(x)[xwdf$probe, ]))
dim(xwdf)
xwdf[, 1:5]

## ----melt----------------------------------------------------------------
xldf = melt(xwdf, id.vars = c("probe", "symbol"),
                  variable.name = "sample")
dim(xldf)
head(xldf)

## ---- leakypipeline, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Sequential data analyses workflows can be leaky. If insufficient information is passed from one stage to the next, the procedure can end up being suboptimal and losing power."----
knitr::include_graphics(c('images/leakypipeline.png'))

## ----vectorization1------------------------------------------------------
a = runif(1e6)
b = runif(length(a))
system.time({
  z1 = numeric(length(a))
  for (i in seq(along = a))
    z1[i] = a[i]^2 * b[i]
})
system.time({
  z2 = a^2 * b
})
identical(z1, z2)

## ----vectorization2, echo = FALSE----------------------------------------
stopifnot(identical(z1, z2))

## ----Rcpp1---------------------------------------------------------------
library("Rcpp")
cppFunction("
  NumericVector myfun(NumericVector x, NumericVector y) {
    int n = x.size();
    NumericVector out(n);
    for(int i = 0; i < n; ++i) {
      out[i] = pow(x[i], 2) * y[i];
    }
    return out;
  }")
z3 = myfun(a, b)
identical(z1, z3)

## ----Rcpp2, echo = FALSE-------------------------------------------------
stopifnot(identical(z1, z3))

