## ----initialize, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
source("../chapter-setup.R"); chaptersetup("/Users/Susan/Courses/CUBook-html/CUBook/Chap2-Stat/FreqandBayesandModels.Rnw", "2")
knitr::opts_chunk$set(dev = 'png', dpi = 100, fig.margin = TRUE, fig.show = 'hold', fig.keep = 'none')

## ---- StatDiagram, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/StatDiagram.png'))


## ---- Parameters, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/Parameters.png'))

## ---- probabilitydiag, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The probabilistic model we obtained in Chapter \\@ref(Chap:Generative). The data are represented as $x$ in green. We can use the observed data to compute the probability if observing $x$ when we know the true value of $\\theta$."----
knitr::include_graphics(c('images/ProbaDiagram.png'))

## ----e99-----------------------------------------------------------------
load("../data/e100.RData")
e99 = e100[-which.max(e100)]

## ----twopoisson, fig.keep = 'high', fig.cap = "The observed distribution of the epitope data without the outlier.", fig.width = 4.5, fig.height = 4----
barplot(table(e99), space = 0.8, col = "chartreuse4")

## ----stat-rooto, fig.keep = 'high', fig.cap = "Rootogram showing the square root of the theoretical values as red dots and the square root of the observed frequencies as drop down rectangles. (We\'ll see a bit below how the `goodfit` function decided which $\\lambda$ to use.)", fig.width= 3, fig.height = 3.5----
library("vcd")
gf1 = goodfit( e99, "poisson")
rootogram(gf1, xlab = "", rect_gp = gpar(fill = "chartreuse4"))

## ----RootogramPoisson, fig.width= 3, fig.height = 3.5--------------------
simp = rpois(100, lambda = 0.05)
gf2 = goodfit(simp, "poisson")
rootogram(gf2, xlab = "")


## ----table100------------------------------------------------------------
table(e100)

## ----table3--------------------------------------------------------------
table(rpois(100, 3))

## ----table101, echo = FALSE----------------------------------------------
counts  =  table(e100)
stopifnot(identical(names(counts), c("0", "1", "2", "7")), all(counts==c(58, 34, 7, 1)))

## ----poism3--------------------------------------------------------------
prod(dpois(c(0, 1, 2, 7), lambda = 3) ^ (c(58, 34, 7, 1)))

## ----anspois-------------------------------------------------------------
prod(dpois(c(0, 1, 2, 7), lambda = 0.4) ^ (c(58, 34, 7, 1)))

## ----functionll----------------------------------------------------------
loglikelihood  =  function(lambda, data = e100) {
  sum(log(dpois(data, lambda)))
}

## ----chap2-r-poislikel-1, fig.keep = 'high', fig.cap = "The red curve is the log-likelihood function. The vertical line shows the value of `m` (the mean) and the horizontal line the log-likelihood of `m`. It looks like `m` maximizes the likelihood.", fig.width=3.5----
lambdas = seq(0.05, 0.95, length = 100)
loglik = vapply(lambdas, loglikelihood, numeric(1))
plot(lambdas, loglik, type = "l", col = "red", ylab = "", lwd = 2,
     xlab = expression(lambda))
m0 = mean(e100)
abline(v = m0, col = "blue", lwd = 2)
abline(h = loglikelihood(m0), col = "purple", lwd = 2)
m0

## ----gfpoisson-----------------------------------------------------------
gf  =  goodfit(e100, "poisson")
names(gf)
gf$par

## ----colorblind, echo = FALSE--------------------------------------------
cb  =  c(rep(0, 110), rep(1, 10))

## ----cb------------------------------------------------------------------
table(cb)


## ----likely1, fig.keep = 'high', fig.cap = "Plot of the likelihood as a function of the probabilities. The likelihood is a function on $[0, 1]$; here we have zoomed into the range of $[(ref:likely1-1), (ref:likely1-2)]$, as the likelihood is practically zero for larger values of $p$.", fig.width = 4----
probs  =  seq(0, 0.3, by = 0.005)
likelihood = dbinom(sum(cb), prob = probs, size = length(cb))
plot(probs, likelihood, pch = 16, xlab = "probability of success",
       ylab = "likelihood", cex=0.6)
probs[which.max(likelihood)]

## ----check, echo = FALSE-------------------------------------------------
stopifnot(abs(probs[which.max(likelihood)]-1/12) < diff(probs[1:2]))


## ----loglike1------------------------------------------------------------
loglikelihood = function(theta, n = 300, k = 40) {
  115 + k * log(theta) + (n - k) * log(1 - theta)
}

## ----chap2-r-loglikelihood-1, fig.keep = 'high', fig.cap = "Plot of the log likelihood function for $n=300$ and $y=40$.", fig.width = 3, fig.height = 3.2----
thetas = seq(0, 1, by = 0.001)
plot(thetas, loglikelihood(thetas), xlab = expression(theta),
  ylab = expression(paste("log f(", theta, " | y)")),type = "l")

## ----staph---------------------------------------------------------------
library("Biostrings")
staph = readDNAStringSet("../data/staphsequence.ffn.txt", "fasta")

## ----firstgenestaph------------------------------------------------------
staph[1]
letterFrequency(staph[[1]], letters = "ACGT", OR = 0)

## ----compareprop---------------------------------------------------------
letterFrq = vapply(staph, letterFrequency, FUN.VALUE = numeric(4),
         letters = "ACGT", OR = 0)
colnames(letterFrq) = paste0("gene", seq(along = staph))
tab10 = letterFrq[, 1:10]
computeProportions = function(x) { x/sum(x) }
prop10 = apply(tab10, 2, computeProportions)
round(prop10, digits = 2)
p0 = rowMeans(prop10)
p0

## ----outerex-------------------------------------------------------------
cs = colSums(tab10)
cs
expectedtab10 = outer(p0, cs, FUN = "*")
round(expectedtab10)

## ----genrandomtabs-------------------------------------------------------
randomtab10 = sapply(cs, function(s) { rmultinom(1, s, p0) } )
all(colSums(randomtab10) == cs)

## ----assertgenrandomtabs, echo = FALSE-----------------------------------
stopifnot(all(colSums(randomtab10) == cs))

## ----chap2-r-quant12-1, fig.keep = 'high', fig.cap = "Histogram of `simulstat`. The value of `S1` is marked by the vertical red line, those of the 0.95 and 0.99 quantiles (see next section) by the dotted lines.", fig.width = 4, fig.height = 3.5----
stat = function(obsvd, exptd = 20 * pvec) {
   sum((obsvd - exptd)^2 / exptd)
}
B = 1000
simulstat = replicate(B, {
  randomtab10 = sapply(cs, function(s) { rmultinom(1, s, p0) })
  stat(randomtab10, expectedtab10)
})
S1 = stat(tab10, expectedtab10)
sum(simulstat >= S1)

hist(simulstat, col = "lavender", breaks = seq(0, 75, length.out=50))
abline(v = S1, col = "red")
abline(v = quantile(simulstat, probs = c(0.95, 0.99)),
       col = c("darkgreen", "blue"), lty = 2)

## ----checksimulstat, echo = FALSE----------------------------------------
stopifnot(max(simulstat)<75, S1<75)

## ----quantiles3, results = "hide"----------------------------------------
qs = ppoints(100)
quantile(simulstat, qs)
quantile(qchisq(qs, df = 30), qs)


## ----chap2-r-qqplot3-1, fig.keep = 'high', fig.cap = "Our simulated statistic\'s distribution compared to $\\chi_{30}^2$ using a QQ-plot, which shows the theoretical **quantiles** for the $\\chi^2_{30}$ distribution on the horizontal axis and the sampled ones on the vertical axis.", fig.width = 3.4, fig.height = 4----
qqplot(qchisq(ppoints(B), df = 30), simulstat, main = "",
  xlab = expression(chi[nu==30]^2), asp = 1, cex = 0.5, pch = 16)
abline(a = 0, b = 1, col = "red")

## ----pvalueBias----------------------------------------------------------
1 - pchisq(S1, df = 30)

## ---- ChargaffColdSpring, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/ChargaffColdSpring.png'))

## ----Chargaff------------------------------------------------------------
load("../data/ChargaffTable.RData")
ChargaffTable

## ----ChargaffBars, fig.keep = 'high', fig.cap = "Barplots for the different rows in `ChargaffTable`. Can you spot the pattern?", fig.margin = FALSE, echo = FALSE, fig.width = 7, fig.height = 3.4----
stopifnot(nrow(ChargaffTable) == 8)
mycolors = c("chocolate", "aquamarine4", "cadetblue4", "coral3",
            "chartreuse4","darkgoldenrod4","darkcyan","brown4")
par(mfrow=c(2, 4), mai = c(0, 0.7, 0.7, 0))
for (i in 1:8) {
  cbp = barplot(ChargaffTable[i, ], horiz = TRUE, axes = FALSE, axisnames = FALSE, col = mycolors[i])
  ax = axis(3, las = 2, labels = FALSE, col = mycolors[i], cex = 0.5, at = c(0, 10, 20))
  mtext(side = 3, at = ax,  text = paste(ax), col = mycolors[i], line = 0, las = 1, cex = 0.9)
  mtext(side = 2, at = cbp, text = colnames(ChargaffTable), col = mycolors[i], line = 0, las = 2, cex = 1)
  title(paste(rownames(ChargaffTable)[i]), col = mycolors[i], cex = 1.1)
}

## ----chap2-r-permstatChf-1, fig.keep = 'high', fig.cap = "Histogram of our statistic `statChf` computed from simulations using per-row permutations of the columns. The value it yields for the observed data is shown by the red line.", fig.width = 3.2, fig.height = 3----
statChf = function(x){
  sum((x[, "C"] - x[, "G"])^2 + (x[, "A"] - x[, "T"])^2)
}
chfstat = statChf(ChargaffTable)
permstat = replicate(100000, {
     permuted = t(apply(ChargaffTable, 1, sample))
     colnames(permuted) = colnames(ChargaffTable)
     statChf(permuted)
})
pChf = mean(permstat <= chfstat)
pChf
hist(permstat, breaks = 100, main = "", col = "lavender")
abline(v = chfstat, lwd = 2, col = "red")

## ----vcdHC---------------------------------------------------------------
HairEyeColor[,, "Female"]

## ----answerHC------------------------------------------------------------
str(HairEyeColor)
? HairEyeColor

## ----Deuto---------------------------------------------------------------
load("../data/Deuteranopia.RData")
Deuteranopia

## ----chisq.test.Deuteranopia---------------------------------------------
chisq.test(Deuteranopia)

## ----chap2-r-HardyWeinberg-1, fig.keep = 'high', fig.cap = "Plot of the log-likelihood for the (ref:chap2-r-HardyWeinberg-1-1) data."----
library("HardyWeinberg")
data("Mourant")
Mourant[214:216,]
nMM = Mourant$MM[216]
nMN = Mourant$MN[216]
nNN = Mourant$NN[216]
loglik = function(p, q = 1 - p) {
  2 * nMM * log(p) + nMN * log(2*p*q) + 2 * nNN * log(q)
}
xv = seq(0.01, 0.99, by = 0.01)
yv = loglik(xv)
plot(x = xv, y = yv, type = "l", lwd = 2,
     xlab = "p", ylab = "log-likelihood")
imax = which.max(yv)
abline(v = xv[imax], h = yv[imax], lwd = 1.5, col = "blue")
abline(h = yv[imax], lwd = 1.5, col = "purple")

## ----phat----------------------------------------------------------------
phat  =  af(c(nMM, nMN, nNN))
phat
pMM   =  phat^2
qhat  =  1 - phat

## ----hweq----------------------------------------------------------------
pHW = c(MM = phat^2, MN = 2*phat*qhat, NN = qhat^2)
sum(c(nMM, nMN, nNN)) * pHW

## ----HWtern, fig.keep = 'high', fig.cap = "This **de Finetti plot** shows the points as barycenters of the three genotypes using the frequencies as weights on each of the corners of the triangle. The Hardy-Weinberg model is the red curve, the acceptance region is between the two purple lines. We see that the US is the furthest from being in HW equilibrium.", fig.margin = FALSE, message = FALSE, warning = FALSE, fig.width = 3.4, fig.height = 3.4, results = FALSE, echo = -1----
par(mai = rep(0.1, 4))
pops = c(1, 69, 128, 148, 192)
genotypeFrequencies = as.matrix(Mourant[, c("MM", "MN", "NN")])
HWTernaryPlot(genotypeFrequencies[pops, ],
        markerlab = Mourant$Country[pops],
        alpha = 0.0001, curvecols = c("red", rep("purple", 4)),
        mcex = 0.75, vertex.cex = 1)

## ----quesTern, echo = -1-------------------------------------------------
HWTernaryPlot(genotypeFrequencies[pops, ],
        markerlab = Mourant$Country[pops],
        alpha = 0.0001, curvecols = c("red", rep("purple", 4)),
        mcex = 0.75, vertex.cex = 1)
HWTernaryPlot(genotypeFrequencies[-pops, ], alpha = 0.0001,
   newframe = FALSE, cex = 0.5)

## ----newMNdata-----------------------------------------------------------
newgf = round(genotypeFrequencies / 50)
HWTernaryPlot(newgf[pops, ],
        markerlab = Mourant$Country[pops],
        alpha = 0.0001, curvecols = c("red", rep("purple", 4)),
        mcex = 0.75, vertex.cex = 1)

## ----chap2-r-seqlogo-1, fig.keep = 'high', fig.cap = "Here is a diagram called a sequence logo for the position dependent multinomial used to model the Kozak motif. It codifies the amount of variation in each of the positions on a log scale. The large letters represent positions where there is no uncertainty about which nucleotide occurs.", fig.margin = FALSE, fig.height=5, fig.width=5----
library("seqLogo")
load("../data/kozak.RData")
kozak
pwm = makePWM(kozak)
seqLogo(pwm, ic.scale = FALSE)

## ----4stateMC, echo = FALSE----------------------------------------------
library("markovchain")
library("igraph")
sequence = toupper(c("a", "c", "a", "c", "g", "t", "t", "t", "t", "c",
"c", "a", "c", "g", "t", "a", "c","c","c","a","a","a","t","a",
"c","g","g","c","a","t","g","t","g","t","g","a","g","c","t","g"))
mcFit   =  markovchainFit(data = sequence)
MCgraph =  markovchain:::.getNet(mcFit$estimate, round = TRUE)
edgelab =  round(E(MCgraph)$weight / 100, 1)

## ----statsfourstateMC, fig.keep = 'high', fig.cap = "Visualisation of a 4-state Markov chain. The probability of each possible digram (e.\\,g., CA) is given by the weight of the edge between the corresponding nodes. So for instance, the probability of CA is given by the edge C$\\to$ A. We\'ll see in Chapter \\@ref(Chap:Images) how to use **R** packages to draw these type of network graphs.", echo = FALSE, fig.width = 4, fig.height = 3.5----
par(mai=c(0,0,0,0))
#plot.igraph(MCgraph, edge.label = edgelab, width = 2, edge.arrow.width = 1.5,
# vertex.size = 40, xlim = c(-1, 1.25))
plot.igraph(MCgraph, edge.label = edgelab,
       #edge.arrow.width = 1.5,
       vertex.size = 40, xlim = c(-1, 1.25))

## ---- FreqBayes-turtles, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Turtles all the way down. Bayesian modeling of the uncertainty of the parameter of a distribution is done by using a random variable whose distribution may depend on parameters whose uncertainty can be modeled as a random variable; these are called hierarchical models."----
knitr::include_graphics(c('images/turtlesalltheway.png'))

## ---- STRDefinition, fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "A short tandem repeat (STR) in DNA occurs when a pattern of two or more nucleotides is repeated and the repeated sequences are directly adjacent to each other. An STR is also known as a microsatellite. The pattern can range in length from 2 to 13 nucleotides, and the number of repeats is highly variable across individuals. STR numbers can be used as genetic signatures, called haplotypes."----
knitr::include_graphics(c('images/STRDefinition.png'))

## ----haplo6--------------------------------------------------------------
haplo6=read.table("../data/haplotype6.txt",header = TRUE)
haplo6

## ----histobeta2, fig.keep = 'high', fig.cap = "Beta distributions with $\\alpha=10,20,50$ and $\\beta=30,60,150$ used as a {prior} for a probability of success. These three distributions have the same mean ($\\frac{\\alpha}{\\alpha +\\beta}$), but different concentrations around the mean.", echo = FALSE, fig.width = 3.5, fig.height = 3.5----
theta = thetas[1:500]
dfbetas = data.frame(theta,
           db1= dbeta(theta,10,30),
           db2 = dbeta(theta, 20, 60),
           db3 = dbeta(theta, 50, 150))
require(reshape2)
datalong  =  melt(dfbetas, id="theta")
ggplot(datalong) +
geom_line(aes(x = theta,y=value,colour=variable)) +
theme(legend.title=element_blank()) +
geom_vline(aes(xintercept=0.25), colour="#990000", linetype="dashed")+
scale_colour_discrete(name  ="Prior",
                          labels=c("B(10,30)", "B(20,60)","B(50,150)"))

## ----chap2-r-histmarginal-1, fig.keep = 'high', fig.cap = "Marginal Distribution of $Y$.", fig.width = 3.5, fig.height = 3.5----
rtheta = rbeta(100000, 50, 350)
y = vapply(rtheta, function(th) {
  rbinom(1, prob = th, size = 300)
}, numeric(1))
hist(y, breaks = 50, col = "orange", main = "", xlab = "")

## ----freqquesvectorize, echo = FALSE-------------------------------------
set.seed(0xbebe)
.w1 = vapply(rtheta, function(th) rbinom(1, prob = th, size = 300), integer(1))
set.seed(0xbebe)
.w2 = rbinom(length(rtheta), rtheta, size = 300)
stopifnot(identical(.w1, .w2))

## ----chap2-r-densityposterior-1, fig.keep = 'high', fig.cap = "Only choosing the values of the distribution with $Y=40$ gives the posterior distribution of $\\theta$. The histogram (green) shows the simulated values for the posteriror distribution, the line the theoretical density of a beta distribution with the theoretical posterior parameters.", fig.width = 3.5, fig.height = 3.5----
thetaPostEmp = rtheta[ y == 40 ]
hist(thetaPostEmp, breaks = 40, col = "chartreuse4", main = "",
  probability = TRUE, xlab = expression("posterior"~theta))
densPostTheory  =  dbeta(thetas, 90, 610)
lines(thetas, densPostTheory, type="l", lwd = 3)

## ----comparetheory1------------------------------------------------------
mean(thetaPostEmp)
dtheta = thetas[2]-thetas[1]
sum(thetas * densPostTheory * dtheta)

## ----comparetheory2, echo = FALSE----------------------------------------
stopifnot(abs(mean(thetaPostEmp) - sum(thetas * densPostTheory * dtheta)) < 1e-3)

## ----mcint---------------------------------------------------------------
thetaPostMC = rbeta(n = 1e6, 90, 610)
mean(thetaPostMC)

## ----chap2-r-qqplotbeta-1, fig.keep = 'high', fig.cap = "QQ-plot of our Monte Carlo sample `thetaPostMC` from the theoretical distribution and our simulation sample `thetaPostEmp`. We could also similarly compare either of these two distributions to the theoretical distribution function `pbeta(., 90, 610)`. If the curve lies on the line $y=x$ this indicates a good agreement. There are some random differences at the tails.", fig.width = 3.5, fig.height = 3.5----
qqplot(thetaPostMC, thetaPostEmp, type = "l", asp = 1)
abline(a = 0, b = 1, col = "blue")

## ----postbeta------------------------------------------------------------
densPost2 = dbeta(thetas, 115, 735)
mcPost2   = rbeta(1e6, 115, 735)

sum(thetas * densPost2 * dtheta)  # mean, by numeric integration
mean(mcPost2)                     # mean, by MC
thetas[which.max(densPost2)]      # MAP estimate

## ---- roulette-chunk-1, eval = TRUE, echo = FALSE, fig.keep = 'high'-----
knitr::include_graphics('images/roulette.png', dpi = 600)

## ----quantilespost-------------------------------------------------------
quantile(mcPost2, c(0.025, 0.975))

## ---- DESeq2-Prediction-Interval, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "An example from  @LoveDESeq2 shows plots of the likelihoods (solid lines, scaled to integrate to 1) and the posteriors (dashed lines) for the green and purple genes and of the prior (solid black line): due to the higher dispersion of the purple gene, its likelihood is wider and less peaked (indicating less information), and the prior has more influence on its posterior than for the green gene. The stronger curvature of the green posterior at its maximum translates to a smaller reported standard error for the MAP logarithmic fold change (LFC) estimate (horizontal error bar)."----
knitr::include_graphics(c('images/DESeq2-Prediction-Interval.png'))

## ----callBios------------------------------------------------------------
library("Biostrings")

## ----BiostringExplore, results = "hide", eval = FALSE--------------------
## GENETIC_CODE
## IUPAC_CODE_MAP
## vignette(package = "Biostrings")
## vignette("BiostringsQuickOverview", package = "Biostrings")

## ----BiostringCheck, echo=FALSE, results = "hide"------------------------
GENETIC_CODE
IUPAC_CODE_MAP

## ----BSgenomes-----------------------------------------------------------
library("BSgenome")
ag = available.genomes()
length(ag)
ag[1:2]

## ----BSGenomeEcoli, results="hide"---------------------------------------
library("BSgenome.Ecoli.NCBI.20080805")
Ecoli
shineDalgarno = "AGGAGGT"
ecoli = Ecoli$NC_010473

## ----window--------------------------------------------------------------
window = 50000
starts = seq(1, length(ecoli) - window, by = window)
ends   = starts + window - 1
numMatches = vapply(seq_along(starts), function(i) {
  countPattern(shineDalgarno, ecoli[starts[i]:ends[i]],
               max.mismatch = 0)
  }, numeric(1))
table(numMatches)

## ----poissonness, fig.keep = 'high', fig.cap = "Evaluation of a Poisson model for motif counts along the sequence \\texttt{Ecoli$NC_010473}.", fig.width = 4, fig.height = 5----
library("vcd")
gf = goodfit(numMatches, "poisson")
summary(gf)
distplot(numMatches, type = "poisson")

## ----pattlocIranges1, results="hide"-------------------------------------
sdMatches = matchPattern(shineDalgarno, ecoli, max.mismatch = 0)

## ----pattlocIranges2-----------------------------------------------------
betweenmotifs = gaps(sdMatches)

## ----chap2-r-expplotdata-1, fig.keep = 'high', fig.cap = "Evaluation of fit to the exponential distribution (black line) of the gaps between the motifs.", fig.width = 3.5, fig.height = 4----
library("Renext")
expplot(width(betweenmotifs), rate = 1/mean(width(betweenmotifs)),
        labels = "fit")

## ----gof, echo = FALSE, eval = FALSE-------------------------------------
## gofExp.test(width(betweenmotifs))

## ----chr8HS--------------------------------------------------------------
library("BSgenome.Hsapiens.UCSC.hg19")
chr8  =  Hsapiens$chr8
CpGtab = read.table("../data/model-based-cpg-islands-hg19.txt",
                    header = TRUE)
nrow(CpGtab)
head(CpGtab)
irCpG = with(dplyr::filter(CpGtab, chr == "chr8"),
         IRanges(start = start, end = end))


## ----grCpG---------------------------------------------------------------
grCpG = GRanges(ranges = irCpG, seqnames = "chr8", strand = "+")
genome(grCpG) = "hg19"

## ----freqandbayes-ideo, fig.keep = 'high', fig.cap = "**[Gviz](https://bioconductor.org/packages/Gviz/)** plot of CpG locations in a selected region of chromosome 8.", fig.height = 4----
library("Gviz")
ideo = IdeogramTrack(genome = "hg19", chromosome = "chr8")
plotTracks(
  list(GenomeAxisTrack(),
    AnnotationTrack(grCpG, name = "CpG"), ideo),
    from = 2200000, to = 5800000,
    shape = "box", fill = "#006400", stacking = "dense")

## ----CGIview-------------------------------------------------------------
CGIview    = Views(unmasked(Hsapiens$chr8), irCpG)
NonCGIview = Views(unmasked(Hsapiens$chr8), gaps(irCpG))

## ----CGIview2------------------------------------------------------------
seqCGI      = as(CGIview, "DNAStringSet")
seqNonCGI   = as(NonCGIview, "DNAStringSet")
dinucCpG    = sapply(seqCGI, dinucleotideFrequency)
dinucNonCpG = sapply(seqNonCGI, dinucleotideFrequency)
dinucNonCpG[, 1]
NonICounts = rowSums(dinucNonCpG)
IslCounts  = rowSums(dinucCpG)

## ----transitions---------------------------------------------------------
TI  = matrix( IslCounts, ncol = 4, byrow = TRUE)
TnI = matrix(NonICounts, ncol = 4, byrow = TRUE)
dimnames(TI) = dimnames(TnI) =
  list(c("A", "C", "G", "T"), c("A", "C", "G", "T"))


## ----MI------------------------------------------------------------------
MI = TI /rowSums(TI)
MI
MN = TnI / rowSums(TnI)
MN

## ----STATI---------------------------------------------------------------
freqIsl = alphabetFrequency(seqCGI, baseOnly = TRUE, collapse = TRUE)[1:4]
freqIsl / sum(freqIsl)
freqNon = alphabetFrequency(seqNonCGI, baseOnly = TRUE, collapse = TRUE)[1:4]
freqNon / sum(freqNon)

## ---- book-chunk-1, eval = TRUE, echo = FALSE, fig.keep = 'high'---------
knitr::include_graphics('images/book_icon.png', dpi = 400)

## ----alphabeta-----------------------------------------------------------
alpha = log((freqIsl/sum(freqIsl)) / (freqNon/sum(freqNon)))
beta  = log(MI / MN)

## ----scorepatt-----------------------------------------------------------
x = "ACGTTATACTACG"
scorefun = function(x) {
  s = unlist(strsplit(x, ""))
  score = alpha[s[1]]
  if (length(s) >= 2)
    for (j in 2:length(s))
      score = score + beta[s[j-1], s[j]]
  score
}
scorefun(x)

## ----scorefun1-----------------------------------------------------------
generateRandomScores = function(s, len = 100, B = 1000) {
  alphFreq = alphabetFrequency(s)
  isGoodSeq = rowSums(alphFreq[, 5:ncol(alphFreq)]) == 0
  s = s[isGoodSeq]
  slen = sapply(s, length)
  prob = pmax(slen - len, 0)
  prob = prob / sum(prob)
  idx  = sample(length(s), B, replace = TRUE, prob = prob)
  ssmp = s[idx]
  start = sapply(ssmp, function(x) sample(length(x) - len, 1))
  scores = sapply(seq_len(B), function(i)
    scorefun(as.character(ssmp[[i]][start[i]+(1:len)]))
  )
  scores / len
}
scoresCGI    = generateRandomScores(seqCGI)
scoresNonCGI = generateRandomScores(seqNonCGI)

## ----chap2-r-ScoreMixture-1, fig.keep = 'high', fig.cap = "Island and non-island scores as generated by the function `generateRandomScores`. This is the first instance of a **mixture** we encounter. We will revisit them in Chapter \\@ref(Chap:Mixtures).", fig.height = 5----
br = seq(-0.6, 0.7, length.out = 50)
h1 = hist(scoresCGI,    breaks = br, plot = FALSE)
h2 = hist(scoresNonCGI, breaks = br, plot = FALSE)
plot(h1, col = rgb(0, 0, 1, 1/4), xlim = c(-0.5, 0.5), ylim=c(0,120))
plot(h2, col = rgb(1, 0, 0, 1/4), add = TRUE)

## ----savescoresforChap4, echo = FALSE, eval=FALSE------------------------
## ###This is for provenance reasons, keep track of how the data
## ###were generated for the EM exercise in Chapter 4.
## Mdata=c(scoresCGI,scoresNonCGI)
## MM1=sample(Mdata[1:1000],800)
## MM2=sample(Mdata[1001:2000],1000)
## Myst=c(MM1,MM2);names(Myst)=NULL
## saveRDS(c(MM1,MM2),"../data/Myst.rds")
## ###True value of m1,m2,s1 and s2
## ###

## ----checkhists, echo = FALSE--------------------------------------------
stopifnot(max(h1$counts) < 120, max(h2$counts) < 120,
          h1$breaks[1] >= br[1], h1$breaks[length(h1$breaks)] <= br[length(br)],
          h2$breaks[1] >= br[1], h2$breaks[length(h2$breaks)] <= br[length(br)])

mtb = read.table("../data/M_tuberculosis.txt", header = TRUE)
head(mtb, n = 4)

## ----ProlMyc-------------------------------------------------------------
pro  =  mtb[ mtb$AmAcid == "Pro", "Number"]
pro/sum(pro)

staph = readDNAStringSet("../data/staphsequence.ffn.txt", "fasta")

## ----staphex-------------------------------------------------------------
staph[1:3, ]
staph

## ----GCfreq--------------------------------------------------------------
letterFrequency(staph[[1]], letters = "ACGT", OR = 0)
GCstaph = data.frame(
  ID = names(staph),
  GC = rowSums(alphabetFrequency(staph)[, 2:3] / width(staph)) * 100
)

## ----chap2-r-SlidingGC-1, fig.keep = 'high', fig.cap = "GC content along sequence 364 of the *Staphylococcus Aureus* genome.", fig.width = 7----
window = 100
gc = rowSums( letterFrequencyInSlidingView(staph[[364]], window,
      c("G","C")))/window
plot(x = seq(along = gc), y = gc, type = "l")

## ----chap2-r-SmoothSlidingGC-1, fig.keep = 'high', fig.cap = "Smoothed GC content along sequence 364 of the *Staphylococcus Aureus* genome.", fig.width = 7----
plot(x = seq(along = gc), y = gc, type = "l")
lines(lowess(x = seq(along = gc), y = gc, f = 0.2), col = 2)

## ----histobeta4, fig.width = 3.5, fig.height = 3.5-----------------------
theta = thetas[1:500]
dfbetas = data.frame(theta,
           db1=  dbeta(theta,0.5,0.5),
           db2= dbeta(theta,1,1),
           db3= dbeta(theta,10,30),
           db4 = dbeta(theta, 20, 60),
           db5 = dbeta(theta, 50, 150))
require(reshape2)
datalong  =  melt(dfbetas, id="theta")
ggplot(datalong) +
geom_line(aes(x = theta,y=value,colour=variable)) +
theme(legend.title=element_blank()) +
geom_vline(aes(xintercept=0.25), colour="#990000", linetype="dashed")+
scale_colour_discrete(name  ="Prior",
                          labels=c("B(0.5,0.5)","U(0,1)=B(1,1)","B(10,30)", "B(20,60)","B(50,150)"))

