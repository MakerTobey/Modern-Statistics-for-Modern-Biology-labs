## ----initialize, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
source("/Users/twenzel/Documents/MSMB20_course/R_files/00-chap.R"); chaptersetup("/Users/twenzel/Documents/MSMB20_course/data/Generative.Rnw", "1")
knitr::opts_chunk$set(dev = 'png', dpi = 100, fig.margin = TRUE, fig.show = 'hold', fig.keep = 'none')

## ---- Pile-ou-face, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/Pile_ou_face.png'))

## ----dpois1--------------------------------------------------------------
dpois(x = 3, lambda = 5)


## ----Poisson5, fig.keep = 'high', fig.cap = "Probabilities of seeing 0,1,2,...,12 mutations, as modeled by the Poisson(5) distribution. The plot shows that we will often see 4 or 5 mutations but rarely as many as 12. The distribution continues to higher numbers ($13,...$), but the probabilities will be successively smaller, and here we don\'t visualize them.", fig.width = 5.9, fig.height = 5, echo = -c(1,5)----
.oldopt = options(digits = 2)
0:12
dpois(x = 0:12, lambda = 5)
barplot(dpois(0:12, 5), names.arg = 0:12, col = "red")
options(.oldopt)



## ----genotype1-----------------------------------------------------------
genotype = c("AA","AO","BB","AO","OO","AO","AA","BO","BO",
             "AO","BB","AO","BO","AB","OO","AB","BB","AO","AO")
table(genotype)

## ----genotype2-----------------------------------------------------------
genotypeF = factor(genotype)
levels(genotypeF)
table(genotypeF)


## ---- twoballs, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Two possible events with unequal probabilities. We model this by a Bernoulli distribution with probability parameter $p=2/3$."----
knitr::include_graphics(c('images/BallsinBoxes2.jpg'))

## ----rbinom1-------------------------------------------------------------
rbinom(15, prob = 0.5, size = 1)

## ----rbinom2-------------------------------------------------------------
rbinom(12, prob = 2/3, size = 1)

## ----rbinom3-------------------------------------------------------------
rbinom(1, prob = 2/3, size = 12)

## ----rbinom4-------------------------------------------------------------
set.seed(235569515)
rbinom(1, prob = 0.3, size = 15)


## ----rbinomhide,echo=FALSE-----------------------------------------------
.freqoutcome <- (0:15)[which.max(dbinom(0:15, prob = 0.3, 15))]

## ----dbinom--------------------------------------------------------------
probabilities = dbinom(0:15, prob = 0.3, size = 15)
round(probabilities, 2)

## ----binombarplot, fig.keep = 'high', fig.cap = "Theoretical distribution of $B(15,0.3)$ . The highest bar is at $x=(ref:binombarplot-1)$. We have chosen to represent theoretical values in red throughout.", fig.width = 7, fig.height = 6----
barplot(probabilities, names.arg = 0:15, col = "red")


## ---- Simeon-Poisson, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Simeon Poisson, after whom the Poisson distribution is named (this is why it always has a capital letter, except in our R code). \\CUP{Source:https://upload.wikimedia.org/wikipedia/commons/b/b7/Simeon_Poisson.jpg}"----
knitr::include_graphics(c('images/Simeon_Poisson.png'))

## ----testbinomvspois, echo = FALSE, fig.width = 4, fig.height = 4--------
plot(dbinom(0:12, prob = 5e-4, size = 1e4),
     dpois(0:12, lambda = 5), asp = 1)
abline(a = 0, b = 1, col = "blue")

## ----dpois2--------------------------------------------------------------
5^3 * exp(-5) / factorial(3)

## ----gen-simpoisson, fig.keep = 'high', fig.cap = "Simulated distribution of B(10000, $10^{-4}$) for (ref:gen-simpoisson-1) simulations.", fig.width = 7, fig.height = 5.5----
rbinom(1, prob = 5e-4, size = 10000)
simulations = rbinom(n = 300000, prob = 5e-4, size = 10000)
barplot(table(simulations), col = "lavender")

## ---- antibody, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "A diagram of an antibody showing several immunoglobulin domains in color."----
knitr::include_graphics(c('images/Antibody_IgG2.png'))

## ----brack100, echo = FALSE----------------------------------------------
`[<-`(rep(0, 100), 22, 1)

## ----typicalP, fig.keep = 'high', fig.cap = "Plot of typical data from our generative model for the background, i.\\,e., for the false positive hits: 100 positions along the protein, at each position the count is drawn from a Poisson(0.5) random variable.", echo=FALSE, fig.width = 6, fig.height = 4----
s100 = rpois(100, lambda=0.5)
barplot(s100, ylim = c(0, 7), width = 0.7, xlim = c(-0.5,100.5),
  names.arg = seq(along = s100), col="lavender")

## ----For_the_record, eval = FALSE, echo = FALSE--------------------------
set.seed(8969311)
e100 = rpois(100,lambda=0.5)
e100[42] = 7
save(e100, file = "../data/e100.RData"

## ----epitopeData, fig.width = 6, fig.height = 4--------------------------
load("../data/e100.RData")
barplot(e100, ylim = c(0, 7), width = 0.7, xlim = c(-0.5, 100.5),
  names.arg = seq(along = e100), col = "darkolivegreen")

## ----epitopedata, fig.keep = 'high', fig.cap = "Output of the ELISA array results for 50 patients in the 100 positions.", echo = FALSE, fig.width = 6, fig.height = 4----
barplot(e100, ylim = c(0, 7), width = 0.7, xlim = c(-0.5, 100.5),
  names.arg = seq(along = e100), col = "darkolivegreen")
text(35, 7, adj = c(-0.05, 0.5), labels = "?", xpd = NA, col = "red",
  cex = 1.25, font = 2)

## ----ppois---------------------------------------------------------------
1 - ppois(6, 0.5)
ppois(6, 0.5, lower.tail = FALSE)

## ----montecarlomethod----------------------------------------------------
maxes = replicate(100000, {
  max(rpois(100, 0.5))
})
table(maxes)

## ----meanmaxes-----------------------------------------------------------
mean( maxes >= 7 )

## ---- Boxes4, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The boxes represent four outcomes or levels of a discrete **categorical** variable. The box on the right represents the more likely outcome."----
knitr::include_graphics(c('images/BallsinBoxes4.jpg'))


## ----dmultinom-----------------------------------------------------------
dmultinom(c(4, 2, 0, 0), prob = rep(1/4, 4))

## ----pvec----------------------------------------------------------------
pvec = rep(1/4, 4)
t(rmultinom(1, prob = pvec, size = 8))

## ---- SampleSize, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/SampleSize.png'))

## ----obsunder0start------------------------------------------------------
obsunder0 = rmultinom(1000, prob = pvec, size = 20)
dim(obsunder0)
obsunder0[, 1:11]

## ----assertpvec, echo = FALSE--------------------------------------------
thep = unique(pvec); stopifnot(length(thep)==1, thep == 0.25)


## ----obsunder0look-------------------------------------------------------
expected0 = pvec * 20
sum((obsunder0[, 1] - expected0)^2 / expected0)
sum((obsunder0[, 2] - expected0)^2 / expected0)
sum((obsunder0[, 3] - expected0)^2 / expected0)

## ----stat----------------------------------------------------------------
stat = function(obsvd, exptd = 20 * pvec) {
  sum((obsvd - exptd)^2 / exptd)
}
stat(obsunder0[, 1])

## ----histS0, fig.keep = 'high', fig.cap = "The histogram of simulated values `S0` of the statistic `stat` under the null (fair) distribution provides an approximation of the **sampling distribution** of the statistic `stat`."----
S0 = apply(obsunder0, 2, stat)
summary(S0)
hist(S0, breaks = 25, col = "lavender", main = "")


## ----quantile------------------------------------------------------------
q95 = quantile(S0, probs = 0.95)
q95

## ----saveS0, echo = FALSE, eval = FALSE----------------------------------
## ## This was done to save this object for its reuse in Chapter 2.
## save(S0, file="../data/S0.RData")

## ---- roulette-chunk-1, eval = TRUE, echo = FALSE, fig.keep = 'high'-----
knitr::include_graphics('images/roulette.png', dpi = 600)

## ----alternativeA--------------------------------------------------------
pvecA = c(3/8, 1/4, 3/12, 1/8)
observed = rmultinom(1000, prob = pvecA, size = 20)
dim(observed)
observed[, 1:7]
apply(observed, 1, mean)
expectedA = pvecA * 20
expectedA

## ----computestat1--------------------------------------------------------
stat(observed[, 1])
S1 = apply(observed, 2, stat)
q95
sum(S1 > q95)
power = mean(S1 > q95)
power

## ----assertS1, echo = FALSE----------------------------------------------
stopifnot(stat(observed[, 1]) < q95)


## ---- ProbaDiagram, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "We have studied how a probability model has a distribution we call this $F$. $F$ often depends on parameters, which are denoted by Greek letters, such as $\\theta$. The observed data are generated via the brown arrow and are represented by Roman letters such as $x$. The vertical bar in the probability computation stands for **supposing that** or **conditional on**"----
knitr::include_graphics(c('images/ProbaDiagram.png'))


## ------------------------------------------------------------------------
dbinom(2, size = 10, prob = 0.3)
pbinom(2, size = 10, prob = 0.3)
sum(sapply(0:2, dbinom, size = 10, prob = 0.3))

## ----poismax, echo = FALSE-----------------------------------------------
poismax = function(lambda, n, m) {
  epsilon = 1 - ppois(m - 1, lambda)
  1 - exp( -n * epsilon)
}

## ----posimaxout----------------------------------------------------------
poismax(lambda = 0.5, n = 100, m = 7)
poismax(lambda = mean(e100), n = 100, m = 7)

## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## BiocManager::install(c("Biostrings", "BSgenome.Celegans.UCSC.ce2"))

## ----BSgenome.Celegans.UCSC.ce2, message=FALSE, warning=FALSE------------
library("BSgenome.Celegans.UCSC.ce2")
Celegans
seqnames(Celegans)
Celegans$chrM
class(Celegans$chrM)
length(Celegans$chrM)

## ----Biostr--------------------------------------------------------------
library("Biostrings")
lfM = letterFrequency(Celegans$chrM, letters=c("A", "C", "G", "T"))
lfM
sum(lfM)
lfM / sum(lfM)

## ----rf------------------------------------------------------------------
t(rmultinom(1, length(Celegans$chrM), p = rep(1/4, 4)))

## ----lengthM-------------------------------------------------------------
length(Celegans$chrM) / 4

## ----oestat--------------------------------------------------------------
oestat = function(o, e) {
  sum((e-o)^2 / e)
}
oe = oestat(o = lfM, e = length(Celegans$chrM) / 4)
oe

## ----oesim---------------------------------------------------------------
B = 10000
n = length(Celegans$chrM)
expected = rep(n / 4, 4)
oenull = replicate(B,
  oestat(e = expected, o = rmultinom(1, n, p = rep(1/4, 4))))

## ----resundernull, echo=FALSE--------------------------------------------
hist(oenull, breaks = 100, col = "skyblue", main="")

## ----assert, echo = FALSE------------------------------------------------
stopifnot( oe/10 > max(oenull) )

