## ----initialize, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
source("../chapter-setup.R"); chaptersetup("/Users/Susan/Courses/CUBook-html/CUBook/Chap5-Mixtures/MixtureModel.Rnw", "5")
knitr::opts_chunk$set(dev = 'png', dpi = 100, fig.margin = TRUE, fig.show = 'hold', fig.keep = 'none')

## ---- t-distribution, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/t_distribution.png'))

## ----twocoins, fig.keep = 'high', fig.cap = "Histogram of 10,000 random draws from a fair mixture of two normals. The left hand part of the histogram is dominated by numbers generated from (A), on the right from (B).", fig.width = 3.5, fig.height = 3.5----
coinflips = (runif(10000) > 0.5)
table(coinflips)
oneFlip = function(fl, mean1 = 1, mean2 = 3, sd1 = 0.5, sd2 = 0.5) {
  if (fl) {
   rnorm(1, mean1, sd1)
  } else {
   rnorm(1, mean2, sd2)
  }
}
fairmix = vapply(coinflips, oneFlip, numeric(1))
library("ggplot2")
library("dplyr")
ggplot(tibble(value = fairmix), aes(x = value)) +
     geom_histogram(fill = "purple", binwidth = 0.1)

## ----efficient-----------------------------------------------------------
means = c(1, 3)
sds   = c(0.25, 0.25)
values = rnorm(length(coinflips),
          mean = ifelse(coinflips, means[1], means[2]),
          sd   = ifelse(coinflips, sds[1],   sds[2]))

## ----limitinghistogram---------------------------------------------------
fair = tibble(
  coinflips = (runif(1e6) > 0.5),
  values = rnorm(length(coinflips),
               mean = ifelse(coinflips, means[1], means[2]),
               sd   = ifelse(coinflips, sds[1],   sds[2])))
ggplot(fair, aes(x = values)) +
     geom_histogram(fill = "purple", bins = 500)

## ----overlaydensity, fig.keep = 'high', fig.cap = "In purple: histogram of half a million draws from the normal distribution $N(\\mu = (ref:overlaydensity-1), \\sigma^2 = (ref:overlaydensity-2)^2)$. The curve is the theoretical density $\\phi(x)$ calculated using the function `dnorm`."----
ggplot(dplyr::filter(fair, coinflips), aes(x = values)) +
   geom_histogram(aes(y = ..density..), fill = "purple",
                  binwidth = 0.01) +
   stat_function(fun = dnorm,
          args = list(mean = means[1], sd = sds[1]), color = "red")

## ----twodensity, fig.keep = 'high', fig.cap = "The theoretical density of the mixture.", fig.width = 3.5, fig.height = 3.5----
fairtheory = tibble(
  x = seq(-1, 5, length.out = 1000),
  f = 0.5 * dnorm(x, mean = means[1], sd = sds[1]) +
      0.5 * dnorm(x, mean = means[2], sd = sds[2]))
ggplot(fairtheory, aes(x = x, y = f)) +
  geom_line(color = "red", size = 1.5) + ylab("mixture density")

## ----histmystery, fig.keep = 'high', fig.cap = "A mixture of two normals that is harder to recognize.", echo = FALSE, fig.width = 3.5, fig.height = 3, results = "hide", message = FALSE, warning = FALSE----
mystery = tibble(
  coinflips = (runif(1e3) > 0.5),
  values = rnorm(length(coinflips),
               mean = ifelse(coinflips, 1, 2),
               sd   = ifelse(coinflips, sqrt(.5), sqrt(.5))))
br2 = with(mystery, seq(min(values), max(values), length.out = 30))
ggplot(mystery, aes(x = values)) +
geom_histogram(fill = "purple", breaks = br2)

## ----chap5-r-betterhistogram-1, fig.keep = 'high', fig.cap = "The mixture from Figure \\@ref(fig:histmystery), but with the two components colored in red and blue.", fig.width = 3.5, fig.height = 3----
head(mystery, 3)
br = with(mystery, seq(min(values), max(values), length.out = 30))
ggplot(mystery, aes(x = values)) +
  geom_histogram(data = dplyr::filter(mystery, coinflips),
     fill = "red", alpha = 0.2, breaks = br) +
  geom_histogram(data = dplyr::filter(mystery, !coinflips),
     fill = "darkblue", alpha = 0.2, breaks = br) 

## ----checkbr, echo = FALSE-----------------------------------------------
stopifnot(identical(br2, br))

## ----chap5-r-comparecomponents-1, fig.keep = 'high', fig.cap = "As Figure \\@ref(fig:chap5-r-betterhistogram-1), with stacked bars for the two mixture components.", fig.width = 3.5, fig.height = 3----
ggplot(mystery, aes(x = values, fill = coinflips)) +
  geom_histogram(data = dplyr::filter(mystery, coinflips),
     fill = "red", alpha = 0.2, breaks = br) +
  geom_histogram(data = dplyr::filter(mystery, !coinflips),
     fill = "darkblue", alpha = 0.2, breaks = br) +
  geom_histogram(fill = "purple", breaks = br, alpha = 0.2)

## ---- book-chunk-1, eval = TRUE, echo = FALSE, fig.keep = 'high'---------
knitr::include_graphics('images/book_icon.png', dpi = 400)

## ----coinmix-------------------------------------------------------------
probHead = c(0.125, 0.25)
for (pi in c(1/8, 1/4)) {
  whCoin = sample(2, 100, replace = TRUE, prob = c(pi, 1-pi))
  K = rbinom(length(whCoin), size = 2, prob = probHead[whCoin])
  print(table(K))
}

## ----mixnorm1------------------------------------------------------------
mus = c(-0.5, 1.5)
u = sample(2, 100, replace = TRUE)
y = rnorm(length(u), mean = mus[u])
duy = tibble(u, y)
head(duy)

## ----mixnorm2------------------------------------------------------------
group_by(duy, u) %>% summarize(mean(y))

## ----mixtools, message = FALSE-------------------------------------------
library("mixtools")
y = c(rnorm(100, mean = -0.2, sd = 0.5),
      rnorm( 50, mean =  0.5, sd =   1))
gm = normalmixEM(y, k = 2, lambda = c(0.5, 0.5),
     mu = c(-0.01, 0.01), sigma = c(1, 1))
gm$lambda
gm$mu
gm$sigma
gm$loglik

## ----mosaics, echo = FALSE, results = "hide", message = FALSE, warning = FALSE----
## PROVENANCE: here's a record of how the data were created
library("mosaics")
library("mosaicsExample")
constructBins(infile = system.file(file.path("extdata", "wgEncodeSydhTfbsGm12878Stat1StdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"),
    fileFormat = "bam", outfileLoc = "../data/",
    byChr = FALSE, useChrfile = FALSE, chrfile = NULL, excludeChr = NULL,
    PET = FALSE, fragLen = 200, binSize = 200, capping = 0)
constructBins(infile = system.file(file.path("extdata", "wgEncodeSydhTfbsGm12878InputStdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"),
    fileFormat = "bam", outfileLoc = "../data/",
    byChr = FALSE, useChrfile = FALSE, chrfile = NULL, excludeChr = NULL,
    PET = FALSE, fragLen = 200, binSize = 200, capping = 0)
datafiles = c("../data/wgEncodeSydhTfbsGm12878Stat1StdAlnRep1_chr22_sorted.bam_fragL200_bin200.txt",
              "../data/wgEncodeSydhTfbsGm12878InputStdAlnRep1_chr22_sorted.bam_fragL200_bin200.txt")
binTFBS = readBins(type = c("chip","input"), fileName = datafiles)

## ----binTFBS-------------------------------------------------------------
binTFBS

## ----chipseqzeros, fig.keep = 'high', fig.cap = "The number of binding sites found in 200nt windows along chromosome 22 in a ChIP-Seq dataset."----
bincts = print(binTFBS)
ggplot(bincts, aes(x = tagCount)) +
  geom_histogram(binwidth = 1, fill = "forestgreen")

## ----ChipseqHistlogY, fig.keep = 'high', fig.cap = "As Figure \\@ref(fig:chipseqzeros), but using a logarithm base 10 scale on the $y$-axis. The fraction of zeros seems elevated compared to that of ones, twos, ...", warning = FALSE----
ggplot(bincts, aes(x = tagCount)) + scale_y_log10() +
   geom_histogram(binwidth = 1, fill = "forestgreen")

## ----chap5-r-nucleotideweights-1, fig.keep = 'high', fig.cap = "Simulation of 7,000 nucleotide mass measurements.", fig.height = 3----
masses = c(A =  331, C =  307, G =  347, T =  322)
probs  = c(A = 0.12, C = 0.38, G = 0.36, T = 0.14)
N  = 7000
sd = 3
nuclt   = sample(length(probs), N, replace = TRUE, prob = probs)
quadwts = rnorm(length(nuclt),
                mean = masses[nuclt],
                sd   = sd)
ggplot(tibble(quadwts = quadwts), aes(x = quadwts)) +
  geom_histogram(bins = 100, fill = "purple")

## ----ecdfZea, fig.keep = 'high', fig.cap = "The observed sample can be seen as a mixture of point masses at each of the values (real point masses would be bars without any width whatsoever).", fig.height = 2, fig.width = 4----
library("HistData")
ZeaMays$diff
ggplot(ZeaMays, aes(x = diff, ymax = 1/15, ymin = 0)) +
  geom_linerange(size = 1, col = "forestgreen") + ylim(0, 0.1)

## ----checkZeaMays, echo = FALSE------------------------------------------
stopifnot(nrow(ZeaMays) == 15)

## ---- samplingdist, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The value of a statistic $\\tau$ is estimated from data (the grey matrices) generated from an underlying distribution $F$. Different samples from $F$ lead to different data, and so to different values of the estimate $\\hat{\\tau}$: this is called **sampling variability**. The distribution of all the $\\hat{\\tau}$\'s is the **sampling distribution**."----
knitr::include_graphics(c('images/BootstrapPrincipleNew.png'))

## ---- bootpple, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The bootstrap simulates the sampling variability by drawing samples not from the underlying true distribution $F$ (as in Figure \\@ref(fig:samplingdist)), but from the empirical distribution function $\\hat{F}_n$."----
knitr::include_graphics(c('images/BootstrapPrinciple2New.png'))

## ----bootmedian, fig.keep = 'high', fig.cap = "The bootstrap sampling distribution of the median of the Zea Mays differences.", fig.height = 2.5, fig.width = 4----
B = 1000
meds = replicate(B, {
  i = sample(15, 15, replace = TRUE)
  median(ZeaMays$diff[i])
})
ggplot(tibble(medians = meds), aes(x = medians)) +
  geom_histogram(bins = 30, fill = "purple")

## ----bootpkg, results = "hide"-------------------------------------------
library("bootstrap")
bootstrap(ZeaMays$diff, B, mean)
bootstrap(ZeaMays$diff, B, median)


## ------------------------------------------------------------------------
c(N3 = choose(5, 3), N15 = choose(29, 15))

## ---- LaplacePortrait, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Laplace knew already that the probability density $$f_Y(y)=\\frac{1}{2\\phi}\\exp\\left(-\\frac{|y-\\theta|}{\\phi}\\right),\\qquad\\phi>0$$ has the median as its location parameter $\\theta$ and the median absolute deviation (MAD) as its scale parameter $\\phi$."----
knitr::include_graphics(c('images/LaplacePortrait.png'))

## ----Laplace1------------------------------------------------------------
w = rexp(10000, rate = 1)

## ----Laplacedistribution, fig.keep = 'high', fig.cap = "Data sampled from a Laplace distribution.", fig.width = 3, fig.height = 2.4----
mu  = 0.3
lps = rnorm(length(w), mean = mu, sd = sqrt(w))
ggplot(data.frame(lps), aes(x = lps)) +
  geom_histogram(fill = "purple", binwidth = 0.1)

## ----ALaplacedistribution, fig.keep = 'high', fig.cap = "Histogram of data generated from an asymmetric Laplace distribution -- a scale mixture of many normals whose means and variances are dependent. We write $X \\sim AL(\\theta, \\mu, \\sigma)$.", fig.width = 3, fig.height = 2.4----
mu = 0.3; sigma = 0.4; theta = -1
w  = rexp(10000, 1)
alps = rnorm(length(w), theta + mu * w, sigma * sqrt(w))
ggplot(tibble(alps), aes(x = alps)) +
  geom_histogram(fill = "purple", binwidth = 0.1)

## ---- promoterlengthsandexpression, out.width = '50%', fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Histogram of real data. On the left the lengths of the promoters shorter than 2000bp from Saccharomyces cerevisiae as studied by  @Kristiansson2009. On the right the log-ratios of microarray gene expression measurements for 20,000 genes  [@Purdom2005], . Both distributions can be modeled by asymmetric Laplace distributions."----
knitr::include_graphics(c('images/LaplaceMixturePromoterLengths.png','images/tcellhist.png'))

## ---- three-worlds, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "How to count the fish in a lake? \\copyright MC Escher."----
knitr::include_graphics(c('images/three-worlds.png'))

## ----gammahist1, fig.keep = 'high', fig.show = "hold", fig.cap = "Histograms of random samples of gamma distributions. The upper histogram shows a gamma$(2,\\frac{1}{3})$, the lower a gamma$(10,\\frac{3}{2})$ distribution. The gamma is a flexible two parameter distribution: \\href{http://en.wikipedia.org/wiki/Gamma_distribution}{see Wikipedia}.", fig.width = 3.5, fig.height = 2.4----
ggplot(tibble(x = rgamma(10000, shape = 2, rate = 1/3)),
   aes(x = x)) + geom_histogram(bins = 100, fill= "purple")
ggplot(tibble(x = rgamma(10000, shape = 10, rate = 3/2)),
   aes(x = x)) + geom_histogram(bins = 100, fill= "purple")

## ----generatepoissongamma, fig.keep = 'high', fig.cap = "Histogram of `gp`, generated via a gamma-Poisson hierachical model.", fig.width = 3.5, fig.height = 2.4----
lambda = rgamma(10000, shape = 10, rate = 3/2)
gp = rpois(length(lambda), lambda = lambda)
ggplot(tibble(x = gp), aes(x = x)) +
  geom_histogram(bins = 100, fill= "purple")

## ----goofy, fig.keep = 'high', fig.cap = "Goodness of fit plot. The **rootogram** shows the theoretical probabilities of the gamma-Poisson distribution (a.k.a.\\ negative binomial) as red dots and the square roots of the observed frequencies as the height of the rectangular bars. The bars all end close to the horizontal axis, which indicates a good fit to the negative binomial distribution.", warning= FALSE, message= FALSE, fig.width = 7, fig.height = 5----
library("vcd")
ofit = goodfit(gp, "nbinomial")
plot(ofit, xlab = "")
ofit$par

## ----setupgammapois, echo = FALSE----------------------------------------
x    = 0:95
mu   = 50
vtot = 80
v1   = vtot - mu
scale = v1/mu    # 0.6
shape = mu^2/v1  # 83.3
p1   = dgamma(x = x, scale = 0.6, shape = 80)
p2   = dpois(x = x, lambda = mu*1.2)
p3   = dnbinom(x = x, mu = mu, size = mu^2/vtot)

## ----mixtures-dgammapois, fig.keep = 'high', fig.cap = "Visualization of the hierarchical model that generates the gamma-Poisson distribution. The top panel shows the density of a gamma distribution with mean (ref:mixtures-dgammapois-1) (vertical black line) and variance (ref:mixtures-dgammapois-2). Assume that in one particular experimental replicate, the value (ref:mixtures-dgammapois-3) is realized. This is our latent variable. The observable outcome is distributed according to the Poisson distribution with that rate parameter, shown in the middle panel. In one particular experiment the outcome may be, say, (ref:mixtures-dgammapois-4), indicated by the dashed green line. Overall, if we repeat these two subsequent random process many times, the outcomes will be distributed as shown in the bottom panel -- the gamma-Poisson distribution.", fig.margin = FALSE, fig.width = 6, fig.height = 8, echo = FALSE----
library("RColorBrewer")
cols = brewer.pal(8, "Paired")
par(mfrow=c(3,1), mai=c(0.5, 0.5, 0.01, 0.01))
xlim = x[c(1, length(x))]
plot(NA, NA, xlim=xlim, ylim=c(0,0.07), type="n", ylab="", xlab="")
polygon(x, p1, col=cols[1])
abline(v=mu, col="black", lwd=3)
abline(v=mu*1.2, col=cols[2], lty=2, lwd=3)
plot(x, p2, col=cols[3], xlim=xlim, ylab="", xlab="", type="h", lwd=2)
abline(v=mu*1.2, col=cols[2], lwd=2)
abline(v=mu*1.1, col=cols[4], lty=2, lwd=3)
plot(x, p3, col=cols[4], xlim=xlim, type="h", lwd=2, ylab="", xlab="")

## ----seriesofpoisson, fig.keep = 'high', fig.show = "hold", fig.cap = "Poisson distributed measurement data, with a range of means from (ref:seriesofpoisson-1) to (ref:seriesofpoisson-2). In the left panel the $y$-axis is proportional to the data, in the right panel it is on a square-root scale. Note how the shapes of the beeswarm clouds change on the left, but less so on the right.", fig.margin = FALSE, out.width = '50%', fig.width = 3.75, fig.height = 3.75----
lambdas = seq(100, 900, by = 100)
simdat = lapply(lambdas, function(l)
    tibble(y = rpois(n = 40, lambda=l), lambda = l)
  ) %>% bind_rows
library("ggbeeswarm")
ggplot(simdat, aes(x = lambda, y = y)) +
  geom_beeswarm(alpha = 0.6, color = "purple")
ggplot(simdat, aes(x = lambda, y = sqrt(y))) +
  geom_beeswarm(alpha = 0.6, color = "purple")

## ----vstpois, echo = -c(1, 3)--------------------------------------------
.o = options(digits = 3)
summarise(group_by(simdat, lambda), sd(y), sd(2*sqrt(y)))
options(.o)

## ----seriesofnb, fig.keep = 'high', fig.cap = "gamma-Poisson distributed measurement data, for a range of $\\mu$ from (ref:seriesofnb-1) to (ref:seriesofnb-2).", fig.width = 3, fig.height = 3----
muvalues = 2^seq(0, 10, by = 1)
simgp = lapply(muvalues, function(mu) {
  u = rnbinom(n = 1e4, mu = mu, size = 4)
  tibble(mean = mean(u), sd = sd(u),
         lower = quantile(u, 0.025),
         upper = quantile(u, 0.975),
         mu = mu)
  } ) %>% bind_rows
head(as.data.frame(simgp), 2)
ggplot(simgp, aes(x = mu, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + geom_errorbar()

## ----chap5-r-pcwlin-1, fig.keep = 'high', fig.cap = "Piecewise linear function that stabilizes the variance of the data in Figure \\@ref(fig:seriesofnb).", fig.width = 3, fig.height = 3----
simgp = mutate(simgp,
  slopes = 1 / sd,
  trsf   = cumsum(slopes * mean))
ggplot(simgp, aes(x = mean, y = trsf)) +
  geom_point() + geom_line() + xlab("")

yvar = readRDS("../data/Myst.rds")$yvar
ggplot(tibble(yvar), aes(x = yvar)) + geom_histogram(binwidth=0.025)
str(yvar)

## ----EM1-----------------------------------------------------------------
pA = runif(length(yvar))
pB = 1 - pA

## ----EM2-----------------------------------------------------------------
iter = 0
loglik = -Inf
delta = +Inf
tolerance = 1e-3
miniter = 50; maxiter = 1000

## ----EM3-----------------------------------------------------------------
while((delta > tolerance) && (iter <= maxiter) || (iter < miniter)) {
  lambda = mean(pA)
  muA = weighted.mean(yvar, pA)
  muB = weighted.mean(yvar, pB)
  sdA = sqrt(weighted.mean((yvar - muA)^2, pA))
  sdB = sqrt(weighted.mean((yvar - muB)^2, pB))

  phiA = dnorm(yvar, mean = muA, sd = sdA)
  phiB = dnorm(yvar, mean = muB, sd = sdB)
  pA   = lambda * phiA
  pB   = (1 - lambda) * phiB
  ptot = pA + pB
  pA   = pA / ptot
  pB   = pB / ptot

  loglikOld = loglik
  loglik = sum(log(pA))
  delta = abs(loglikOld - loglik)
  iter = iter + 1
}
param = tibble(group = c("A","B"), mean = c(muA,muB), sd = c(sdA,sdB))
param
iter

## ----SusansPreviousCode, echo = FALSE, results = "hide"------------------
expectMaximize=function(datavec=yvar,lambda=0.5,tol=10^(-3))
{####function to decompose a normal mixture using Exp-Maximization
datavec = data.frame(datavec)
nr = nrow(datavec)
colnames(datavec)= "yvar"
ngr2=round(lambda*nr)
gri=c(rep(1,nr-ngr2),rep(2,ngr2))
##Compute initial subgroups, proportion of group 2 is lambda
#
## Initialize by making a random cut on the data into two pieces.
summari = datavec  %>%
     mutate(ug = sample(gri,nr,replace = TRUE)) %>%
     group_by(ug) %>%
     summarise(mns=mean(yvar),sds=sd(yvar))
params= tibble(m1 = summari$mns[1], m2 = summari$mns[2],
               s1 = summari$sds[1], s2 = summari$sds[2])
eps=1 ; iter=0; logli = 0;
while(eps>tol){
phi1=dnorm(datavec$yvar,params$m1,params$s1)
phi2=dnorm(datavec$yvar,params$m2,params$s2)
wts=(lambda*phi2) / ( (1-lambda)*phi1 + lambda*phi2)
params$m2=weighted.mean(datavec$yvar,wts)
params$m1=weighted.mean(datavec$yvar,1-wts)
params$s2=sqrt(weighted.mean((datavec$yvar-params$m2)^2,wts))
params$s1=sqrt(weighted.mean((datavec$yvar-params$m1)^2,1-wts))
lambda=sum(wts)/nr
logliOld=logli
logli=sum(log(wts))
eps=max(abs(logliOld- logli))
iter=iter+1
}
results=list("params"=params,"lambda"=lambda,"wts"=wts,"iter"=iter,"logli"=logli)
return(results)
}

# We can try this out on our data:
res = expectMaximize(tol=10^(-5),lambda=0.5)
res$iter
res$params
res$logli

## ----compareEM-----------------------------------------------------------
gm = mixtools::normalmixEM(yvar, k = 2)
gm$lambda
gm$mu
gm$sigma
gm$loglik

library("flexmix")
data("NPreg")

## ----flexmix2------------------------------------------------------------
m1 = flexmix(yn ~ x + I(x^2), data = NPreg, k = 2)

## ----npreg, fig.keep = 'high', fig.cap = "The points seem to come from two different generative processes, one is linear; the other quadratic."----
ggplot(NPreg, aes(x = x, y = yn)) + geom_point()

## ----parm1c1-------------------------------------------------------------
parameters(m1, component = 1)
parameters(m1, component = 2)

## ----tableNP-------------------------------------------------------------
table(NPreg$class, clusters(m1))

## ----summ1, results=FALSE------------------------------------------------
summary(m1)

## ----npregC, fig.keep = 'high', fig.cap = "Regression example using \x60flexmix\x60 with the points colored according to their estimated class. You can see that at the intersection we have an `identifiability\' problem: we cannot distinguish points that belong to the straight line from ones that belong to the parabole."----
NPreg = mutate(NPreg, gr = factor(class))
ggplot(NPreg, aes(x = x, y = yn, group = gr)) +
   geom_point(aes(colour = gr, shape = gr)) +
   scale_colour_hue(l = 40, c = 180)

