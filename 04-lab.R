pkgs_needed = c("tidyverse", "ggplot2", "mixtools", "HistData",
                "bootstrap", "ggbeeswarm", "pasilla", "matrixStats", "DESeq2")
BiocManager::install(setdiff(pkgs_needed, installed.packages()))

library("tidyverse")
library("ggplot2")
library("mixtools")
library("HistData")
library("bootstrap")
library("ggbeeswarm")
library("pasilla")
library("matrixStats")
library("DESeq2")


coinflips = (runif(10000) > 0.5)
table(coinflips)

sds   = c(0.5, 0.5)
means = c(  1,   3)

fairmix = rnorm(length(coinflips),
                mean = ifelse(coinflips, means[1], means[2]),
                sd   = ifelse(coinflips, sds[1],   sds[2]))
fairdf = data.frame(fairmix)
ggplot(fairdf, aes(x = fairmix)) +
  geom_histogram(fill = "skyblue", binwidth = 0.2)


B = 1e+6
coinflips = (runif(B) > 0.5)
fairmix = rnorm(length(coinflips),
                mean = ifelse(coinflips, means[1], means[2]),
                sd   = ifelse(coinflips, sds[1],   sds[2]))
fairdf = data.frame(fairmix)
ggplot(fairdf, aes(x = fairmix)) +
  geom_histogram(fill = "sienna", bins = 500)


faird = data.frame(fairdf,coinflips)
setone = dplyr::filter(faird, coinflips)
fnpoints = faird[sample(nrow(faird), 1000), ]
ggplot(setone, aes(x= fairmix)) +
  geom_histogram(aes(y = ..density..), fill = "thistle1", binwidth = 0.01) +
  stat_function(fun = dnorm, data = fnpoints,
                args = list(mean = means[1], sd = sds[1]), color = "springgreen")


xs = seq(-1, 5, length = 1000)
dens2 = 0.5 * dnorm(xs, mean = means[1], sd = sds[1])+
  0.5 * dnorm(xs, mean = means[2], sd = sds[2])
fairtheory = data.frame(xs,dens2)
ggplot(fairtheory) + aes(x = xs, y = dens2) +
  geom_line(color = "springgreen") + ylab("mixture density")



set.seed(1233341)
sds = rep(sqrt(0.5), 2)
means = c(1, 2)
coinflips = (runif(1000) > 0.5)

output = rnorm(length(coinflips),
               mean = ifelse(coinflips, means[1], means[2]),
               sd   = ifelse(coinflips, sds[1],   sds[2]))

ht = hist(output, nclass = 30, plot = FALSE)
maxcount = max(ht$count)
mysterydata = data.frame(x = output, group = ifelse(coinflips, "A", "B"))
xmin = min(output); xmax = max(output)

ggplot(mysterydata, aes(x = x)) +
  geom_histogram(fill = "orange", bins = 30)+
  xlim(c(xmin, xmax)) + ylim(c(0, maxcount))

head(mysterydata, 5)

ggplot(mysterydata, aes(x = x, group= group)) +
  geom_histogram(data = dplyr::filter(mysterydata, group == "A"),
                 fill = "red",  alpha = 0.3, bins = 30) +
  geom_histogram(data = dplyr::filter(mysterydata, group == "B"),
                 fill = "darkblue", alpha = 0.3, bins = 30) +
  xlim(c(xmin,xmax))+ ylim(c(0, maxcount))

ggplot(mysterydata, aes(x = x)) +
  geom_histogram(fill = "orange", alpha = 0.4, bins=30) +
  geom_histogram(data = dplyr::filter(mysterydata, group == "A"),
                 fill = "red",  alpha = 0.4 , bins=30) +
  geom_histogram(data = dplyr::filter(mysterydata, group == "B"),
                 fill = "darkblue", alpha = 0.4, bins=30) +
  xlim(c(xmin,xmax))+ ylim(c(0, maxcount)) 


## demixing distributions

set.seed(198435)
mus = c(-0.5,1.5)
u = sample(2, 20, replace = TRUE)
y = rnorm(length(u), mean = mus[u])
duy = data.frame(u, y)
group_by(duy, u)[1:6, ]

group_by(duy, u) %>% summarize(mean(y))

library("mixtools")
n     = c( 100,  50)
mu    = c(-0.2, 0.5)
sigma = c( 0.5,   1)

y = c(rnorm(n[1], mu[1], sd = sigma[1]), rnorm(n[2], mu[2], sigma[2]))
gm = normalmixEM(y, k = 2, lambda = c(0.5, 0.5),
                 mu = c(-0.02, 0.02), sigma = c(1, 1))


## bootstrap

library("HistData")
ZeaMays$diff

ggplot(data.frame(ZeaMays, y = 1/15),
       aes(x = diff, ymax = 1/15, ymin = 0)) +
  geom_linerange(size=1, col= "forestgreen") +
  ylim(0, 0.25)

set.seed(1)
B = 10000
difference = ZeaMays$diff
samplesIdx = replicate(B, sample(15, 15, replace = TRUE))
samplingDist = apply(samplesIdx, 2, function(x) median(difference[x]) )

ggplot(data.frame(samplingDist), aes(x = samplingDist)) +
  geom_histogram(bins = 30, fill = "skyblue")

ggplot(data.frame(samplingDist), aes(x = samplingDist)) +
  geom_histogram(bins = 30, fill = "skyblue") +
  geom_vline(xintercept = median(difference))

quantile(samplingDist, probs = .975)


## infinite mixtures

ggplot(data.frame(x = rgamma(10000, shape = 2, rate = 1/3)),
       aes(x = x)) + geom_histogram(bins = 100, fill= "purple")

ggplot(data.frame(x = rgamma(10000, shape = 10, rate = 3/2)),
       aes(x = x)) + geom_histogram(bins = 100, fill= "purple")

lambda = rgamma(100000, shape = 10, rate = 3/2)
gp = rpois(length(lambda), lambda = lambda)

ggplot(data.frame(x = gp), aes(x = x)) +
  geom_histogram(bins = 100, fill= "purple")

## Variance Stabilization

lambdas = seq(100, 500, by = 100)

simdat = lapply(lambdas, function(l)
  data.frame(y = rpois(n = 100, lambda = l), lambda = l)) %>% 
  bind_rows()

library("ggbeeswarm")
ggplot(simdat, aes( x=lambda, y=y)) +
  geom_beeswarm(alpha = 0.6, color="purple")

ggplot(simdat, aes( x=lambda, y=sqrt(y))) +
  geom_beeswarm(alpha = 0.6, color="purple")

summarise(group_by(simdat, lambda), sd(y), sd(2*sqrt(y)))


muvalues = 2^seq(0, 10, by = 1)
simgp = lapply(muvalues, function(mu) {
  u = rnbinom(n = 1e4, mu = mu, size = 4)
  data.frame(mean = mean(u), sd = sd(u),
             lower = quantile(u, 0.025),
             upper = quantile(u, 0.975),
             mu = mu)
} ) %>% bind_rows
head(as.data.frame(simgp), 2)

ggplot(simgp, aes(x = mu, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + geom_errorbar()


## Intro to RNA-Seq

fn = system.file("extdata", "pasilla_gene_counts.tsv",
                 package = "pasilla", mustWork = TRUE)
counts = as.matrix(read.csv(fn, sep = "\t", row.names = "gene_id"))
dim(counts)

head(counts)

sf = estimateSizeFactorsForMatrix(counts)
ncounts  = counts / matrix(sf, byrow = TRUE, 
                           ncol = ncol(counts), nrow = nrow(counts))
head(ncounts)

# only untreated samples
uncounts = ncounts[, grep("^untreated", colnames(ncounts)), drop = FALSE]
head(uncounts)

ggplot(data.frame( mean = rowMeans(uncounts), var  = rowVars( uncounts)),
       aes(x = log(mean), y = log(var))) +
  geom_hex() + 
  geom_abline(slope = c(1, 2), intercept = c(0, -4),
              color = c("forestgreen", "red")) +
  coord_fixed() + theme(legend.position = "none") 
