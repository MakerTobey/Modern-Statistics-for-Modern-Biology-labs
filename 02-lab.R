# Lab2 - simluations

sample_from_poisson = rpois(n = 10, lambda = 5)
sample_from_poisson[1:3]
mean(sample_from_poisson)

set.seed(19720505)
sample_from_poisson = rpois(10, 5)
sample_from_poisson[1:3]
mean(sample_from_poisson)

set.seed(1)
sample_from_poisson = rpois(10000, 5)
mean(sample_from_poisson==5)
dpois(5,5)
mean(sample_from_poisson <= 5)
ppois(5,5)

num_replicates = 50000
nexps = 5
rate = 0.1
set.seed(0xdada)
x1 = replicate(num_replicates, {
  sum(rexp(n = nexps, rate = rate))
}) # end of replicate
head(x1)

hist(x1, freq = FALSE, ylim = c(0, 0.02))
lines(sort(x1), dgamma(sort(x1), shape = nexps, scale = 1/rate), 
      col = "blue", lwd = 2)

set.seed(0xdada)
x1 = sapply(seq_len(num_replicates), function(i) {
  sum(rexp(n = nexps, rate = rate))
}
) # end of sapply
head(x1) 

set.seed(0xdada)
x1 = vapply(seq_len(num_replicates), function(i) {
  sum(rexp(n = nexps, rate = rate))
}, # end of anonymous function
FUN.VALUE = numeric(1)
) # end of vapply
head(x1) 

hist(x1, freq = FALSE, ylim = c(0, 0.02))
lines(sort(x1), dgamma(sort(x1), shape = nexps, scale = 1/rate), 
      col = "blue", lwd = 2)

#pkgs_needed = c("Biostrings", "BSgenome.Celegans.UCSC.ce2")
letsinstall = setdiff(pkgs_needed, installed.packages()) 
if (length(letsinstall) > 0) {
  BiocManager::install(letsinstall)
}
library("BSgenome.Celegans.UCSC.ce2")
library("Biostrings")
