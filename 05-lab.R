pkgs_needed = c("dbscan","tidyverse","GGally", "pheatmap",
                "flowCore","flowViz","flowPeaks", "ggcyto") # packages for cytometry data analysis)
letsinstall = setdiff(pkgs_needed, installed.packages()) 
if (length(letsinstall) > 0) {
  BiocManager::install(letsinstall)
}
library("tidyverse")
library("flowCore")
library("flowViz")
library("flowPeaks")
library("ggcyto")
download.file(url = "http://web.stanford.edu/class/bios221/data/Bendall_2011.fcs",
              destfile = "Bendall_2011.fcs", mode = "wb")
fcsB = read.FCS("Bendall_2011.fcs")
slotNames(fcsB)

colnames(fcsB) #41 variables

exprs(fcsB) #i don't know how many events

markersB = read_csv(url("http://web.stanford.edu/class/bios221/data/Bendall_2011_markers.csv"))

mt = match(markersB$isotope, colnames(fcsB))
stopifnot(!any(is.na(mt)))
colnames(fcsB)[mt] = markersB$marker

flowPlot(fcsB, plotParameters = c("Cell_length", "DNA191"), logy=TRUE)

# `densityplotvis()` is from `densityplot` package
densityplot(~`CD3all`, fcsB)

# variance stabilising transformation
asinhT = arcsinhTransform(a = 0.1, b = 1)
cols_to_transform <- setdiff(colnames(fcsB), c("Time", "Cell_length", "absoluteEventNumber"))
trans1 = transformList(cols_to_transform, asinhT)
fcsBT = transform(fcsB, trans1)
densityplot(~`CD3all`, fcsBT)

kf = kmeansFilter("CD3all"=c("Pop1","Pop2"), filterId="myKmFilter")
fres = filter(fcsBT, kf)
summary(fres)

fcsBT1 = split(fcsBT, fres, population="Pop1")
fcsBT2 = split(fcsBT, fres, population="Pop2")

dat = data.frame(exprs(fcsBT)[ , c("CD56", "CD3all")])
fp = flowPeaks(dat)
plot(fp)

?flowPeaks

library("ggcyto")
# Untransformed data
ggcyto(fcsB,aes(x = CD4)) + geom_histogram(bins = 60) 
# Transformed data
ggcyto(fcsBT, aes(x=CD4)) + geom_histogram(bins=90) 
# ggcyto automatic plotting
autoplot(fcsBT, "CD4")

ggcyto(fcsBT, aes(x = CD4, y = CD8)) + geom_density2d(colour="black") 

ggcyto(fcsBT, aes(x = CD4, y = CD8)) + geom_hex(bins = 50) 
# ggcyto automatic plotting
autoplot(fcsBT, "CD4", "CD8", bins = 50)


## density based clustering 
library("dbscan")
# Select a small subset of 5 protein markers
mc5 = exprs(fcsBT)[, c("CD4", "CD8", "CD20", "CD3all", "CD56")]
res5 = dbscan::dbscan(mc5, eps = .65, minPts = 30)
mc5df = data.frame(mc5)
mc5df$cluster = as.factor(res5$cluster)

res5 #7 clusters, 4616 cells in cl3

ggplot(mc5df, aes(x = CD4, y = CD8, colour = cluster)) + geom_density2d()

ggplot(mc5df,aes(x = CD3all, y = CD20, colour = cluster))+ geom_density2d()


# validating
simul4 = lapply(c(0,8), function(x){
  lapply(c(0,8), function(y){
    data.frame(x = rnorm(100, x, 2),
               y = rnorm(100, y, 2), 
               class = paste(x, y, sep = "")
    )
  }) %>% do.call(rbind,.)
}) %>% do.call(rbind,.)

ggplot(simul4, aes(x = x, y = y)) +
  geom_point(aes(color = class), size = 2)

# Compute the kmeans within group wss for k=1 to 12
wss = rep(0,8)
# for a single cluster the WSS statistic is just sum of squares of centered data
wss[1] = sum(apply(scale(simul4[,1:2], scale = F), 2, function(x){ x^2 }))
# for k = 2, 3, ... we perform kmeans clustering and compute the associated WSS statistic
for (k in 2:8) {
  km4 <- kmeans(simul4[,1:2],k)
  wss[k] =  sum(km4$withinss)
}
# Now, we are ready to plot the computed statistic:
ggplot(data.frame(k = 1:length(wss), wss = wss)) +
  geom_point(aes(x = k, y = wss), color = "blue", size= 3) +
  xlab('k') + ylab('WSS(k)')


# hierarchical clustering
load(url("http://web.stanford.edu/class/bios221/data/Morder.RData"))
dim(Morder)

D <- dist(t(Morder))
gene_clust <- hclust(d = D)
plot(gene_clust)

# ?t() the matrix needs to be transposed to compare bases in columns (variables)
?hclust
gene_clust2 <- hclust(d = D, method = "ward.D2")
plot(gene_clust2)

# we don't transpose the matrix now (samples are rows)
D_samples <- dist(Morder)
sample_clust <- hclust(d = D_samples)
plot(sample_clust)

?abline()
plot(sample_clust)
abline(h =12, untf = FALSE)
# 3

library(pheatmap)
pheatmap(Morder, fontsize_col=12, fontsize_row = 15) 
