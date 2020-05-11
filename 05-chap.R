## ----initialize, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
source("../chapter-setup.R"); chaptersetup("/Users/Susan/Courses/CUBook-html/CUBook/Chap6-Clustering/Cluster.Rnw", "5")
knitr::opts_chunk$set(dev = 'png', dpi = 100, fig.margin = TRUE, fig.show = 'hold', fig.keep = 'none')

## ---- starlings-copyrightfree, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/starlings_copyrightfree.jpg'))

## ---- SnowMapSmallest, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "John Snow\'s map of cholera cases: small barcharts at each house indicate a clustering of diagnosed cases."----
knitr::include_graphics(c('images/SnowMapSmallest.png'))

## ---- book-chunk-1, eval = TRUE, echo = FALSE, fig.keep = 'high'---------
knitr::include_graphics('images/book_icon.png', dpi = 400)

## ---- bombings, fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Here is a map of the location of the bombs that were dropped on London on September, 7th, 1940 as depicted by the website of the British National Archives [http://bombsight.org](http://bombsight.org)."----
knitr::include_graphics(c('images/RedBombsLondon.png'))

## ---- cancerHC, fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The breast cancer samples (shown from The Cancer Genome Atlas (TCGA) and the Molecular Taxonomy of Breast Cancer International Consortium (METABRIC)) can be split into groups using their miRNA expression  [@Aure2017], . The authors show in the lower plots that the survival times in different % groups were different. Thus these clusters were biologically and clinically relevant. The promise of such analyses is that the groups can be used to provide more specific, optimized treatments."----
knitr::include_graphics(c('images/BreastCancerSubType_Biomed.png'))

## ---- ClusteringA, fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "We decompose the choices made in a clustering algorithm according to the steps taken: starting from an observations-by-features rectangular table $X$, we choose an observations-to-observations distance measure and compute the distance matrix, here schematized by the triangle. The distances are used to construct the clusters. On the left, we schematize agglomerative methods, that build a hierarchical clustering tree; on the right, partitioning methods that separate the data into subsets. Both types of methods require a choice to be made: the number $k$ of clusters. For partitionning approaches such as $k$-means this choice has to be made at the outset; for hierarchical clustering this can be deferred to the end of the analysis."----
knitr::include_graphics(c('images/ClusteringA.png'))


## ---- fourdistances, out.width = '25%', eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Equal-distance contour plots according to four different distances: points on any one curve are all the same distance from the center point."----
knitr::include_graphics(c('images/FourDistances_a.png','images/FourDistances_b.png','images/FourDistances_c.png','images/FourDistances_d.png'))

## ----Mahalanobis, fig.keep = 'high', fig.cap = "An example for the use of Mahalanobis distances to measure the distance of a new data point (red) from two cluster centers.", echo=FALSE,eval=TRUE,messages=FALSE,warnings=FALSE,fig.width=4.7,fig.height=4.5----
library(MASS)
set.seed(101)
n <- 60000
S1=matrix(c(1,.72,.72,1), ncol=2)
S2=matrix(c(1.5,-0.6,-0.6,1.5),ncol=2)
mu1=c(.5,2.5)
mu2=c(6.5,4)
X1 = mvrnorm(n, mu=c(.5,2.5), Sigma=matrix(c(1,.72,.72,1), ncol=2))
X2 = mvrnorm(n,mu=c(6.5,4), Sigma=matrix(c(1.5,-0.6,-0.6,1.5),ncol=2))
# A color palette from blue to yellow to red
library(RColorBrewer)
k = 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))
plot(X1, xlim=c(-4,12),ylim=c(-2,9), xlab="Orange", ylab="Red", pch='.', cex=1)
points(X2, pch='.', cex=1)
 # Draw the colored contour lines
## compute 2D kernel density, see MASS book, pp. 130-131
z1 = kde2d(X1[,1], X1[,2], n=50)
z2 = kde2d(X2[,1], X2[,2], n=50)
contour(z1, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
contour(z2, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
points(3.2,2,pch=20,cex=2.2,col="red")
#text(3.1,1.4,"(3.2,2)",col="red")
lines(c(3.2,6.5),c(2,4),col="red",lwd=3)
lines(c(3.2,.5),c(2,2.5),col="red",lwd=3)
#text(-1.5,-1.5,expression(paste(D[1]^{2},"(p,",m[1],")=19.7",sep="")),cex=1.1)
#text(6.3,-1.5,expression(paste(D[2]^{2},"(p,",m[2],")=16",sep="")),cex=1.1)

## ---- DistanceTriangle, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The lower triangle of distances can be computed by any of a hundred different functions in various **R** packages (`vegdist` in **[vegan](https://cran.r-project.org/web/packages/vegan/)**, `daisy` in **[cluster](https://cran.r-project.org/web/packages/cluster/)**, `genetic_distance` in **[gstudio](https://cran.r-project.org/web/packages/gstudio/)**, `dist.dna` in **[ape](https://cran.r-project.org/web/packages/ape/)**, `Dist` in **[amap](https://cran.r-project.org/web/packages/amap/)**, `distance` in **[ecodist](https://cran.r-project.org/web/packages/ecodist/)**, `dist.multiPhylo` in **[distory](https://cran.r-project.org/web/packages/distory/)**, `shortestPath` in **[gdistance](https://cran.r-project.org/web/packages/gdistance/)**, % `dudi.dist` and `dist.genet` in **[ade4](https://cran.r-project.org/web/packages/ade4/)**). "----
knitr::include_graphics(c('images/DistanceTriangle.png'))

## ----xz------------------------------------------------------------------
mx  = c(0, 0, 0, 1, 1, 1)
my  = c(1, 0, 1, 1, 0, 1)
mz  = c(1, 1, 1, 0, 1, 1)
mat = rbind(mx, my, mz)
dist(mat)
dist(mat, method = "binary")

## ----morder--------------------------------------------------------------
load("../data/Morder.RData")
sqrt(sum((Morder[1, ] - Morder[2, ])^2))
as.matrix(dist(Morder))[2, 1]

## ----HIVmut--------------------------------------------------------------
mut = read.csv("../data/HIVmutations.csv")
mut[1:3, 10:16]

## ----answer-hiv, echo = 2:6----------------------------------------------
.o = options(digits = 3)
library("vegan")
mutJ = vegdist(mut, "jaccard")
mutC = sqrt(2 * (1 - cor(t(mut))))
mutJ
as.dist(mutC)
options(.o)

## ---- xkcd-birds, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "An example of computing the cophenetic distance (xkcd)."----
knitr::include_graphics(c('images/birds_and_dinosaurs.png'))


## ----clust-kmeansastep1, fig.keep = 'high', fig.show = "hold", fig.cap = "An example run of the $k$-means algorithm. The initial, randomly chosen centers (black circles) and groups (colors) are shown in the top panel. The group memberships are assigned based on their distance to centers. At each iteration, the group centers are redefined, and the points reassigned to the cluster centers.", fig.margin = FALSE, out.width = '33%', echo = FALSE, fig.width = 3.2, fig.height = 3, warning = FALSE----
set.seed(248811)
Xmat = matrix(runif(100), ncol = 2)
nk = 3
cents = Xmat[sample(nrow(Xmat), nk, replace = FALSE), ]
# default distance: Euclidean
dist1 = function(vec){dist(rbind(vec, cents[1,]))}
dist2 = function(vec){dist(rbind(vec, cents[2,]))}
dist3 = function(vec){dist(rbind(vec, cents[3,]))}
dists123 = cbind(apply(Xmat, 1, dist1),
                 apply(Xmat, 1, dist2),
                 apply(Xmat, 1, dist3))
clust0 = apply(dists123, 1, which.min)
out1 = kmeans(Xmat, cents, iter.max=1)
out2 = kmeans(Xmat, cents, iter.max=3)
data0 = data.frame(x = Xmat[,1],
                   y = Xmat[,2],
                   cluster = as.factor(clust0))
data1 = data.frame(x = Xmat[,1],
                   y = Xmat[,2],
                   cluster = as.factor(out1$cluster))
data2 = data.frame(x = Xmat[,1],
                   y = Xmat[,2],
                   cluster = as.factor(out2$cluster))
library("ggplot2")
.mp = function(v, cdg) {
  ggplot(data = v, aes(x = x, y = y)) +
  geom_point(aes(col = cluster, shape = cluster), size = 5) + xlab("") + ylab("") +
  geom_point(data = cdg, fill = "black", size = 7, shape = 1) +
  scale_shape_discrete(solid = TRUE, guide = FALSE) + guides(col = FALSE) + coord_fixed()
}

# centers of clusters:
cdg = data.frame(x = cents[,1],y = cents[,2])
.mp(data0, cdg)

cents = out1$centers
cdg1 = data.frame(x=cents[,1],y=cents[,2])
.mp(data1, cdg1)

cents = out2$centers
cdg2 = data.frame(x=cents[,1],y=cents[,2])
.mp(data2, cdg2)

## ----chap5-r-quiltclust-1, fig.keep = 'high', fig.cap = "Comparison of clustering results (rows), for different numbers of included genes and for varying numbers of clusters, $k$. Each column of the heatmap corresponds to a cell, and the colors represent the cluster assignments.", fig.height = 6, fig.width = 6, echo = 3:9----
.oldMar = par("mar")
par(mar = c(1.1, 6, 4.1, 1.1))
library("clusterExperiment")
data("fluidigm", package = "scRNAseq")
se = fluidigm[, fluidigm$Coverage_Type == "High"]
assays(se) = list(normalized_counts = 
   round(limma::normalizeQuantiles(assay(se))))
ce = clusterMany(se, clusterFunction = "pam", ks = 5:10, run = TRUE,
  isCount = TRUE, reduceMethod = "var", nFilterDims = c(60, 100, 150))
clusterLabels(ce) = sub("FilterDims", "", clusterLabels(ce))
plotClusters(ce, whichClusters = "workflow", axisLine = -1)
par(mar = .oldMar)

## ----quiltclust2-check, echo = FALSE-------------------------------------
stopifnot(length(clusterLabels(ce)) == 18)

## ---- book-chunk-2, eval = TRUE, echo = FALSE, fig.keep = 'high'---------
knitr::include_graphics('images/book_icon.png', dpi = 400)

## ----flowCore------------------------------------------------------------
library("flowCore")
library("flowViz")
fcsB = read.FCS("../data/Bendall_2011.fcs")
slotNames(fcsB)

## ----RenameCols, message=FALSE, warning=FALSE----------------------------
markersB = readr::read_csv("../data/Bendall_2011_markers.csv")
mt = match(markersB$isotope, colnames(fcsB))
stopifnot(!any(is.na(mt)))
colnames(fcsB)[mt] = markersB$marker

## ----ObviousClusters, fig.keep = 'high', fig.cap = "Cell measurements that show clear clustering in two dimensions."----
flowPlot(fcsB, plotParameters = colnames(fcsB)[2:3], logy = TRUE)

## ----v1v3, warning = FALSE-----------------------------------------------
v1 = seq(0, 1, length.out = 100)
plot(log(v1), asinh(v1), type = 'l')
plot(v1, asinh(v1), type = 'l')
v3 = seq(30, 3000, length = 100)
plot(log(v3), asinh(v3), type= 'l')

## ----plotTransformations, fig.keep = 'high', fig.show = "hold", fig.cap = "The left plot shows the histogram of the CD3all variable: the cells are clustered around 0 with a few large values. On the right, we see that after an asinh transformation, the cells cluster and fall into two groups or types.", fig.margin = FALSE, out.width = '50%', fig.width = 3, fig.height = 3----
asinhtrsf = arcsinhTransform(a = 0.1, b = 1)
fcsBT = transform(fcsB,
  transformList(colnames(fcsB)[-c(1, 2, 41)], asinhtrsf))
densityplot( ~`CD3all`, fcsB)
densityplot( ~`CD3all`, fcsBT)

## ----fcsBT---------------------------------------------------------------
kf = kmeansFilter("CD3all" = c("Pop1","Pop2"), filterId="myKmFilter")
fres = flowCore::filter(fcsBT, kf)
summary(fres)
fcsBT1 = flowCore::split(fcsBT, fres, population = "Pop1")
fcsBT2 = flowCore::split(fcsBT, fres, population = "Pop2")

## ----chap5-r-flowCD3CD56-1, fig.keep = 'high', fig.cap = "After transformation these cells were clustered using `kmeans`.", message = FALSE, results = "hide"----
library("flowPeaks")
fp = flowPeaks(Biobase::exprs(fcsBT)[, c("CD3all", "CD56")])
plot(fp)

## ----chap5-r-groupcontourCD3CD56-1, fig.keep = 'high', fig.cap = "Like Figure \\@ref(fig:chap5-r-flowCD3CD56-1), using contours.", fig.width = 5, fig.height = 5----
flowPlot(fcsBT, plotParameters = c("CD3all", "CD56"), logy = FALSE)
contour(fcsBT[, c(40, 19)], add = TRUE)

## ----SamSPECTRAL, eval=FALSE, echo = FALSE-------------------------------
## ## This was too slow, didn't work for flow cytometry data well.
## library("SamSPECTRAL")
## mc2 = fcsBT@exprs[,c(40,33)]

## ----ggcytoCD4CD8, fig.width = 7,fig.height=7----------------------------
library("ggcyto")
library("labeling")
ggcd4cd8=ggcyto(fcsB,aes(x=CD4,y=CD8))
ggcd4=ggcyto(fcsB,aes(x=CD4))
ggcd8=ggcyto(fcsB,aes(x=CD8))
p1=ggcd4+geom_histogram(bins=60)
p1b=ggcd8+geom_histogram(bins=60)
asinhT = arcsinhTransform(a=0,b=1)
transl = transformList(colnames(fcsB)[-c(1,2,41)], asinhT)
fcsBT = transform(fcsB, transl)
p1t=ggcyto(fcsBT,aes(x=CD4))+geom_histogram(bins=90)
p2t=ggcyto(fcsBT,aes(x=CD4,y=CD8))+geom_density2d(colour="black")
p3t=ggcyto(fcsBT,aes(x=CD45RA,y=CD20))+geom_density2d(colour="black")

## ----dbscanfcs5, fig.keep = 'high', fig.show = "hold", fig.cap = "These two plots show the results of clustering with `dbscan` using five markers. Here we only show the projections of the data into the CD4-CD8 and C3all-CD20 planes.", fig.margin = FALSE, out.width = '50%'----
library("dbscan")
mc5 = Biobase::exprs(fcsBT)[, c(15,16,19,40,33)]
res5 = dbscan::dbscan(mc5, eps = 0.65, minPts = 30)
mc5df = data.frame(mc5, cluster = as.factor(res5$cluster))
table(mc5df$cluster)
ggplot(mc5df, aes(x=CD4,    y=CD8,  col=cluster))+geom_density2d()
ggplot(mc5df, aes(x=CD3all, y=CD20, col=cluster))+geom_density2d()

## ----mc6-----------------------------------------------------------------
mc6 = Biobase::exprs(fcsBT)[, c(15, 16, 19, 33, 25, 40)]
res = dbscan::dbscan(mc6, eps = 0.65, minPts = 20)
mc6df = data.frame(mc6, cluster = as.factor(res$cluster))
table(mc6df$cluster)

## ----mc7-----------------------------------------------------------------
mc7 = Biobase::exprs(fcsBT)[, c(11, 15, 16, 19, 25, 33, 40)]
res = dbscan::dbscan(mc7, eps = 0.95, minPts = 20)
mc7df = data.frame(mc7, cluster = as.factor(res$cluster))
table(mc7df$cluster)


## ---- LinnaeusClass, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "A snippet of Linn\\ae{}us\' taxonomy that clusters organisms according to feature similarities."----
knitr::include_graphics(c('images/LinnaeusClass-01.png'))

## ---- Cluster-dendrogramorder, fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Three representations of the **same** hierarchical clustering tree."----
knitr::include_graphics(c('images/SameTree-01.png'))

## ---- single, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "In the single linkage method, the distance between groups $C_1$ and $C_2$ is defined as the distance between the closest two points from the groups."----
knitr::include_graphics(c('images/ClusterStepChoiceSingle1b.png'))

## ---- complete, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "In the complete linkage method, the distance between groups $C_1$ and $C_2$ is defined as the maximum distance between pairs of points from the two groups."----
knitr::include_graphics(c('images/ClusterStepChoiceComplete1b.png'))

## ---- betweenwithin, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The Ward method maximizes the between group sum of squares (red edges), while minimizing the sums of squares within groups (black edges)."----
knitr::include_graphics(c('images/BetweenWithinb.png'))

## ---- mobile, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Hierarchical clustering output has similar properties to a mobile: the branches can rotate freely from their suspension points."----
knitr::include_graphics(c('images/CalderHand.png'))

## ----threetreeshapes, eval = FALSE, echo = FALSE-------------------------
## ####This is for documentation purposes, had to paste
## ####trees together par(mfrow)) not working for pheatmap
## library("pheatmap")
## load("../data/d14.RData")
## pheatmap(d14,clustering_distance_rows=d14,treeheight_col =200,
## cellwidth=20,cellheight=10,lwd=5,treeheight_row=0,clustering_method = "single",
## labels_col=1:11,main="single")
## pheatmap(d14,clustering_distance_rows=d14,treeheight_col =200,cellwidth=20,
## cellheight=10,lwd=5,treeheight_row=0,clustering_method = "average",
## labels_col=1:11,main="average")
## pheatmap(d14,clustering_distance_rows=d14,treeheight_col =200,cellwidth=20,
## cellheight=10,lwd=5,treeheight_row=0,clustering_method = "complete",labels_col=1:11,
## main="complete")

## ----fortherecord-2, echo=FALSE, eval=FALSE------------------------------
## ####### For the eecord: this is what was done to the data
## ####### Melanoma/Tcell Data: Peter Lee, Susan Holmes, PNAS.
## load("../data/Msig3transp.RData")
## celltypes=factor(substr(rownames(Msig3transp),7,9))
## status=factor(substr(rownames(Msig3transp),1,3))
## Msig2=as.matrix(Msig3transp)
## rownames(Msig2)=substr(rownames(Msig2),1,9)
## hm1=heatmap(as.matrix(dist(Msig2)))
## Morder=Msig2[hm1$rowInd,]
## save(Morder,file="../data/Morder.RData")
## write.table(Morder,"../data/Morder.txt")

## ----hclust30Tcells, eval=FALSE, echo=FALSE------------------------------
## library("gplots")
## library("pheatmap")
## library("RColorBrewer")
## load("../data/Morder.RData")
## celltypes=factor(substr(rownames(Morder),7,9))
## status=factor(substr(rownames(Morder),1,3))
## ##Just the Euclidean distance
## pheatmap(as.matrix(dist(Morder)),cluster_rows=FALSE,
##         cluster_cols=FALSE,cellwidth=10,cellheight=10)
## ###Manhattan
## pheatmap(as.matrix(dist(Morder,"manhattan")),cluster_rows=FALSE,
##         cluster_cols=FALSE,cellwidth=10,cellheight=10)

## ----corT,eval=FALSE,echo=FALSE------------------------------------------
## pheatmap(corT,clustering_distance_rows=distcor,
##        annotation_row=samplesdata[,c("celltypes","status")],
##        show_rownames = FALSE, show_colnames = FALSE)
## pheatmap(corT,clustering_distance_rows=distcor,treeheight_row =150,
##        annotation_row=samplesdata[,c("celltypes","status")],
##        show_rownames = FALSE, show_colnames = FALSE)
## pheatmap(corT,clustering_distance_rows=distcor,treeheight_row =150,
##        annotation_row=samplesdata[,c("celltypes","status")],
##        treeheight_col =150,
##        show_rownames = FALSE, show_colnames = FALSE)
## pheatmap(corT,clustering_distance_rows=distcor,treeheight_row =150,
##           annotation_col=samplesdata[,c("celltypes","status")],
##        annotation_row=samplesdata[,c("celltypes","status")],
##        treeheight_col =150,
##        show_rownames = FALSE, show_colnames = FALSE)

## ---- treeshapes, out.width = '33%', fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Three hierarchical clustering plots made with different agglomeration choices. Note the comb-like structure for single linkage in the left. The average and complete linkage trees only differ by the lengths of their inner branches."----
knitr::include_graphics(c('images/single14heatmap.png','images/average14heatmap.png','images/complete14heatmap.png'))

## ---- apeclust14, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "This tree can be drawn in many different ways. The ordering of the leaves as it is appears here is $(8,11,9,10,7,5,6,1,4,2,3)$."----
knitr::include_graphics(c('images/apeclust14.png'))

## ----Fake4---------------------------------------------------------------
library("dplyr")
simdat = lapply(c(0, 8), function(mx) {
  lapply(c(0,8), function(my) {
    tibble(x = rnorm(100, mean = mx, sd = 2),
           y = rnorm(100, mean = my, sd = 2),
           class = paste(mx, my, sep = ":"))
   }) %>% bind_rows
}) %>% bind_rows
simdat
simdatxy = simdat[, c("x", "y")] # without class label

## ----chap5-r-simdat-1, fig.keep = 'high', fig.cap = "The `simdat` data colored by the class labels. Here, we know the labels since we generated the data -- usually we do not know them.", fig.width = 3, fig.height = 2.3----
ggplot(simdat, aes(x = x, y = y, col = class)) + geom_point() +
  coord_fixed()

## ----WSS, fig.keep = 'high', fig.cap = "The barchart of the WSS statistic as a function of $k$ shows that the last substantial jump is just before $k=4$. This indicates that the best choice for these data is $k=4$.", fig.width = 3, fig.height = 3----
wss = tibble(k = 1:8, value = NA_real_)
wss$value[1] = sum(scale(simdatxy, scale = FALSE)^2)
for (i in 2:nrow(wss)) {
  km  = kmeans(simdatxy, centers = wss$k[i])
  wss$value[i] = sum(km$withinss)
}
ggplot(wss, aes(x = k, y = value)) + geom_col()

## ----chap5-r-CHIndex-1, fig.keep = 'high', fig.cap = "The Calinski-Harabasz index, i.\\,e., the ratio of the between and within group variances for different choices of $k$, computed on the `simdat` data.", fig.width = 3, fig.height = 3----
library("fpc")
library("cluster")
CH = tibble(
  k = 2:8,
  value = sapply(k, function(i) {
    p = pam(simdatxy, i)
    calinhara(simdatxy, p$cluster)
  })
)
ggplot(CH, aes(x = k, y = value)) + geom_line() + geom_point() +
  ylab("CH index")

## ---- roulette-chunk-3, eval = TRUE, echo = FALSE, fig.keep = 'high'-----
knitr::include_graphics('images/roulette.png', dpi = 600)

## ----chap5-r-GapStat-1, fig.keep = 'high', fig.cap = "The gap statistic, see Question \\@ref(ques:cluster-ques-gapstat).", messages = FALSE, warnings = FALSE----
library("cluster")
library("ggplot2")
pamfun = function(x, k)
  list(cluster = pam(x, k, cluster.only = TRUE))

gss = clusGap(simdatxy, FUN = pamfun, K.max = 8, B = 50,
              verbose = FALSE)
plot_gap = function(x) {
  gstab = data.frame(x$Tab, k = seq_len(nrow(x$Tab)))
  ggplot(gstab, aes(k, gap)) + geom_line() +
    geom_errorbar(aes(ymax = gap + SE.sim,
                      ymin = gap - SE.sim), width=0.1) +
    geom_point(size = 3, col=  "red")
}
plot_gap(gss)

## ----Hiiragi-------------------------------------------------------------
library("Hiiragi2013")
data("x")

## ----submatH-------------------------------------------------------------
selFeats = order(rowVars(Biobase::exprs(x)), decreasing = TRUE)[1:50]
embmat = t(Biobase::exprs(x)[selFeats, ])
embgap = clusGap(embmat, FUN = pamfun, K.max = 24, verbose = FALSE)
k1 = maxSE(embgap$Tab[, "gap"], embgap$Tab[, "SE.sim"])
k2 = maxSE(embgap$Tab[, "gap"], embgap$Tab[, "SE.sim"],
           method = "Tibs2001SEmax")
c(k1, k2)

## ----checkassertion, echo = FALSE----------------------------------------
stopifnot("firstSEmax" == eval(formals(maxSE)$method)[1])

## ----chap5-r-gapHiiragi-1, fig.keep = 'high', fig.cap = "The gap statistic for the **[Hiiragi2013](https://bioconductor.org/packages/Hiiragi2013/)** data.", fig.width = 4.5, fig.height = 4.5----
plot(embgap, main = "")
cl = pamfun(embmat, k = k1)$cluster
table(pData(x)[names(cl), "sampleGroup"], cl)

## ---- BootstrapClusterNew, out.width = '50%', fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Different samples from the same distribution $F$ lead to different clusterings. On the left, we see the true sampling variability. The bootstrap simulates this sampling variability by drawing subsamples using the empirical distribution function $\\hat{F}_n$ as shown on the right."----
knitr::include_graphics(c('images/BootstrapClusterNew.png','images/BootstrapCluster2New.png'))

## ----BootstrapCluster----------------------------------------------------
clusterResampling = function(x, ngenes = 50, k = 2, B = 250,
                             prob = 0.67) {
  mat = Biobase::exprs(x)
  ce = cl_ensemble(list = lapply(seq_len(B), function(b) {
    selSamps = sample(ncol(mat), size = round(prob * ncol(mat)),
                      replace = FALSE)
    submat = mat[, selSamps, drop = FALSE]
    sel = order(rowVars(submat), decreasing = TRUE)[seq_len(ngenes)]
    submat = submat[sel,, drop = FALSE]
    pamres = pam(t(submat), k = k)
    pred = cl_predict(pamres, t(mat[sel, ]), "memberships")
    as.cl_partition(pred)
  }))
  cons = cl_consensus(ce)
  ag = sapply(ce, cl_agreement, y = cons)
  list(agreements = ag, consensus = cons)
}

## ----ce1-----------------------------------------------------------------
iswt = (x$genotype == "WT")
cr1 = clusterResampling(x[, x$Embryonic.day == "E3.25" & iswt])
cr2 = clusterResampling(x[, x$Embryonic.day == "E3.5"  & iswt])

## ----figClue1, fig.keep = 'high', fig.show = "hold", fig.cap = "Cluster stability analysis with E3.25 and E3.5 samples. Left: beeswarm plots of the cluster agreements with the consensus, for the `B` clusterings; $1$ indicates perfect agreement, lower values indicate lower degrees of agreement. Right: membership probabilities of the consensus clustering. For E3.25, the probabilities are diffuse, indicating that the individual clusterings often disagree, whereas for E3.5, the distribution is bimodal, with only one ambiguous sample.", fig.margin = FALSE, out.width = '50%', fig.width = 2.4, fig.height = 4----
ag1 = tibble(agreements = cr1$agreements, day = "E3.25")
ag2 = tibble(agreements = cr2$agreements, day = "E3.5")
ggplot(bind_rows(ag1, ag2), aes(x = day, y = agreements)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(cex = 1.5, col = "#0000ff40")
mem1 = tibble(y = sort(cl_membership(cr1$consensus)[, 1]),
              x = seq(along = y), day = "E3.25")
mem2 = tibble(y = sort(cl_membership(cr2$consensus)[, 1]),
              x = seq(along = y), day = "E3.5")
ggplot(bind_rows(mem1, mem2), aes(x = x, y = y, col = day)) +
  geom_point() + facet_grid(~ day, scales = "free_x")


## ----seqradius, fig.keep = 'high', fig.cap = "Although both groups have noise distributions with the same variances, the apparent radii of the groups are very different. The $10^{5}$ instances in `seq2` have many more opportunities for errors than what we see in `seq1`, of which there are only $10^{3}$. Thus we see that frequencies are important in clustering the data.", fig.width=3.1, fig.height=2.6----
library("mixtools")
library("ggplot2")
seq1 = rmvnorm(n = 1e3, mu = -c(1, 1), sigma = 0.5 * diag(c(1, 1)))
seq2 = rmvnorm(n = 1e5, mu =  c(1, 1), sigma = 0.5 * diag(c(1, 1)))
twogr = data.frame(
  rbind(seq1, seq2),
  seq = factor(c(rep(1, nrow(seq1)),
                 rep(2, nrow(seq2))))
)
colnames(twogr)[1:2] = c("x", "y")
ggplot(twogr, aes(x = x, y = y, colour = seq,fill = seq)) +
  geom_hex(alpha = 0.5, bins = 50) + coord_fixed()

## ---- book-chunk-3, eval = TRUE, echo = FALSE, fig.keep = 'high'---------
knitr::include_graphics('images/book_icon.png', dpi = 400)

## ----seqradius2----------------------------------------------------------
n    = 2000
len  = 200
perr = 0.001
seqs = matrix(runif(n * len) >= perr, nrow = n, ncol = len)

## ----dists---------------------------------------------------------------
dists = as.matrix(dist(seqs, method = "manhattan"))

## ----diameter, fig.keep = 'high', fig.cap = "The diameter of a set of sequences as a function of the number of sequences.", fig.width = 3.1, fig.height = 2.6----
library("tibble")
dfseqs = tibble(
  k = 10 ^ seq(log10(2), log10(n), length.out = 20),
  diameter = vapply(k, function(i) {
    s = sample(n, i)
    max(dists[s, s])
    }, numeric(1)))
ggplot(dfseqs, aes(x = k, y = diameter)) + geom_point()+geom_smooth()

## ----seqradiusex, fig.keep = 'high', fig.cap = "`distplot` for the `simseq10K` data."----
simseq10K = replicate(1e5, sum(rpois(200, 0.0005)))
mean(simseq10K)
vcd::distplot(simseq10K, "poisson")

## ----rates, results = FALSE, message = FALSE, warning = FALSE------------
derepFs = readRDS(file="../data/derepFs.rds")
derepRs = readRDS(file="../data/derepRs.rds")
library("dada2")
ddF = dada(derepFs, err = NULL, selfConsist = TRUE)
ddR = dada(derepRs, err = NULL, selfConsist = TRUE)

## ----rerrorprofile1, fig.keep = 'high', fig.cap = "Forward transition error rates as provided by `plotErrors(ddF)`. This shows the frequencies of each type of nucleotide transition as a function of quality.", fig.margin = FALSE, fig.width = 6, warning = FALSE----
plotErrors(ddF)

## ----dada, results = FALSE, message = FALSE, eval = TRUE-----------------
dadaFs = dada(derepFs, err=ddF[[1]]$err_out, pool = TRUE)
dadaRs = dada(derepRs, err=ddR[[1]]$err_out, pool = TRUE)

## ----merge, eval=TRUE----------------------------------------------------
mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs)

## ----seqtab--------------------------------------------------------------
seqtab.all = makeSequenceTable(mergers[!grepl("Mock",names(mergers))])

## ----checkdada, echo = FALSE---------------------------------------------
dadada = unique(vapply(dadaRs, class, character(1)))
stopifnot(is.list(dadaRs), identical("dada", dadada))

## ----answer-mergers, echo = FALSE----------------------------------------
length(dadaRs)
length(dadaFs)
class(dadaRs)
names(dadaRs)
mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs)
class(mergers)
length(mergers)

## ----chimeras------------------------------------------------------------
seqtab = removeBimeraDenovo(seqtab.all)

library("cluster")
pam4 = pam(simdatxy, 4)
sil = silhouette(pam4, 4)
plot(sil, col=c("red","green","blue","purple"), main="Silhouette")

## ----ggplotdistheatmap,include=FALSE,warning=FALSE,eval=FALSE------------
## jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "BuPu"))
## paletteSize <- 256
## jBuPuPalette <- jBuPuFun(paletteSize)
## dd=as.matrix(dist.dune)
## prdune <- data.frame(sample = colnames(dd),
##                         probe = rownames(dd),
##                         dist = dd)
## ggplot(prdune, aes(x = probe, y = sample, fill = dist)) +
##   geom_tile() +
##   scale_fill_gradient2(low = jBuPuPalette[1],
##                        mid = jBuPuPalette[paletteSize/2],
##                        high = jBuPuPalette[paletteSize],
##                        midpoint = (max(prdune$dist) + min(prdune$dist)) / 2,
##                        name = "Distance")

## ----heatmapDist,include=FALSE,eval=FALSE--------------------------------
## ##solution should use pheatmap
## library("graphics")
## library("gplots")
## #rc <- rainbow(40, start=0, end=0.99)
## rc= heat.colors(21, alpha = 1)
## dr=round(as.matrix(dist.dune),1)
## heatmap.2(1-as.matrix(dist.dune),symm = TRUE, margins = c(3,3),Rowv = NA, Colv = NA,col=rc,
## distfun=function(c) as.dist(c), cellnote=dr,key=FALSE)

library("kernlab")
data("spirals")
clusts = kmeans(spirals,2)$cluster
plot(spirals, col = c("blue", "red")[clusts])
data("spirals", package = "kernlab")
res.dbscan = dbscan::dbscan(spirals, eps = 0.16, minPts = 3)
plot(spirals,col=c("blue","red","forestgreen")[res.dbscan$cluster])

## ----checkdbscan, echo = FALSE-------------------------------------------
stopifnot(identical(range(res.dbscan$cluster), c(1L, 3L)))

## ----specc, eval = FALSE, echo = FALSE-----------------------------------
## sc = specc(spirals, centers=2)
## plot(spirals, col=sc)

base_dir = "../data"
miseq_path = file.path(base_dir, "MiSeq_SOP")
filt_path = file.path(miseq_path, "filtered")
fnFs = sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs = sort(list.files(miseq_path, pattern="_R2_001.fastq"))
sampleNames = sapply(strsplit(fnFs, "_"), `[`, 1)
if (!file_test("-d", filt_path)) dir.create(filt_path)
filtFs = file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs = file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
fnFs = file.path(miseq_path, fnFs)
fnRs = file.path(miseq_path, fnRs)
print(length(fnFs))

## ----chap5-r-profile-1, fig.keep = 'high', fig.show = "hold", fig.cap = "Quality scores. The lines show positional summary statistics: green is the mean, orange is the median, and the dashed orange lines are the 25th and 75th quantiles.", fig.margin = FALSE, out.width = '50%', fig.width = 5.8, fig.height=5----
plotQualityProfile(fnFs[1:2]) + ggtitle("Forward")
plotQualityProfile(fnRs[1:2]) + ggtitle("Reverse")

## ----profilerev----------------------------------------------------------
  ii = sample(length(fnFs), 4)
  plotQualityProfile(fnFs[ii]) + ggtitle("Forward")
  plotQualityProfile(fnRs[ii]) + ggtitle("Reverse")

## ----filter--------------------------------------------------------------
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
        maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,  trimLeft=10,
        compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

## ----derep---------------------------------------------------------------
derepFs = derepFastq(filtFs, verbose = FALSE)
derepRs = derepFastq(filtRs, verbose = FALSE)
names(derepFs) = sampleNames
names(derepRs) = sampleNames

