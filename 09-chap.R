## ----initialize, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
source("../chapter-setup.R"); chaptersetup("/Users/Susan/Courses/CUBook-html/CUBook/Chap9-MVA2/MDSCA.Rnw", "9")
knitr::opts_chunk$set(dev = 'png', dpi = 100, fig.margin = TRUE, fig.show = 'hold', fig.keep = 'none')

## ---- Brighton-West-Pier-20090214-sunset, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/Brighton_West_Pier_20090214_sunset.jpg'))

## ----chap9-r-HeatEuro-1, fig.keep = 'high', fig.cap = "A heatmap of the distances between some of the cities. The function has re-arranged the order of the cities, grouping the closest ones.", height=4, width=7.5----
library("pheatmap")
load("../data/distEuroN.RData")
seteuro = as.matrix(distEuroN)[1:12, 1:12]
pheatmap(seteuro, cluster_rows = TRUE,
  treeheight_row = 0.0001, treeheight_col = 0.8,
  fontsize_col = 8, cellwidth = 13, cellheight = 13)

## ----cmdscaleDistE-------------------------------------------------------
MDSEuro = cmdscale(distEuroN, eig = TRUE)

## ----plotbareig, fig.keep = 'high', fig.cap = "Screeplot of the first 5 eigenvalues. The drop after the first two eigenvalues is very visible.", fig.height = 4, fig.width = 4----
library("tibble")
plotbar = function(res, m = 9) {
  tibble(eig = res$eig[seq_len(m)], k = seq(along = eig)) %>%
  ggplot(aes(x = k, y = eig)) +
    scale_x_discrete("k", limits = seq_len(m)) + theme_minimal() +
    geom_bar(stat="identity", width=0.5, color="orange", fill="pink")
}
plotbar(MDSEuro, m = 5)

## ----Europe, fig.width = 6,fig.height = 4, warning = FALSE---------------
plotbar(MDSEuro, m = length(MDSEuro$eig))

## ----chap9-r-Cities-1, fig.keep = 'high', fig.cap = "MDS map of European cities based on their distances.", fig.width = 5, fig.height = 4,  warning = FALSE----
MDSeur = tibble(
  PCo1 = MDSEuro$points[, 1],
  PCo2 = MDSEuro$points[, 2],
  labs = rownames(MDSEuro$points))
g = ggplot(MDSeur, aes(x = PCo1, y = PCo2, label = labs)) +
  geom_point(color = "red") + xlim(-1950, 2000) + ylim(-1150, 1150) +
  coord_fixed() + geom_text(size = 4, hjust = 0.3, vjust = -0.5)
g

## ----chap9-r-EuroLL-1, fig.keep = 'high', fig.show = "hold", fig.cap = "Left: same as Figure \\@ref(fig:chap9-r-Cities-1), but with axes flipped. Right: true latitudes and longitudes.", fig.margin = FALSE, out.width = '50%', fig.width = 5, fig.height = 4, warning = FALSE----
g %+% mutate(MDSeur, PCo1 = -PCo1, PCo2 = -PCo2)
Eurodf = readRDS("../data/Eurodf.rds")
ggplot(Eurodf, aes(x = Long,y = Lat, label = rownames(Eurodf))) +
   geom_point(color = "blue") + geom_text(hjust = 0.5, vjust = -0.5)

## ----EarthRadius, echo = FALSE-------------------------------------------
earthradius = 6371

## ----D2------------------------------------------------------------------
Eurodf = readRDS("../data/Eurodf.rds")
X = as.matrix(Eurodf)
DdotD = as.matrix(dist(X)^2)

## ----CheckMDS------------------------------------------------------------
n = nrow(X)
H = diag(rep(1,n))-(1/n) * matrix(1, nrow = n, ncol = n)
Xc = sweep(X,2,apply(X,2,mean))
Xc[1:2, ]
HX = H %*% X
HX[1:2, ]
apply(HX, 2, mean)

## ------------------------------------------------------------------------
B0 = H  %*% DdotD %*% H
B2 = HX %*% t(HX)
B2[1:3, 1:3] / B0[1:3, 1:3]
max(abs(-0.5 * B0 - B2))

## ----ekmandis------------------------------------------------------------
ekm = read.table("../data/ekman.txt", header=TRUE)
rownames(ekm) = colnames(ekm)
disekm = 1 - ekm - diag(1, ncol(ekm))
disekm[1:5, 1:5]
disekm = as.dist(disekm)

## ----ekmanMDSeig, fig.keep = 'high', fig.cap = "The screeplot shows us that the phenomenon is two dimensional, giving a clean answer to Ekman\'s question.", fig.width=3, fig.height=3----
mdsekm = cmdscale(disekm, eig = TRUE)
plotbar(mdsekm)

## ----chap9-r-ekmanMDS-1, fig.keep = 'high', fig.cap = "The layout of the scatterpoints in the first two dimensions has a horseshoe shape. The labels and colors show that the arch corresponds to the wavelengths.", fig.width = 4.6, fig.height = 4.5, message = FALSE----
dfekm = as_tibble(mdsekm$points[,1:2])%>%
  setNames(paste0("MDS", 1:2)) %>%
  mutate(
    name = rownames(ekm),
    rgb = photobiology::w_length2rgb(
          as.numeric(sub("w", "", name))))
library("ggrepel")
ggplot(dfekm, aes(x = MDS1, y = MDS2)) +
  geom_point(col = dfekm$rgb, size = 4) +
  geom_text_repel(aes(label = name)) + coord_fixed()

## ----simstress, output = FALSE, results = "hide", message = FALSE, warning = FALSE----
library("vegan")
nmds.stress = function(x, sim = 100, kmax = 4) {
  sapply(seq_len(kmax), function(k)
    replicate(sim, metaMDS(x, k = k, autotransform = FALSE)$stress))
}
stress = nmds.stress(disekm, sim = 100)
dim(stress)

## ----chap9-r-NMDSscreeplot-1, fig.keep = 'high', fig.cap = "Several replicates at each dimension were run to evaluate the stability of the <span style=\"font-variant:small-caps;\">stress</span>. We see that the <span style=\"font-variant:small-caps;\">stress</span> drops dramatically with two or more dimensions, thus indicating that a two dimensional solution is appropriate here.", fig.height = 4, fig.width = 4----
dfstr = reshape2::melt(stress, varnames = c("replicate","dimensions"))
ggplot(dfstr, aes(y = value, x = dimensions, group = dimensions)) +
  geom_boxplot()

## ----chap9-r-Shepardsplot-1, fig.keep = 'high', fig.cap = "The Shepard\'s plot compares the original distances or dissimilarities (along the horizonal axis) to the reconstructed distances, in this case for $k=2$ (vertical axis).", message = FALSE, results = 'hide', warning=FALSE----
nmdsk2 = metaMDS(disekm, k = 2, autotransform = FALSE)
stressplot(nmdsk2, pch = 20)

## ----ekmannonMDS-code----------------------------------------------------

## ----ekmannonMDS-plot, fig.width = 4.6, fig.height = 4.5-----------------
dfnmdsk2 = as_tibble(nmdsk2$points[,1:2]) %>%
           setNames(paste0("NmMDS", 1:2)) %>%
           bind_cols(select(dfekm, rgb, name))
ggplot(dfnmdsk2, aes(x = NmMDS1, y = NmMDS2)) +
  geom_point(col = dfekm$rgb, size = 4) +
  geom_text_repel(aes(label = name))

## ---- chap9-r-ekmannonMDS-1, out.width = '50%', fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Comparison of the output from the classical multidimensional scaling on the left (same as Figure \\@ref(fig:chap9-r-ekmanMDS-1)) and the nonmetric version on the right."----
knitr::include_graphics(c('images/chap9-r_ekmanMDS-plot-1.png','images/chap9-r_ekmannonMDS-plot-1.png'))



## ----vsn28Exprs----------------------------------------------------------
IBDchip = readRDS("../data/vsn28Exprd.rds")
library("ade4")
library("factoextra")
library("sva")

## ----IBDchip-------------------------------------------------------------
class(IBDchip)
dim(IBDchip)
tail(IBDchip[,1:3])
summary(IBDchip[nrow(IBDchip),])

## ----assayIBD------------------------------------------------------------
assayIBD = IBDchip[-nrow(IBDchip), ]
day = factor(IBDchip[nrow(IBDchip), ])

## ----screepc12, fig.keep = 'high', fig.cap = "The screeplot shows us that the phenomenon can be usefully represented in two dimensions.", echo = FALSE, fig.width = 4, fig.height = 4----
rankthreshPCA = function(x, threshold = 3000) {
  ranksM = apply(x, 2, rank)
  ranksM[ranksM < threshold] = threshold
  ranksM = threshold - ranksM
  dudi.pca(t(ranksM), scannf = FALSE, nf = 2)
}
pcaDay12 = rankthreshPCA(assayIBD[,day!=3])
fviz_eig(pcaDay12, bar_width = 0.6) + ggtitle("")

## ----rankthreshPCA-code, warning = FALSE---------------------------------
rankthreshPCA = function(x, threshold = 3000) {
  ranksM = apply(x, 2, rank)
  ranksM[ranksM < threshold] = threshold
  ranksM = threshold - ranksM
  dudi.pca(t(ranksM), scannf = FALSE, nf = 2)
}
pcaDay12 = rankthreshPCA(assayIBD[,day!=3])
day12 = day[day!=3]
fviz(pcaDay12, element="ind", axes=c(1,2), geom=c("point","text"),
  habillage = day12, repel = TRUE, palette = "Dark2",
  addEllipses = TRUE, ellipse.type = "convex") + ggtitle("") +
  coord_fixed()

## ----rankthreshPCA, fig.keep = 'high', fig.cap = "We have used colors to identify the different days and have kept the sample labels as well. We have also added convex hulls for each day. The group mean is identified as the point with the larger symbol (circle, triangle or square).", fig.width = 5.2, fig.height = 3.5, warning = FALSE, echo = FALSE----
rtPCA1 <- fviz(pcaDay12, element="ind", axes=c(1,2), geom=c("point","text"),
  habillage = day12, repel = TRUE, palette = "Dark2",
  addEllipses = TRUE, ellipse.type = "convex") + ggtitle("") +
  coord_fixed()
rtPCA1

## ----Threesetspca123-code, warning=FALSE---------------------------------
pcaDay123 = rankthreshPCA(assayIBD)
fviz(pcaDay123, element="ind", axes=c(1,2), geom=c("point","text"),
  habillage = day, repel=TRUE, palette = "Dark2",
  addEllipses = TRUE, ellipse.type = "convex") + ggtitle("") +
  coord_fixed()

## ----Threesetspca123, fig.keep = 'high', fig.show = "hold", fig.cap = "When comparing the three day analysis to that of the first two days, we notice the inversion of signs in the coordinates on the second axis: this has no biological relevance. The important finding is that group 3 overlaps heavily with group 1 indicating that it was the protocol change on Day 2 which created the variability.", fig.margin = FALSE, out.width = '50%',fig.width=5.2,fig.height=3.5,warning=FALSE,echo = FALSE----
rtPCA1
fviz(pcaDay123, element="ind", axes=c(1,2), geom=c("point","text"),
  habillage = day, repel=TRUE, palette = "Dark2",
  addEllipses = TRUE, ellipse.type = "convex") + ggtitle("") +
  coord_fixed()

## ----exCodeEllipse-------------------------------------------------------
fviz_pca_ind(pcaDay123, habillage = day, labelsize = 3,
  palette = "Dark2", addEllipses = TRUE, ellipse.level = 0.69)

## ----screepc123, fig.keep = 'high', fig.cap = "The eigenvalue screeplot the case of 3 groups is extremely similar to that with two groups shown in Figure \\@ref(fig:screepc12).", echo=FALSE, fig.width=4, fig.height=4----
fviz_eig(pcaDay123,bar_width=0.6) + ggtitle("")

## ----CombatIBD, fig.keep = 'high', fig.cap = "The modified data with the batch effects removed now show three batch-groups heavily overlapping and centered almost at the origin.",message = FALSE, warning = FALSE,fig.width = 5.2, fig.height = 3.5----
model0 = model.matrix(~1, day)
combatIBD = ComBat(dat = assayIBD, batch = day, mod = model0)
pcaDayBatRM = rankthreshPCA(combatIBD)
fviz(pcaDayBatRM, element = "ind", geom = c("point", "text"),
  habillage = day, repel=TRUE, palette = "Dark2", addEllipses = TRUE,
  ellipse.type = "convex", axes =c(1,2)) + coord_fixed() + ggtitle("")


## ------------------------------------------------------------------------
library("SummarizedExperiment")
sampletypes = c("IBS","CTL")
status = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
colData = DataFrame(day = day, treatment = factor(sampletypes[status]))
chipse = SummarizedExperiment(assays = list(abund = assayIBD),
                              colData = colData)


## ----SEfilter------------------------------------------------------------
chipse[,day==2]

## ----SingleCell----------------------------------------------------------
corese = readRDS("../data/normse.rds")
norm = assays(corese)$normalizedValues

## ----coldatacore---------------------------------------------------------
length(unique(colData(corese)$Batch))

## ----screeplotnorm, fig.keep = 'high', fig.cap = "Screeplot of the PCA of the normalized data.", fig.width = 4, fig.height = 4----
respca = dudi.pca(t(norm), nf = 3, scannf = FALSE)
plotbar(respca, 15)
PCS = respca$li[, 1:3]


## ----setupcolors, echo=FALSE---------------------------------------------
library("RColorBrewer")
publishedClusters = colData(corese)[, "publishedClusters"]
batch = colData(corese)$Batch
col_clus = c("transparent", "#1B9E77", "antiquewhite2", "cyan", "#E7298A",
      "#A6CEE3", "#666666", "#E6AB02", "#FFED6F", "darkorchid2",
          "#B3DE69", "#FF7F00", "#A6761D", "#1F78B4")
#col_batch
names(col_clus) = sort(unique(publishedClusters))

## ----screeplotnorm-2-----------------------------------------------------
library("rgl")
batch = colData(corese)$Batch
plot3d(PCS,aspect=sqrt(c(84,24,20)),col=col_clus[batch])
plot3d(PCS,aspect=sqrt(c(84,24,20)),
col = col_clus[as.character(publishedClusters)])

## ---- 3dplotsnorm, out.width = '50%', fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Two-dimensional screenshots of three-dimensional **[rgl](https://cran.r-project.org/web/packages/rgl/)** plots. On the left the points are colored according to batch numbers and on the right according to the original clustering. We can see that the batch effect has been effectively removed and that the cells show the original clustering."----
knitr::include_graphics(c('images/plotnormpcabatch1.png','images/plotnormpclust1.png'))

## ----HIV, echo = FALSE---------------------------------------------------
HIV <- data.frame(Patient = c('AHX112', 'AHX717', 'AHX543'), 
                           Mut1 = c(0, 1, 1),
                           Mut2 = c(0, 0, 0),
                           Mut3 = c(0, 1, 0),
                           '...' = rep(' ', 3))
        knitr::kable(HIV, format = 'html', table.attr = 'class=\"margintab marginnote\"',
        caption = 'Sample by mutation matrix.')

## ----crossHIV, echo = FALSE----------------------------------------------
crossHIV <- data.frame(Patient = c('Mut1', 'Mut2', 'Mut3'), 
                               Mut1 = c(853, 29, 10),
                               Mut2 = c(29, 853, 52),
                               Mut3 = c(10, 52, 853),
                               '...' = rep(' ', 3))
        knitr::kable(crossHIV, format = 'html', table.attr = 'class=\"margintab marginnote\"',
                     caption = 'Cross-tabulation of the HIV mutations showing two-way co-occurrences.')

## ----HIVnnrti, fig.keep = 'high', fig.cap = "The dependencies between HIV mutations is clearly a three dimensional phenomenon, the three first eigenvalues show a clear signal in the data."----
cooc = read.delim2("../data/coccurHIV.txt", header = TRUE, sep = ',')
cooc[1:4,1:11]
HIVca=dudi.coa(cooc,nf=4,scannf=FALSE)
fviz_eig(HIVca,geom="bar",bar_width=0.6)+ggtitle("")

## ---- HIV3d, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "A screenshot of the output from an interactive 3d plotting function (`plot3d`)."----
knitr::include_graphics(c('images/scatter3d-HIV.png'))

## ----CA123-HIV-----------------------------------------------------------
library("rgl")
CA1=HIVca$li[,1];CA2=HIVca$li[,2];CA3=HIVca$li[,3]
plot3d(CA1,CA2,CA3,aspect=FALSE,col="purple")

## ----HIVca, fig.keep = 'high', fig.show = "hold", fig.cap = "Two planar maps of the mutations defined with the horizontal axis corresponding to the first eigenvector of the CA and the vertical axis being the second axis on the left and the third on the right, notice the difference in heights.", fig.margin = FALSE, out.width = '50%', fig.width=4.5----
fviz_ca_row(HIVca,axes = c(1, 2),geom="text", col.row="purple",
  labelsize=3)+ggtitle("") + xlim(-0.55, 1.7) + ylim(-0.53,1.1) +
  theme_bw() +  coord_fixed()
fviz_ca_row(HIVca,axes = c(1, 3), geom="text",col.row="purple",
    labelsize=3)+ggtitle("")+ xlim(-0.55, 1.7)+ylim(-0.5,0.6) +
    theme_bw() + coord_fixed()

## ----HIVnnrtiMut13a, fig.width=4.5---------------------------------------
fviz_ca_row(HIVca,axes = c(1, 3), geom="text", col.row="purple",
      labelsize=3)+ ggtitle("")+ theme_minimal() +
      coord_fixed()

## ----chisquaredtest,warning=FALSE----------------------------------------
HairColor = HairEyeColor[,,2]
chisq.test(HairColor)

## ---- HairEye, echo=FALSE------------------------------------------------
knitr::kable(HairColor, format = 'html', table.attr = 'class=\"margintab marginnote\"',
                     caption = 'Cross tabulation of students hair and eye color')

## ----ExpectedEyes,fig.width=5,fig.height=6-------------------------------
rowsums=as.matrix(apply(HairColor,1,sum))
rowsums
colsums=as.matrix(apply(HairColor,2,sum))
t(colsums)
HCexp=rowsums%*%t(colsums)/sum(colsums)

## ----HCexp, fig.keep = 'high', fig.cap = "Here is a schematic representation of the expected table \x60HCexp\x60. We see that it has the \'rectangular\' property chracteristic of rank one matrices we saw in chapter \\@ref(Chap:Multivariate). The boxes are all white.", echo=FALSE,fig.show="hide",fig.width=4,fig.height=3.8----
mosaicplot(HCexp, shade = TRUE, las = 1, type="pearson",
           cex.axis = 0.7, main="")

## ------------------------------------------------------------------------
sum((HairColor  - HCexp)^2/HCexp)

## ----MosaicHair, fig.keep = 'high', fig.cap = "Visualization of the departure from independence. Now, the boxes are proportional in size to the actual observed counts and we no longer have a \'rectangular\' property. The departure from independence is measured in Chisquared distance for each of the boxes and colored according to whether the residuals are large and positive. Dark blue indicates a positive association, for instance between blue eyes and blonde hair, red indicates a negative association such as in the case of blond hair and brown eyes.", fig.margin = FALSE,fig.width=4.5,fig.height=4.5----
round(t(HairColor-HCexp))
library("vcd")
mosaicplot(HairColor,shade=TRUE,las=1,type="pearson",cex.axis=0.7,main="")

## ----HCscatter, fig.keep = 'high', fig.cap = "The CA plot gives a representation of a large proportion of the chisquare distance between the data and the values expected under independence. The first axis shows a contrast between black haired and blonde haired students, mirrored by the brown eye, blue eye opposition. In CA the two categories play symmetric roles and we can interpret the proximity of Blue eyes and Blond hair has meaning that there is strong co-occurence of these categories.", fig.margin = FALSE,fig.width=6,fig.height= 3----
HC=as.data.frame.matrix(HairColor)
coaHC=dudi.coa(HC,scannf=FALSE,nf=2)
round(coaHC$eig[1:3]/sum(coaHC$eig)*100)
fviz_ca_biplot(coaHC,repel=TRUE,col.col="brown", col.row="purple") +
ggtitle("") + ylim(c(-0.5,0.5))

## ----VeganCCA,fig.width=4,fig.height=3-----------------------------------
  library("vegan")
  res.ca=vegan::cca(HairColor)
  plot(res.ca,scaling=3)

## ---- ProustProxy, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/ProustProxy.png'))

## ----LakeCAr, fig.keep = 'high', fig.show = "hold", fig.cap = "The locations near the lake are ordered along an arch as shown in the top plot. In the biplot underneath, we can see which plants are most frequent at which locations by looking at the red triangles closest to the blue points. ", fig.margin = FALSE, out.width = '50%', fig.width=4,fig.height=2.5----
load("../data/lakes.RData")
lakelike[1:3,1:8]
reslake=dudi.coa(lakelike,scannf=FALSE,nf=2)
round(reslake$eig[1:8]/sum(reslake$eig),2)
fviz_ca_row(reslake,repel=TRUE)+ggtitle("")+ylim(c(-0.55,1.7))
fviz_ca_biplot(reslake,repel=TRUE)+ggtitle("")+ylim(c(-0.55,1.7))

## ---- CellTree, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The four cell populations studied here are representative of three sequential states (PS,NP,HF) and two possible final branches (4SG and 4SFG$^{-}$)."----
knitr::include_graphics(c('images/CellTree.png'))

## ----NormData, echo=FALSE,results=FALSE,eval=FALSE-----------------------
## ###For the record only, make our own data from these
## ###Provenance tracking.....
## Nor=read.csv("../data/nbt.3154-S3.csv",row.names=1)
## dim(Nor)
## blom=as.matrix(Nor)

## ----Descriptors,results='hide', echo=FALSE, eval=FALSE------------------
## desc1=unlist(strsplit(rownames(blom),"_"))
## desc=desc1[seq(1,7867,2)]
## gr4sfg=which(substr(rownames(blom),1,5)=="4SFGA")
## gr4sf=which(substr(rownames(blom),1,4)=="4SGA")
## gr1=which(substr(rownames(blom),1,2)=="PS")
## gr2=which(substr(rownames(blom),1,2)=="NP")
## gr3=which(substr(rownames(blom),1,2)=="HF")
## colscells=c("blue","green","orange","red","purple")
## colnb=rep(0,3934)
## colnb[gr1]=1
## colnb[gr2]=2
## colnb[gr3]=3
## colnb[gr4sf]=4
## colnb[gr4sfg]=5
## typesort=rep(0,3934)
## typesort[which(nchar(desc)<5 & substr(rownames(blom),3,3)=="A")]="sortA"
## typesort[which(nchar(desc)<5 & substr(rownames(blom),3,3)=="B")]="sortB"
## typesort[which(nchar(desc)>4)]="sortA"
## ftable(typesort)
## celltypes=as.factor(c("PS","NP","HF","4SG","4SGF-")[colnb])
## cellcol = colscells[colnb]
## colCells = DataFrame(celltypes=celltypes, cellcol=colscells[colnb])
## Moignard= SummarizedExperiment(assays=list(assayCells = blom),
##                    rowData=colCells)
## ## saveRDS(Moignard,file="../data/Moignard.rds")

## ----Distances-----------------------------------------------------------
Moignard = readRDS("../data/Moignard.rds")
cellt = rowData(Moignard)$celltypes
colsn = c("red", "purple", "orange", "green", "blue")
blom = assay(Moignard)
dist2n.euclid = dist(blom)
dist1n.l1     = dist(blom, "manhattan")

## ----colorforexs, echo = FALSE-------------------------------------------
## Set up for 3d work, not needed for pdf
colc = rowData(Moignard)$cellcol

## ----CMDS----------------------------------------------------------------
ce1Mds = cmdscale(dist1n.l1,     k = 20, eig = TRUE)
ce2Mds = cmdscale(dist2n.euclid, k = 20, eig = TRUE)
perc1  = round(100*sum(ce1Mds$eig[1:2])/sum(ce1Mds$eig))
perc2  = round(100*sum(ce2Mds$eig[1:2])/sum(ce2Mds$eig))

## ----CMDSplotscree, fig.keep = 'high', fig.show = "hold", fig.cap = "Screeplots from MDS on $\\ell_1$ (left) and $L_2$ (right) distances. We see that the eigenvalues are extremely similar and both point to a $2$ dimensional phenomenon.", fig.height = 2.7, fig.width = 2.2----
plotbar(ce1Mds, m = 4)
plotbar(ce2Mds, m = 4)

## ----CMDSplotL2, fig.keep = 'high', fig.show = "hold", fig.cap = "Moignard cell data colored according to the cell types (blue: PS, green: NP, yellow: HF, red: 4SG, purple: 4SFG$^-$) in the two dimensional MDS plots created.On the left (A) using $\\ell_1$ distances and on the right (B) using the L2 distances.", fig.margin = FALSE, out.width = '50%', fig.height=5,fig.width=5----
c1mds = as_tibble(ce1Mds$points[, 1:2]) %>%
            setNames(paste0("L1_PCo", 1:2))
ggplot(c1mds,aes(x = L1_PCo1,y = L1_PCo2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) +
  scale_colour_manual(values = colsn) + guides(color = FALSE)
c2mds = as_tibble(ce1Mds$points[, 1:2]) %>%
            setNames(paste0("L2_PCo", 1:2))
ggplot(c2mds,aes(x=L2_PCo1,y=L2_PCo2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) +
   scale_colour_manual(values = colsn) + guides(color = FALSE)

## ----colorlegend, echo = FALSE, fig.width = 2, fig.height = 3------------
## Hack from https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
ggpcolor = ggplot(c1mds,aes(x=L1_PCo1,y=L1_PCo2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) +
  scale_colour_manual(values=colsn, name = "cell type")
g_legend = function(a) {
     gt = ggplot_gtable(ggplot_build(a))
     leg = which(sapply(gt$grobs, function(x) x$name) == "guide-box")
     gt$grobs[[leg]]
}
grid.draw(g_legend(ggpcolor))

## ----tsnecells, fig.keep = 'high', fig.show = "hold", fig.cap = "The four cell populations studied here are representative of three sequential states (PS,NP,HF) and two possible final branches (4SG and 4SFG$^{-}$). The plot on the left was obtained by choosing 2 dimensions for t-sne at a perplexity of 30. The lower plot has obtained by choosing 3 dimensions, we can see that this third t-SNE axis represented here as the horizontal axis.", fig.margin = FALSE, out.width = '50%', message=FALSE, warning=FALSE, fig.width=5, fig.height=4----
library("Rtsne")
restsne = Rtsne(blom, dims = 2, perplexity = 30, verbose = FALSE,
                max_iter = 900)
dftsne = as_tibble(restsne$Y[, 1:2]) %>%
                            setNames(paste0("taxis", 1:2))
ggplot(dftsne,aes(x = taxis1, y = taxis2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) +
   scale_color_manual(values = colsn) + guides(color = FALSE)
restsne3 = Rtsne(blom, dims = 3, perplexity = 30, verbose = FALSE,
                 max_iter = 900)
dftsne3 = as_tibble(restsne3$Y[, 1:3]) %>%
            setNames(paste0("taxis", 1:3))
ggplot(dftsne3,aes(x = taxis3, y = taxis2, group = cellt)) +
      geom_point(aes(color = cellt), alpha = 0.6) +
      scale_colour_manual(values = colsn) + guides(color = FALSE)

## ---- tsne3d, out.width = '50%', fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Moignard cell data colored according to the cell types (blue: PS, green: NP, yellow: HF, red: 4SG, purple: 4SFG$^-$) in the three-dimensional t-SNE layouts. We can see that the purple cells (4SFG$^-$) segregate at the outer shell on the top of the point cloud."----
knitr::include_graphics(c('images/tsnemoignard3scrop.png','images/tsnemoignard3crop.png'))



## ----multitable-setup----------------------------------------------------
mb_path = "../data/metabolites.csv"
library("genefilter")
load("../data/microbe.rda")
metab = read.csv(mb_path, row.names = 1) %>% as.matrix

## ----multitable-filtering------------------------------------------------
metab   = metab[rowSums(metab == 0) <= 3, ]
microbe = prune_taxa(taxa_sums(microbe) > 4, microbe)
microbe = filter_taxa(microbe, filterfun(kOverA(3, 2)), TRUE)
metab = log(1 + metab, base = 10)
X = as.matrix(otu_table(microbe))
X = log(1 + X, base=10)

## ----RVtest--------------------------------------------------------------
colnames(metab)=colnames(X)
pca1 = dudi.pca(t(metab), scal = TRUE, scann = FALSE)
pca2 = dudi.pca(t(X), scal = TRUE, scann = FALSE)
rv1 = RV.rtest(pca1$tab, pca2$tab, 999)
rv1

## ----multitable-sparse-cca-----------------------------------------------
library("PMA")
ccaRes = CCA(t(X), t(metab), penaltyx = 0.15, penaltyz = 0.15)
ccaRes

## ----multitablepluginpca, echo=FALSE-------------------------------------
combined = cbind(t(X[ccaRes$u != 0, ]),
                 t(metab[ccaRes$v != 0, ]))
pcaRes = dudi.pca(combined, scannf = FALSE, nf = 3)
# annotation
genotype     = substr(rownames(pcaRes$li), 1, 2)
sampleType  = substr(rownames(pcaRes$l1), 3, 4)
featureType = grepl("\\.", colnames(combined))
featureType = ifelse(featureType, "Metabolite", "OTU")
sampleInfo  = data.frame(pcaRes$li, genotype, diet=sampleType)
featureInfo = data.frame(pcaRes$c1,
                          feature = substr(colnames(combined), 1, 6))

## ----multitableinterpretpca, fig.keep = 'high', fig.cap = "A PCA triplot produced from the CCA selected features from muliple data types (metabolites and OTUs).", fig.margin = FALSE, fig.height = 5, fig.width = 6, echo=FALSE----
ggplot() +  geom_point(data = sampleInfo,
  aes(x = Axis1, y = Axis2, col = diet, shape = genotype), size = 3) +
  geom_label_repel(data = featureInfo,
  aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = featureType),
      size = 2, segment.size = 0.3,
      label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = featureInfo,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = featureType),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed()+
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pcaRes$eig[1] / sum(pcaRes$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pcaRes$eig[2] / sum(pcaRes$eig), 2)),
       fill = "Feature Type", col = "Sample Type")


## ----ccpna-correspondence-analysis---------------------------------------
ps1=readRDS("../data/ps1.rds")
ps1p=filter_taxa(ps1, function(x) sum(x) > 0, TRUE)
psCCpnA = ordinate(ps1p, "CCA",
                 formula = ps1p ~ ageBin + family_relationship)

## ----ccpna-join-data, echo=FALSE-----------------------------------------
library("dplyr")
tax = data.frame(tax_table(ps1p),stringsAsFactors = FALSE)
tax$seq = rownames(tax)
mainOrders = c("Clostridiales", "Bacteroidales",
               "Lactobacillales", "Coriobacteriales")
tax$Order[!(tax$Order %in% mainOrders)] = "Other"
tax$Order = factor(tax$Order, levels = c(mainOrders, "Other"))
tax$otu_id = seq_len(ncol(otu_table(ps1p)))
scoresCCpnA = vegan::scores(psCCpnA)
sites = data.frame(scoresCCpnA$sites)
sites$SampleID = rownames(sites)
sites = left_join(sites, sample_data(ps1p))
species = data.frame(scoresCCpnA$species)
species$otu_id = seq_along(colnames(otu_table(ps1p)))
species = left_join(species, tax)

## ----ccpnaplotage, fig.keep = 'high', fig.cap = "The mouse and taxa scores generated by CCpnA. The sites (mice samples) are triangles; species are circles, respectively. The separate panels indicate different age groups.", fig.margin = FALSE, fig.height = 4, fig.width = 7----
evalProp = 100 * psCCpnA$CCA$eig[1:2] / sum(psCCpnA$CA$eig)
ggplot() +
 geom_point(data = sites,aes(x =CCA2, y =CCA1),shape =2,alpha=0.5) +
 geom_point(data = species,aes(x =CCA2,y =CCA1,col = Order),size=1)+
 geom_text_repel(data = species %>% filter(CCA2 < -2),
                   aes(x = CCA2, y = CCA1, label = otu_id),
                   size = 2, segment.size = 0.1) +
 facet_grid(. ~ ageBin) +
 guides(col = guide_legend(override.aes = list(size = 2))) +
 labs(x = sprintf("Axis2 [%s%% variance]", round(evalProp[2])),
        y = sprintf("Axis1 [%s%% variance]", round(evalProp[1]))) +
 scale_color_brewer(palette = "Set1") + theme(legend.position="bottom")

## ----ccpnaplotlitter, fig.keep = 'high', fig.cap = "The analogue to Figure \\@ref(fig:ccpnaplotage), faceting by litter membership rather than age bin.", fig.margin = FALSE, fig.height = 4, fig.width = 7, echo=FALSE----
ggplot() +
  geom_point(data = sites, aes(x = CCA2, y = CCA1), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA2, y = CCA1, col = Order), size = 1) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                  aes(x = CCA2, y = CCA1, label = otu_id),
                  size = 2, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 2))) +
  labs(x = sprintf("Axis2 [%s%% variance]", round(evalProp[2])),
       y = sprintf("Axis1 [%s%% variance]", round(evalProp[1]))) +
  scale_color_brewer(palette = "Set1") + theme(legend.position="bottom")

ibd.pres = ifelse(assayIBD[, 1:28] > 8.633, 1, 0)

## ----chap9-r-Threesetscoa, fig.keep = 'high', fig.show = "hold", fig.cap = "Correspondence analysis on binary data.", warning = FALSE, message = FALSE, fig.hold=TRUE----
IBDca = dudi.coa(ibd.pres, scannf = FALSE, nf = 4)
fviz_eig(IBDca, geom = "bar", bar_width = 0.7) +
    ylab("Percentage of chisquare") + ggtitle("")
fviz(IBDca, element = "col", axes =c(1, 2), geom = "point",
     habillage = day, palette = "Dark2", addEllipses = TRUE, color = day,
     ellipse.type = "convex", alpha = 1, col.row.sup =  "blue",
     select = list(name = NULL, cos2 = NULL, contrib = NULL),
     repel = TRUE)

d1 <- t(data.frame(
                              quiet = c(2770, 2150, 2140, 875, 1220, 821, 2510),
                              angry = c(2970, 1530, 1740, 752, 1040, 710, 1730),
                              clever = c(1650, 1270, 1320, 495, 693, 416, 1420),
                              depressed = c(1480, 957, 983, 147, 330, 102, 1270),
                              happy = c(19300, 8310, 8730, 1920, 4220, 2610, 9150),
                              lively = c(1840, 1250, 1350, 659, 621, 488, 1480),
                              perplexed = c(110,  71,  80,  19,  23,  15, 109),
                              virtuous = c(179,  80, 102,  20,  25,  17, 165)))
               colnames(d1) <- c('black','blue','green','grey','orange','purple','white')
knitr::kable(d1, caption = 'Contingency table of co-occurring terms from search engine results.')

## ----chap9-r-ColorBiplot-1, fig.keep = 'high', fig.cap = "Correspondence Analysis allows for a symmetrical graphical representation of two categorical variables, in this case colors and emotions for a contingency table of co-occurrences such as Table \\@ref(tab:colors).", echo = FALSE, fig.height = 5, fig.width = 4.5----
colorsentiment = read.csv("../data/colorsentiment.csv")
colsent = xtabs(colorsentiment[,3] ~ colorsentiment[,2] + colorsentiment[,1])
coldf = data.frame(unclass(colsent))
coldf = round(coldf / 1000)
# xtable::xtable(round(coldf),display=rep("d", 8))
colorfactor = names(coldf)
veganout = vegan::cca(coldf)
colorfactor[c(4,7)] = c("darkgrey", "grey")
ordiplot(veganout, scaling = 3, type = "none", xlim =c(-1.2, 0.75), ylim =c(-0.7, 1))
text(veganout, "sites", pch = 21, col = "red", bg = "yellow", scaling = 3)
text(veganout, "species", pch = 21, col = colorfactor, bg = "black", cex=1.2, scaling = 3)

## ---- PlatoTableImage, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/PlatoTableImage.png'))

platof = read.table("../data/platof.txt", header = TRUE)
platof[1:4, ]
resPlato = dudi.coa(platof, scannf = FALSE, nf = 2)
fviz_ca_biplot(resPlato, axes=c(2, 1)) + ggtitle("")
fviz_eig(resPlato, geom = "bar", width = 0.6) + ggtitle("")

## ----PercentageInertia---------------------------------------------------
names(resPlato)
sum(resPlato$eig)
percentageInertia=round(100*cumsum(resPlato$eig)/sum(resPlato$eig))
percentageInertia
percentageInertia[2]

## ----MakeMatrices--------------------------------------------------------
load("../data/lakes.RData")
lakelike[ 1:3, 1:8]
lakelikeh[1:3, 1:8]
e_coa  = dudi.coa(lakelike,  scannf = FALSE, nf = 2)
e_pca  = dudi.pca(lakelike,  scannf = FALSE, nf = 2)
eh_coa = dudi.coa(lakelikeh, scannf = FALSE, nf = 2)
eh_pca = dudi.pca(lakelikeh, scannf = FALSE, nf = 2)

## ----adescatter, results = "hide"----------------------------------------
scatter(e_pca)
scatter(e_coa)
s.label(e_pca$li)
s.label(e_coa$li)

s.label(eh_pca$co)
s.label(eh_pca$li)
s.label(eh_coa$li)
s.label(eh_coa$co)

## ----RawData-------------------------------------------------------------
moignard_raw = as.matrix(read.csv("../data/nbt.3154-S3-raw.csv",
                                  row.names = 1))
dist2r.euclid = dist(moignard_raw)
dist1r.l1 = dist(moignard_raw, "manhattan")
cells.cmds  = cmdscale(dist1r.l1,     k = 20, eig = TRUE)
cells2.cmds = cmdscale(dist2r.euclid, k = 20, eig = TRUE)
sum(cells.cmds$eig[1:2]) / sum(cells.cmds$eig)
sum(cells2.cmds$eig[1:2]) / sum(cells2.cmds$eig)

## ----KernelD, results = "hide"-------------------------------------------
library("kernlab")
# Kernelized Distances
laplacedot1=laplacedot(sigma=1/3934)
rbfdot1=rbfdot(sigma =(1/3934)^2 )
Klaplace_cellsn <- kernelMatrix(laplacedot1, blom)
KGauss.cellsn <- kernelMatrix(rbfdot1, blom)
Klaplace_rawcells <- kernelMatrix(laplacedot1, moignard_raw)
KGauss.rawcells <- kernelMatrix(rbfdot1, moignard_raw)

# Use kernelized distances to protect against outliers
# and allows discovery of non linear  components
dist1kr=1-Klaplace_rawcells
dist2kr=1-KGauss.rawcells
dist1kn=1-Klaplace_cellsn
dist2kn=1-KGauss.cellsn

cells.kcmds = cmdscale(dist1kr,k=20,eig=TRUE)
cells2.kcmds =cmdscale(dist2kr,k=20,eig=TRUE)
kperc1=round(100*sum(cells.kcmds$eig[1:4])/
       sum(cells.kcmds$eig[which(cells.kcmds$eig>0)]))
kperc2=round(100*sum(cells2.kcmds$eig[1:4])/
       sum(cells2.kcmds$eig[which(cells2.kcmds$eig>0)]))
cellsn.kcmds=cmdscale(dist1kn,k=20,eig=TRUE)
cellsn2.kcmds=cmdscale(dist2kn,k=20,eig=TRUE)

## ----KernelMDSplots, warning = FALSE, results = "hide"-------------------
colc = rowData(Moignard)$cellcol
library("scatterplot3d")
scatterplot3d(cellsn2.kcmds$points[,1:3],color=colc,pch=20,size=0.8,
        xlab="Axis k1",ylab="Axis k2",zlab="Axis k3",angle=15)
scatterplot3d(cellsn2.kcmds$points[,1:3],color=colc,pch=20,size=0.8,
        xlab="Axis k1",ylab="Axis k2",zlab="Axis k3",angle=-70)

## ----CodeB, results = "hide"---------------------------------------------
library("rgl")
plot3d(cellsn2.kcmds$points[, 1:3], col = colc, size = 3,
       xlab = "Axis1", ylab = "Axis2", zlab = "Axis3")
plot3d(cellsn2.kcmds$points[, c(1,2,4)], col = colc, size = 3,
       xlab = "Axis1", ylab = "Axis2", zlab = "Axis4")
# Using an L1 distance instead.
plot3d(cellsn.kcmds$points[, 1:3], col = colc, size = 3,
       xlab = "Axis1", ylab = "Axis2", zlab = "Axis3")
plot3d(cellsn.kcmds$points[, c(1,2,4)], col = colc, size = 3,
       xlab = "Axis1", ylab = "Axis2", zlab = "Axis4")

## ----lpc3d---------------------------------------------------------------
library("LPCM")
library("diffusionMap")
dmap1 = diffuse(dist1n.l1, neigen = 10)
combs = combn(4, 3)
lpcplots = apply(combs, 2, function(j) lpc(dmap1$X[, j], scale = FALSE))

## ----PLOT3DLPC-----------------------------------------------------------
library("rgl")
for (i in seq_along(lpcplots))
  plot(lpcplots[[i]], type = "l", lwd = 3,
  xlab = paste("Axis", combs[1, i]),
  ylab = paste("Axis", combs[2, i]),
  zlab = paste("Axis", combs[3, i]))

## ---- diffusionmap3, out.width = '50%', eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Diffusion map projection for Axes 1, 3 and 4. The lower figure shows the smoothed path followed by the cells in their development."----
knitr::include_graphics(c('images/TripleArm.png','images/SmoothLineP134h7.png'))

## ----scatter3smoothedline------------------------------------------------
outlpce134=lpc(dmap1$X[,c(1,3,4)],scale=FALSE,h=0.5)
plot3d(dmap1$X[,c(1,3,4)], col=colc,xlab="Axis1",ylab="Axis3",
        zlab="Axis4",pch=20)
plot3d(outlpce134$LPC,type="l",lwd=7,add=TRUE)
outlpce134=lpc(dmap1$X[,c(1,3,4)],scale=FALSE,h=0.7)
plot3d(outlpce134$LPC,type="l",lwd=7,xlab="Axis1",ylab="Axis3",
zlab="Axis4")
plot3d(dmap1$X[,c(1,3,4)], col=colc,xlab="",ylab="",
zlab="",add=TRUE)

knitr::include_graphics(c('images/dmap134.png'))

## ----Diffuse, messages = FALSE, warnings = FALSE-------------------------
library("diffusionMap")
dmap2 = diffuse(dist2n.euclid, neigen = 11)
dmap1 = diffuse(dist1n.l1, neigen = 11)
plot(dmap2)

## ----scp3d---------------------------------------------------------------
library("scatterplot3d")
scp3d = function(axestop = 1:3, dmapRes = dmap1, color = colc,
           anglea = 20, pch = 20)
scatterplot3d(dmapRes$X[, axestop], color = colc,
    xlab = paste("Axis",axestop[1]), ylab = paste("Axis", axestop[2]),
    zlab = paste("Axis",axestop[3]), pch = pch, angle = anglea)

## ----dmap3dplots---------------------------------------------------------
scp3d()
scp3d(anglea=310)
scp3d(anglea=210)
scp3d(anglea=150)

## ----CodeD---------------------------------------------------------------
# interactive plot
library("rgl")
plot3d(dmap1$X[,1:3], col=colc, size=3)
plot3d(dmap1$X[,2:4], col=colc, size=3)

