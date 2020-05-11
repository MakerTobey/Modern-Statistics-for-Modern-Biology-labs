## ----initialize, cache = FALSE, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
source("../chapter-setup.R"); chaptersetup("/Users/Susan/Courses/CUBook-html/CUBook/Chap3-Graphics/rgraphics.Rnw", "3")
knitr::opts_chunk$set(dev = 'png', dpi = 100, fig.margin = TRUE, fig.show = 'hold', fig.keep = 'none')

## ----introplot, echo = FALSE, warning = FALSE, fig.width = 4.2, fig.height = 4.2, dev = setdiff(knitr::opts_chunk$get()$dev, "postscript")----
library("xkcd")
library("showtext")
library("sysfonts")
library("tibble")

introplotdata = tibble(
  y = c(seq(-8, 1, length=25)^2, rep(1, 5), seq(1, 5,length=25)^2)^2,
  x = seq(1, 55, length.out = length(y)))

dataman = tibble(
  x = 30,
  y = 400,
  scale = 100,
  ratioxy = 0.1,
  angleofspine =  -pi/2 ,
  anglerighthumerus = -pi/6,
  anglelefthumerus  = pi * 7/6,
  anglerightradius = 0,
  angleleftradius = 0,
  angleleftleg  = 19*pi/12,
  anglerightleg = 17*pi/12,
  angleofneck   = 1.4*pi)

mapping = do.call(aes_string, colnames(dataman) %>% (function(x) setNames(as.list(x), x)))

ggplot(introplotdata) + geom_line(aes(x = x, y = y), size = 2) +
   xkcdaxis(c(0, 50), c(0, 1000)) + xlab("Time to make plot in minutes") +
   ylab("Time to understand plot in minutes") + xkcdman(mapping, dataman) +
   theme(axis.title.x = element_text(margin = margin(15, 0, 0, 0)))

## ---- rgraphics-plotter, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The ZUSE Plotter Z64 (presented in 1961). Source: [https://en.wikipedia.org/wiki/Plotter](https://en.wikipedia.org/wiki/Plotter)."----
knitr::include_graphics(c('images/Automatisches_Zeichengeraet_ZUSE_Z64_ubt.JPG'))

## ----rgraphics-basicplotting1, fig.keep = 'high', fig.cap = "Plot of concentration vs.\\ density for an ELISA assay of DNase."----
head(DNase)
plot(DNase$conc, DNase$density)

## ----rgraphics-basicplotting2, fig.keep = 'high', fig.cap = "Same data as in Figure \\@ref(fig:rgraphics-basicplotting1) but with better axis labels and a different plot symbol."----
plot(DNase$conc, DNase$density,
  ylab = attr(DNase, "labels")$y,
  xlab = paste(attr(DNase, "labels")$x, attr(DNase, "units")$x),
  pch = 3,
  col = "blue")

## ----rgraphics-basicplotting3, fig.keep = 'high', fig.show = "hold", fig.cap = "Histogram of the density from the ELISA assay, and boxplots of these values stratified by the assay run. The boxes are ordered along the axis in lexicographical order because the runs were stored as text strings. We could use R\'s type conversion functions to achieve numerical ordering.", fig.width=3, fig.height=3----
hist(DNase$density, breaks=25, main = "")
boxplot(density ~ Run, data = DNase)

## ---- rgraphics-cells, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Single-section immunofluorescence image of the E3.5 mouse blastocyst stained for Serpinh1, a marker of primitive endoderm (blue), Gata6 (red) and Nanog (green)."----
knitr::include_graphics(c('images/Yusukecells-2.jpg'))


## ----loadHiiragi---------------------------------------------------------
library("Hiiragi2013")
data("x")
dim(Biobase::exprs(x))

## ----xpData--------------------------------------------------------------
head(pData(x), n = 2)

## ----groupSize-----------------------------------------------------------
library("dplyr")
groups = group_by(pData(x), sampleGroup) %>%
  summarise(n = n(), color = unique(sampleColour))
groups


## ----explainpipe, eval = FALSE-------------------------------------------
## f(x) %>% g(y) %>% h
## h(g(f(x), y))

## ----rgraphics-figredobasicplottingwithggplot, fig.keep = 'high', fig.cap = "Our first **[ggplot2](https://cran.r-project.org/web/packages/ggplot2/)** figure, similar to the base graphics Figure \\@ref(fig:rgraphics-basicplotting1).", fig.width = 3.5, fig.height = 3----
library("ggplot2")
ggplot(DNase, aes(x = conc, y = density)) + geom_point()

## ----rgraphics-qplot1, fig.keep = 'high', fig.cap = "A barplot, produced with the `ggplot` function from the table of group sizes in the mouse single cell data.", fig.width=5, fig.height=4----
ggplot(groups, aes(x = sampleGroup, y = n)) +
  geom_bar(stat = "identity")

## ----checkgeombar, echo = FALSE------------------------------------------
## check an assertion made in the text above
stopifnot(formals(ggplot2::geom_bar)$stat=="count")

## ----groupColor----------------------------------------------------------
groupColor = setNames(groups$color, groups$sampleGroup)

## ----rgraphics-qplot2, fig.keep = 'high', fig.cap = "Similar to Figure \\@ref(fig:rgraphics-qplot1), but with colored bars and better bar labels.\\label{rgraphics-qplot2}", fig.width=5, fig.height=4----
ggplot(groups, aes(x = sampleGroup, y = n, fill = sampleGroup)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = groupColor, name = "Groups") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


## ----ggplotobject--------------------------------------------------------
gg = ggplot(DNase, aes(x = conc, y = density)) + geom_point()

## ----ggpprintobject------------------------------------------------------
gg
print(gg)

## ----plotsave1-----------------------------------------------------------
ggsave("DNAse-histogram-demo.pdf", plot = gg)

## ----plotsave2, echo = FALSE, results = "hide"---------------------------
file.remove("DNAse-histogram-demo.pdf")

## ----loadlib, echo = FALSE, message = FALSE------------------------------
library("mouse4302.db")

## ----findprobepairs, echo = FALSE, eval = FALSE--------------------------
## # I used this code to find the below two probes
## idx = order(rowVars(Biobase::exprs(x)), decreasing=TRUE)[seq_len(2000)]
## cc  = cor(t(Biobase::exprs(x)[idx,]))
## cco = order(cc)[seq(1, 1001, by=2) ]
## jj2 = rownames(Biobase::exprs(x))[ idx[ (cco-1) %/% length(idx) + 1 ] ]
## jj1 = rownames(Biobase::exprs(x))[ idx[ (cco-1) %%  length(idx) + 1 ] ]
## dftx = as.data.frame(t(Biobase::exprs(x)))
## par(ask=TRUE)
## for(i in seq(along = cco)) {
##   df = AnnotationDbi::select(mouse4302.db,
##    keys = c(jj1[i], jj2[i]), keytype = "PROBEID",
##    columns = c("SYMBOL", "GENENAME"))
##   print(ggplot(dftx, aes( x = get(jj1[i]), y = get(jj2[i]))) +
##   geom_point(shape = 1) +
##   xlab(paste(jj1[i], df$SYMBOL[1])) +
##   ylab(paste(jj2[i], df$SYMBOL[2])) +
##   ggtitle(round(cc[jj1[i], jj2[i]], 3)) +
##   geom_smooth(method = "loess"))
## }

## ----rgraphics-scp2layers1, fig.keep = 'high', fig.cap = "A scatterplot with three layers that show different statistics of the same data: points (`geom_point`), a smooth regression line and a confidence band (the latter two from `geom_smooth`).", fig.width = 3.5, fig.height = 3.5----
dftx = data.frame(t(Biobase::exprs(x)), pData(x))
ggplot( dftx, aes( x = X1426642_at, y = X1418765_at)) +
  geom_point( shape = 1 ) +
  geom_smooth( method = "loess" )

## ----checkclassdftx, echo=FALSE------------------------------------------
stopifnot(is(dftx, "data.frame"))

## ----rgraphics-scp2layers2, fig.keep = 'high', fig.cap = "As Figure \\@ref(fig:rgraphics-scp2layers1), but in addition with points colored by the time point and cell lineage (as defined in Figure \\@ref(fig:rgraphics-qplot2)). We can now see that the expression values of the gene (ref:rgraphics-scp2layers2-1) (targeted by the probe 1418765_at) are consistently high in the early time points, whereas its expression goes down in the EPI samples at days 3.5 and 4.5. In the FGF4-KO, this decrease is delayed - at E3.5, its expression is still high. Conversely, the gene (ref:rgraphics-scp2layers2-2) (1426642_at) is off in the early timepoints and then goes up at days 3.5 and 4.5. The PE samples (green) show a high degree of cell-to-cell variability.", fig.width = 3.5, fig.height = 3.5----
ggplot( dftx, aes( x = X1426642_at, y = X1418765_at ))  +
  geom_point( aes( color = sampleColour), shape = 19 ) +
  geom_smooth( method = "loess" ) +
  scale_color_discrete( guide = FALSE )


## ----mouse4302.db, results="hide", message=FALSE-------------------------
library("mouse4302.db")

## ----select, warning = FALSE---------------------------------------------
AnnotationDbi::select(mouse4302.db,
   keys = c("1426642_at", "1418765_at"), keytype = "PROBEID",
   columns = c("SYMBOL", "GENENAME"))

## ----rgraphics-hists, fig.keep = 'high', fig.cap = "Histogram of probe intensities for one particular sample, cell number 20, which was from day E3.25.", fig.width = 3.5, fig.height = 2.5----
dfx = as.data.frame(Biobase::exprs(x))
ggplot(dfx, aes(x = `20 E3.25`)) + geom_histogram(binwidth = 0.2)

## ----figbpgg1------------------------------------------------------------
pb = ggplot(groups, aes(x = sampleGroup, y = n))

## ----rgraphics-figbpempty, fig.keep = 'high', fig.cap = "`pb`: without a geometric object, the plot remains empty.", fig.width = 3.2, fig.height = 2.5----
class(pb)
pb

## ----rgraphics-bpgg3, fig.keep = 'high', fig.cap = "The graphics object `bp` in its full glory.", fig.width = 5, fig.height = 4----
pb = pb + geom_bar(stat = "identity")
pb = pb + aes(fill = sampleGroup)
pb = pb + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pb = pb + scale_fill_manual(values = groupColor, name = "Groups")
pb

## ----rgraphics-bpgg7, fig.keep = 'high', fig.cap = "A barplot in a polar coordinate system.", fig.width = 5, fig.height = 4----
pb.polar = pb + coord_polar() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") + ylab("")
pb.polar

## ----genes2ps1-----------------------------------------------------------
selectedProbes = c( Fgf4 = "1420085_at", Gata4 = "1418863_at",
                   Gata6 = "1425463_at",  Sox2 = "1416967_at")

## ----genes2ps2, echo = FALSE, eval = FALSE-------------------------------
## # How I found the selectedProbes:
## AnnotationDbi::select(mouse4302.db,
##    keys = c("Fgf4", "Sox2", "Gata6", "Gata4"), keytype = "SYMBOL",
##    columns = c("PROBEID"))

## ----genes2ps3, echo = FALSE, warning = FALSE----------------------------
selectedProbes2 = AnnotationDbi::select(mouse4302.db,
   keys = selectedProbes, keytype = "PROBEID", columns = c("SYMBOL"))
stopifnot(identical(sort(selectedProbes2$SYMBOL), sort(names(selectedProbes))),
          all(selectedProbes[selectedProbes2$SYMBOL] == selectedProbes2$PROBEID))

## ----melt----------------------------------------------------------------
library("reshape2")
genes = melt(Biobase::exprs(x)[selectedProbes, ],
             varnames = c("probe", "sample"))
head(genes)

## ----symbol--------------------------------------------------------------
genes$gene =
  names(selectedProbes)[match(genes$probe, selectedProbes)]

## ----rgraphics-onedbp1, fig.keep = 'high', fig.cap = "Barplots showing the means of the distributions of expression measurements from four probes.", fig.width = 3, fig.height = 3.75----
ggplot(genes, aes( x = gene, y = value)) +
  stat_summary(fun.y = mean, geom = "bar")

## ----rgraphics-onedbp2, fig.keep = 'high', fig.cap = "Barplots with error bars indicating standard error of the mean.", fig.width = 3.75, fig.height = 3.75----
library("Hmisc")
ggplot(genes, aes( x = gene, y = value, fill = gene)) +
  stat_summary(fun.y = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
               width = 0.25)

## ----rgraphics-onedboxpl, fig.keep = 'high', fig.cap = "Boxplots.", fig.width = 3.75, fig.height = 3.75----
p = ggplot(genes, aes( x = gene, y = value, fill = gene))
p + geom_boxplot()

## ----rgraphics-onedviolin, fig.keep = 'high', fig.cap = "Violin plots.", fig.width = 3.75, fig.height = 3.75----
p + geom_violin()

## ----rgraphics-oneddot, fig.keep = 'high', fig.show = "hold", fig.cap = "Left: dot plots, made using `geom_dotplot` from **[ggplot2](https://cran.r-project.org/web/packages/ggplot2/)**. Right: beeswarm plots, made using `geom_beeswarm` from **[ggbeeswarm](https://cran.r-project.org/web/packages/ggbeeswarm/)**.", fig.margin = FALSE, out.width = '50%', fig.width = 5, fig.height = 5----
p + geom_dotplot(binaxis = "y", binwidth = 1/6,
       stackdir = "center", stackratio = 0.75,
       aes(color = gene))
library("ggbeeswarm")
p + geom_beeswarm(aes(color = gene))

## ----rgraphics-oneddens, fig.keep = 'high', fig.cap = "Density plots.", fig.width = 3.75, fig.height = 3.75----
ggplot(genes, aes( x = value, color = gene)) + geom_density()

## ----rgraphics-ecdfexample, fig.keep = 'high', fig.cap = "Sorted values of `simdata` versus their index. This is the empirical cumulative distribution function of `simdata`.", fig.width = 3.75, fig.height = 3.75----
simdata = rnorm(70)
tibble(index = seq(along = simdata),
          sx = sort(simdata)) %>%
ggplot(aes(x = sx, y = index)) + geom_step()

## ----rgraphics-onedecdf, fig.keep = 'high', fig.cap = "Empirical cumulative distribution functions (ECDF).", fig.width = 3.75, fig.height = 3.75----
ggplot(genes, aes( x = value, color = gene)) + stat_ecdf()

## ---- rgraphics-Lawrence-TCGA-Nature-2013-Fig1, fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Part of Figure 1 from  @Lawrence:Nature:2013. Each dot corresponds to a tumour-normal pair, with vertical position indicating the total frequency of somatic mutations in the exome. The resulting curves are, in essence, ECDF plots, and conceptually this plot is similar to Figure \\@ref(fig:rgraphics-onedecdf), just that the graphs are rotated by 90 degrees (i.e., the roles of $x$- and $y$-axis are exchanged) and the curves for the individual tumor types are horizontally displaced to keep them better apart."----
knitr::include_graphics(c('images/Lawrence-TCGA-Nature-2013-Fig1.png'))

## ----modes, eval = FALSE, echo = FALSE-----------------------------------
## # I used the functon bimodality_coefficient from the modes package to identify the most
## # bimodal looking array, number 64
## library("modes")
## j0 = which.max(vapply(seq_len(ncol(x)), function(j){
##       bimodality_coefficient(Biobase::exprs(x)[, j])
##     }, numeric(1)))

## ----rgraphics-onedtrsf, fig.keep = 'high', fig.show = "hold", fig.cap = "Histograms of the same data, with and without logarithm transform. On the top, the data are shown on the scale on which they are stored in the data object `x`, which resulted from logarithm (base 2) transformation of the microarray fluorescence intensities [@Irizarry:Biostat:2003], ; on the bottom, after re-exponentiating them back to the fluorescence scale. For better use of space, we capped the $x$-axis range at 1500.", fig.margin = FALSE, out.width = '50%', fig.width = 3.75, fig.height = 3.5----
ggplot(dfx, aes(x = `64 E4.5 (EPI)`)) + geom_histogram(bins = 100)
ggplot(dfx, aes(x = 2 ^ `64 E4.5 (EPI)`)) + 
  geom_histogram(binwidth = 20) + xlim(0, 1500)

## ----rgraphics-twodsp1, fig.keep = 'high', fig.cap = "Scatterplot of (ref:rgraphics-twodsp1-1) expression measurements for two of the samples.", fig.width = 3.75, fig.height = 3.75, dev = "png"----
scp = ggplot(dfx, aes(x = `59 E4.5 (PE)` ,
                      y = `92 E4.5 (FGF4-KO)`))
scp + geom_point()

## ----rgraphics-twodsp2, fig.keep = 'high', fig.cap = "As Figure \\@ref(fig:rgraphics-twodsp1), but with semi-transparent points to resolve some of the overplotting.", fig.width = 3.75, fig.height = 3.75, dev = "png"----
scp  + geom_point(alpha = 0.1)

## ----rgraphics-twodsp3, fig.keep = 'high', fig.cap = "As Figure \\@ref(fig:rgraphics-twodsp1), but rendered as a contour plot of the 2D density estimate.", fig.width = 3.75, fig.height = 3.75----
scp + geom_density2d()

## ----rgraphics-twodsp4, fig.keep = 'high', fig.show = "hold", fig.cap = "Left: as Figure \\@ref(fig:rgraphics-twodsp3), but with smaller smoothing bandwidth and tighter binning for the contour lines. Right: with color filling.", fig.margin = FALSE, out.width = '50%', fig.width = 3.75, fig.height = 3.75----
scp + geom_density2d(h = 0.5, bins = 60)
library("RColorBrewer")
colorscale = scale_fill_gradientn(
    colors = rev(brewer.pal(9, "YlGnBu")),
    values = c(0, exp(seq(-5, 0, length.out = 100))))

scp + stat_density2d(h = 0.5, bins = 60,
          aes( fill = ..level..), geom = "polygon") +
  colorscale + coord_fixed()

## ----rgraphics-twodsp6, fig.keep = 'high', fig.show = "hold", fig.cap = "Hexagonal binning. Left: default parameters. Right: finer bin sizes and customized color scale.", fig.margin = FALSE, out.width = '50%', fig.width = 5.25, fig.height = 3.75----
scp + geom_hex() + coord_fixed()
scp + geom_hex(binwidth = c(0.2, 0.2)) + colorscale +
  coord_fixed()

## ----rgraphics-banking, fig.keep = 'high', fig.show = "hold", fig.cap = "The sunspot data. In the upper panel, the plot shape is roughly quadratic, a frequent default choice. In the lower panel, a technique called **banking** was used to choose the plot shape. (Note: the placement of the tick labels is not great in this plot and would benefit from customization.)", fig.width = 3.75, fig.height = 3.75----
library("ggthemes")
sunsp = tibble(year   = time(sunspot.year),
               number = as.numeric(sunspot.year))
sp = ggplot(sunsp, aes(x = year, y = number)) + geom_line()
sp
ratio = with(sunsp, bank_slopes(year, number))
sp + coord_fixed(ratio = ratio)



## ----rgraphics-facet1, fig.keep = 'high', fig.cap = "An example of **faceting**: the same data as in Figure \\@ref(fig:rgraphics-scp2layers1), but now split by the categorical variable `lineage`.", fig.margin = FALSE, fig.width = 8, fig.height = 2----
library("magrittr")
dftx$lineage %<>% sub("^$", "no", .)
dftx$lineage %<>% factor(levels = c("no", "EPI", "PE", "FGF4-KO"))

ggplot(dftx, aes( x = X1426642_at, y = X1418765_at)) +
  geom_point() + facet_grid( . ~ lineage )

## ----rgraphics-facet2, fig.keep = 'high', fig.cap = "**Faceting**: the same data as in Figure \\@ref(fig:rgraphics-scp2layers1), split by the categorical variables `Embryonic.day` (rows) and `lineage` (columns).", fig.margin = FALSE, fig.width = 8, fig.height = 6----
ggplot( dftx,
  aes( x = X1426642_at, y = X1418765_at)) + geom_point() +
   facet_grid( Embryonic.day ~ lineage )

## ----rgraphics-facet3, fig.keep = 'high', fig.cap = "**Faceting**: the same data as in Figure \\@ref(fig:rgraphics-scp2layers1), split by the continuous variable `X1450989_at` and arranged by `facet_wrap`.", fig.width = 4, fig.height = 4----
ggplot(mutate(dftx, Tdgf1 = cut(X1450989_at, breaks = 4)),
   aes( x = X1426642_at, y = X1418765_at)) + geom_point() +
   facet_wrap( ~ Tdgf1, ncol = 2 )

## ----plotly, eval = FALSE------------------------------------------------
## library("plotly")
## plot_ly(economics, x = ~ date, y = ~ unemploy / pop)

## ---- rgraphics-rglvolcano, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "**[rgl](https://cran.r-project.org/web/packages/rgl/)** rendering of the `volcano` data, the topographic information for Maunga Whau (Mt Eden), one of about 50 volcanos in the Auckland volcanic field."----
knitr::include_graphics(c('images/chap3-rglvolcano.png'))

## ----volcano1------------------------------------------------------------
data("volcano")
volcanoData = list(
  x = 10 * seq_len(nrow(volcano)),
  y = 10 * seq_len(ncol(volcano)),
  z = volcano,
  col = terrain.colors(500)[cut(volcano, breaks = 500)]
)
library("rgl")
with(volcanoData, persp3d(x, y, z, color = col))

## ----volcano2, echo = FALSE----------------------------------------------
rgl.snapshot(filename = file.path("figure", "chap3-r_rglvolcano.png"))

## ----volcano3, echo = FALSE----------------------------------------------
.volcanocut = cut(volcano, breaks = 500)
stopifnot(!any(is.na(.volcanocut)), all(as.integer(.volcanocut) %in% 1:500))

## ----rgraphics-simplecolorpie, fig.keep = 'high', fig.cap = "Basic R colors. ", fig.width = exp(1), fig.height = exp(1), eval = TRUE, echo = FALSE----
par(mai = rep(0,4))
pie(rep(1, 8), col=1:8)

## ----simplecolorpieShow, eval = FALSE, echo = TRUE-----------------------
## pie(rep(1, 8), col=1:8)

## ----rgraphics-RColorBrewer, fig.keep = 'high', fig.cap = "RColorBrewer palettes.", fig.width = 4, fig.height = 8----
display.brewer.all()

## ----color3--------------------------------------------------------------
head(brewer.pal.info)
table(brewer.pal.info$category)

## ----color4--------------------------------------------------------------
brewer.pal(4, "RdYlGn")

## ----rgraphics-colorRampPalette, fig.keep = 'high', fig.cap = "A quasi-continuous color palette derived by interpolating between the colors `darkorange3`, `white` and `darkblue`.", fig.width = 3, fig.height = .7----
mypalette  = colorRampPalette(
    c("darkorange3", "white","darkblue")
  )(100)
head(mypalette)
par(mai = rep(0.1, 4))
image(matrix(1:100, nrow = 100, ncol = 10), col = mypalette,
        xaxt = "n", yaxt = "n", useRaster = TRUE)

## ----rgraphics-heatmap, fig.keep = 'high', fig.cap = "A heatmap of relative expression values, i.\\,e., logarithmic fold change compared to the average expression of that gene (row) across all samples (columns). The color scale uses a diverging palette whose midpoint is at 0.", fig.margin = FALSE, fig.width = 7, fig.height = 7----
library("pheatmap")
topGenes = order(rowVars(Biobase::exprs(x)), decreasing = TRUE)[1:500]
rowCenter = function(x) { x - rowMeans(x) }
pheatmap( rowCenter(Biobase::exprs(x)[ topGenes, ] ),
  show_rownames = FALSE, show_colnames = FALSE,
  breaks = seq(-5, +5, length = 101),
  annotation_col =
    pData(x)[, c("sampleGroup", "Embryonic.day", "ScanDate") ],
  annotation_colors = list(
    sampleGroup = groupColor,
    genotype = c(`FGF4-KO` = "chocolate1", `WT` = "azure2"),
    Embryonic.day = setNames(brewer.pal(9, "Blues")[c(3, 6, 9)],
                             c("E3.25", "E3.5", "E4.5")),
    ScanDate = setNames(brewer.pal(nlevels(x$ScanDate), "YlGn"),
                        levels(x$ScanDate))
  ),
  cutree_rows = 4
)

## ----groupColor2---------------------------------------------------------
groupColor[1]

## ----hexvals, echo = FALSE-----------------------------------------------
hexvals = sapply(1:3, function(i) substr(groupColor[1], i*2, i*2+1))
decvals = strtoi(paste0("0x", hexvals))

## ----somecolors, echo = FALSE, results = "hide"--------------------------
library("colorspace")
library("grid")

plothcl = function(h, c, l, what, x0 = 0.5, y0 = 0.5, default.units = "npc", ...) {
  switch(what,
         "c" = {
           stopifnot(length(l)==1)
           n = length(c)
         },
         "l" = {
           stopifnot(length(c)==1)
           n = length(l)
         },
         stop("Sapperlot"))

  cr = seq(0.1, 0.5, length = n+1)
  dr = 0.05 / n

  for (j in seq_len(n)) {
    r = c(cr[j]+dr, cr[j+1]-dr)
    for(i in 1:(length(h)-1)){
      phi = seq(h[i], h[i+1], by=1)/180*pi
      px = x0 + c(r[1]*cos(phi), r[2]*rev(cos(phi)))
      py = y0 + c(r[1]*sin(phi), r[2]*rev(sin(phi)))
      mycol = switch(what,
        "c" = hcl(h=mean(h[i+(0:1)]), c=c[j], l=l),
        "l" = hcl(h=mean(h[i+(0:1)]), c=c, l=l[j]))
      grid.polygon(px, py, gp=gpar(col=mycol, fill=mycol),
                   default.units=default.units,...)
    }
  }
}

## ----rgraphics-hcl, fig.keep = 'high', fig.show = "hold", fig.cap = "\\label{rgraphics-hcl}Circles in HCL colorspace. Upper panel: The luminance $L$ is fixed at $75$, while the angular coordinate $H$ (hue) varies from 0 to 360 and the radial coordinate $C=0, 10, ..., 60$. Lower panel: constant chroma $C=50$, $H$ as above, and varying luminance $L=10, 20, ..., 90$.", fig.width = 3.6, fig.height = 3.6, echo = FALSE----
plothcl( h = seq(0, 360, by=3), c = seq(5, 75, by=10), l = 75,   what="c")
grid.newpage()
plothcl( h = seq(0, 360, by=3), c = 55, l = seq(20, 100, by=10), what="l")

## ----rgraphics-MA, fig.keep = 'high', fig.show = "hold", fig.cap = "The effect of rank transformation on the visual perception of dependency.", fig.width = 3.75, fig.height = 3.75, dev = "png"----
gg = ggplot(tibble(A = Biobase::exprs(x)[, 1], M = rnorm(length(A))),
            aes(y = M))
gg + geom_point(aes(x = A), size = 0.2)
gg + geom_point(aes(x = rank(A)), size = 0.2)

## ---- rgraphics-otherfont, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "As Figure \\@ref(fig:rgraphics-onedecdf), with font \"Bauhaus 93\"."----
knitr::include_graphics(c('images/chap3-ecdfs-bauhaus93font.png'))

## ----rgraphics-mathnot, fig.keep = 'high', fig.cap = "Volume $\\Omega$ of the $\\nu$-dimensional sphere with radius $\\rho=1$, for $\\nu=1,...,15$.", fig.width = 3, fig.height = 2.5----
volume = function(rho, nu)
            pi^(nu/2) * rho^nu / gamma(nu/2+1)

ggplot(tibble(nu    = 1:15,
  Omega = volume(1, nu)), aes(x = nu, y = Omega)) +
geom_line() +
xlab(expression(nu)) + ylab(expression(Omega)) +
geom_text(label =
"Omega(rho,nu)==frac(pi^frac(nu,2)~rho^nu, Gamma(frac(nu,2)+1))",
  parse = TRUE, x = 6, y = 1.5)

## ----rgraphics-timesfont, fig.keep = 'high', fig.cap = "As Figure \\@ref(fig:rgraphics-onedecdf), with a different font.", fig.width = 3.75, fig.height = 3----
ggplot(genes, aes( x = value, color = gene)) + stat_ecdf() +
  theme(text = element_text(family = "Times"))

## ----otherfont, eval = FALSE, echo = FALSE-------------------------------
## ggplot(genes, aes( x = value, color = gene)) + stat_ecdf() + theme(text = element_text(family = "Bauhaus 93"))

## ---- rgraphics-EBI-genomebrowser-rnaseq, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Screenshot from Ensembl genome browser, showing gene annotation of a genomic region as well as a read pile-up visualization of an RNA-Seq experiment."----
knitr::include_graphics(c('images/EBI-genomebrowser-rnaseq.png'))

## ----rgraphics-ideogram-1, fig.keep = 'high', fig.cap = "Chromosome 1 of the human genome: ideogram plot.", fig.width=4, fig.height=1.5, warning = FALSE----
library("ggbio")
data("hg19IdeogramCyto", package = "biovizBase")
plotIdeogram(hg19IdeogramCyto, subchr = "chr1")

## ----rgraphics-darned1, fig.keep = 'high', fig.cap = "Karyogram with RNA editing sites. `exReg` indicates whether a site is in the coding region (C), 3\'- or 5\'-UTR.", warning = FALSE----
library("GenomicRanges")
data("darned_hg19_subset500", package = "biovizBase")
autoplot(darned_hg19_subset500, layout = "karyogram",
         aes(color = exReg, fill = exReg))

## ----rgraphics-darned2, fig.keep = 'high', fig.cap = "Improved version of Figure \\@ref(fig:rgraphics-darned1)."----
data("ideoCyto", package = "biovizBase")
dn = darned_hg19_subset500
seqlengths(dn) = seqlengths(ideoCyto$hg19)[names(seqlengths(dn))]
dn = keepSeqlevels(dn, paste0("chr", c(1:22, "X")))
autoplot(dn, layout = "karyogram", aes(color = exReg, fill = exReg))

## ----whatisdarned1-------------------------------------------------------
darned_hg19_subset500[1:2,]

## ----whatisdarned2, echo = FALSE-----------------------------------------
stopifnot(is(darned_hg19_subset500, "GRanges"), identical(start(darned_hg19_subset500),end(darned_hg19_subset500)))

ggcars = ggplot(mtcars, aes(x = hp, y = mpg)) + geom_point()
ggcars
ggcars + theme_bw()
ggcars + theme_minimal()

knitr::include_graphics(c('images/xkcdgraph.png'))

