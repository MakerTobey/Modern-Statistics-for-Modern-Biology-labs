## ----initialize, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
source("../chapter-setup.R"); chaptersetup("/Users/Susan/Courses/CUBook-html/CUBook/Chap13-ImageAnalysis/ImageAnalysis.Rnw", "13")
knitr::opts_chunk$set(dev = 'png', dpi = 100, fig.margin = TRUE, fig.show = 'hold', fig.keep = 'none')
library("magrittr")
.h = list(width = 4, height = 3.3)

## ---- chap13-colorLabelscellbodies, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/chap13-colorLabelscellbodies.png'))

## ----data1, fig.cap=""---------------------------------------------------
library("EBImage")
imagefile = system.file("images", "mosquito.png",
                        package = "MSMB")
mosq = readImage(imagefile)

## ----vis1, eval = FALSE--------------------------------------------------
## display(mosq)

## ----mosquito, fig.keep = 'high', fig.cap = "Mosquito discovered deceased in the suburbs of Decatur, Georgia (credit: CDC / Janice Haney Carr).", dev = "png", dpi = 300, fig.width=dim(mosq)[1]/300, fig.height=dim(mosq)[2]/300----
display(mosq, method = "raster")
text(x = 85, y = 800, label = "A mosquito",
     adj = 0, col = "orange", cex = 1.5)

## ----vis3a, eval = FALSE-------------------------------------------------
## imagefile = system.file("images", "hiv.png",
##                         package = "MSMB")
## hivc = readImage(imagefile)
## display(hivc)

## ----vis3b, eval = TRUE, echo = FALSE------------------------------------
imagefile = system.file("images", "hiv.png",
                        package = "MSMB")
hivc = readImage(imagefile)

## ---- hiv, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Scanning electron micrograph of HIV-1 virions budding from a cultured lymphocyte (credit: CDC / C. Goldsmith, P. Feorino, E.L. Palmer, W.R. McManus)."----
knitr::include_graphics(c('images/chap13-hivc.png'))

## ----image-oneminus, fig.keep = 'high', fig.cap = "Tiled display of four images of cell nuclei from the **[EBImage](https://bioconductor.org/packages/EBImage/)** package.", fig.margin = FALSE, dev = "png"----
nuc = readImage(system.file("images", "nuclei.tif",
                            package = "EBImage"))
display(1 - nuc, method = "raster", all = TRUE)

## ----nucleitiledoneframe, dev = "png"------------------------------------
display(1 - nuc, method = "raster", frame = 2)

## ----how1----------------------------------------------------------------
class(mosq)

## ----how2, echo = FALSE, eval = FALSE------------------------------------
## showClass("Image")

## ----dim-----------------------------------------------------------------
dim(mosq)

## ----mosqhist, fig.keep = 'high', fig.cap = "Histogram of the pixel intensities in `mosq`. Note that the range is between 0 and 1.", fig.width = .h$width, fig.height = .h$height----
hist(mosq)

## ----check, echo=FALSE---------------------------------------------------
stopifnot(all(mosq>=0 & mosq<=1), isTRUE(all.equal(max(mosq), 1)), isTRUE(all.equal(min(mosq), 0)))

## ----how3----------------------------------------------------------------
imageData(mosq)[1:3, 1:6]

## ----show1---------------------------------------------------------------
mosq

## ----show2---------------------------------------------------------------
hivc

## ----show3, echo=FALSE---------------------------------------------------
stopifnot(colorMode(mosq)==Grayscale, colorMode(hivc)==Color, dim(nuc)[3]==4)

## ----how6, results = "hide"----------------------------------------------
nuc
dim(imageData(nuc))

## ----checkassertion, echo = FALSE----------------------------------------
stopifnot(all(c("  frames.total : 4 ", "  frames.render: 4 ") %in%
              capture.output(EBImage:::showImage(nuc))))

## ----write1--------------------------------------------------------------
writeImage(hivc, "hivc.jpeg", quality = 85)

## ----objectsize----------------------------------------------------------
object.size(hivc) %>% format(units = "Mb")
(object.size(hivc) / prod(dim(hivc))) %>% format %>% paste("per pixel")
file.info("hivc.jpeg")$size
16 * 3 * 8

## ---- manip1, out.width = '25%', fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The original mosquito image (top left) and three different image transformations (subtraction, multiplication, power transformation).\\label{manip1}"----
knitr::include_graphics(c('images/chap13-mosq.png','images/chap13-mosqinv.png','images/chap13-mosqcont.png','images/chap13-mosqexp.png'))

## ---- mosqcrop, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Cropping: `mosqcrop`"----
knitr::include_graphics(c('images/chap13-mosqcrop.png'))

## ----manip1a-------------------------------------------------------------
mosqinv = normalize(-mosq)

## ----manip3a-------------------------------------------------------------
mosqcont = mosq * 3
mosqexp = mosq ^ (1/3)

## ---- mosqthresh, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Threshold: `mosqthresh`"----
knitr::include_graphics(c('images/chap13-mosqthresh.png'))

## ----manip4a-------------------------------------------------------------
mosqcrop   = mosq[100:438, 112:550]
mosqthresh = mosq > 0.5
mosqtransp = transpose(mosq)

## ---- mosqtransp, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Transposition: `mosqtransp`"----
knitr::include_graphics(c('images/chap13-mosqtransp.png'))

## ----checkassertionont, echo = FALSE-------------------------------------
stopifnot(identical(t(mosq), transpose(mosq)))

## ---- flipflop, out.width = '25%', fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Spatial transformations: rotation (top left), translation (top right), reflection about the central horizontal axis (`flip`, bottom left), reflection about the central vertical axis (`flop`, bottom right). \\label{flipflop}"----
knitr::include_graphics(c('images/chap13-mosqrot.png','images/chap13-mosqshift.png','images/chap13-mosqflip.png','images/chap13-mosqflop.png'))

## ----spattrans1----------------------------------------------------------
mosqrot   = EBImage::rotate(mosq, angle = 30)
mosqshift = translate(mosq, v = c(40, 70))
mosqflip  = flip(mosq)
mosqflop  = flop(mosq)

## ----MSMB, results="hide"------------------------------------------------
imagefiles = system.file("images", c("image-DAPI.tif",
  "image-FITC.tif", "image-Cy3.tif"), package="MSMB")
cells = readImage(imagefiles)

## ----checkdim, echo=FALSE------------------------------------------------
stopifnot(dim(cells)[3]==3)

## ---- LauferCells, fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Human colon cancer cells (HCT116). The four images show the same cells: the leftmost image corresponds to DAPI staining of the cells\' DNA, the second to immunostaining against $\\alpha$-tubulin, the third to actin. They are displayed as gray-scale images. The rightmost image is obtained by overlaying the three images as color channels of an RGB image (red: actin, green: $\\alpha$-tubulin, blue: DNA)."----
knitr::include_graphics(c('images/chap13-LauferCells.png'))

## ----range---------------------------------------------------------------
apply(cells, 3, range)

## ----fixrange------------------------------------------------------------
cells[,,1]   = 32 * cells[,,1]
cells[,,2:3] = 16 * cells[,,2:3]
apply(cells, 3, range)

## ----writeCells, echo = FALSE--------------------------------------------
# combined is also defined in dplyr, tile seems also popular....
writeImage(EBImage::tile(EBImage::combine(
  toRGB(getFrame(cells, 1)),
  toRGB(getFrame(cells, 2)),
  toRGB(getFrame(cells, 3)),
  rgbImage(red   = getFrame(cells, 3),
           green = getFrame(cells, 2),
           blue  = getFrame(cells, 1))), nx = 4, lwd = 5, fg.col = "white"),
  files = paste0(knitr::opts_chunk$get("fig.path"), "LauferCells.png")) #$

## ----defw----------------------------------------------------------------
w = makeBrush(size = 51, shape = "gaussian", sigma = 7)
nucSmooth = filter2(getFrame(cells, 1), w)

## ----image-filter2, fig.keep = 'high', fig.cap = "The middle row of the weight matrix, \\texttt{w[(ref:image-filter2-1), ]}. \\label{image-filter2}", fig.width = 3, fig.height = 2.6, dev = "pdf"----
library("tibble")
library("ggplot2")
tibble(w = w[(nrow(w)+1)/2, ]) %>%
  ggplot(aes(y = w, x = seq(along = w))) + geom_point()

## ---- nucSmooth, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "`nucSmooth`, a smoothed version of the DNA channel in the image object `cells` (the original version is shown in the leftmost panel of Figure \\@ref(fig:LauferCells))"----
knitr::include_graphics(c('images/chap13-nucSmooth.png'))

## ----smooth--------------------------------------------------------------
cellsSmooth = Image(dim = dim(cells))
sigma = c(1, 3, 3)
for(i in seq_along(sigma))
  cellsSmooth[,,i] = filter2( cells[,,i],
         filter = makeBrush(size = 51, shape = "gaussian",
                            sigma = sigma[i]) )

## ----illuminationartifact1-----------------------------------------------
py = seq(-1, +1, length.out = dim(cellsSmooth)[1])
px = seq(-1, +1, length.out = dim(cellsSmooth)[2])
illuminationGradient = Image(
     outer(py, px, function(x, y) exp(-(x^2+y^2))))
nucBadlyIlluminated = cellsSmooth[,,1] * illuminationGradient

## ----illuminationartifact2-----------------------------------------------
disc = makeBrush(21, "disc")
disc = disc / sum(disc)
localBackground = filter2(nucBadlyIlluminated, disc)
offset = 0.02
nucBadThresh = (nucBadlyIlluminated - localBackground > offset)

## ---- illumination, out.width = '25%', fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "From left to right: 1. `illuminationGradient`, a function that has its maximum at the center and falls off towards the sides, and which simulates uneven illumination sometimes seen in images. 2. `nucBadlyIlluminated`, the image that results from multiplying the DNA channel in `cellsSmooth` with `illuminationGradient`. 3. `localBackground`, the result of applying a linear filter with a bandwidth that is larger than the objects to be detected. 4. `nucBadThresh`, the result of adaptive thresholding. The nuclei at the periphery of the image are reasonably well identified, despite the drop off in signal strength."----
knitr::include_graphics(c('images/chap13-illuminationGradient.png','images/chap13-nucBadlyIlluminated.png','images/chap13-localBackground.png','images/chap13-nucBadThresh.png'))

## ----adathresh-----------------------------------------------------------
nucThresh =
  (cellsSmooth[,,1] - filter2(cellsSmooth[,,1], disc) > offset)

## ----morphopen1----------------------------------------------------------
nucOpened = EBImage::opening(nucThresh,
                  kern = makeBrush(5, shape = "disc"))

## ---- morphop, out.width = '20%', fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Different steps in the segmentation of the nuclei."----
knitr::include_graphics(c('images/chap13-nucThresh.png','images/chap13-nucOpened.png','images/chap13-colorLabelsnucSeed.png','images/chap13-nucMask.png','images/chap13-colorLabelsnuclei.png'))

## ----imageProcessing14---------------------------------------------------
nucSeed = bwlabel(nucOpened)
table(nucSeed)

## ----imageProcessing17, eval = FALSE-------------------------------------
## display(colorLabels(nucSeed))

## ----imageProcessing15a--------------------------------------------------
nucMask = cellsSmooth[,,1] - filter2(cellsSmooth[,,1], disc) > 0

## ----imageProcessing15b--------------------------------------------------
nucMask = fillHull(nucMask)

## ----imageProcessing16---------------------------------------------------
nuclei = propagate(cellsSmooth[,,1], nucSeed, mask = nucMask)

## ---- voronoiPaint, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Example of a Voronoi segmentation, indicated by the gray lines, using the nuclei (indicated by black regions) as seeds.\\label{voronoiPaint}."----
knitr::include_graphics(c('images/chap13-voronoiPaint.png'))

## ----voronoiExample------------------------------------------------------
zeros        = Image(dim = dim(nuclei))
voronoiExamp = propagate(seeds = nuclei, x = zeros, lambda = 100)
voronoiPaint = paintObjects(voronoiExamp, 1 - nucOpened)

## ----voronoiEx-----------------------------------------------------------
head(table(voronoiExamp))
ind = which(voronoiExamp == 13, arr.ind = TRUE)
head(ind, 3)

## ---- histcellbody-2, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Zoom into Figure \\@ref(fig:histcellbody-1)."----
knitr::include_graphics(c('images/chap13-r_histcellbody-2.png'))

## ----histcellbody-1, fig.keep = 'high', fig.cap = "Histogram of the actin channel in `cellsSmooth`, after taking the logarithm.", fig.width = .h$width, fig.height = .h$height----
hist(log(cellsSmooth[,,3]) )
hist(log(cellsSmooth[,,3]), xlim = -c(3.6, 3.1), breaks = 300)

## ----checkhistcellsSmooth, echo=FALSE------------------------------------
stopifnot(mean(cellsSmooth[,,3]>=exp(-3.6) & cellsSmooth[,,3]<=exp(-3.1)) > 0.68)

## ----musigmaEstimator, message=FALSE-------------------------------------
library("genefilter")
bgPars = function(x) {
  x    = log(x)
  loc  = half.range.mode( x )
  left = (x - loc)[ x < loc ]
  wid  = sqrt( mean(left^2) )
  c(loc = loc, wid = wid, thr = loc + 6*wid)
}
cellBg = apply(cellsSmooth, MARGIN = 3, FUN = bgPars)
cellBg

## ----histcellbody-3, fig.keep = 'high', fig.cap = "As in Figure \\@ref(fig:histcellbody-2), but with `loc` and `thr` shown by vertical lines.", fig.width = .h$width, fig.height = .h$height----
hist(log(cellsSmooth[,,3]), xlim = -c(3.6, 3.1), breaks = 300)
abline(v = cellBg[c("loc", "thr"), 3], col = c("brown", "red"))

## ----cytoplasmMask-------------------------------------------------------
cytoplasmMask = (cellsSmooth[,,2] > exp(cellBg["thr", 2])) |
       nuclei | (cellsSmooth[,,3] > exp(cellBg["thr", 3]))

## ----imageProcessing22---------------------------------------------------
cellbodies = propagate(x = cellsSmooth[,,3], seeds = nuclei,
                       lambda = 1.0e-2, mask = cytoplasmMask)

## ---- cellbodies, out.width = '20%', fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Steps in the segmentation of the cell bodies.\\label{cellbodies}"----
knitr::include_graphics(c('images/chap13-cytoplasmMask.png','images/chap13-colorLabelscellbodies.png','images/chap13-nucSegOnNuc.png','images/chap13-nucSegOnAll.png','images/chap13-cellSegOnAll.png'))

## ----imageProcessing25---------------------------------------------------
cellsColor = rgbImage(red   = cells[,,3],
                      green = cells[,,2],
                      blue  = cells[,,1])

nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(cells[,,1]),
                            col = "#ffff00")
nucSegOnAll  = paintObjects(nuclei, tgt = cellsColor,
                            col = "#ffff00")
cellSegOnAll = paintObjects(cellbodies, tgt = nucSegOnAll,
                            col = "#ff0080")

## ----writeImages, echo = FALSE-------------------------------------------
for(v in c("illuminationGradient", "nucBadlyIlluminated", "localBackground", "nucBadThresh",
           "nucThresh", "nucOpened", "colorLabels(nucSeed)", "nucMask", "colorLabels(nuclei)",
           "cytoplasmMask", "colorLabels(cellbodies)", "nucSegOnNuc", "nucSegOnAll", "nucSmooth", "cellSegOnAll",
           "voronoiPaint", "mosq", "mosqinv", "mosqcont", "mosqexp", "mosqcrop", "mosqthresh", "mosqtransp",
           "mosqrot", "mosqshift", "mosqflip", "mosqflop", "hivc"))
  writeImage(eval(parse(text = v)), files =
    paste0(knitr::opts_chunk$get("fig.path"), gsub("[[:punct:]]", "", v), ".png")) #$
  ## 'eval(parse())' since we also want to evaluate expressions "colorLabels(nucSeed)" etc.

## ----baserfeats----------------------------------------------------------
meanNucInt       = tapply(cells[,,1], nuclei, mean)
meanActIntInNuc  = tapply(cells[,,3], nuclei, mean)
meanActIntInCell = tapply(cells[,,3], cellbodies, mean)

## ----pairsint, fig.keep = 'high', fig.cap = "Pairwise scatterplots of per-cell intensity descriptors. \\label{pairsint}", fig.margin = FALSE, fig.width = 5.3, fig.height = 5.3----
library("GGally")
ggpairs(tibble(meanNucInt, meanActIntInNuc, meanActIntInCell))

## ----imageProcessing27---------------------------------------------------
F1 = computeFeatures(nuclei,     cells[,,1], xname = "nuc",
                                             refnames = "nuc")
F2 = computeFeatures(cellbodies, cells[,,2], xname = "cell",
                                             refnames = "tub")
F3 = computeFeatures(cellbodies, cells[,,3], xname = "cell",
                                             refnames = "act")
dim(F1)

## ----showF1--------------------------------------------------------------
 F1[1:3, 1:5]

## ---- sixpanelslymphsmall, fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Biopsy of an enlarged lymph node revealed an intact capsule and obliterated sinuses (upper left panel, stained with hematoxylin and eosin, original magnification $\\times$ 100). The infiltrate was composed of an admixture of small lymphocytes, macrophages, and plasma cells (upper right panel, hematoxylin and eosin, original magnification $\\times$ 400). The infiltrate was composed of a mixture of CD3 positive T-cells (including both CD4 and CD8 positive cells) and CD20 positive B-cells. Numerous macrophages were also CD4 positive. (From: Hurley et al., Diagnostic Pathology (2008) 3:13)\\label{sixpanelslymphsmall}"----
knitr::include_graphics(c('images/SixPanelsLymphsmall.png'))

## ---- stainedlymphnode, fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "A stained lymph node; this image is the basis for the spatial data in `brcalymphnode`."----
knitr::include_graphics(c('images/testscan_1_2_RGB99-4525D.jpg'))

## ----readlymphnodedata---------------------------------------------------
library("readr")
library("dplyr")
cellclasses = c("T_cells", "Tumor", "DCs", "other_cells")
brcalymphnode = lapply(cellclasses, function(k) {
    read_csv(file.path("..", "data",
             sprintf("99_4525D-%s.txt", k))) %>%
    transmute(x = globalX,
              y = globalY,
              class = k)
}) %>% bind_rows %>% mutate(class = factor(class))

brcalymphnode
table(brcalymphnode$class)

## ----checkcellnumber, echo = FALSE---------------------------------------
tabtab = table(brcalymphnode$class)
within = function(x, a, b) (x>a & x<b)
stopifnot(all(within(tabtab[c("T_cells", "Tumor", "DCs")], c(100000, 27000, 800), c(110000, 28000, 1000))))

## ----brcalntcells, fig.keep = 'high', fig.cap = "Scatterplot of the $x$ and $y$ positions of the T- and tumor cells in `brcalymphnode`. The locations were obtained by a segmentation algorithm from a high resolution version of Figure \\@ref(fig:stainedlymphnode). Some rectangular areas in the T-cells plot are suspiciously empty, this could be because the corresponding image tiles within the overall composite image went missing, or were not analyzed.", fig.margin = FALSE, dev = "png", dpi = 300, fig.width=9, fig.height=4.5, pointsize=24----
ggplot(filter(brcalymphnode, class %in% c("T_cells", "Tumor")),
   aes(x = x, y = y, col = class)) + geom_point(shape = ".") +
   facet_grid( . ~ class) + guides(col = FALSE)

## ----spatstat1, results = "hide"-----------------------------------------
library("spatstat")

## ----spatstat2-----------------------------------------------------------
ln = with(brcalymphnode,
  ppp(x = x, y = y, marks = class, xrange = range(x), yrange = range(y)))
ln

## ----checkclassln, echo = FALSE------------------------------------------
stopifnot(identical(class(ln), "ppp"))

## ----convhull1, results="hide"-------------------------------------------
library("geometry")
coords = cbind(ln$x, ln$y)
chull = convhulln( coords )

## ----convhull2-----------------------------------------------------------
pidx = integer(nrow(chull) + 1)
pidx[1:2] = chull[1, ]
  wh = which(chull == pidx[j-1], arr.ind = TRUE)
  stopifnot(nrow(wh )== 1)
  wh[, "col"] = 3 - wh[, "col"] ## 2->1, 1->2
  pidx[j] = chull[wh]
pidx = rev(pidx)

## ----convhull, fig.keep = 'high', fig.cap = "Polygon describing the convex hull of the points in `ln`.", fig.width=4, fig.height=4----
ggplot(tibble(x = ln$x, y = ln$y)[pidx, ], aes(x = x, y = y)) +
  geom_point() + geom_path() + coord_fixed()

## ----spatstat3-----------------------------------------------------------
ln = with(brcalymphnode,
   ppp(x = x, y = y, marks = class, poly = coords[ pidx, ],
       check = FALSE))
ln

## ---- raindrops, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Rain drops falling on the floor are modelled by a Poisson process. The number of drops falling on a particular spot only depends on the rate $\\lambda$ (and on the size of the spot), but not on what happens at other spots.\\label{raindrops}"----
knitr::include_graphics(c('images/Rain-Drops-small.png'))

## ----densityppz1---------------------------------------------------------
d = density(subset(ln, marks == "Tumor"), edge=TRUE, diggle=TRUE)
plot(d)

## ----densityppp1, fig.keep = 'high', fig.cap = "Intensity estimate for the cells marked `Tumor` in `ppp`. The support of the estimate is the polygon that we specified earlier on (Figure \\@ref(fig:convhull)).", fig.width = 3.75, fig.height = 3.5, echo = FALSE----
par(mai = c(0, 0, 0.2, 0.7))
plot(d)

## ----densityppz0---------------------------------------------------------
d0 = density(subset(ln, marks == "Tumor"), edge = FALSE)
plot(d0)

## ----densityppp0, fig.keep = 'high', fig.cap = "As Figure \\@ref(fig:densityppp1), but without edge correction \\label{densityppp0}", fig.width = 3.75, fig.height = 3.5, echo = FALSE----
par(mai = c(0, 0, 0.2, 0.7))
plot(d0)

## ----relrisk-calc--------------------------------------------------------
rr = relrisk(ln, sigma = 250)

## ----relrisk, fig.keep = 'high', fig.cap = "Estimates of the spatially varying probability of each of the cell claases, conditional on there being cells.", fig.margin = FALSE, fig.width=7, fig.height=7----
plot(rr)

## ----checkConditional, echo=FALSE----------------------------------------
m = rr[[1]]$v
for(i in 2:length(rr)) m = m + rr[[i]]$v
stopifnot(all(is.na(m) | abs(m-1)<1e-6))

## ----Gestshow------------------------------------------------------------
gln = Gest(ln)
gln

## ----Gest, fig.keep = 'high', fig.cap = "Estimates of $G$, using three different edge effect corrections --which here happen to essentially lie on top of each other-- and the theoretical distribution for a homogenous Poisson process. \\label{Gest}", fig.width = 4.25, fig.height = 4.25, results="hide"----
library("RColorBrewer")
plot(gln, xlim = c(0, 10), lty = 1, col = brewer.pal(4, "Set1"))

## ----Linhom--------------------------------------------------------------
Lln = Linhom(subset(ln, marks == "T_cells"))
Lln

## ----Images-Lln, fig.keep = 'high', fig.cap = "Estimate of $L_{\\scriptsize \\mbox{inhom}}$, Equations (\\@ref(eq:kinhom)) and (\\@ref(eq:Lest)), of the T cell pattern. \\label{Images-Lln}", fig.width = 4, fig.height = 5.9, results = "hide"----
plot(Lln, lty = 1, col = brewer.pal(3, "Set1"))

## ----pcfdo---------------------------------------------------------------
pcfln = pcf(Kinhom(subset(ln, marks=="T_cells")))

## ----Images-pcf, fig.keep = 'high', fig.show = "hold", fig.cap = "Estimate of the pair correlation function, Equation (\\@ref(eq:pcf)), of the T cell pattern. \\label{Images-pcf}", fig.width=5, fig.height=5, results="hide"----
plot(pcfln, lty = 1)
plot(pcfln, lty = 1, xlim = c(0,10))

