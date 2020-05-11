pkgs_needed = c("EBImage","magrittr","tibble","ggplot2","genefilter",
                "GGally", "MSMB")
letsinstall = setdiff(pkgs_needed, installed.packages()) 
if (length(letsinstall) > 0) {
  BiocManager::install(letsinstall)
}

library("EBImage")
library("magrittr")
library("tibble")
library("ggplot2")
library("genefilter")
library("GGally")
knitr::opts_chunk$set(echo = TRUE)

imagefile = system.file("images", "mosquito.png",
                        package = "MSMB")
mosq = readImage(imagefile)

EBImage::display(mosq)

EBImage::display(mosq, method = "raster")
text(x = 85, y = 800, label = "A mosquito",
     adj = 0, col = "orange", cex = 1.5)

imagefile = system.file("images", "hiv.png", package = "MSMB")
hivc = readImage(imagefile)
EBImage::display(hivc, method = "raster")

nuc = readImage(system.file("images", "nuclei.tif", package = "EBImage"))
EBImage::display(1 - nuc, method = "raster", all = TRUE)

EBImage::display(1 - nuc, method = "raster", frame = 2)

writeImage(hivc, "hivc.jpeg", quality = 85)

class(mosq)
dim(mosq)
hist(mosq)
imageData(mosq)[1:3, 1:6]
dim(imageData(nuc))

mosqinv = 1 - mosq
EBImage::display(mosqinv, method = "raster")
EBImage::display(normalize(-mosq), method = "raster")
mosqcont = mosq * 3
EBImage::display(mosqcont, method = "raster")
mosqexp = mosq ^ (1/3)
EBImage::display(mosqexp, method = "raster")
mosqcrop   = mosq[100:438, 112:550]
EBImage::display(mosqcrop, method = "raster")
mosqthresh = mosq > 0.5
EBImage::display(mosqthresh, method = "raster")
mosqtransp = transpose(mosq)
EBImage::display(mosqtransp, method = "raster")

mosqrot   = EBImage::rotate(mosq, angle = 30)
mosqshift = translate(mosq, v = c(40, 70))
mosqflip  = flip(mosq)
mosqflop  = flop(mosq)
EBImage::display(mosqflip, method = "raster")

#linear filters
imagefiles = system.file("images", 
                         c("image-DAPI.tif", "image-FITC.tif", "image-Cy3.tif"),
                         package="MSMB")
cells = readImage(imagefiles)
cells

apply(cells, 3, range)

cells[,,1] = 32 * cells[,,1]
cells[,,2:3] = 16 * cells[,,2:3]
apply(cells, 3, range)

w = makeBrush(size = 51, shape = "gaussian", sigma = 7)
tibble(w = w[(nrow(w)+1)/2, ]) %>%
  ggplot(aes(y = w, x = seq(along = w))) + geom_point()

nucSmooth = filter2(getFrame(cells, 1), w)
EBImage::display(nucSmooth, method = "raster")

cellsSmooth = Image(dim = dim(cells))
sigma = c(1, 3, 3)
for(i in seq_along(sigma)) {
  cellsSmooth[, ,i] = filter2( 
    cells[,,i],
    filter = makeBrush(size = 51, shape = "gaussian",
                       sigma = sigma[i])
  )
}
EBImage::display(cellsSmooth, method = "raster", all = TRUE)

disc = makeBrush(21, "disc")
disc = disc / sum(disc)
offset = 0.02
nucThresh = (cellsSmooth[,,1] - filter2( cellsSmooth[,,1], disc ) > offset)
EBImage::display(nucThresh, method = "raster")

#morphological operations 
nucOpened = EBImage::opening(nucThresh, kern = makeBrush(3, shape = "disc"))
EBImage::display(nucOpened, method = "raster")

#segmentation
nucSeed = bwlabel(nucOpened)

EBImage::display(colorLabels(nucSeed), method = "raster")

nucMask = cellsSmooth[,, 1] - filter2(cellsSmooth[ , , 1], disc) > 0

nucMask = fillHull(nucMask)

nuclei = propagate(cellsSmooth[,,1], nucSeed, mask = nucMask)
EBImage::display(nuclei,method = "raster")

#Segmenting the cell bodies
hist( log(cellsSmooth[,,3]) )
hist( log(cellsSmooth[,,3]), xlim = -c(3.6, 3.1), breaks = 300)

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

hist(log(cellsSmooth[,, 3]), xlim = -c(3.6, 3.1), breaks = 300)
abline(v = cellBg[c("loc", "thr"), 3], col = c("brown", "red"))

cytoplasmMask = (cellsSmooth[,,2] > exp(cellBg["thr", 2])) |
  nuclei | (cellsSmooth[,,3] > exp(cellBg["thr", 3]))
EBImage::display(cytoplasmMask, method = "raster")

cellbodies = propagate(x = cellsSmooth[ , , 3], seeds = nuclei,
                       lambda = 1.0e-2, mask = cytoplasmMask)
EBImage::display(colorLabels(cellbodies), method = "raster")

nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(cells[, , 1]), col = "#ffff00")
EBImage::display(nucSegOnNuc, method = "raster")

cellsColor = rgbImage(red   = cells[,, 3], 
                      green = cells[,, 2], 
                      blue  = cells[,, 1])
nucSegOnAll  = paintObjects(nuclei, tgt = cellsColor, col = "#ffff00")
EBImage::display(nucSegOnAll, method = "raster")

cellSegOnAll = paintObjects(cellbodies, tgt = nucSegOnAll, col = "#ff0080")
EBImage::display(cellSegOnAll, method = "raster")

#Feature extraction
meanNucInt       = tapply(cells[,,1], nuclei, mean)
meanActIntInNuc  = tapply(cells[,,3], nuclei, mean)
meanActIntInCell = tapply(cells[,,3], cellbodies, mean)

library("GGally")
ggpairs(tibble(meanNucInt, meanActIntInNuc, meanActIntInCell))

F1 = computeFeatures(nuclei,     cells[,,1], xname = "nuc",
                     refnames = "nuc")
F2 = computeFeatures(cellbodies, cells[,,2], xname = "cell",
                     refnames = "tub")
F3 = computeFeatures(cellbodies, cells[,,3], xname = "cell",
                     refnames = "act")
dim(F1)
F1[1:3, 1:5]
