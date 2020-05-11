## ----initialize, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
source("../chapter-setup.R"); chaptersetup("/Users/Susan/Courses/CUBook-html/CUBook/Chap16-Learning/Supervised.Rnw", "16")
knitr::opts_chunk$set(dev = 'png', dpi = 100, fig.margin = TRUE, fig.show = 'hold', fig.keep = 'none')

## ---- BuildWall, out.width = '50%', eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/BuildWall.png','images/EWall.png'))

## ----chap16-r-overfitting-1, fig.keep = 'high', fig.cap = "An example for **overfitting**: two regression lines are fit to data in the $(x, y)$-plane (black points). We can think of such a line as a rule that predicts the $y$-value, given an $x$-value. Both lines are smooth, but the fits differ in what is called their **bandwidth**, which intuitively can be interpreted their stiffness. The blue line seems overly keen to follow minor wiggles in the data, while the orange line captures the general trend but is less detailed. The effective number of parameters needed to describe the blue line is much higher than for the orange line. Also, if we were to obtain additional data, it is likely that the blue line would do a **worse** job than the orange line in modeling the new data. We\'ll formalize these concepts --training error and test set error-- later in this chapter. Although exemplified here with line fitting, the concept applies more generally to prediction models.", fig.width = 3, fig.height = 3, echo = FALSE----
ov = tibble(
  x = seq(0, 30, by = 1),
  y = 2 + 0.01 * x^2 + 0.1 * x + 2 * rnorm(length(x)))
ggplot(ov, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(span = 0.2, col = "dodgerblue3", se = FALSE) +
  geom_smooth(span = 0.8, col = "darkorange1", se = FALSE)

## ---- Supervised-fourtypes, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "In supervised learning, we assign two different roles to our variables. We have labeled the explanatory variables $X$ and the response variable(s) $Y$. There are also two different sets of observations: the training set $X_\\ell$ and $Y_\\ell$ and the test set $X_v$ and $Y_v$. (The subscripts refer to alternative names for the two sets: \"learning\" and \"validation\".)"----
knitr::include_graphics(c('images/fourquad.png'))

## ----diabetes------------------------------------------------------------
library("readr")
library("magrittr")
diabetes = read_csv("../data/diabetes.csv", col_names = TRUE)
diabetes
diabetes$group %<>% factor

## ----chap16-r-ldagroups-1, fig.keep = 'high', fig.cap = "We see already from the one-dimensional distributions that some of the individual variables could potentially predict which group a patient is more likely to belong to. Our goal will be to combine variables to improve over such one-dimensional prediction models.", fig.width = 3.5, fig.height = 7.5----
library("ggplot2")
library("reshape2")
ggplot(melt(diabetes, id.vars = c("id", "group")),
       aes(x = value, col = group)) +
 geom_density() + facet_wrap( ~variable, ncol = 1, scales = "free") +
 theme(legend.position = "bottom")

## ---- Supervised-cellshape, fig.margin = FALSE, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The data were images of $2\\times10^9$ nuclei from movies. The images were segmented to identify the nuclei, and numeric features were computed for each nucleus, corresponding to size, shape, brightness and lots of other more or less abstract quantitative summaries of the joint distribution of pixel intensities. From the features, the cells were classified into 16 different nuclei morphology classes, represented by the rows of the barplot. Representative images for each class are shown in black and white in the center column. The class frequencies, which are very unbalanced, are shown by the lengths of the bars."----
knitr::include_graphics(c('images/Neumann2010Fig1b.png'))

## ----chap16-r-scatterdiabetes-1, fig.keep = 'high', fig.cap = "Scatterplot of two of the variables in the `diabetes` data. Each point is a sample, and the color indicates the diabetes type as encoded in the `group` variable.", fig.width = 3.5, fig.height = 3----
ggdb = ggplot(mapping = aes(x = insulin, y = glutest)) +
  geom_point(aes(colour = group), data = diabetes)
ggdb

## ----ldaresults----------------------------------------------------------
library("MASS")
diabetes_lda = lda(group ~ insulin + glutest, data = diabetes)
diabetes_lda
ghat = predict(diabetes_lda)$class
table(ghat, diabetes$group)
mean(ghat != diabetes$group)

## ----make1Dgrid----------------------------------------------------------
make1Dgrid = function(x) {
  rg = grDevices::extendrange(x)
  seq(from = rg[1], to = rg[2], length.out = 100)
}

## ----diabetes_grid-1-----------------------------------------------------
diabetes_grid = with(diabetes,
  expand.grid(insulin = make1Dgrid(insulin),
              glutest = make1Dgrid(glutest)))

## ----diabetes_grid-2-----------------------------------------------------
diabetes_grid$ghat =
  predict(diabetes_lda, newdata = diabetes_grid)$class

## ----centers-------------------------------------------------------------
centers = diabetes_lda$means

## ----unitcircle----------------------------------------------------------
unitcircle = exp(1i * seq(0, 2*pi, length.out = 90)) %>%
          {cbind(Re(.), Im(.))}
ellipse = unitcircle %*% solve(diabetes_lda$scaling)

## ----ellipses------------------------------------------------------------
ellipses = lapply(seq_len(nrow(centers)), function(i) {
  (ellipse +
   matrix(centers[i, ], byrow = TRUE,
          ncol = ncol(centers), nrow = nrow(ellipse))) %>%
     cbind(group = i)
}) %>% do.call(rbind, .) %>% data.frame
ellipses$group %<>% factor

## ----chap16-r-modeldiabetes-1, fig.keep = 'high', fig.cap = "As Figure \\@ref(fig:chap16-r-scatterdiabetes-1), with the classification regions from the LDA model shown. The three ellipses represent the class centers and the covariance matrix of the LDA model; note that there is only one covariance matrix, which is the same for all three classes. Therefore also the sizes and orientations of the ellipses are the same for the three classes, only their centers differ. They represent contours of equal class membership probability.", fig.width = 5, fig.height = 4----
ggdb + geom_raster(aes(fill = ghat),
            data = diabetes_grid, alpha = 0.25, interpolate = TRUE) +
    geom_point(data = as_tibble(centers), pch = "+", size = 8) +
    geom_path(aes(colour = group), data = ellipses) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))

## ----chap16-r-diabetes-lda-uniform-prior-1, fig.keep = 'high', fig.cap = "As Figure \\@ref(fig:chap16-r-modeldiabetes-1), but with uniform class priors.", fig.width = 5, fig.height = 4----
diabetes_up = lda(group ~ insulin + glutest, data = diabetes,
  prior = with(diabetes, rep(1/nlevels(group), nlevels(group))))

diabetes_grid$ghat_up =
  predict(diabetes_up, newdata = diabetes_grid)$class

stopifnot(all.equal(diabetes_up$means, diabetes_lda$means))

ellipse_up  = unitcircle %*% solve(diabetes_up$scaling)
ellipses_up = lapply(seq_len(nrow(centers)), function(i) {
  (ellipse_up +
   matrix(centers[i, ], byrow = TRUE,
          ncol = ncol(centers), nrow = nrow(ellipse_up))) %>%
     cbind(group = i)
}) %>% do.call(rbind, .) %>% data.frame
ellipses_up$group %<>% factor

ggdb + geom_raster(aes(fill = ghat_up),
            data = diabetes_grid, alpha = 0.4, interpolate = TRUE) +
    geom_point(data = data.frame(centers), pch = "+", size = 8) +
    geom_path(aes(colour = group), data = ellipses_up) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))

## ----all5diab------------------------------------------------------------
diabetes_lda5 = lda(group ~ relwt + glufast + glutest +
           steady + insulin, data = diabetes)
diabetes_lda5
ghat5 = predict(diabetes_lda5)$class
table(ghat5, diabetes$group)
mean(ghat5 != diabetes$group)

## ----loadHiiragi2--------------------------------------------------------
library("Hiiragi2013")
data("x")
probes = c("1426642_at", "1418765_at", "1418864_at", "1416564_at")
embryoCells = t(Biobase::exprs(x)[probes, ]) %>% as_tibble %>%
  mutate(Embryonic.day = x$Embryonic.day) %>%
  dplyr::filter(x$genotype == "WT")

## ----annoHiiragi, warning = FALSE----------------------------------------
annotation(x)
library("mouse4302.db")
anno = AnnotationDbi::select(mouse4302.db, keys = probes,
         columns = c("SYMBOL", "GENENAME"))
anno
mt = match(anno$PROBEID, colnames(embryoCells))
colnames(embryoCells)[mt] = anno$SYMBOL

## ----assertprobeid, echo = FALSE-----------------------------------------
stopifnot(!any(is.na(mt)))

## ----chap16-HiiragiFourGenesPairs-1, fig.keep = 'high', fig.cap = "Expression values of the discriminating genes, with the prediction target Embryonic.day shown by color.", fig.margin = FALSE, fig.width = 6, fig.height = 6----
library("GGally")
ggpairs(embryoCells, mapping = aes(col = Embryonic.day),
  columns = anno$SYMBOL, upper = list(continuous = "points"))

## ----ldacells, fig.width=8, fig.height=4---------------------------------
ec_lda = lda(Embryonic.day ~ Fn1 + Timd2 + Gata4 + Sox7,
             data = embryoCells)
round(ec_lda$scaling, 1)

## ----chap16-r-edcontour-1, fig.keep = 'high', fig.cap = "LDA classification regions for Embryonic.day.", fig.width = 4.5, fig.height = 3.5----
ec_rot = predict(ec_lda)$x %>% as_tibble %>%
           mutate(ed = embryoCells$Embryonic.day)
ec_lda2 = lda(ec_rot[, 1:2], predict(ec_lda)$class)
ec_grid = with(ec_rot, expand.grid(
  LD1 = make1Dgrid(LD1),
  LD2 = make1Dgrid(LD2)))
ec_grid$edhat = predict(ec_lda2, newdata = ec_grid)$class
ggplot() +
  geom_point(aes(x = LD1, y = LD2, colour = ed), data = ec_rot) +
  geom_raster(aes(x = LD1, y = LD2, fill = edhat),
            data = ec_grid, alpha = 0.4, interpolate = TRUE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed()

## ----chap16-r-qdamouse-1, fig.keep = 'high', fig.cap = "QDA for the mouse cell data. Shown are all pairwise plots of the four features. In each plot, the other two features are set to the median.", fig.margin = FALSE, fig.width = 9, fig.height = 9----
library("gridExtra")

ec_qda = qda(Embryonic.day ~ Fn1 + Timd2 + Gata4 + Sox7,
             data = embryoCells)

variables = colnames(ec_qda$means)
pairs = combn(variables, 2)
lapply(seq_len(ncol(pairs)), function(i) {
  grid = with(embryoCells,
    expand.grid(x = make1Dgrid(get(pairs[1, i])),
                y = make1Dgrid(get(pairs[2, i])))) %>%
    `colnames<-`(pairs[, i])

  for (v in setdiff(variables, pairs[, i]))
    grid[[v]] = median(embryoCells[[v]])

  grid$edhat = predict(ec_qda, newdata = grid)$class

  ggplot() + geom_point(
      aes_string(x = pairs[1, i], y = pairs[2, i],
      colour = "Embryonic.day"), data = embryoCells) +
    geom_raster(
      aes_string(x = pairs[1, i], y = pairs[2, i], fill = "edhat"),
      data = grid, alpha = 0.4, interpolate = TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    if (i != ncol(pairs)) theme(legend.position = "none")
}) %>% grid.arrange(grobs = ., ncol = 2)

## ----ladallvariables, warning = TRUE, error = TRUE, results = "hide"-----
lda(t(Biobase::exprs(x))[, 1:1000], x$Embryonic.day)
qda(t(Biobase::exprs(x))[, 1:1000], x$Embryonic.day)

## ----chap16-r-learnbyheart-1, fig.keep = 'high', fig.cap = "Misclassification rate of LDA applied to random data. While the number of observations `n` is held constant (at 20), we are increasing the number of features `p` starting from 2 up to 21. The misclassification rate becomes almost zero as `p` approaches 20. The LDA model becomes so elaborate and over-parameterized that it manages to learn the random labels \"by heart\". (As `p` becomes even larger, the \"performance\" degrades again somewhat, apparently due to numerical properties of the `lda` implementation used here.)", warning = FALSE, fig.width = 3, fig.height = 3----
library("dplyr")
p = 2:21
n = 20

mcl = lapply(p, function(pp) {
  replicate(100, {
    xmat = matrix(rnorm(n * pp), nrow = n)
    resp = sample(c("apple", "orange"), n, replace = TRUE)
    fit  = lda(xmat, resp)
    pred = predict(fit)$class
    mean(pred != resp)
  }) %>% mean %>% tibble(mcl = ., p = pp)
}) %>% bind_rows

ggplot(mcl, aes(x = p, y = mcl)) + geom_line() + geom_point() +
  ylab("Misclassification rate")

## ---- book-chunk-1, eval = TRUE, echo = FALSE, fig.keep = 'high'---------
knitr::include_graphics('images/book_icon.png', dpi = 400)

## ----chap16-r-mclcv-1, fig.keep = 'high', fig.cap = "Cross-validation: the misclassification rate of LDA applied to random data, when evaluated on test data that were not used for learning, hovers around 0.5 independent of `p`. The misclassification rate on the training data is also shown. It behaves similar to what we already saw in Figure \\@ref(fig:chap16-r-learnbyheart-1).", warning = FALSE----
estimate_mcl_loocv = function(x, resp) {
  vapply(seq_len(nrow(x)), function(i) {
    fit  = lda(x[-i, ], resp[-i])
    ptrn = predict(fit, newdata = x[-i,, drop = FALSE])$class
    ptst = predict(fit, newdata = x[ i,, drop = FALSE])$class
    c(train = mean(ptrn != resp[-i]), test = (ptst != resp[i]))
  }, FUN.VALUE = numeric(2)) %>% rowMeans %>% t %>% as_tibble
}

xmat = matrix(rnorm(n * last(p)), nrow = n)
resp = sample(c("apple", "orange"), n, replace = TRUE)

mcl = lapply(p, function(k) {
  estimate_mcl_loocv(xmat[, 1:k], resp)
}) %>% bind_rows %>% data.frame(p) %>% melt(id.var = "p")

ggplot(mcl, aes(x = p, y = value, col = variable)) + geom_line() +
  geom_point() + ylab("Misclassification rate")

## ----chap16-r-curseofdim, fig.keep = 'high', fig.cap = "As we increase the number of features included in the model, the misclassification rate initially improves; as we start including more and more irrelevant features, it increases again, as we are fitting noise.", warning = FALSE, fig.width = 3.5, fig.height = 3----
p   = 2:20
mcl = replicate(100, {
  xmat = matrix(rnorm(n * last(p)), nrow = n)
  resp = sample(c("apple", "orange"), n, replace = TRUE)
  xmat[, 1:6] = xmat[, 1:6] + as.integer(factor(resp))

  lapply(p, function(k) {
    estimate_mcl_loocv(xmat[, 1:k], resp)
  }) %>% bind_rows %>% cbind(p = p) %>% melt(id.var = "p")
}, simplify = FALSE) %>% bind_rows

mcl =  group_by(mcl, p, variable) %>%
   summarise(value = mean(value))

ggplot(mcl, aes(x = p, y = value, col = variable)) + geom_line() +
   geom_point() + ylab("Misclassification rate")

## ---- chap16-BiasVarianceTradeoff, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Idealized version of Figure \\@ref(fig:chap16-r-curseofdim), from @HastieTibshiraniFriedman. A recurrent goal in machine learning is finding the sweet spot in the variance--bias trade-off."----
knitr::include_graphics(c('images/BiasVarianceTradeoff.png'))

## ----chap16-r-cursedimans1-1, fig.keep = 'high', fig.cap = "Side length of a $p$-dimensional hybercube expected to contain 10 points out of 1 million uniformly distributed ones, as a function of the $p$. While for $p=1$, this length is conveniently small, namely $10/10^6=10^{-5}$, for larger $p$ it approaches 1, i.\\,e., becomes the same as the range of each the features. This means that a \"local neighborhood\" of 10 points encompasses almost the same data range as the whole dataset.", fig.width = 3.5, fig.height = 2.5----
sideLength = function(p, pointDensity = 1e6, pointsNeeded = 10)
   (pointsNeeded / pointDensity) ^ (1 / p)
ggplot(tibble(p = 1:400, sideLength = sideLength(p)),
       aes(x = p, y = sideLength)) + geom_line(col = "red") +
  geom_hline(aes(yintercept = 1), linetype = 2)

## ----chap16-r-cursedimans2-1, fig.keep = 'high', fig.cap = "Fraction of a unit cube\'s total volume that is in its \"shell\" (here operationalised as those points that are closer than 0.01 to its surface) as a function of the dimension $p$.", fig.width = 3.5, fig.height = 2.5----
tibble(
  p = 1:400,
  volOuterCube = 1 ^ p,
  volInnerCube = 0.98 ^ p,  # 0.98 = 1 - 2 * 0.01
  `V(shell)` = volOuterCube - volInnerCube) %>%
ggplot(aes(x = p, y =`V(shell)`)) + geom_line(col = "blue")

## ----chap16-r-cursedimans3-1, fig.keep = 'high', fig.cap = "Coefficient of variation (CV) of the distance between randomly picked points in the unit hypercube, as a function of the dimension. As the dimension increases, everybody is equally far away from everyone else: there is almost no variation in the distances any more.", fig.width = 3.5, fig.height = 2.5----
n = 1000
df = tibble(
  p = round(10 ^ seq(0, 4, by = 0.25)),
  cv = vapply(p, function(k) {
    x1 = matrix(runif(k * n), nrow = n)
    x2 = matrix(runif(k * n), nrow = n)
    d = sqrt(rowSums((x1 - x2)^2))
    sd(d) / mean(d)
  }, FUN.VALUE = numeric(1)))
ggplot(df, aes(x = log10(p), y = cv)) + geom_line(col = "orange") +
  geom_point()

## ----confusiontable, eval=FALSE------------------------------------------
## table(truth, response)


## ---- Supervised-bullseye, out.width = '50%', eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "In the upper bull\'s eye, the estimates are systematically off target, but in a quite reproducible manner. The green segment represents the bias. In the lower bull\'s eye, the estimates are not biased, as they are centered in the right place, however they have high variance. We can distinguish the two scenarios since we see the result from many shots. If we only had one shot and missed the bull\'s eye, we could not easily know whether that\'s because of bias or variance."----
knitr::include_graphics(c('images/TargetBias.png','images/TargetVariance.png'))


## ----colon1, results = "hide"--------------------------------------------
library("ExperimentHub")
eh = ExperimentHub()
zeller = eh[["EH361"]]

## ----colon1b-------------------------------------------------------------
table(zeller$disease)

## ----colon2--------------------------------------------------------------
zellerNC = zeller[, zeller$disease %in% c("n", "cancer")]

## ----ehzellertest, echo = FALSE------------------------------------------
stopifnot(is.numeric(Biobase::exprs(zellerNC)), !any(is.na(Biobase::exprs(zellerNC))))

## ----zellerpData-1-------------------------------------------------------
pData(zellerNC)[ sample(ncol(zellerNC), 3), ]


## ----zellerpData-2-------------------------------------------------------
formatfn = function(x)
   gsub("|", "| ", x, fixed = TRUE) %>% lapply(strwrap)

rownames(zellerNC)[1:4]
rownames(zellerNC)[nrow(zellerNC) + (-2:0)] %>% formatfn

## ----chap16-r-zellerHist-1, fig.keep = 'high', fig.cap = "Histograms of the distributions for two randomly selected features. The distributions are highly skewed, with many zero values and a thin, long tail of non-zero values.", fig.width = 3, fig.height = 4----
ggplot(melt(Biobase::exprs(zellerNC)[c(510, 527), ]), aes(x = value)) +
    geom_histogram(bins = 25) +
    facet_wrap( ~ Var1, ncol = 1, scales = "free")

## ----glmnet--------------------------------------------------------------
library("glmnet")
glmfit = glmnet(x = t(Biobase::exprs(zellerNC)),
                y = factor(zellerNC$disease),
                family = "binomial")

## ----colonPred-----------------------------------------------------------
predTrsf = predict(glmfit, newx = t(Biobase::exprs(zellerNC)),
                   type = "class", s = 0.04)
table(predTrsf, zellerNC$disease)

## ----chap16-r-plotglmfit-1, fig.keep = 'high', fig.cap = "Regularization paths for `glmfit`.", fig.width = 3.6, fig.height = 3.2, echo = -1----
par(mai = c(0.5, 0.5, 0.575, 0.05))
plot(glmfit, col = brewer.pal(8, "Dark2"), lwd = sqrt(3), ylab = "")


## ----chap16-r-colonCV-1, fig.keep = 'high', fig.cap = "Diagnostic plot for `cv.glmnet`: shown is a measure of cross-validated prediction performance, the deviance, as a function of $\\lambda$. The dashed vertical lines show `lambda.min` and `lambda.1se`.", fig.width = 4, fig.height = 4----
cvglmfit = cv.glmnet(x = t(Biobase::exprs(zellerNC)),
                     y = factor(zellerNC$disease),
                     family = "binomial")
plot(cvglmfit)

## ----lambda.min----------------------------------------------------------
cvglmfit$lambda.min

## ----lambda.1se----------------------------------------------------------
cvglmfit$lambda.1se

## ----predictwithlambda1se------------------------------------------------
s0 = cvglmfit$lambda.1se
predict(glmfit, newx = t(Biobase::exprs(zellerNC)),type = "class", s = s0) %>%
    table(zellerNC$disease)

## ----zellercoef----------------------------------------------------------
coefs = coef(glmfit)[, which.min(abs(glmfit$lambda - s0))]
topthree = order(abs(coefs), decreasing = TRUE)[1:3]
as.vector(coefs[topthree])
formatfn(names(coefs)[topthree])

## ----chap16-r-colonCVTrsf-1, fig.keep = 'high', fig.cap = "like Figure \\@ref(fig:chap16-r-colonCV-1), but using an $\\text{asinh}$ transformation of the data.", fig.width = 4, fig.height = 4----
cv.glmnet(x = t(asinh(Biobase::exprs(zellerNC))),
          y = factor(zellerNC$disease),
          family = "binomial") %>% plot

## ----chap16-r-mousecvglmfit-1, fig.keep = 'high', fig.cap = "Cross-validated misclassification error versus penalty parameter for the mouse cells data.", fig.width = 4, fig.height = 4----
sx = x[, x$Embryonic.day == "E3.25"]
embryoCellsClassifier = cv.glmnet(t(Biobase::exprs(sx)), sx$genotype,
                family = "binomial", type.measure = "class")
plot(embryoCellsClassifier)

## ----checkclaimMouseCellsClassifier, echo = FALSE------------------------
stopifnot(sum((diff(embryoCellsClassifier$cvm) * diff(embryoCellsClassifier$lambda)) < 0) <= 2)

## ----chap16-r-mousecellsrowttst-1, fig.keep = 'high', fig.cap = "Histogram of p-values for the per-feature $t$-tests between genotypes in the E3.25 cells.", fig.width = 4, fig.height = 2.5----
mouse_de = rowttests(sx, "genotype")
ggplot(mouse_de, aes(x = p.value)) +
  geom_histogram(boundary = 0, breaks = seq(0, 1, by = 0.01))

## ----mousecellsnn1-------------------------------------------------------
dists = as.matrix(dist(scale(t(Biobase::exprs(x)))))
diag(dists) = +Inf

## ----mousecellsnn2-------------------------------------------------------
nn = sapply(seq_len(ncol(dists)), function(i) which.min(dists[, i]))
table(x$sampleGroup, x$sampleGroup[nn]) %>% `colnames<-`(NULL)

## ----caret1, message = FALSE---------------------------------------------
library("caret")
caretMethods = names(getModelInfo())
head(caretMethods, 8)
length(caretMethods)

## ----caret2--------------------------------------------------------------
getModelInfo("nnet", regex = FALSE)[[1]]$parameter

## ----caret3, results = "hide", message  = FALSE--------------------------
trnCtrl = trainControl(
  method = "repeatedcv",
  repeats = 3,
  classProbs = TRUE)
tuneGrid = expand.grid(
  size = c(2, 4, 8),
  decay = c(0, 1e-2, 1e-1))
nnfit = train(
  Embryonic.day ~ Fn1 + Timd2 + Gata4 + Sox7,
  data = embryoCells,
  method = "nnet",
  tuneGrid  = tuneGrid,
  trControl = trnCtrl,
  metric = "Accuracy")

## ----ML-nnfit, fig.keep = 'high', fig.cap = "Parameter tuning of the neural net by cross-validation.", fig.width = 3.75, fig.height = 4.25----
nnfit
plot(nnfit)
predict(nnfit) %>% head(10)

## ----kernelsvm, eval = FALSE, echo = FALSE-------------------------------
## library("kernlab")
## kfunction= function(linear =0, quadratic=0)
## {  k = function (v,w){ linear*sum((v)*(w)) + quadratic*sum((v^2)*(w^2))}
##   class(k) = "kernel"
##   return(k) }
## subx=subx[,2:3]
## svp = ksvm(subx,dftxy$tg,type="C-svc",C = 100, kernel=kfunction(1,0),scaled=c())
## plot(c(min(subx[,1]), max(subx[,1])),c(min(subx[,2]), max(subx[,2])),
##             type='n',xlab='x1',ylab='x2')
## ymat = ymatrix(svp)
## points(subx[-SVindex(svp),1], subx[-SVindex(svp),2],
##          pch = ifelse(ymat[-SVindex(svp)] < 0, 2, 1))
## points(subx[SVindex(svp),1], subx[SVindex(svp),2],
##          pch = ifelse(ymat[SVindex(svp)] < 0, 17, 16))
## 
## # Extract w and b from the model
## w = colSums(coef(svp)[[1]] * subx[SVindex(svp),])
## b = b(svp)
## # Draw the lines
## abline(b/w[2],-w[1]/w[2])
## abline((b+1)/w[2],-w[1]/w[2],lty=2)
## abline((b-1)/w[2],-w[1]/w[2],lty=2)

