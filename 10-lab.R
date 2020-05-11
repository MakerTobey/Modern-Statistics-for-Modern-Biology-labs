pkgs_needed = c("MASS","ExperimentHub", "tidyverse","glmnet",
                "RColorBrewer","caret", "magrittr")
letsinstall = setdiff(pkgs_needed, installed.packages()) 
if (length(letsinstall) > 0) {
  BiocManager::install(setdiff(pkgs_needed, installed.packages()))
}

library("tidyverse")
library("MASS")
library("ExperimentHub")
library("glmnet")
library("RColorBrewer")
library("caret")
library("magrittr")

# Linear discrimination
# Diabetes data set
diabetes = read_csv(url("http://web.stanford.edu/class/bios221/data/diabetes.csv"), 
                    col_names = TRUE)
diabetes
diabetes$group %<>% factor

diabetes.long = gather(diabetes, variable, value, -c(id,group))
ggplot(diabetes.long, aes(x = value, col = group)) +
  geom_density() + facet_wrap( ~ variable, ncol = 2, scales = "free") 

# LDA

ggdb = ggplot(mapping = aes(x = insulin, y = glutest)) +
  geom_point(aes(colour = group), data = diabetes)
ggdb

diabetes_lda = lda(group ~ insulin + glutest, data = diabetes) # from MASS package
diabetes_lda

ghat = predict(diabetes_lda)$class
table(predicted = ghat, truth = diabetes$group)

mean(ghat != diabetes$group)

make1Dgrid = function(x) {
  rg = grDevices::extendrange(x)
  seq(from = rg[1], to = rg[2], length.out = 100)
}

diabetes_grid = with(diabetes,
                     expand.grid(insulin = make1Dgrid(insulin),
                                 glutest = make1Dgrid(glutest)))

diabetes_grid$ghat = predict(diabetes_lda, newdata = diabetes_grid)$class

centers = diabetes_lda$means

ggdb + geom_raster(
  data = diabetes_grid, alpha = 0.25,
  aes(fill = ghat),
  interpolate = TRUE ) +
  geom_point(data = as_tibble(centers), pch = "+", size = 8) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))

diabetes_lda5 = lda(group ~ relwt + glufast + glutest +
                      steady + insulin, data = diabetes)
diabetes_lda5

ghat5 = predict(diabetes_lda5)$class
table(ghat5, diabetes$group)

mean(ghat5 != diabetes$group)


# Cross-validation

estimate_mcl_loocv = function(x, resp) {
  vapply(seq_len(nrow(x)), function(i) {
    fit  = lda(x[-i, ], resp[-i])
    ptrn = predict(fit, newdata = x[-i,, drop = FALSE])$class
    ptst = predict(fit, newdata = x[ i,, drop = FALSE])$class
    c(train = mean(ptrn != resp[-i]), test = (ptst != resp[i]))
  }, FUN.VALUE = numeric(2)) %>% rowMeans %>% t %>% as_tibble
}

n = 20
p   = 2:20
# This will take a while
mcl = replicate(100, {
  xmat = matrix(rnorm(n * last(p)), nrow = n)
  resp = sample(c("apple", "orange"), n, replace = TRUE)
  xmat[, 1:6] = xmat[, 1:6] + as.integer(factor(resp))
  
  lapply(p, function(k) {
    estimate_mcl_loocv(xmat[, 1:k], resp)
  }) %>% bind_rows %>% cbind(p = p) %>% gather(variable, value, -p)
}, simplify = FALSE) %>% bind_rows 

mcl =  group_by(mcl, p, variable) %>%
  summarise(value = mean(value))

ggplot(mcl, aes(x = p, y = value, col = variable)) + geom_line() +
  geom_point() + ylab("Misclassification rate")


# Variance-bias trade-off

library("ExperimentHub")
eh = ExperimentHub()

zeller = eh[["EH361"]]

table(zeller$disease)

zellerNC = zeller[, zeller$disease %in% c("n", "cancer")]

pData(zellerNC)[ sample(ncol(zellerNC), 3), ]

rownames(zellerNC)[1:4]

rownames(zellerNC)[nrow(zellerNC) + (-2:0)]

tidy_zeller_subset <- as.data.frame(t(exprs(zellerNC)[c(510, 527), ])) %>% 
  mutate(Var2 = colnames(exprs(zellerNC))) %>% 
  gather(Var1, value, - Var2)

ggplot(tidy_zeller_subset, aes(x = value)) +
  geom_histogram(bins = 25) +
  facet_wrap( ~ Var1, ncol = 1, scales = "free")

library("glmnet")
glmfit = glmnet(x = t(exprs(zellerNC)),
                y = factor(zellerNC$disease),
                family = "binomial")

pred = predict(glmfit, newx = t(exprs(zellerNC)), type = "class", s = 0.04)
confusion_table = table(predicted = pred, truth = zellerNC$disease)
confusion_table

plot(glmfit, xvar = "norm", col = RColorBrewer::brewer.pal(12, "Set3"), lwd = sqrt(3))

fitted_beta = coef(glmfit, s=0.04)
sum(abs(fitted_beta))

sum(abs(coef(glmfit, s=0.1)))

set.seed(0xdada2)
cvglmfit = cv.glmnet(x = t(exprs(zellerNC)),
                     y = factor(zellerNC$disease),
                     family = "binomial")
plot(cvglmfit)

cvglmfit$lambda.min

s0 = cvglmfit$lambda.1se
s0

coefs = coef(glmfit)[, which.min(abs(glmfit$lambda - s0))]

