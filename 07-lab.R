pkgs_needed = c("tidyverse","GGally", "factoextra", "ade4")
BiocManager::install(setdiff(pkgs_needed, installed.packages()))
library("tidyverse")

turtles = read.table(url("https://web.stanford.edu/class/bios221/data/PaintedTurtles.txt"),
                     header=TRUE)
head(turtles)

download.file(url = "https://web.stanford.edu/class/bios221/data/athletes.RData",
              destfile = "athletes.RData",mode = "wb")
load("athletes.RData")
athletes[1:3,]

library("GGally")
nggpairs(dplyr::select(turtles, -matches("sex")), 
        axisLabels = "none")

library("pheatmap")
pheatmap(cor(athletes), cell.width = 10, cell.height = 10)

scaledTurtles = cbind(scale(dplyr::select(turtles, -matches("sex"))), 
                      dplyr::select(turtles, matches("sex")))
ggplot(scaledTurtles, aes(x = width, y = height, group = sex)) +
  geom_point(aes(color = sex)) + coord_fixed()

test <- dplyr::select(turtles, -matches("sex"))
test  

library("ggplot2")
athletes_sc = scale(athletes)
n = nrow(athletes_sc)
athletes_sc = data.frame(athletes_sc)

p = ggplot(athletes_sc, aes(x = weight,y = disc)) +
  geom_point(size = 2, shape = 21)
p + geom_point(aes(y = rep(0, n)), colour = "red") +
  geom_segment(aes(xend = weight, yend = rep(0, n)), linetype = "dashed") +
  coord_fixed()

#q = ggplot(athletes_sc, aes(x = weight,y = disc)) +
#  geom_point(aes(y = rep(0, n)), colour = "red") +
#  geom_errorbar(aes(ymin=len-var, ymax=len+(var(len)/2), width=.2, position=position_dodge(.9))

reg1 = lm(disc ~ weight, data = athletes_sc)
a = reg1$coefficients[1] # Intercept
b = reg1$coefficients[2] # slope
p + geom_abline(intercept = a, slope = b, col = "blue", lwd = 1.5) +
  geom_segment(aes(xend = weight, yend = reg1$fitted),
               colour = "red", arrow = arrow(length = unit(0.15,"cm"))) +  
  coord_fixed()

X = cbind(athletes_sc$disc, athletes_sc$weight)
svda = svd(X)
pc = X %*% svda$v[, 1] %*% t(svda$v[, 1])
bp = svda$v[2, 1] / svda$v[1, 1]
ap = mean(pc[, 2]) - bp * mean(pc[, 1])

p + geom_segment(xend = pc[,1], yend = pc[,2]) + 
  geom_abline(intercept = ap, slope = bp, col = "purple", lwd = 1.5) + 
  coord_fixed()

ppdf = data.frame(
  PC1n = -svda$u[,1] * svda$d[1], 
  PC2n =  svda$u[,2] * svda$d[2])

ggplot(ppdf,aes(x = PC1n, y=PC2n)) + geom_point() + ylab("PC2 ") +
  geom_hline(yintercept = 0, color = "purple", lwd = 1.5, alpha = 0.5) +
  geom_point(aes(x = PC1n, y = 0),color = "red") + xlab("PC1 ")+
  geom_segment(aes(xend = PC1n,yend = 0), color = "red") +
  coord_fixed()

pca_athletes = princomp(X)
svda$v
pca_athletes$loadings
head(pca_athletes$scores)
head(svda$u %*% diag(svda$d))
c(sd(ppdf$PC1n), sd(ppdf$PC2n))

turtles3var = dplyr::select(turtles, -matches("sex"))
apply(turtles3var, 2, mean)
apply(turtles3var, 2, var)
cor(turtles3var)
library("factoextra")
pca1 = princomp(scale(turtles3var))
# or alternatively:
# pca1 = ade4::dudi.pca(scale(turtles3var), scannf = FALSE)

fviz_eig(pca1, geom = "bar", width = 0.4)

fviz_pca_biplot(pca1, label = "var", col.ind = turtles[,1]) 

library("ade4")
pca.ath = dudi.pca(athletes_sc, scannf = FALSE)
pca.ath$eig

fviz_eig(pca.ath, geom = "bar", bar_width = 0.3) + ggtitle("")

fviz_pca_var(pca.ath, col.circle = "black") + ggtitle("")

athletes_2 = athletes
m_columns = grep("^m[[:digit:]]+$", colnames(athletes_2))
athletes_2[, m_columns] = athletes_2[, m_columns] ^(-1)
cor(athletes_2) %>% round(1)

pcan.ath = dudi.pca(scale(athletes_2), nf = 2, scannf = FALSE)
pcan.ath$eig

fviz_pca_var(pcan.ath, col.circle="black") + ggtitle("")

fviz_pca_ind(pcan.ath) + ggtitle("") + ylim(c(-2.5,5.7)) + coord_fixed()
