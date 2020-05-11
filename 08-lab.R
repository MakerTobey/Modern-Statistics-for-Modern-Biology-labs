pkgs_needed = c("tidyverse","genefilter","ggrepel",
                "igraph","statnet", "ggnetwork", "rworldmap", "intergraph",
                "PMA", "phyloseq","vegan", "ade4", "impute")
BiocManager::install(setdiff(pkgs_needed, installed.packages()))

knitr::opts_chunk$set(echo = TRUE, fig.dim = c((1+sqrt(5))/2, 1) * 5, cache = TRUE, autodep = TRUE)

library("tidyverse")
library("igraph")
library("statnet")
library("ggnetwork")
library("rworldmap")

v1 = letters[1:5] # vertex names
A1 = rbind(c(0,1,0,0,1),
           c(1,0,0,1,0),
           c(0,0,0,0,1),
           c(0,0,1,0,1),
           c(1,0,0,0,0))
dimnames(A1) = list(v1, v1) # assign the names to A
A1

P1 = A1 / rowSums(A1)  
P1 %>% cbind(., Σ = rowSums(P1)) %>% rbind(., Σ = c(colSums(P1), sum(P1)))

A1net = network(A1, directed = TRUE) 
A1net

set.seed(0xdada)
plot(A1net, label = v1) 

A1df = ggnetwork(A1net)
head(A1df)

ggf = ggplot(A1df, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges() + 
  geom_nodes(aes(x = x, y = y), size = 6, color = "#8856a7") +
  geom_nodetext(aes(label = vertex.names), size = 4, color = "white") +
  theme_blank() + 
  theme(legend.position = "none")
ggf

#vertex visualisation
heatmap(A1, col = grey(c(0.9, 0.1)), symm = TRUE)



#STRING DB example
ccnb1 = read.table(
  url("http://web.stanford.edu/class/bios221/data/ccnb1datsmall.txt"), 
  header = TRUE, comment.char = "", check.names = FALSE)
head(ccnb1)

edge_list = as.matrix(ccnb1[, c("#node1", "node2")])
edge_list[1:3, ]

vertex_names = unique(as.vector(edge_list))  
n_edges = length(vertex_names)                  

M = matrix(0, 
           nrow = n_edges, ncol = n_edges,
           dimnames= list(vertex_names, vertex_names))
M[edge_list] = ccnb1$coexpression     # fill in the edge weights
M = M + t(M)                          # make symmetric
A = ifelse(M > 0, 1, 0)

ccnb1_net = network(A, directed = FALSE)
ccnb1_net

network(edge_list, directed = FALSE)

set.seed(0xbeef)
plot(ccnb1_net, label = vertex_names)    

set.seed(0xfish)
plot(ccnb1_net, label = vertex_names)    

# colors: white below 0.9, grey from 0.9 to 1 with darker closer to 1
breaks = c(0, seq(0.9, 1, length = 11))               # breaks in the color bins
cols = grey(1 - c(0, seq(0.5, 1, length = 10)))       # colors

# color the vertex for CCNB1 blue, its neighbors red, and the others white
ccnb1ind = which(vertex_names == "CCNB1")                  # index for ccnb1
vcols = rep("white", n_edges)
vcols[ccnb1ind] = "blue"
vcols[which(M[, ccnb1ind] > 0 | M[ccnb1ind, ] > 0 )] = "red"  

# now actually make the heat map
par(mar = rep(0,4))                                    # make plot margins 0
heatmap(M, symm = TRUE, ColSideColors = vcols, RowSideColors = vcols,
        col = cols, breaks = breaks,  frame = T)
legend("topleft", c("Neighbors(CCNB1)", "CCNB1"), fill = c("red","blue"),  
       bty = "n", inset = 0, xpd = T,  border = F)


#Minimum spanning trees 

load(url("http://web.stanford.edu/class/bios221/data/dist2009.RData"))
# clean up names a bit
country09 = sapply(attr(dist2009, "Labels"), function(x) {
  strsplit(x,"_") %>% unlist %>% `[`(length(.))
})  %>%
  unname %>%
  ifelse(. ==  "U.S.", "United States", .)
attr(dist2009, "Labels") = country09
str(dist2009)

# this takes a few minutes
mstree2009 = ape::mst(dist2009)
gr09 = graph.adjacency(mstree2009, mode = "undirected")
gg = ggnetwork(gr09, arrow.gap = 0, layout_with_fr(gr09))

ggplot(gg, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black",alpha = 0.5,curvature = 0.1) +
  geom_nodes(aes(color = name), size = 2) +
  geom_nodetext(aes(label = name), color = "black",size = 2) +
  theme_blank() +
  guides(color = guide_legend(keyheight = 0.1,keywidth = 0.1,
                              title = "Countries"))

?ggnetwork
#gedgelist <-network.edgelist(gg,network.initialize(4),ignore.eval=FALSE)

mat = match(country09, countriesLow$NAME)
lat2009 = countriesLow$LAT[mat]
lon2009 = countriesLow$LON[mat]
coords2009 = data.frame(lat = lat2009, lon = lon2009, country = country09)
x = jitter(coords2009$lon, amount = 15)
y = jitter(coords2009$lat, amount = 8)
layoutCoordinates = cbind(x, y)
labc = names(table(country09)[which(table(country09) > 1)])
idx = match(labc, countriesLow$NAME)
latc = countriesLow$LAT[idx]
lonc = countriesLow$LON[idx]
dfc = data.frame(latc, lonc, labc)
dfctrans = dfc
dfctrans[, 1] = (dfc[, 1]+31)/(93)
dfctrans[, 2] = (dfc[, 2]+105)/(238)

ggeo09 = ggnetwork(gr09, arrow.gap = 0, layout= layoutCoordinates)
ggplot(ggeo09,  aes(x = x,  y = y, xend = xend,  yend = yend)) +
  geom_edges(color  = "black", alpha = 0.5, curvature = 0.1) +
  geom_nodes(aes(color = name), size = 2) +
  theme_blank() +
  geom_label(
    data = dfctrans, aes(x = lonc, xend = lonc,
                         y = latc, yend = latc, label = labc,
                         fill = labc),colour = "white", alpha = 0.5, size = 3)+
  theme(legend.position = "none")



#MICROBIOME analysis
library("ggplot2")
library("PMA")
library("dplyr")
library("ade4")
library("genefilter")
library("ggrepel")
library("phyloseq")
library("vegan")

download.file("https://cdn.rawgit.com/spholmes/F1000_workflow/891463f6/data/ps.rds",
              "ps.rds", mode = "wb")
ps = readRDS("ps.rds")

class(ps)

?tax_table
tax1 = tax_table(ps)
summary(tax1)
length(!is.na(tax1))
sum(length(!is.na(tax1)))

?sample_data
head(sample_data(ps))

sample_data(ps) %>%
  filter(sex == "F") %>%
  pull("sex") %>%
  length()
# 173

sample_data(ps) %>%
  filter(sex == "F") %>%
  distinct(host_subject_id)
# 6

ps = subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

ggplot(sample_data(ps), aes(x = age)) + geom_histogram(bins = 40)

sample_data(ps)$age_binned = cut(
  sample_data(ps)$age, breaks = c(0, 100, 200, 400))

pstrans = transform_sample_counts(ps, asinh)

plot_tree(pstrans, "treeonly")

#PCoA
set.seed(10)     
out_wuf_asinh = ordinate(pstrans, method= "MDS", distance = "wunifrac")
evals = out_wuf_asinh$values$Eigenvalues

plot_ordination(pstrans, out_wuf_asinh, color = "age_binned", shape = "sex") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2]/evals[1]))

head(evals)
sume = sum(evals)
evals[1]/sume
evals[3]/sume
# 7.6%

outlier_idx = rownames(sample_data(pstrans)) %in% c("M3D149","M2D19","M1D9", "M1D149", "F5D165", "F6D165")
pstrans2 = prune_samples(!outlier_idx, pstrans)

#repeat for data without outliers
set.seed(10)     
out_wuf_asinh2 = ordinate(pstrans2, method= "MDS", distance = "wunifrac")
evals = out_wuf_asinh2$values$Eigenvalues
plot_ordination(pstrans2, out_wuf_asinh, color = "age_binned", shape = "sex") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2]/evals[1]))


out_bc_asinh = ordinate(pstrans2, method = "MDS", distance = "bray")

evals_bc = out_bc_asinh$values$Eigenvalues
plot_ordination(pstrans2, out_bc_asinh, color = "age_binned") +
  coord_fixed(sqrt(evals_bc[2] / evals_bc[1])) +
  labs(col = "Binned Age")

out_dpcoa_asinh = ordinate(pstrans2, method = "DPCoA")

evals_dpcoa = out_dpcoa_asinh$eig
plot_ordination(pstrans2, out_dpcoa_asinh, color = "age_binned",
                shape = "family_relationship") +
  coord_fixed(sqrt(evals_dpcoa[2] / evals_dpcoa[1])) +
  labs(col = "Binned Age", shape = "Litter")

plot_ordination(pstrans2, out_dpcoa_asinh, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals_dpcoa[2] / evals_dpcoa[1]))

