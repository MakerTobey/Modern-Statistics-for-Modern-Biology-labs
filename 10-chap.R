## ----initialize, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
source("../chapter-setup.R"); chaptersetup("/Users/Susan/Courses/CUBook-html/CUBook/Chap11-NetworksTrees/Graphs.Rnw", "11")
knitr::opts_chunk$set(dev = 'png', dpi = 100, fig.margin = TRUE, fig.show = 'hold', fig.keep = 'none')

## ---- Darwin-Tree-1837, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/Darwin_Tree_1837.png'))

## ----chemo, fig.keep = 'high', fig.cap = "This graph represents genetic interactions that appeared to be modified in cancer patients where perturbations appeared in a chemokine signaling network.", echo = FALSE, fig.width = 4, fig.height = 4, warning = FALSE----
dats = read.table("../data/small_chemokine.txt", header = TRUE)
library("igraph")
library("tibble")
library("ggplot2")
library("reshape")
library("ggnetwork")
library("ggrepel")
library("dplyr")
gr = graph_from_data_frame(dats[,c("node1", "node2")], directed = FALSE)
E(gr)$weight = 1
V(gr)$size = centralization.degree(gr)$res
ggd=ggnetwork(gr,layout= "fruchtermanreingold",cell.jitter = 0)
gg=ggplot(ggd,aes(x=x, y=y, xend=xend, yend=yend))+
geom_edges(color="black", curvature=0.1, size=0.85, alpha=1/2)+
geom_nodes(aes(x=x, y=y),size=3,alpha=1/2, color="orange") +
geom_nodelabel_repel(aes(label = vertex.names), size=4, color="#8856a7") +
theme_blank() +theme(legend.position="none")
print(gg)

## ----igraphplot1, echo = FALSE, eval = 1:3-------------------------------
library("igraph")
edges1 = matrix(c(1,3,2,3,3,4,4,5,4,6),byrow = TRUE, ncol = 2)
g1 = graph_from_edgelist(edges1, directed = FALSE)
plot(g1, vertex.size = 25, edge.width = 5, vertex.color = "coral")

## ----adjacencyplot1, echo = FALSE----------------------------------------
ggplotadjacency = function(a){
  n = nrow(a)
  p = ncol(a)
  melted_a  =  melt(a)
  melted_a$value = as.factor(melted_a$value)
  cols = c("white", "darkblue")
  ggplot(data = melted_a, aes(x = X1, y = X2, fill=value)) +
    geom_tile(colour="black") +
    coord_fixed(ratio = 1,ylim=c(0.5,n+0.5),xlim=c(0.5,p+0.5))+
    scale_fill_manual(values=cols)+scale_x_discrete(name="",breaks=1:p)+
    scale_y_reverse(name="",breaks=n:1)+theme_bw() +
    theme(axis.text = element_text(size = 10),
      legend.key = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      panel.border=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line =element_line(color="white")
    )
}

## ----adjmatrix, fig.keep = 'high', fig.cap = "Here is the adjacency matrix of the small undirected graph represented below in Figure \\@ref(fig:igraphplot). We see that ${\\mathbf A}$ is symmetric ${n \\times n}$ matrix of $0$\'s and $1$\'s. ", fig.width=4, fig.height=4, echo = FALSE----
ggplotadjacency(as_adj(g1, sparse = FALSE))

## ----igraphplot, fig.keep = 'high', fig.cap = "A small undirected graph with numbered nodes.", echo = TRUE, eval = TRUE, ref.label = 'igraphplot1', fig.width=4, fig.height=4----

## ----graphex1wh----------------------------------------------------------
edges = "1,3\n2,3\n3,4\n4,6\n4,5"
sg  = graph_from_data_frame(read.csv(textConnection(edges),
                      header = FALSE), directed = FALSE)
sg

## ---- bipartite, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "This **bipartite** graph connects each taxon to the sites where it was observed."----
knitr::include_graphics(c('images/BipartiteTaxa.png'))

## ----chap11-r-finches-1, fig.keep = 'high', fig.cap = "The **bipartite** `finches` graph.", fig.width = 7, fig.height = 7----
library("networksis")
data("finch")
plot(finch, vertex.col =
       ifelse(nchar(network.vertex.names(finch)) == 1,
         "forestgreen", "gold3"),
     vertex.cex = 2.5, displaylabels = TRUE)

## ----chap11-r-ggnetworkong1-1, fig.keep = 'high', fig.cap = "A **[ggnetwork](https://cran.r-project.org/web/packages/ggnetwork/)** example.", message = FALSE, fig.width = 3, fig.height = 3----
library("ggnetwork")
g1df = ggnetwork(g1)
ggplot(g1df, aes(x = x, y = y, xend = xend, yend = yend)) +
 geom_edges() + geom_nodes(size = 6,color = "#8856a7") +
 geom_nodetext(aes(label = vertex.names), size = 4, color = "white") +
 theme_blank() + theme(legend.position = "none")

## ----T1net, cache = TRUE, fig.height = 9, fig.width = 9, eval = TRUE, results = FALSE, warnings = FALSE, messages = FALSE----
library("markovchain")
statesNames = c("A", "C", "G","T")
T1MC = new("markovchain", states = statesNames, transitionMatrix =
  matrix(c(0.2,0.1,0.4,0.3,0,1,0,0,0.1,0.2,0.2,0.5,0.1,0.1,0.8,0.0),
    nrow = 4,byrow = TRUE, dimnames=list(statesNames,statesNames)))
plot(T1MC, edge.arrow.size = 0.4, vertex.color="purple",
      edge.arrow.width = 2.2, edge.width = 5, edge.color = "blue",
      edge.curved = TRUE, edge.label.cex = 2.5, vertex.size= 32,
      vertex.label.cex = 3.5, edge.loop.angle = 3,
      vertex.label.family="sans", vertex.label.color = "white")

## ---- fourstateMC, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "A four state Markov chain with arrows representing possible transitions between states."----
knitr::include_graphics(c('images/MarkovACGT.png'))

## ----completechemokine, fig.keep = 'high', fig.cap = "Perturbed chemokine subnetwork uncovered in  @YuGXNA using differential gene expression patterns in sorted T-cells. Notice the clique-like structure of the genes CXCR3, CXCL13, CCL19, CSCR5 and CCR7 in the right hand corner.", fig.margin = FALSE, fig.width = 5, fig.height = 4----
library("ggnetwork")
oldpar=par(mar=c(5.1,4.7,4.1,2.6))
datf = read.table("../data/string_graph.txt", header = TRUE)
grs = graph_from_data_frame(datf[, c("node1", "node2")], directed = FALSE)
E(grs)$weight = 1
V(grs)$size = centralization.degree(grs)$res
ggdf = ggnetwork(grs, layout = "fruchtermanreingold", cell.jitter = 0)
ggplot(ggdf, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black", curvature = 0.1, size = 0.95, alpha = 0.8)+
  geom_nodes(aes(x = x, y = y), size = 3, alpha = 0.5, color = "orange") +
geom_nodelabel_repel(aes(label = vertex.names), size=4, color="#8856a7") +
#  geom_nodetext(aes(label = vertex.names), size = 4, color = "#8856a7") +
  theme_blank() + theme(legend.position = "none")
par(oldpar)


## ----testGSEABase, eval=FALSE--------------------------------------------
## library("GSEABase")
## ## Not run, this requires a login to the website.
## fl  =  "/path/to/msigdb_v5.1.xml"
## gss  =  getBroadSets(fl) # read entire msigdb
## organism(gss[[1]])
## table(sapply(gss, organism))


## ----simulationFishertest------------------------------------------------
universe = c(rep("Yellow", 500), rep("Blue", 100), rep("Red", 400))
countblue = replicate(20000, {
  pick75 = sample(universe, 75, replace = FALSE)
  sum(pick75 == "Blue")
})
summary(countblue)

## ----histblue, fig.keep = 'high', fig.cap = "We can see that even in (ref:histblue-1) simulations, no blue count comes close to being 25. We can reject such an event as having happened by chance and conclude that the blue are **enriched**.", echo = FALSE----
ggplot(data.frame(countblue), aes(x = countblue)) +
  geom_histogram(binwidth = 1, colour = "white", fill = "purple", center = 0.5, lwd = 1)

## ----checkcountblue, echo = FALSE----------------------------------------
stopifnot(all(countblue<22))

## ----GOplotEC, fig.keep = 'high', fig.cap = "This graph shows the correspondence between GO terms and significantly changed genes in a study on differential expression in endothelial cells from two steady state tissues (brain and heart, see  @Nolan:2013). After normalization a differential expression analysis was performed giving a list of genes. A gene-annotation enrichment analysis of the set of differentially expressed genes (adjusted p-value < 0.05) was then performed with the **[GOplot](https://cran.r-project.org/web/packages/GOplot/)** package.", fig.margin = FALSE, fig.width = 10,fig.height = 9, warning = FALSE, message = FALSE----
library("GOplot")
data("EC")
circ  =  circle_dat(EC$david, EC$genelist)
chord =  chord_dat(circ, EC$genes, EC$process)
GOChord(chord, limit = c(0, 5))

## ----BioNetsimple--------------------------------------------------------
library("BioNet")
library("DLBCL")
data("dataLym")
data("interactome")
interactome
pval = dataLym$t.pval
names(pval)  =  dataLym$label
subnet = subNetwork(dataLym$label, interactome)
subnet = rmSelfLoops(subnet)
subnet

## ----auxiliaryggplotfit, include = FALSE---------------------------------
## Function to qqplot the output from fitBumModel
ggplotqqbu = function(x) {
  n = length(x$pvalues)
  probs = (rank(sort(x$pvalues)) - 0.5)/n
  quantiles = unlist(sapply(probs, uniroot, f = BioNet:::.pbum.solve,
        interval = c(0, 1), lambda = x$lambda, a = x$a)[1, ])
  df = data.frame(fitted = quantiles, observed = sort(x$pvalues))
  ggplot(df, aes(x = fitted, y = observed)) +
      xlab("Theoretical quantiles") + ylab("Observed quantiles") +
      geom_point(size = 0.3, alpha = 0.3, color = "red")+
      geom_segment(aes(x = 0, y = 0, xend = 1, yend= 1 ), color = "blue") +
      coord_fixed(ratio = 1)
}

hist1.bum = function(x, breaks = 50){
  n = length(x$pvalues)
  bumdata = seq(from = 0.006, to = 1, length=n)
  ys = x$lambda + (1 - x$lambda) * x$a * bumdata^(x$a -1)
  dat = data.frame(pvalues = x$pvalues, xxs = bumdata, y = ys)
  y0 = BioNet:::piUpper(x)
  ggplot(dat, aes(x = pvalues)) +
    geom_histogram(aes( y= ..density..), binwidth = .02, fill = "orange", alpha = 0.75)+
    geom_hline(yintercept = y0, color = "blue3", alpha = 0.5, size = 1.2)+
    geom_line(aes(x = xxs, y = y), col = "red3", size = 1.3, alpha = 0.5)+
    xlab("p-values") +
    annotate("text", x = -0.03, y = y0 + 0.5, label = "pi[0]", parse = TRUE, size = 8)
}

## ----originalBioNet, eval = FALSE, include = FALSE-----------------------
## ## Original analysis as done by the authors using both p-values.
## library(BioNet)
## library(DLBCL)
## data("dataLym")
## data("interactome")
## interactome
## pvals  =  cbind(t = dataLym$t.pval, s = dataLym$s.pval)
## rownames(pvals)  =  dataLym$label
## pval  =  aggrPvals(pvals, order = 2, plot = FALSE)
## subnet  =  subNetwork(dataLym$label, interactome)
## subnet  =  rmSelfLoops(subnet)
## subnet

## ----fitBUM--------------------------------------------------------------
fb = fitBumModel(pval, plot = FALSE)
fb
scores=scoreNodes(subnet, fb, fdr = 0.001)

## ----plotFITBum, fig.keep = 'high', fig.cap = "The qqplot shows the quality of the fit of beta-uniform mixture model to the data. The red points have the theoretical quantiles from the beta distribution as the x coordinates the observed quantiles and the y coordinates. The blue line shows that this model fits nicely.", fig.width = 4, fig.height = 4, echo = FALSE----
ggp=ggplotqqbu(fb)
print(ggp)

## ----histFITBum, fig.keep = 'high', fig.cap = "A histogram of the mixture components for the p-values, the beta in red and the uniform in blue, $\\pi_0$ is the mixing proportion assigned to the null component whose distribution should be uniform. ",fig.width=4,fig.height=4,echo=FALSE----
ggh=hist1.bum(fb)
print(ggh)

## ----subgraphHeinz-------------------------------------------------------
hotSub  =  runFastHeinz(subnet, scores)
hotSub
logFC=dataLym$diff
names(logFC)=dataLym$label

## ----plotBioNet, fig.keep = 'high', fig.cap = "The subgraph found as maximally enriched for differential expression between ABC and GCB B-cell lymphoma. The nodes are colored in red and green: green shows an upregulation in ACB and red an upregulation in GBC. The shape of the nodes depicts the score: rectangles indicate a negative score, circles a positive score.", fig.margin = FALSE,fig.width=7.2,fig.height=7----
plotModule(hotSub, layout = layout.davidson.harel, scores = scores,
                  diff.expr = logFC)

## ---- IntroTree, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "As mathematical objects, the hierarchical clustering trees (studied in Chapter \\@ref(Chap:Clustering)) are the same as phylogenetic trees. They are **rooted binary** trees with labels at the tips."----
knitr::include_graphics(c('images/SameTreePhyl.png'))

## ---- HIVtree, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "This phylogenetic tree describes the history of different HIV/SIV strains in Africa  [@Wertheim:2009], [Figure from]."----
knitr::include_graphics(c('images/hivtreeplos.png'))



## ---- JCKimura, echo = FALSE---------------------------------------------
knitr::kable(data.frame("$Q=\\begin{array}{lcccc} & A &T &C & G\\\\A&-3\\alpha & \\alpha &\\alpha&\\alpha\\\\T&\\alpha & -3\\alpha&\\alpha &\\alpha\\\\C&\\alpha & \\alpha &-3\\alpha&\\alpha\\\\G&\\alpha &\\alpha &\\alpha &-3\\alpha\\\\\\end{array}$", "$Q=\\begin{array}{lcccc} & A &T &C & G\\\\A&-\\alpha-2\\beta & \\beta &\\beta&\\alpha\\\\T&\\beta & -\\alpha-2\\beta &\\alpha &\\beta\\\\C&\\beta & \\alpha &-\\alpha-2\\beta&\\beta\\\\G&\\alpha &\\beta &\\beta &-\\alpha-2\\beta\\\\\\end{array}$"), caption = "Two examples of rate matrices, on the left:the Jukes-Cantor (`JC69`) model, on the right is shown the Kimura (`K80`) two parameter model.", col.names = c("",""), align = "lr")


## ----Atree---------------------------------------------------------------
library("phangorn")
library("ggtree")
load("../data/tree1.RData")

## ----Atree1b, fig.keep = 'high', fig.cap = "This is the tree we use as our *true* parameter. We generate nucleotides one at a time from the root and `dropping\' them down the tree. With some probability proportional to the edge lengths, mutations occur down the branches.", fig.width = 5, fig.height = 5----
ggtree(tree1, lwd = 2, color = "darkgreen", alpha = 0.8, right = TRUE) +
  geom_tiplab(size = 7, angle = 90, offset = 0.05) +
  geom_point(aes(shape = isTip, color = isTip), size = 5, alpha = 0.6)

## ----ggtreeAlignment, fig.keep = 'high', fig.cap = "The tree on the left was used to generate the sequences on the right according to a Jukes Cantor model the nucleotide frequencies generated at the root were quite unequal, with \x60A\x60 and \x60C\x60 being generated more rarely. As the sequences percolate down the tree, mutations occur, they are more likely to occur on the longer branches. ", fig.margin = FALSE, fig.width = 7, fig.height = 4.5, warning = FALSE----
seqs6 = simSeq(tree1, l = 60, type = "DNA", bf = c(1, 1, 3, 3)/8, rate = 0.1)
seqs6
mat6df = data.frame(as.character(seqs6))
p = ggtree(tree1, lwd = 1.2) + geom_tiplab(aes(x = branch), size = 5, vjust = 2)
gheatmap(p, mat6df[, 1:60], offset = 0.01, colnames = FALSE)

## ----igraphsteiner, fig.keep = 'high', fig.cap = "A Steiner tree, the inner points are represented as squares. The method for creating the shortest tree thatpasses through all outer 1,2,5,6 is to create two inside (\"ancester\") points 3 and 4.",echo=FALSE----
edges1 = matrix(c(1,3,2,3,3,4,4,5,4,6),byrow = TRUE, ncol = 2)
g1 = graph_from_edgelist(edges1,directed = FALSE)
plot(g1, vertex.size=20,edge.width=5, vertex.color="coral",
    vertex.shape=c("circle","square")[c(1,1,2,2,1,1)])

## ----njtree1, fig.keep = 'high', fig.cap = "Trees built with a neighbor joining algorithm are very fast to compute and are often used as initial values for more expensive estimation procedures such as the maximum likelihood or parsimony."----
tree.nj = nj(dist.ml(seqs6, "JC69"))
ggtree(tree.nj) + geom_tiplab(size = 7) + ggplot2::xlim(0, 0.8)

## ----pmltree1------------------------------------------------------------
fit = pml(tree1, seqs6, k = 4)

## ----readseqs------------------------------------------------------------
library("dada2")
seqtab = readRDS("../data/seqtab.rds")
seqs = getSequences(seqtab)
names(seqs) = seqs

## ----tax, eval=FALSE-----------------------------------------------------
## fastaRef = "../tmp/rdp_train_set_16.fa.gz"
## taxtab = assignTaxonomy(seqtab, refFasta = fastaRef)

## ----taxtabload----------------------------------------------------------
taxtab = readRDS(file= "../data/taxtab16.rds")
dim(taxtab)

## ----taxtabanswer--------------------------------------------------------
taxtab %>% `rownames<-`(NULL) %>% head

## ----readDNAalign--------------------------------------------------------
readLines("../data/mal2.dna.txt") %>% head(12) %>% cat(sep="\n")

## ----alignwithdecipher, message=FALSE, warning=FALSE---------------------
library("DECIPHER")
alignment = AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = FALSE)

## ----treeFIT, eval = TRUE------------------------------------------------
phangAlign = phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
dm = phangorn::dist.ml(phangAlign)
treeNJ = phangorn::NJ(dm) # Note, tip order != sequence order
fit = phangorn::pml(treeNJ, data = phangAlign)
fitGTR = update(fit, k = 4, inv = 0.2)
fitGTR = phangorn::optim.pml(fitGTR, model = "GTR", optInv = TRUE,
         optGamma = TRUE,  rearrangement = "stochastic",
         control = phangorn::pml.control(trace = 0))

## ----samdatprepare, eval = FALSE, echo = FALSE, results=FALSE------------
## ###For the record
## ###mothur standard data, cleaned up here
## mimarksPath = "../data/MIMARKS_Data_combined.csv"
## samdf = read.csv(mimarksPath, header = TRUE)
## samdf$SampleID = paste0(gsub("00","",samdf$host_subject_id),
##                         "D", samdf$age-21)
## samdf = samdf[!duplicated(samdf$SampleID),]
## # Fixing an odd discrepancy
## rownames(seqtab) = gsub("124", "125", rownames(seqtab))
## all(rownames(seqtab) %in% samdf$SampleID)
## rownames(samdf) = samdf$SampleID
## keepCols = c("collection_date","biome","target_gene","target_subfragment",
## "host_common_name", "host_subject_id", "age", "sex", "body_product",
##  "tot_mass","diet", "family_relationship", "genotype", "SampleID")
## samdf = samdf[rownames(seqtab), keepCols]
## write.csv(samdf,"../data/MIMARKS_Data_clean.csv")

## ----samdat, eval = TRUE-------------------------------------------------
mimarksPathClean = "../data/MIMARKS_Data_clean.csv"
samdf = read.csv(mimarksPathClean, header = TRUE)

## ----phyloseqOBJ, eval = FALSE-------------------------------------------
## library("phyloseq")
## physeq = phyloseq(tax_table(taxtab), sample_data(samdf),
##         otu_table(seqtab, taxa_are_rows = FALSE), phy_tree(fitGTR$tree))

## ----onelineprune,eval=FALSE---------------------------------------------
## ps1 =  prune_samples(rowSums(otu_table(ps0)) > 5000, ps0)

## ----what-phyla,eval=FALSE-----------------------------------------------
## prev0 = apply(X = otu_table(ps0),
##               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
##               FUN = function(x) {sum(x > 0)})
## prevdf = data.frame(Prevalence = prev0,
##                      TotalAbundance = taxa_sums(ps0),
##                      tax_table(ps1))
## keepPhyla = table(prevdf$Phylum)[table(prevdf$Phylum) >  5]
## prevdf1   = subset(prevdf, Phylum %in% names(keepPhyla))
## ps2v      = subset_taxa(ps1v, Phylum %in% names(keepPhyla))

## ----deseq-transform, warning = FALSE------------------------------------
library("DESeq2")
library("phyloseq")
ps1 = readRDS("../data/ps1.rds")
ps_dds = phyloseq_to_deseq2(ps1, design = ~ ageBin + family_relationship)
geometricmean = function(x)
   if (all(x == 0)) { 0 } else { exp(mean(log(x[x != 0]))) }
geoMeans = apply(counts(ps_dds), 1, geometricmean)
ps_dds = estimateSizeFactors(ps_dds, geoMeans = geoMeans)
ps_dds = estimateDispersions(ps_dds)
abund = getVarianceStabilizedData(ps_dds)

## ----shorten-------------------------------------------------------------
rownames(abund) = substr(rownames(abund), 1, 5) %>% make.names(unique = TRUE)

## ----structssi-unadjp----------------------------------------------------
library("structSSI")
el = phy_tree(ps1)$edge
el0 = el
el0 = el0[rev(seq_len(nrow(el))), ]
el_names = c(rownames(abund), seq_len(phy_tree(ps1)$Nnode))
el[, 1] = el_names[el0[, 1]]
el[, 2] = el_names[el0[, 2]]
unadj_p = treePValues(el, abund, sample_data(ps1)$ageBin)

## ----structssi-test, results=FALSE---------------------------------------
hfdr_res = hFDR.adjust(unadj_p, el, 0.75)
summary(hfdr_res)
#plot(hfdr_res, height = 5000) # not run: opens in a browser

## ---- structssi-hfdr, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "A screenshot of a subtree with many differentially abundant microbes, as determined by the hierarchical testing procedure. Currently the user is hovering over the node associated with microbe GCGAG.33; this causes the adjusted p-value (0.029) to appear."----
knitr::include_graphics(c('images/structssi-screenshot.png'))

## ----structssi-tax-------------------------------------------------------
library("dplyr")
options(digits = 3)
tax = tax_table(ps1)[, c("Family", "Genus")] %>% data.frame
tax$seq = rownames(abund)
hfdr_res@p.vals$seq = rownames(hfdr_res@p.vals)
tax %>%  left_join(hfdr_res@p.vals[,-3]) %>%
  arrange(adjp) %>% head(9) %>% dplyr::select(1,2,4,5)

## ----ST-setup, echo=FALSE------------------------------------------------
library("igraph")
library("ggnetwork")
library("ggplot2")
pts = structure(c(0, 0, 1, 1, 1.5, 2, 0, 1, 1, 0, 0.5, 0.5),
                .Dim = c(6L,2L))
matxy = pts
distxy = stats::dist(matxy)
g = graph.adjacency(as.matrix(distxy), weighted = TRUE)
mst1 = igraph::mst(g)

## ----MST, fig.keep = 'high', fig.show = "hold", fig.cap = "A **spanning tree** is a tree that goes through all points at least once, the top graph with red edges is such a tree.", echo=FALSE,fig.width=4,fig.height=2.2----
gred=igraph::make_ring(6) - igraph::edges(6)
ggred=ggnetwork(gred, arrow.gap=0, layout = matxy)
ggplot(ggred, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(size=1.5, color = "red",alpha=0.8) +
  geom_nodes(size = 6) +
  theme_blank()
ggmst1=ggnetwork(mst1,arrow.gap=0,layout=matxy)
#gg=ggnetwork(gr1, arrow.gap=0, layout = matxy)
ggplot(ggmst1, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(size=1.5, color = "steelblue",alpha=0.8) +
  geom_nodes(size = 6) +
  theme_blank()

## ----HIVMSTi, fig.keep = 'high', fig.cap = "The minimum spanning tree computed from DNA distances between HIV sequences from samples taken in 2009 and whose country of origin was known, data as published in the \x60HIVdb\x60 database  [@HIVdb], .", fig.margin = FALSE,fig.width=7.5,fig.height=6----
load("../data/dist2009c.RData")
country09 = attr(dist2009c, "Label")
mstree2009 = ape::mst(dist2009c)
gr09 = graph.adjacency(mstree2009, mode = "undirected")
gg = ggnetwork(gr09, arrow.gap = 0, layout = "fruchtermanreingold")
ggplot(gg, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black", alpha = 0.5, curvature = 0.1) +
  geom_nodes(aes(color = vertex.names), size = 2) +  theme_blank() +
  geom_nodetext(aes(label = vertex.names), color = "black", size = 2.5) +
  theme(plot.margin = unit(c(0, 1, 1, 6), "cm"))+
  guides(color = guide_legend(keyheight = 0.09, keywidth = 0.09,
      title = "Countries")) + theme(legend.position = c(0, 0.14),
      legend.background = element_blank(),
      legend.text = element_text(size = 7))

## ----QuesAnsw, echo = TRUE-----------------------------------------------
ggplot(gg, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black", alpha = 0.5, curvature = 0.1) +
  geom_nodes(aes(color = vertex.names), size = 2) +
  geom_nodetext_repel(aes(label = vertex.names), color="black", size = 2) +
  theme_blank() +
  guides(color = guide_legend(keyheight = 0.3, keywidth = 0.3,
         override.aes = list(size = 6), title = "Countries"))

## ----HIVmap, fig.keep = 'high', fig.cap = "A minimum spanning tree between HIV cases using **a jitter** of the geographic location of each case was made to reduce overlapping of points; the DNA distances between patient strains were used as the input to an undirected minimum spanning tree algorithm.", fig.margin = FALSE, fig.width = 8.5, fig.height = 4,message = FALSE, warning =  FALSE----
library("rworldmap")
mat = match(country09, countriesLow$NAME)
coords2009 = data.frame(
  lat = countriesLow$LAT[mat],
  lon = countriesLow$LON[mat],
  country = country09)
layoutCoordinates = cbind(
  x = jitter(coords2009$lon, amount = 15),
  y = jitter(coords2009$lat, amount = 8))
labc = names(table(country09)[which(table(country09) > 1)])
matc = match(labc, countriesLow$NAME)
dfc = data.frame(
  latc = countriesLow$LAT[matc],
  lonc = countriesLow$LON[matc],
  labc)
dfctrans = dfc
dfctrans[, 1] = (dfc[,1] + 31) / 93
dfctrans[, 2] = (dfc[,2] + 105) / 238
ggeo09 = ggnetwork(gr09, arrow.gap = 0, layout = layoutCoordinates)
ggplot(ggeo09, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black", alpha = 0.5, curvature = 0.1) +
  geom_nodes(aes(color = vertex.names), size = 2) +
  theme_blank() +
  geom_label(data = dfctrans, aes(x = lonc, xend = lonc, y = latc, yend = latc,
       label = labc, fill = labc), colour = "white", alpha = 0.5, size = 3) +
   theme(legend.position = "none")

## ----WWtest, fig.keep = 'high', fig.cap = "Seeing the number of runs in a one-dimensional, two-sample, nonparametric Wald-Wolfowitz test can indicate whether the two groups have the same distributions.", fig.margin = FALSE, echo=FALSE, fig.width=10, fig.height=2----
dfbr=data.frame(measure=c(rnorm(15,0.9),rnorm(15,1.8)),
  group=as.factor(c(rep("men",15),rep("women",15))))
ggplot(dfbr,aes(x=measure,group=group,y=0)) + ylim(-0.25,+0.25) +
  geom_point(aes(col=group,x=measure,y=0,shape=group),size=5,alpha=0.6)+
  scale_color_manual(values=c("blue","red"))+
  theme_bw() + geom_hline(yintercept = 0) +
  theme(panel.border = element_blank(),
  axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.major.y = element_blank() ,
        panel.grid.minor.y = element_blank() )+ coord_fixed()


## ----preparegraph--------------------------------------------------------
ps1  = readRDS("../data/ps1.rds")
sampledata = data.frame( sample_data(ps1))
d1 = as.matrix(phyloseq::distance(ps1, method="jaccard"))
gr = graph.adjacency(d1,  mode = "undirected", weighted = TRUE)
net = igraph::mst(gr)
V(net)$id = sampledata[names(V(net)), "host_subject_id"]
V(net)$litter = sampledata[names(V(net)), "family_relationship"]

## ----mstplot, fig.keep = 'high', fig.cap = "The minimum spanning tree based on Jaccard dissimilarity and annotated with the mice ID and litter factors", fig.width=4, fig.height=4.5----
gnet=ggnetwork(net)
ggplot(gnet, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "darkgray") +
  geom_nodes(aes(color = id, shape = litter)) + theme_blank()+
  theme(legend.position="bottom")

## ----MSTJaccardplain-----------------------------------------------------
library("phyloseqGraphTest")
gt = graph_perm_test(ps1, "host_subject_id", distance="jaccard",
                     type="mst",  nperm=1000)
gt$pval

## ----mstJaccard, fig.keep = 'high', fig.cap = "The permutation histogram of the number of pure edges in the network obtained from the minimal spanning tree with Jaccard similarity.", fig.width=4, fig.height=2.1----
plot_permutations(gt)

## ----ggnetworkphyl-------------------------------------------------------
net = make_network(ps1, max.dist = 0.35)
sampledata = data.frame(sample_data(ps1))
V(net)$id = sampledata[names(V(net)), "host_subject_id"]
V(net)$litter = sampledata[names(V(net)), "family_relationship"]
netg = ggnetwork(net)

## ----ggnetworkplotJ, fig.keep = 'high', fig.cap = "A co-occurrence network created by using a threshold on the Jaccard dissimilarity matrix. The colors represent which mouse the sample came from; the shape represents which litter the mouse was in.", fig.margin = FALSE, fig.width=6, fig.height=5----
ggplot(netg, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "darkgray") +
  geom_nodes(aes(color = id, shape = litter)) + theme_blank()+
    theme(plot.margin = unit(c(0, 5, 2, 0), "cm"))+
    theme(legend.position = c(1.4, 0.3),legend.background = element_blank(),
          legend.margin=margin(0, 3, 0, 0, "cm"))+
         guides(color=guide_legend(ncol=2))

## ----mst-----------------------------------------------------------------
gt = graph_perm_test(ps1, "family_relationship",
        grouping = "host_subject_id",
        distance = "jaccard", type = "mst", nperm= 1000)
gt$pval

## ----mstpermplotNest, fig.keep = 'high', fig.cap = "The permutation histogram obtained from the minimal spanning tree with Jaccard similarity.", fig.width=4, fig.height=2.1----
plot_permutations(gt)

## ----knn1test------------------------------------------------------------
gtnn1 = graph_perm_test(ps1, "family_relationship",
                      grouping = "host_subject_id",
                      distance = "jaccard", type = "knn", knn = 1)
gtnn1$pval

## ----knn1plot, fig.keep = 'high', fig.cap = "The graph obtained from a nearest-neighbor graph with Jaccard similarity.", fig.width=4, fig.height = 2.7----
plot_test_network(gtnn1)

## ----adjacencyplot2, echo = TRUE, eval = FALSE, ref.label = 'adjacencyplot1'----
## ----adjmatrix2, echo = TRUE, eval = FALSE, ref.label = 'adjmatrix1'-----
## ---- PSB-MC-s, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "This figure was originally created for the study done in  @DigiulioCallahan:2015 where the nodes designate states of the vaginal microbiome and arrows represent transitions between states with different probabilities."----
knitr::include_graphics(c('images/PSB_MC_s.png'))

## ----plot-MC-plus,eval=FALSE---------------------------------------------
## library("markovchain")
## # Make Markov chain object
## mcPreg  =  new("markovchain", states=CSTs,
##               transitionMatrix = trans, name="PregCST")
## mcPreg
## # Set up igraph of the markov chain
## netMC  =  markovchain:::.getNet(mcPreg, round = TRUE)

## ----CSTMarkov, eval=FALSE-----------------------------------------------
## wts  =  E(netMC)$weight/100
## edgel  =  get.edgelist(netMC)
## elcat  =  paste(edgel[,1], edgel[,2])
## elrev  =  paste(edgel[,2], edgel[,1])
## edge.curved  =  sapply(elcat, function(x) x %in% elrev)
## samdf_def  =  data.frame(sample_data(ps))
## samdf_def  =  samdf_def[samdf$Preterm | samdf$Term,] # Only those definitely assigned, i.e. not marginal
## premat  =  table(samdf_def$CST, samdf_def$Preterm)
## rownames(premat)  =  markovchain::states(mcPreg)
## colnames(premat)  =  c("Term", "Preterm")
## premat
## premat  =  premat/rowSums(premat)
## vert.CSTclrs  =  CSTColors

## ----CSTMC, eval=FALSE---------------------------------------------------
## default.par  =  par(no.readonly = TRUE)
## # Define color scale
## # Plotting function for markov chain
## plotMC  =  function(object, ...) {
##     netMC  =  markovchain:::.getNet(object, round = TRUE)
##     plot.igraph(x = netMC, ...)
## }
## # Color bar for the markov chain visualization, gradient in strength of preterm association
## color.bar  =  function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title=NULL) {
##     scale = (length(lut)-1)/(max-min)
##     cur.par = par(no.readonly = TRUE)
##     par(mar = c(0, 4, 1, 4) + 0.1, oma = c(0, 0, 0, 0) + 0.1)
##     par(ps = 10, cex = 0.8)
##     par(tcl=-0.2, cex.axis=0.8, cex.lab = 0.8)
##     plot(c(min,max), c(0,10), type='n', bty='n', xaxt='n', xlab=", yaxt='n', ylab=", main=title)
##     axis(1, c(0, 0.5, 1))
##     for (i in 1:(length(lut)-1)) {
##       x = (i-1)/scale + min
##       rect(x,0,x+1/scale,10, col=lut[i], border=NA)
##     }
## }
## 
## 
## pal  =  colorRampPalette(c("grey50", "maroon", "magenta2"))(101)
## vert.clrs  =  sapply(states(mcPreg), function(x) pal[1+round(100*premat[x,"Preterm"])])
## vert.sz  =  4 + 2*sapply(states(mcPreg),
##               function(x) nrow(unique(sample_data(ps)[sample_data(ps)$CST==x,"SubjectID"])))
## vert.sz  =  vert.sz * 0.85
## vert.font.clrs  =  c("white", "white", "white", "white", "white")
## # E(netMC) to see edge list, have to define loop angles individually by the # in edge list, not vertex
## edge.loop.angle = c(0, 0, 0, 0, 3.14, 3.14, 0, 0, 0, 0, 3.14, 0, 0, 0, 0, 0)-0.45
## 
## layout  =  matrix(c(0.6,0.95, 0.43,1, 0.3,0.66, 0.55,0.3, 0.75,0.65), nrow = 5, ncol = 2, byrow = TRUE)
## 
## # Colored by association with preterm birth
## layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(1,10))
## color.bar(pal, min=0, max=1, nticks=6, title="Fraction preterm")
## par(mar=c(0,1,1,1)+0.1)
## edge.arrow.size=0.8
## edge.arrow.width=1.4
## edge.width = (15*wts + 0.1)*0.6
## edge.labels  =  as.character(E(netMC)$weight/100)
## edge.labels[edge.labels<0.4]  =  NA  # labels only for self-loops
## plotMC(mcPreg, edge.arrow.size=edge.arrow.size, edge.arrow.width = edge.arrow.width,
##        edge.width=edge.width, edge.curved=edge.curved,
##        vertex.color=vert.clrs, vertex.size=(vert.sz),
##        vertex.label.font = 2, vertex.label.cex = 1,
##        vertex.label.color = vert.font.clrs, vertex.frame.color = NA,
##        layout=layout, edge.loop.angle = edge.loop.angle)
## par(default.par)

## ---- ccnb1img, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "This network was created with the STRING website by setting a 2 step neighborhood around the CCNB1 gene for co-expression levels $\\geq$ 0.900. "----
knitr::include_graphics(c('images/ccnb1img.png'))

## ----ccnb1data-----------------------------------------------------------
dat = read.table("../data/ccnb1datsmall.txt", header = TRUE, comment.char = "")
v = levels(unlist(dat[,1:2]))        # vertex names
n = length(v)                        # number of vertices
e = matrix(match(as.character(unlist(dat[,1:2])), v),ncol=2) # edge list
w = dat$coexpression                 # edge weights

## ----ccnbet--------------------------------------------------------------
M = matrix(0, n, n)
M[e] = w
M = M + t(M)
dimnames(M) = list(v, v)
A = 1*(M > 0)


## ----plotgraph-----------------------------------------------------------
library(igraph)
net = network(e, directed=FALSE)
par(mar=rep(0,4))
plot(net, label=v)

## ----heatmapCCNB1, fig.keep = 'high', fig.cap = "This represents the adjacency of the CCNB1 network -- 2 step neighborhood with co-expression levels $\\geq$ 0.900, generated from R (darker is closer to 1, we ignore values < 0.9)."----
breaks  =  c(0, seq(0.9, 1, length=11))
cols  =  grey(1-c(0,seq(0.5,1,length=10)))
ccnb1ind  =  which(v == "CCNB1")
vcols  =  rep("white",n)
vcols[ccnb1ind]  =  "blue"
vcols[which(M[,ccnb1ind]>0 | M[ccnb1ind,])]  =  "red"
par(mar = rep(0, 4))
heatmap(M, symm = TRUE, ColSideColors = vcols, RowSideColors = vcols,
        col = cols, breaks = breaks,  frame = TRUE)
legend("topleft", c("Neighbors(CCNB1)", "CCNB1"),
       fill = c("red","blue"),
       bty = "n", inset = 0, xpd = TRUE,  border = FALSE)

## library("ape")
## library("phangorn")
## GAG=read.dna("../data/DNA_GAG_20.txt")

## ----knn-2-test, fig.width=4.5, fig.height=4.5---------------------------
gt = graph_perm_test(ps1,"family_relationship",
       distance = "bray", grouping = "host_subject_id",
       type = "knn", knn = 2)
gt$pval

## ----knn-2-plot, fig.keep = 'high', fig.show = "hold", fig.cap = "The graph and permutation histogram obtained from a two nearest-neighbor graph with Jaccard similarity.", fig.width=4.5, fig.height = 2.7----
plot_test_network(gt)
permdf = data.frame(perm=gt$perm)
obs = gt$observed
ymax = max(gt$perm)
ggplot(permdf, aes(x = perm)) + geom_histogram(bins = 20) +
  geom_segment(aes(x = obs, y = 0, xend = obs, yend = ymax/10), color = "red") +
  geom_point(aes(x = obs, y = ymax/10), color = "red") + xlab("Number of pure edges")
#plot_permutations(gt)

