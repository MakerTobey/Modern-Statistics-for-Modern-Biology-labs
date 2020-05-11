pkgs_needed = c("tidyverse", "Biostrings", 
                "parathyroidSE", "EnsDb.Hsapiens.v86")
letsinstall = setdiff(pkgs_needed, installed.packages()) 
if (length(letsinstall) > 0) {
  BiocManager::install(letsinstall)
}

library("dplyr")
library("ggplot2")
library("Biostrings")
library("readr")

births2006 = readr::read_csv(file = "https://www.huber.embl.de/msmb/course_spring_2020/data/births2006.csv.gz")

births = mutate(births2006,
                WEEKEND = ifelse(DOB_WK %in% c(1, 7), "Weekend", "Weekday"),
                HEALTH  = c("CritLow", "Low", "Normal")[ 1+findInterval(APGAR5, c(3.5, 6.5)) ],
                ESTGEST = replace(ESTGEST, ESTGEST==99, NA))

set.seed(1)
births_small = births[ sample(nrow(births), 40000), ]
head(births_small)

ppp = ggplot(births_small) 
ppp + geom_bar(aes(x = factor(DOB_WK)))
ppp + geom_bar(aes(x = factor(DOB_WK), fill = DMETH_REC))

ppp = ggplot(births_small) 
ppp + geom_histogram(aes(x = MAGER), binwidth = 1)
ppph = ggplot(births_small) + 
  geom_histogram(aes(x = MAGER, fill = DMETH_REC), binwidth = 1)
ppph + facet_wrap( ~ DOB_WK)

ppph + facet_grid(DOB_WK ~ SEX)

ggplot(dplyr::filter(births, !(DMETH_REC %in% "Unknown"))) +
  geom_histogram(aes(x = MAGER, fill = factor(TBO_REC)), binwidth = 1) +
  facet_grid(WEEKEND ~ DMETH_REC, scale="free_y", drop = TRUE) +
  geom_vline(xintercept = seq(15, 45, by=5), alpha=0.2, color="white") +
  labs(title = "Births in USA 2006", fill="Birth\nOrder")

ggplot(births_small, aes(x = MAGER)) + 
  geom_histogram(aes(y = ..density..), binwidth = 1, fill = "grey", col = "black") +
  stat_density(col = "red", fill = NA)

ggplot(births_small, aes(x = MAGER)) + 
  geom_histogram(aes(y = ..density..), binwidth = 1, fill="grey", col="black") +
  stat_density(col="red", fill=NA) +
  stat_function(fun = dlnorm, col = "blue", 
                args = list(mean = mean(log(births_small$MAGER)),
                            sd =    sd(log(births_small$MAGER))))

ppp2 = ggplot(dplyr::filter(births_small, !is.na(WTGAIN) & !is.na(DBWT)), 
              aes(x = WTGAIN, y = DBWT)) +
  labs(x = "Weight Gain by Mother", y = "Birth Weight in Grams")
ppp2 + geom_point()
ppp2 + stat_summary_2d(aes(z = ESTGEST), fun = median) + 
  labs(title = "median number of weeks of gestation")
ppp2 + geom_point() + stat_smooth(method = lm) 
ppp2 + geom_point() + stat_smooth(aes(col = SEX))
ppp2 + geom_hex(bins = 30) + stat_smooth(aes(col = SEX))

# we repeat again the code above so you don't gave to scroll to look at it.
ppp2 = ggplot(dplyr::filter(births_small, !is.na(WTGAIN) & !is.na(DBWT)), 
              aes(x = WTGAIN, y = DBWT)) +
  labs(x = "Weight Gain by Mother", y = "Birth Weight in Grams")

ppp2 + stat_bin2d() + scale_fill_gradient(trans = "sqrt")

ppp2 + stat_summary_2d(aes(z = ESTGEST), fun = mean) + 
  scale_y_log10(limits = 10^c(2, 4)) +
  scale_fill_gradient2(midpoint = 24) +
  labs(title="mean number of weeks of gestation")

ppp3 = ggplot(dplyr::filter(births, DPLURAL == "4 Quadruplet"), 
              aes(x = UPREVIS, y = MAGER)) + 
  geom_point(aes(size = ESTGEST, shape = DMETH_REC, col = DMEDUC)) +
  stat_smooth(aes(col = DMETH_REC), method = "lm")
ppp3

ppp3 + scale_size(range=c(3, 6)) + scale_color_brewer(palette = "Set1") +
  scale_shape(solid = FALSE)

####
oestat = function(o, e){
  sum( (e-o)^2/e )
}

set.seed(1)
B = 10000
# here we pick an arbitrary length / not the same as for Celegans
n = 2847
expected = rep(n/4 ,4)
oenull = replicate(
  B, oestat(e=expected, o=rmultinom(1,size = n, prob = rep(1/4,4))))

ggplot(data.frame(null_stats = oenull)) +
  geom_histogram(aes(x = null_stats), bins = 100, boundary=0)

ggplot(data.frame(stat = oenull), aes(sample = stat)) +
  stat_qq(distribution = stats::qchisq, dparams = list(df = 3)) +
  stat_qq_line(distribution = stats::qchisq, dparams = list(df = 3)) 

####
library("parathyroidSE")
library("EnsDb.Hsapiens.v86")

data("parathyroidGenesSE", package = "parathyroidSE")
metadata(parathyroidGenesSE)$MIAME 

abstract(metadata(parathyroidGenesSE)$MIAME)

genes = read.csv(textConnection(
  "name, group
   ESR1,  estrogen
   ESR2,  estrogen
   CASR,  parathyroid
   VDR,   parathyroid
   JUN,   parathyroid
   CALR,  parathyroid
   ORAI2, parathyroid"), 
  stringsAsFactors = FALSE, strip.white = TRUE)

ens = ensembldb::select(EnsDb.Hsapiens.v86,
                        keys = list(GenenameFilter(genes$name), 
                                    TxBiotypeFilter("protein_coding")),
                        columns = c("GENEID", "GENENAME"))

ens = 
  dplyr::filter(ens, GENEID %in% rownames(parathyroidGenesSE)) %>%
  mutate(group = genes$group[match(GENENAME, genes$name)])

ens

countData = assay( parathyroidGenesSE ) 
gene.counts = t(countData[ens$GENEID, ])
colnames(gene.counts) = ens$GENENAME
dat = cbind(data.frame(colData( parathyroidGenesSE)), data.frame(gene.counts))
head(dat)

ggplot(dat, aes(col = patient, x = treatment, y = ESR1)) +
  geom_point(size = 3) + 
  facet_grid( . ~ time)

#### Questions
# Try to answer the following questions to check your understanding of the topics covered in this lab.
# From the plot of the parathyroid data, answer the following.
#Quiz question 1 : How many patient sample are there?
head(countData)
test = data.frame(colData( parathyroidGenesSE))
head (test)
test %>% group_by(patient) %>% summarise(Count = n())
#answer is 4  

#Quiz question 2 : How many time points are there?
test %>% group_by(time) %>% summarise(Count = n())
#answer is 2 

#Quiz question 3 : There were 3 treatments: “Control”, “DPN”, and “OHT”. How many measurements were taken from patient sample 2 under the DPN treatment?
test %>% group_by(patient)
#answer is 
#Make your own plot of VDR versus CASR. (That is CASR, not CALR).
ggplot(test, aes(time = patient, x = treatment, y = CASR)) +
  geom_point(size = 3)

#Quiz question 4 : Which patient sample has the highest recorded level of CASR?
# 2

#Quiz question 5 : Which of the following pairs of patient samples seem to be well separated in this plot (i.e., for which two patient samples can you draw a line on the plot that perfectly separates them)?
# 2

#Quiz question 6 : NA
#--

#Quiz question 7 : Which patient sample looks different from the other three when you plot VDR versus ORAI2?
ggplot(dat, aes(col = patient, x = treatment, y = ORAI2)) +
  geom_point(size = 3) + 
  facet_grid( . ~ time)
# 3 noise

#Quiz question 8 : Plot ORAI2 versus JUN. Which can you separate best?
ggplot(dat, aes(col = patient, x = treatment, y = ORAI2)) +
  geom_point(size = 3) + 
  facet_grid( . ~ time)

#Quiz question 9 : Plot CASR versus (ESR1 + ESR2). Fit a separate linear model for each treatment (Control, DPN, OHT). Which linear models are increasing?
  
#Quiz question 10 : What is the maximum number of shapes that you are allowed to use in ggplot2 by default?
  #26
#Quiz question 11 : Write the name of the function that you could use to make more than the maximum number of default shapes allowed. Hint: this function has “values” as one of the arguments ____(…, values = (…)).
#scale_shape_manual
#Quiz question 12 : What do Themes do in ggplot2?
  #styles
  