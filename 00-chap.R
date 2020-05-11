## ----initialize, echo = FALSE, message = FALSE, error = FALSE, warning = FALSE----
knitr::opts_chunk$set(dev = "png", dpi = 100, fig.margin = TRUE, fig.show = "hold", fig.keep = "none")

## ---- WringingFlood, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high'----
knitr::include_graphics(c('images/WringingFlood.png'))

## ---- intro-Fisher, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "The hypothesis testing paradigm recommended by R.A.\\ Fisher starts with the formulation of a null hypothesis and the design of an experiment before the collection of any data. We could think in a similarly schematic way about model fitting -- just replace *Hypothesis H0* by *Parametric Model* and \\textit{Compute p-value} by *Fit Parameters*."----
knitr::include_graphics(c('images/FisherParadigm.png'))

## ---- intro-iterative, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "JW Tukey recommended starting any analysis with the data and wrote: \"No catalogue of techniques can convey a willingness to look for what can be seen, whether or not anticipated.\"  [@Holmes1985], [{\\tiny Image source}]."----
knitr::include_graphics(c('images/iterativeparadigm.png'))

## ---- roulette-chunk-1, eval = TRUE, echo = FALSE, fig.keep = 'high'-----
knitr::include_graphics('images/roulette.png', dpi = 600)

## ---- structuralescal, eval = TRUE, echo = FALSE, fig.show = 'hold', fig.keep = 'high', fig.cap = "Analyzing data is not a one step process. Each step involves visualizing and decomposing some of the complexity in the data. Tukey\'s iterative data structuration can be conceptualized as $Total=V_1+V_2+V_3$"----
knitr::include_graphics(c('images/structuralescal.png'))


