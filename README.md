Julia code implementing an lmer extension to the rERP method (Smith and Kutas, 2015a, 2015b. See also Brouwer, Delogu, Crocker, 2021). 

This code accompanies the publication by Aurnhammer, Delogu, Schulz, Brouwer, Crocker (subm.).

The original ERP data as well as the fitted datasets described in the article can be downloaded from https://osf.io/nrb4t/ . Datasets are to be placed in the ./data/ directory.

PACKAGE and CODE information

The Julia code implements data pre-processing and fitting of lmerEPR models.
Developed on Julia version 1.5.1.
Using Packages: 

[336ed68f] CSV v0.7.7

[a93c6f00] DataFrames v0.21.7

[1313f7d8] DataFramesMeta v0.5.1

[ff71e718] MixedModels v3.0.0-DEV `https://github.com/JuliaStats/MixedModels.jl.git#Master`

[2dfb63ee] PooledArrays v0.5.3

R code used for plotting estimated waveforms, residuals, coefficients, z-values. 
Developed on R version 3.6.1. Using Packages:

data.table 1.12.8

ggplot2 3.2.1

grid 3.6.1

gridExtra 2.3


ABBREVIATIONS

cl : Cloze

srp : log(Cloze) (Suprisal)

rcnoun : association rating of the adverbial clause noun

rcverb : association rating of the adverbial clause verb

spr : self-paced reading

ERP : event-related potential

est : estimated data

res : residual error

##############
# References # 
##############

Brouwer, H., Delogu, F., & Crocker, M. W. (2021). Splitting Event‐Related Potentials:
  Modeling Latent Components using Regression‐based Waveform Estimation. 
  European Journal of Neuroscience.

Julia: A Fast Dynamic Language for Technical Computing. 
  Jeff Bezanson, Stefan Karpinski, Viral B. Shah, Alan Edelman.
  (2012). arXiv: 1209.5145.

R Core Team (2019). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna, Austria.
  URL https://www.R-project.org/.

Smith, N. J., & Kutas, M. (2015a). Regression‐based estimation of ERP waveforms:
  I. The rERP framework. Psychophysiology, 52(2), 157-168.

Smith, N. J., & Kutas, M. (2015b). Regression‐based estimation of ERP waveforms: 
  II. Nonlinear effects, overlap correction, and practical considerations. 
  Psychophysiology, 52(2), 169-181.

