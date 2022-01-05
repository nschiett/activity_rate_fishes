-   [Stereo-video monitoring and physiological trials reveal metabolic
    demands of reef fishes in the
    wild](#stereo-video-monitoring-and-physiological-trials-reveal-metabolic-demands-of-reef-fishes-in-the-wild)
    -   [Instructions](#instructions)
    -   [Details](#details)
    -   [Datasets](#datasets)

Stereo-video monitoring and physiological trials reveal metabolic demands of reef fishes in the wild
====================================================================================================

This repository contains the code and data to reproduce all tables and
figures presented in Schiettekatte, Conte et al. “Stereo-video
monitoring and physiological trials reveal metabolic demands of reef
fishes in the wild.”

Instructions
------------

All analyses for the current project were done in R, using the drake
pipeline [drake](https://github.com/ropensci/drake). You can use drake
to compile all models, figures, and tables. To do so, first install
`drake` from CRAN:

``` r
install.packages("drake")
```

Next you need to open an R session with working directory set to the
root of the project.

We use a number of packages, listed in `r/packages.R`. If needed,
missing packages should first be installed.

Then, to generate all figures, analyses, and tables, simply run:

``` r
drake::r_make()
```

All output will be automatically rendered inside the folder called
output.

Details
-------

The elements of this project include:

1.  data: this folder includes all the raw data necessary to run the
    script.
2.  output: this folder includes all outputs from the R-script,
    including figures and the manuscript in docx.
3.  text: this folder includes the raw manuscript, and cover letter.
4.  R: this folder contains three .R files, `packages.R`, `functions.R`,
    and `plan.R`.  
    `packages.R` contains all packages needed, `functions.R` contains
    all functions to reproduce the analyses, and `plan.R` provides the
    code that binds each step of the workflow.

Datasets
--------

1.  data\_fulton\_2007 = data from Fulton et al. (2007) containing
    maximum swimming speeds of fishes.

2.  data\_respirometry = Respirometry: these data are derived from
    intermittent respirometry, performed in the lab at CRIOBE. Each row
    is an individual fish and its associated metadata. SMR = standard
    metabolic rate in mg O2 per h. MaxMR = maximum metabolic rate in mg
    O2 per h. MeanTemp…C. = temperature in ºCelsius. Weight..kg. = mass
    in g.

3.  video\_swimming\_speeds = Measurements of size and swimming speed of
    our seven model species, obtained through video analysis using the
    software VidSync.

4.  data\_uvc\_moo\_2016\_fsp = Underwater visual census data, snapshot
    of a long-term monitoring program, kindly shared by CRIOBE, Moorea.

This paper was produced using the following software and associated
packages:

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value                       
    ##  version  R version 3.6.3 (2020-02-29)
    ##  os       Ubuntu 16.04.6 LTS          
    ##  system   x86_64, linux-gnu           
    ##  ui       X11                         
    ##  language en_GB                       
    ##  collate  en_GB.UTF-8                 
    ##  ctype    en_GB.UTF-8                 
    ##  tz       Europe/Paris                
    ##  date     2020-10-21                  
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package     * version date       lib source        
    ##  assertthat    0.2.1   2019-03-21 [1] CRAN (R 3.6.0)
    ##  backports     1.1.8   2020-06-17 [1] CRAN (R 3.6.3)
    ##  callr         3.4.3   2020-03-28 [1] CRAN (R 3.6.2)
    ##  cli           2.0.2   2020-02-28 [1] CRAN (R 3.6.2)
    ##  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.6.0)
    ##  desc          1.2.0   2018-05-01 [1] CRAN (R 3.6.0)
    ##  devtools      2.3.0   2020-04-10 [1] CRAN (R 3.6.2)
    ##  digest        0.6.25  2020-02-23 [1] CRAN (R 3.6.1)
    ##  ellipsis      0.3.1   2020-05-15 [1] CRAN (R 3.6.3)
    ##  evaluate      0.14    2019-05-28 [1] CRAN (R 3.6.0)
    ##  fansi         0.4.1   2020-01-08 [1] CRAN (R 3.6.1)
    ##  fs            1.4.1   2020-04-04 [1] CRAN (R 3.6.2)
    ##  glue          1.4.1   2020-05-13 [1] CRAN (R 3.6.3)
    ##  htmltools     0.4.0   2019-10-04 [1] CRAN (R 3.6.1)
    ##  knitr         1.28    2020-02-06 [1] CRAN (R 3.6.1)
    ##  magrittr      1.5     2014-11-22 [1] CRAN (R 3.6.0)
    ##  memoise       1.1.0   2017-04-21 [1] CRAN (R 3.6.0)
    ##  pkgbuild      1.0.8   2020-05-07 [1] CRAN (R 3.6.3)
    ##  pkgload       1.1.0   2020-05-29 [1] CRAN (R 3.6.3)
    ##  prettyunits   1.1.1   2020-01-24 [1] CRAN (R 3.6.1)
    ##  processx      3.4.2   2020-02-09 [1] CRAN (R 3.6.1)
    ##  ps            1.3.3   2020-05-08 [1] CRAN (R 3.6.3)
    ##  R6            2.4.1   2019-11-12 [1] CRAN (R 3.6.1)
    ##  Rcpp          1.0.4.6 2020-04-09 [1] CRAN (R 3.6.2)
    ##  remotes       2.1.1   2020-02-15 [1] CRAN (R 3.6.1)
    ##  rlang         0.4.6   2020-05-02 [1] CRAN (R 3.6.2)
    ##  rmarkdown     2.1     2020-01-20 [1] CRAN (R 3.6.1)
    ##  rprojroot     1.3-2   2018-01-03 [1] CRAN (R 3.6.0)
    ##  sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.6.0)
    ##  stringi       1.4.6   2020-02-17 [1] CRAN (R 3.6.1)
    ##  stringr       1.4.0   2019-02-10 [1] CRAN (R 3.6.0)
    ##  testthat      2.3.2   2020-03-02 [1] CRAN (R 3.6.2)
    ##  usethis       1.6.1   2020-04-29 [1] CRAN (R 3.6.3)
    ##  withr         2.2.0   2020-04-20 [1] CRAN (R 3.6.2)
    ##  xfun          0.14    2020-05-20 [1] CRAN (R 3.6.3)
    ##  yaml          2.2.1   2020-02-01 [1] CRAN (R 3.6.1)
    ## 
    ## [1] /home/nschiettekatte/R/x86_64-pc-linux-gnu-library/3.6
    ## [2] /usr/local/lib/R/site-library
    ## [3] /usr/lib/R/site-library
    ## [4] /usr/lib/R/library

All code written by Nina M. D. Schiettekatte
(<nina.schiettekatte@gmail.com>) and Francesca Conte. Please contact us
for any issue or question.
