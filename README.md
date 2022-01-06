# SThe use of stereo-video monitoring and physiological trials to estimate metabolic demands of reef fishes in the wild

This repository contains the code and data to reproduce all tables and
figures presented in Schiettekatte, Conte et al. “The use of
stereo-video monitoring and physiological trials to estimate metabolic
demands of reef fishes in the wild”.

## Instructions

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

## Details

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

## Datasets

1.  data_fulton_2007 = data from Fulton et al. (2007) containing maximum
    swimming speeds of fishes.

2.  data_respirometry = Respirometry: these data are derived from
    intermittent respirometry, performed in the lab at CRIOBE. Each row
    is an individual fish and its associated metadata. SMR = standard
    metabolic rate in mg O2 per h. MaxMR = maximum metabolic rate in mg
    O2 per h. MeanTemp…C. = temperature in ºCelsius. Weight..kg. = mass
    in g.

3.  video_swimming_speeds = Measurements of size and swimming speed of
    our seven model species, obtained through video analysis using the
    software VidSync.

4.  data_uvc_moo_2016_fsp = Underwater visual census data, snapshot of a
    long-term monitoring program, kindly shared by CRIOBE, Moorea.

This paper was produced using the following software and associated
packages:

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value                       
    ##  version  R version 4.1.2 (2021-11-01)
    ##  os       Ubuntu 20.04.3 LTS          
    ##  system   x86_64, linux-gnu           
    ##  ui       X11                         
    ##  language en_US                       
    ##  collate  en_US.UTF-8                 
    ##  ctype    en_US.UTF-8                 
    ##  tz       Pacific/Honolulu            
    ##  date     2022-01-06                  
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package     * version date       lib source        
    ##  cachem        1.0.6   2021-08-19 [3] CRAN (R 4.1.1)
    ##  callr         3.7.0   2021-04-20 [1] CRAN (R 4.1.0)
    ##  cli           3.0.1   2021-07-17 [1] CRAN (R 4.1.0)
    ##  crayon        1.4.1   2021-02-08 [1] CRAN (R 4.1.0)
    ##  desc          1.4.0   2021-09-28 [3] CRAN (R 4.1.1)
    ##  devtools      2.4.3   2021-11-30 [3] CRAN (R 4.1.2)
    ##  digest        0.6.27  2020-10-24 [1] CRAN (R 4.1.0)
    ##  ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.1.0)
    ##  evaluate      0.14    2019-05-28 [1] CRAN (R 4.1.0)
    ##  fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.1.0)
    ##  fs            1.5.0   2020-07-31 [1] CRAN (R 4.1.0)
    ##  glue          1.4.2   2020-08-27 [1] CRAN (R 4.1.0)
    ##  htmltools     0.5.2   2021-08-25 [1] CRAN (R 4.1.1)
    ##  knitr         1.37    2021-12-16 [1] CRAN (R 4.1.2)
    ##  lifecycle     1.0.0   2021-02-15 [1] CRAN (R 4.1.0)
    ##  magrittr      2.0.1   2020-11-17 [1] CRAN (R 4.1.0)
    ##  memoise       2.0.1   2021-11-26 [3] CRAN (R 4.1.2)
    ##  pkgbuild      1.3.1   2021-12-20 [3] CRAN (R 4.1.2)
    ##  pkgload       1.2.4   2021-11-30 [3] CRAN (R 4.1.2)
    ##  prettyunits   1.1.1   2020-01-24 [1] CRAN (R 4.1.0)
    ##  processx      3.5.2   2021-04-30 [1] CRAN (R 4.1.0)
    ##  ps            1.6.0   2021-02-28 [1] CRAN (R 4.1.0)
    ##  purrr         0.3.4   2020-04-17 [1] CRAN (R 4.1.0)
    ##  R6            2.5.1   2021-08-19 [1] CRAN (R 4.1.0)
    ##  remotes       2.4.0   2021-06-02 [1] CRAN (R 4.1.0)
    ##  rlang         0.4.11  2021-04-30 [1] CRAN (R 4.1.0)
    ##  rmarkdown     2.11    2021-09-14 [1] CRAN (R 4.1.2)
    ##  rprojroot     2.0.2   2020-11-15 [1] CRAN (R 4.1.0)
    ##  rstudioapi    0.13    2020-11-12 [1] CRAN (R 4.1.0)
    ##  sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 4.1.0)
    ##  stringi       1.7.3   2021-07-16 [1] CRAN (R 4.1.0)
    ##  stringr       1.4.0   2019-02-10 [1] CRAN (R 4.1.0)
    ##  testthat      3.1.1   2021-12-03 [3] CRAN (R 4.1.2)
    ##  usethis       2.0.1   2021-02-10 [1] CRAN (R 4.1.0)
    ##  withr         2.4.2   2021-04-18 [1] CRAN (R 4.1.0)
    ##  xfun          0.29    2021-12-14 [1] CRAN (R 4.1.2)
    ##  yaml          2.2.1   2020-02-01 [1] CRAN (R 4.1.0)
    ## 
    ## [1] /home/nina/R/x86_64-pc-linux-gnu-library/4.1
    ## [2] /usr/local/lib/R/site-library
    ## [3] /usr/lib/R/site-library
    ## [4] /usr/lib/R/library

All code written by Francesca Conte and Nina M. D. Schiettekatte
(<nina.schiettekatte@gmail.com>). Please contact us for any issue or
question.
