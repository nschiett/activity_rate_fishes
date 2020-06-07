source("R/packages.R")
source("R/functions.R")
source("R/plan.R")

if(!is.na(
  tryCatch((readd(main_text_pdf)), error = function(cond){return(NA)})
)){
  clean(main_text_pdf)
}
if(!is.na(
  tryCatch((readd(main_text_doc)), error = function(cond){return(NA)})
)){
  clean(main_text_doc)
}

config <- drake_config(plan = plan, lock_envir = FALSE)

