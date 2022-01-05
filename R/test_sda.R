
gC_to_kJ  <- 39

loadd(field_summary)
par <- read_csv("data/params_sst_glob.csv") %>%
  filter(Species %in% unique(field_summary$species), v_m == 27)

par$size = 10


cnp_out <- function(x){
  
  d <- data[x,]
  size <- purrr::simplify(d$size)
  p <- d %>%
    select(k_m, Qc_m ,Qn_m, Qp_m, Dc_m, Dn_m, Dp_m, alpha_m, f0_m,          
           theta_m, lwa_m,  lwb_m, r_m, h_m,  v_m, linf_m,        
           F0nz_m, F0pz_m, ac_m, an_m, ap_m) %>% purrr::simplify() %>% as.list()
  
  fit <- fishflux::cnp_model_mcmc(TL = size, param = p)
  
  out <- fishflux::extract(fit, c("Ic", "Sc", "w1", "F0c"))
  
  out
}

results <- lapply(1:nrow(par), cnp_out) %>%
  plyr::ldply() 
results <- cbind(data, results)

test <- results %>%
  # meal mass
  mutate(mm = Ic_mean * 100/Dc_m) %>%
  # SDA 
  mutate(sda = exp((0.33 * log(w1_mean)) + (0.67 * log(mm) ) - 0.45)/gC_to_kJ) %>%
  mutate(sda_ratio = (sda )/F0c_mean)


0.33 log bm + 0.67 log mm − 0.45
