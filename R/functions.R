
#### 1) metabolism model ####

run_mod_mr <- function(respiro){
  # tall dataframe
  respiro <- filter(respiro, !Species == "Chlorurus spirulus")
  respiro_tall <- tidyr::gather(respiro, "MR_type", "MR_gd", SMR_gd, MMR_gd)
  # fit model
  priors <- c(set_prior("normal(0.75,0.1)", class = "b", coef = "log10Weight_g", resp = "log10SMRgd"),
              set_prior("normal(0.75,0.1)", class = "b", coef = "log10Weight_g", resp = "log10MMRgd"))
  fit <- brm(
    mvbind(log10(SMR_gd), log10(MMR_gd)) ~ 1 + log10(Weight_g) + (1|Species) + (0 + log10(Weight_g)|Species),
             data = respiro, prior = priors, control = list(adapt_delta = 0.9))

    return(fit)
}

get_table_mr <- function(mod_mr){
  mod_mr %>%
    spread_draws(b_log10SMRgd_Intercept,
                 b_log10SMRgd_log10Weight_g,
                 r_Species__log10SMRgd[species, smr_type],
                 b_log10MMRgd_Intercept,
                 b_log10MMRgd_log10Weight_g,
                 r_Species__log10MMRgd[species, mmr_type]) %>% 
    pivot_wider(names_from = smr_type, values_from = r_Species__log10SMRgd) %>%
    dplyr::rename(int_smr_sp = Intercept, b_smr_sp = log10Weight_g) %>%
    pivot_wider(names_from = mmr_type, values_from = r_Species__log10MMRgd) %>%
    dplyr::rename(int_mmr_sp = Intercept, b_mmr_sp = log10Weight_g) %>%
    dplyr::mutate(slope_smr_sp = b_log10SMRgd_log10Weight_g + b_smr_sp,
                  slope_mmr_sp = b_log10MMRgd_log10Weight_g + b_mmr_sp,
                  intercept_smr_sp = exp10(b_log10SMRgd_Intercept + int_smr_sp),
                  intercept_mmr_sp = exp10(b_log10MMRgd_Intercept + int_mmr_sp)) %>%
    dplyr::select(species, slope_smr_sp, slope_mmr_sp, intercept_smr_sp, intercept_mmr_sp) %>%
    group_by(species) %>%
    median_hdci() %>%
    mutate_if(is.numeric, round, 4) %>%
    mutate_at(.vars = vars(starts_with("slope")), round, 2) %>%
    dplyr::mutate(
       species = gsub("\\.", " ", species), 
       `SMR slope` = paste0(slope_smr_sp, " (", slope_smr_sp.lower, ";", slope_smr_sp.upper, ")"),
       `SMR (weight = 1g)` = paste0(intercept_smr_sp, " (", intercept_smr_sp.lower, ";", intercept_smr_sp.upper, ")"),
       `MMR slope` = paste0(slope_mmr_sp, " (", slope_mmr_sp.lower, ";", slope_mmr_sp.upper, ")"),
       `MMR (weight = 1g)` = paste0(intercept_mmr_sp, " (", intercept_mmr_sp.lower, ";", intercept_mmr_sp.upper, ")")
  ) %>%
    select(species, `SMR slope`,`SMR (weight = 1g)`, `MMR slope`, `MMR (weight = 1g)` )
}


#### 2) speed model ####

run_mod_speed <- function(speed){
  fit <- brm(log10(Speed_cms) ~ 1 + log10(Length_cm) + (log10(Length_cm) | Species),
      data = speed, control = list(adapt_delta = 0.999, max_treedepth = 20), family = "student")
  return(fit)
}

get_table_speed <- function(mod_speed){
  mod_speed %>%
    spread_draws(b_Intercept,
                 b_log10Length_cm,
                 r_Species[species, type]) %>% 
    pivot_wider(names_from = type, values_from = r_Species) %>%
    dplyr::rename(int_sp = Intercept, b_sp = log10Length_cm) %>%
    dplyr::mutate(slope_sp = b_log10Length_cm + b_sp,
                  intercept_sp = b_Intercept + int_sp) %>%
    dplyr::select(species, slope_sp, intercept_sp) %>%
    group_by(species) %>%
    median_hdci() %>%
    mutate_if(is.numeric, round, 2) %>%
    dplyr::mutate(
      species = gsub("\\.", " ", species), 
      slope = paste0(slope_sp, " (", slope_sp.lower, ";", slope_sp.upper, ")"),
      intercept = paste0(intercept_sp, " (", intercept_sp.lower, ";", intercept_sp.upper, ")")
    ) %>%
    dplyr::select(species, slope, intercept)
}


#### 3) speedmax  ####

### wrangle speedmax
wrangle_speedmax <- function(speedmax){
  
  # check for name errors
  errors <- fishflux::name_errors(unique(speedmax$Species))
  
  # corrected names using taxize
  corr <- errors %>%
    gnr_resolve(data_source_ids = c(3, 4), 
                with_canonical_ranks = T, 
                highestscore = T) %>%
    select(submitted_name, matched_name2) %>% 
    unique()
  
  speedmax2 <- speedmax %>% 
    # join operation:
    left_join(corr, by = c("Species" = "submitted_name")) %>%
    mutate(Species = case_when(!is.na(matched_name2) ~ matched_name2,
                               TRUE ~ Species)) %>%
    select(-matched_name2)
  
  # second check errors
  fishflux::name_errors(unique(speedmax2$Species))
  
  # manually replace remaining errors
  speedmax2 <- speedmax2 %>%
    mutate(Species = recode(Species,
                            "Apogon nigrofasciatus"  = "Ostorhinchus nigrofasciatus",
                            "Cheilinus chlorurus"    = "Cheilinus chlorourus" ,
                            "Coris schroederii"      = "Coris batuensis",
                            "Oxycheilinus digrammus" = "Oxycheilinus digramma",
                            "Scolopsis bilineatus"   = "Scolopsis bilineata"
    ))
  
  body <- species(unique(speedmax2$Species), fields=c("Species","BodyShapeI"))
  
  # get aspect ratio's 
  ar <- lapply(unique(speedmax2$Species),
               function(x){
                 print(x)
                 fishflux::aspect_ratio(x)}) %>%
    bind_rows() %>%
    select(Species = species, aspect_ratio)
  
  # Combine plus filter out the five families
  speedmax2 <- speedmax2 %>%
    left_join(ar) %>%
    left_join(body) %>%
    filter(Family %in% c(
      "Chaetodontidae", "Acanthuridae", "Pomacentridae", "Serranidae",
      "Balistidae"
    ))
  return(speedmax2)
}

# run model
run_mod_speedmax <- function(speedmax){
  fit <- brm(log10(Ucrit_cms) ~ 1 + log10(Length_cm) + aspect_ratio + 
               (log10(Length_cm) | Family:BodyShapeI), 
             data = speedmax, control = list(adapt_delta = 0.95), 
             family = "student", iter = 4000)
  return(fit)
}

get_table_speedmax <- function(mod_speedmax){
  mod_speedmax %>%
    spread_draws(b_Intercept,
                 b_log10Length_cm,
                 r_Family[family, type]) %>% 
    pivot_wider(names_from = type, values_from = r_Family) %>%
    dplyr::rename(int_fam = Intercept, b_fam = log10Length_cm) %>%
    dplyr::mutate(slope_fam = b_log10Length_cm + b_fam,
                  intercept_fam = b_Intercept + int_fam) %>%
    dplyr::select(family, slope_fam, intercept_fam) %>%
    group_by(family) %>%
    median_hdci() %>%
    mutate_if(is.numeric, round, 2) %>%
    dplyr::mutate(
      slope = paste0(slope_fam, " (", slope_fam.lower, ";", slope_fam.upper, ")"),
      intercept = paste0(intercept_fam, " (", intercept_fam.lower, ";", intercept_fam.upper, ")")
    ) %>%
    dplyr::select(family, slope, intercept)
}

##### FMR ####

create_reference <- function(moorea){
  
  # min and max size found in Moorea
  minmax <- moorea %>% filter(!species %in% c("Chaetodon pelewensis", "Chlorurus spilurus")) %>% 
    group_by(species) %>% 
    dplyr::summarize(min = min(Size), max = max(Size))
  
  # get aspect ratios 
  ar <- lapply(c("Zebrasoma scopas", 
                 "Chromis iomelas",
                 "Naso lituratus",
                 "Odonus niger",
                 "Chaetodon ornatissimus",
                 "Ctenochaetus striatus",
                 "Cephalopholis argus"),
               function(x){
                 print(x)
                 fishflux::aspect_ratio(x)}) %>%
    bind_rows() %>%
    dplyr::select(species, aspect_ratio)
  
  # body shape
  body <- species(unique(ar$species), fields=c("Species","BodyShapeI")) %>%
    select(species = Species, BodyShapeI)
  
  
  ref <-
    lapply(1:6, function(i){
      res <- data.frame(
        species = minmax$species[i],
        Size = minmax$min[i]:minmax$max[i]
      )
    }) %>% bind_rows() %>%   
    # add Chromis iomelas that is not in moorea dataframe
    rbind(data.frame(   
      species = "Chromis iomelas",
      Size = 2:8
    )) %>%
    # add aspect ratio
    left_join(ar) %>%
    # add body
    left_join(body)
  
  fams <- data.frame(species = unique(ref$species), 
                     Family = c("Serranidae", "Chaetodontidae", "Acanthuridae", "Acanthuridae", "Balistidae", "Acanthuridae", "Pomacentridae"))
  ref <- left_join(ref, fams)
  
  # length weight relationship 
  lw <- lapply(unique(ref$species), fishflux::find_lw, mirror = "se") %>% bind_rows() %>%
    dplyr::select(species, lwa = lwa_m, lwb = lwb_m)
  
  ref <- ref %>% 
    left_join(lw) %>%
    mutate(weight = lwa * (Size^lwb))
  
  return(ref)
}

exp10 <- function(x){
  10^x
}


make_plot_mr <- function(mod_mr, respiro){
  
  pred_smr <- fitted(mod_mr, summary = FALSE, nsamples = 1000)[,,1] %>% as.vector()
  pred_mmr <- fitted(mod_mr, summary = FALSE, nsamples = 1000)[,,2] %>% as.vector()
  
  fit_mr <- mod_mr$data %>% 
    slice(rep(1:n(), each = 1000)) %>%
    mutate(smr = exp10(pred_smr), mmr = exp10(pred_mmr))

  plot <-
    ggplot(fit_mr) +
    geom_point(aes(x = Weight_g, y = SMR_gd, color = Species), alpha = 0.7, size = 2,
               data = mod_mr$data) +
    geom_point(aes(x = Weight_g, y = MMR_gd, color = Species), alpha = 0.7, size = 2,
               shape = 17, data = mod_mr$data) +
    stat_lineribbon(aes(y = smr, x = Weight_g, color = Species, fill = Species), 
                    .width = c(.8), alpha = 0.2)  +
    stat_lineribbon(aes(y = smr, x = Weight_g, color = Species, fill = Species), 
                    .width = c(0), alpha = 0.8)  +
    stat_lineribbon(aes(y = mmr, x = Weight_g, color = Species, fill = Species), 
                    .width = c(.8), alpha = 0.2, linetype = 2)  +
    stat_lineribbon(aes(y = mmr, x = Weight_g, color = Species, fill = Species), 
                    .width = c(0), alpha = 0.8, linetype = 2)  +
    theme_bw() +
    labs(x = "Weight (g)", y = expression(paste("Metabolic rate (g ", O[2], ".", d^{-1} ,")")),
         fill = "Species", color = "Species", linetype = "", shape = "") +
    scale_color_fish_d(option = "Chaetodon_ephippium") +
    scale_fill_fish_d(option = "Chaetodon_ephippium") +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    scale_linetype_discrete(labels = c("MMR", "SMR")) +
    scale_shape_discrete(labels = c("MMR", "SMR")) +
    theme(legend.text = element_text(face = "italic"), 
          text = element_text( size = 14))
  #plot
  ggsave(plot, filename = "output/plots/plot_mr.png" ,width = 10, height = 6)
  
  return(plot)  
  
}


make_plot_speed <- function(reference, mod_speed, mod_speedmax, speed){

newdata <- reference %>%
  select(Species = species, Length_cm = Size, aspect_ratio, Family, BodyShapeI)

pred_speed <- fitted(mod_speed, newdata, summary = FALSE, nsamples = 1000) %>%
  as.vector()

fit_speed <- reference %>% slice(rep(1:n(), each = 1000)) %>%
  mutate(v = exp10(pred_speed))

pred_speedmax <- fitted(mod_speedmax, newdata, summary = FALSE, nsamples = 1000) %>%
  as.vector()

fit_speedmax <- reference %>% slice(rep(1:n(), each = 1000)) %>%
  mutate(vmax = exp10(pred_speedmax)) 

plot <-
ggplot(fit_speed) +
  geom_point(aes(x = Length_cm, y = Speed_cms, color = Species), alpha = 0.5, size = 1,
             data = filter(speed, ! Species %in%c("Chlorurus spilurus", "Chaetodon pelewensis"), Mode == "Swimming")) +
  stat_lineribbon(aes(y = (v), x = Size, color = species, fill = species), 
                 .width = c(.8), alpha = 0.2)  +
  stat_lineribbon(aes(y = (v), x = Size, color = species, fill = species), 
                  .width = c(0), alpha = 0.8)  +
  stat_lineribbon(aes(y = (vmax), x = Size, color = species, fill = species), 
                  .width = c(0), alpha = 0.8, data = fit_speedmax, linetype = 2, size = 1)  +
  theme_bw() +
  labs(x = "Length (cm)", y = expression(paste("Swimming speed (cm ", s^{-1} ,")")), fill = "Species", color = "Species") +
  scale_color_fish_d(option = "Chaetodon_ephippium") +
  scale_fill_fish_d(option = "Chaetodon_ephippium") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  theme(legend.text = element_text(face = "italic"), 
        text = element_text(size = 14))

ggsave(plot, filename = "output/plots/plot_speed.png", width = 10, height = 6)

return(plot)  

}


get_fmr <- function(reference, mod_speed, mod_speedmax, mod_mr){
  
  newdata <- reference %>%
    select(Species = species, Length_cm = Size, aspect_ratio, Family, BodyShapeI)
  
  pred_speed <- fitted(mod_speed, newdata, summary = FALSE, nsamples = 4000) %>%
    as.vector()
  
  fit_speed <- reference %>% slice(rep(1:n(), each = 4000)) %>%
    mutate(v = exp10(pred_speed))
  
  pred_speedmax <- fitted(mod_speedmax, newdata, summary = TRUE)
  
  fit_speedmax <- reference %>% 
    mutate(vmax = exp10(pred_speedmax[,1])) 
  
  newdata <- reference %>%
    select(Species = species, Weight_g = weight)
  
  pred_mr <- fitted(mod_mr, newdata, robust = TRUE) 
  
  fit_mr <- reference %>% 
    mutate(smr = exp10(pred_mr[,1,1]), mmr = exp10(pred_mr[,1,2]))
  
  fit_speedmax$id <- 1:nrow(fit_speedmax)
  
  fit <- left_join(fit_speed, fit_speedmax) %>% left_join(fit_mr) %>%
    mutate(
      b = (log10(mmr) - log10(smr))/vmax
    ) %>%
    mutate(
      FMR = exp10(log10(smr) + (b * v)) 
    ) %>%
    mutate(
      FSA = (FMR + smr)/(2*smr)
    ) %>%
    group_by(Family, species, Size, weight) %>%
    dplyr::summarise(FMR_m = median(FMR),
                     logFMR_m = mean(log10(FMR)),
                     logFMR_sd = sd(log10(FMR)),
                     FMR_lb = quantile(FMR, 0.025),
                     FMR_ub = quantile(FMR, 0.975),
                     FSA_m = median(FSA),
                     FSA_lb = quantile(FSA, 0.025),
                     FSA_ub = quantile(FSA, 0.975)
    )
    
    
   
    ## get FAS ##
    newdata <- reference %>%
      select(Species = species, Weight_g = weight)
    
    pred_smr <- fitted(mod_mr, newdata, summary = FALSE, nsamples = 4000) %>%
      as.vector()
    
    pred_smr <- fitted(mod_mr, newdata, summary = FALSE, nsamples = 4000)[,,1] %>% as.vector()
    pred_mmr <- fitted(mod_mr, newdata, summary = FALSE, nsamples = 4000)[,,2] %>% as.vector()
    
    fit_mr <- reference %>% 
      slice(rep(1:n(), each = 4000)) %>%
      mutate(smr = exp10(pred_smr), mmr = exp10(pred_mmr))
    
    fas <- fit_mr %>%  
      dplyr::mutate(FAS = mmr/smr) %>%
      group_by(species, Size) %>%
      dplyr::summarize(FAS_m = median(FAS),
                      FAS_lb = quantile(FAS, 0.025),
                      FAS_ub = quantile(FAS, 0.975),
                      smr_m = median(smr),
                      smr_lb = quantile(smr, 0.025),
                      smr_ub = quantile(smr, 0.975),
                      mmr_m = median(mmr),
                      mmr_lb = quantile(mmr, 0.025),
                      mmr_ub = quantile(mmr, 0.975))
    
    result <- fit %>%
      left_join(fas)
    

    return(result)
}

get_slopes <- function(field_summary, mod_mr){
  
  fit_fmr <- brm(logFMR_m|se(logFMR_sd) ~ 1 + log10(weight) 
                 + (1|species) + ( 0 + log10(weight)|species),
                 data = field_summary)
  
  summary(fit_fmr)
  
  slope_fmr <- fit_fmr %>% 
    spread_draws(b_log10weight, r_species[species, log10weight]) %>%
    mutate(slope = b_log10weight + r_species) %>%
    filter(log10weight == "log10weight") %>%
    median_qi() %>%
    select(species, fmr_alpha_m = slope, fmr_alpha_lb = slope.lower, fmr_alpha_ub = slope.upper) %>%
    mutate(species = gsub("\\.", " ", species))
  
  ## slope SMR 
  slope_mr <- mod_mr %>% 
    spread_draws(b_log10SMRgd_log10Weight_g, b_log10MMRgd_log10Weight_g, 
                 r_Species__log10SMRgd[species, log10Weight_g], r_Species__log10MMRgd[species, log10Weight_g]) %>%
    mutate(slope_smr = b_log10SMRgd_log10Weight_g + r_Species__log10SMRgd,
           slope_mmr = b_log10MMRgd_log10Weight_g + r_Species__log10MMRgd) %>%
    filter(log10Weight_g == "log10Weight_g") %>%
    median_qi()  %>%
    select(species, smr_alpha_m = slope_smr, smr_alpha_lb = slope_smr.lower, smr_alpha_ub = slope_smr.upper,
           mmr_alpha_m = slope_mmr, mmr_alpha_lb = slope_mmr.lower, mmr_alpha_ub = slope_mmr.upper) %>%
    mutate(species = gsub("\\.", " ", species)) 
  
  slopes <- slope_mr %>% left_join(slope_fmr) 
  
  return(slopes)
}

make_plot_fsa <- function(field_summary){

result <- field_summary %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(size_stand = (Size - min(Size)) /(max(Size)-min(Size)))

p_fsa <-
  ggplot(result) +
  geom_line(aes( x = Size, y = FSA_m, color = species), size = 1)+
  geom_point(aes( x = Size, y = FSA_m, color = species), size = 2, alpha = 1)+
  scale_color_fish_d(option = "Chaetodon_ephippium", guide = "none") +
  theme_bw() +
  labs(x = "Length (cm)", y = "FSA",
       fill = "Species", color = "Species", size = "Standardized length") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"), 
        text = element_text( size = 14)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.text = element_text(face = "italic"),
        legend.spacing = unit(3, "cm")) 

ggsave( "output/plots/FSA_FAS.png", p_fsa, width = 10, height = 6)

return(p_fsa)
}


make_plot_slopes <- function(slopes){
  
  tall <- slopes %>%
    pivot_longer(2:10) %>%
    separate(name, c("mr", "alpha", "q")) %>%
    pivot_wider(names_from = q, values_from = value)
  
  p_scale <-
  ggplot(tall) +
    geom_point(aes(x = m, color = species, y = mr), position = position_dodgev(height =  0.7), size = 2) +
    geom_errorbarh(aes(xmax = ub, xmin =  lb, color = species, y = mr), 
                   height = 0, position = position_dodgev(height =  0.7), size = 1, linetype = 1) +
    scale_color_fish_d(option = "Chaetodon_ephippium") +
    geom_vline(aes(xintercept = m, color = species), data = filter(tall, mr == "smr"), linetype = 2) +
    scale_y_discrete(labels = c("FMR", "MMR", "SMR")) +
    labs(x = "Scaling coefficient", y = "", color = "") +
    theme_bw() +
    theme(legend.text = element_text(face = "italic"),
          text = element_text(size = 14),
          legend.position = "bottom")
  
  ggsave("output/plots/scaling_sp.png", p_scale, width = 10, height = 6)
  
  return(p_scale)
}

make_plot_combined <- function(plot_fsa, plot_slopes){
  
  p <-
   plot_slopes + labs(tag = 'a)') + plot_fsa + labs(tag = 'b)')  + 
    plot_layout(ncol = 1,  widths = c(2, 2), heights = c(3,3))
  
  ggsave("output/plots/combined.png", p, height = 10, width = 8)
  
  return(p)
}

make_plot_community_mr <- function(field_summary, moorea){

com <- moorea %>% left_join(field_summary) %>%
  drop_na() %>%
  ungroup() %>%
  filter(!Site_name %in% c("Moorea Entre 2 Baies", "Moorea Haapiti")) %>%
  dplyr::mutate(FMR = (FMR_m + smr_m)/(2)) %>%
  mutate(fmr = FMR * Abundance, smr = smr_m * Abundance) %>%
  mutate(fmr = fmr/Method_area,
         smr = smr/Method_area) %>%
  group_by(Site_name, Lat, Long) %>%
  dplyr::summarise(smr = sum(smr),
            fmr = sum(fmr)) %>%
  group_by(Site_name) %>%
  dplyr::summarise(SMR = mean(smr),
                   FMR = mean(fmr)) %>%
  mutate(ratio = FMR/SMR, l = SMR * 2) %>%
  group_by(Site_name) %>%
  pivot_longer(c(SMR, FMR)) 

p_com <-
ggplot(com) +
  geom_bar(aes(x = name, y = value, fill = name), stat = "identity", alpha = 0.9) +
  geom_hline(aes(yintercept = 1.7 * value), data = filter(com, name == "SMR"), linetype = 2) +
  facet_wrap(~Site_name, scales = "free") +
  scale_fill_fish_d(option = "Callanthias_australis", 
                     labels = c("FMR", "SMR"), begin = 0.4, end = 1) +  
  theme_bw()+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        text = element_text( size = 12))  +
  labs(fill = "")

ggsave("output/plots/community_plots.png", p_com, width = 10, height = 8)

return(p_com)

}

make_plot_community_ab <- function(moorea){
  moorea_s <- moorea %>%
    filter(!Site_name %in% c("Moorea Entre 2 Baies", "Moorea Haapiti")) %>%
    group_by(Site_name, Method_area, species) %>%
    dplyr::summarise(ab = sum(Abundance)) %>%
    mutate(ab_m2 = ab/Method_area) %>%
    filter(!species %in% c("Chmorurus spilurus", "Chaetodon pelewensis"))
  
  plot <-
  ggplot(moorea_s)+ 
    geom_bar(aes(x = as.factor(Site_name), y = ab_m2, fill = species, color = species), stat = "identity")+
    xlab("")+
    labs(y = expression(paste("Fish abundance ", (m^{-2}))))+
    scale_fill_fish_d(option = "Chaetodon_ephippium") +
    scale_color_fish_d(option = "Chaetodon_ephippium") +
    theme_bw()+
    coord_flip() +
    theme(legend.text = element_text(face="italic"), text = element_text(family = "Calibri", size = 14))
  ggsave("output/plots/community_plot_abundance.png", plot, width = 10, height = 8)
  return(plot)
}


load_fig1 <- function(){
  fig <- rasterGrob(as.raster(readPNG("text/figures/figure1.png")), interpolate = FALSE)
  return(fig)
}




