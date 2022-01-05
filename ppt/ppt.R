
theme_ppt <- function(){
  theme_classic() +
    theme(rect = element_rect(fill = cols()[1], color = cols()[1]), 
          panel.background = element_rect(fill = cols()[1], color = cols()[1]), 
          plot.background = element_rect(fill = cols()[1], color = cols()[1]), 
          
          text = element_text(color = cols()[3], size = 10),
          axis.text = element_text(color = cols()[3], size = 9),
          axis.title = element_text(color = cols()[3], size = 10),
          axis.ticks = element_line(color = cols()[3]), 
          axis.line = element_line(colour = cols()[3], size = 1) 
    ) 
}

cols <- function(){
  c("#00183aff","#c2c2c2ff","#e6e6e6ff","#6ac2e5ff",
    "#85ffbcff","#fffc5cff","#f6ab13ff","#df3416ff")
}

save_plot<- function(plot, name, width = 27, height = 12){
  ggsave(filename = paste0("ppt/plots/", name, ".png"), plot = plot, 
         width = width, height = height,
         units = "cm")
}

source("R/packages.R")
drake::loadd(field_summary)

field_summary <-
field_summary %>%
  filter(!species == "Chromis iomelas")

a <- 
ggplot(field_summary) +
  geom_line(aes(x = weight, y = FMR_m), color = cols()[7], size = 1) +
  geom_line(aes(x = weight, y = smr_m), color = cols()[6], size = 1) +
  geom_line(aes(x = weight, y = mmr_m, group = species), color = cols()[8], size = 1) +
  facet_wrap(~species, scales = "free_x", nrow = 1) +
  theme_ppt() + theme(strip.text = element_blank()) +
  labs(x = "", y = expression(paste("Metabolic rate (g ", O[2], d^{-1} ,")")))
a
b <-
ggplot(field_summary) +
  geom_line(aes(x = weight, y = FMR_m/weight), color = cols()[7], size = 1) +
  geom_line(aes(x = weight, y = smr_m/weight), color = cols()[6], size = 1) +
  geom_line(aes(x = weight, y = mmr_m/weight, group = species), color = cols()[8], size = 1) +
  facet_wrap(~species, scales = "free_x", nrow = 1) +
  theme_ppt() + theme(strip.text = element_blank()) +
  labs(x = "", y = expression(paste("Mass-specific MR (g ", O[2], d^{-1} , g^{-1}, ")")))

c <-
  ggplot(field_summary) +
  geom_line(aes(x = weight, y = FSA_m), color = cols()[7], size = 1) +
  facet_wrap(~species, scales = "free_x", nrow = 1) +
  theme_ppt() + theme(strip.text = element_blank()) +
  labs(y = "Activity scope", x = "Weight (g)")

loadd(slopes)

tall <- slopes %>%
  pivot_longer(2:10) %>%
  separate(name, c("mr", "alpha", "q")) %>%
  pivot_wider(names_from = q, values_from = value) %>%
  filter(!species == "Chromis iomelas")

d <-
  ggplot(tall) +
  geom_vline(xintercept = 0.75, linetype = 2, color = "white") +
  geom_point(aes(x = m, color = mr, y = "mr"), size = 2,
             position = position_dodgev(height =  0.7)) +
  geom_errorbarh(aes(xmax = ub, xmin =  lb, color = mr, y = "mr"), 
                 height = 0, position = position_dodgev(height =  0.7), size = 1, linetype = 1) +
  facet_wrap(~species, nrow = 1) +
  theme_ppt() + theme(strip.text = element_blank()) +
  labs(y = "Activity scope", x = "Metabolic scaling coefficient") +
  scale_color_manual(values = c(cols()[c(7, 8,6)])) +
  theme(legend.position = "none", axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank())
d


plot <- 
d +b  + c + plot_layout(nrow = 3) &
  theme(panel.background = element_rect(fill = cols()[1], color = cols()[1]), 
        plot.background = element_rect(fill = cols()[1], color = cols()[1]))

save_plot(plot, "c7_p1")

####### swimming ############
colnames(field_summary)

a <- 
  ggplot(field_summary) +
  geom_line(aes(x = weight, y = FMR_m), color = cols()[7], size = 1) +
  geom_line(aes(x = weight, y = smr_m), color = cols()[6], size = 1) +
  geom_line(aes(x = weight, y = mmr_m, group = species), color = cols()[8], size = 1) +
  facet_wrap(~species, scales = "free_x", nrow = 1) +
  theme_ppt() + theme(strip.text = element_blank()) +
  labs(x = "Weight (g)", y = expression(paste("Metabolic rate (g ", O[2], d^{-1} ,")")))
a

loadd(reference, mod_speed, mod_speedmax, speed)
  
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
    ggplot(fit_speed[!fit_speed$species == "Chromis iomelas",]) +
    stat_lineribbon(aes(y = (v), x = Size), color = cols()[7],
                    .width = c(0), alpha = 0.8)  +
    stat_lineribbon(aes(y = (vmax), x = Size),  color = cols()[8],
                    .width = c(0), alpha = 0.8, data = fit_speedmax[!fit_speedmax$species == "Chromis iomelas",], 
                    linetype = 1, size = 1)  +
    theme_bw() +
    labs(x = "Length (cm)", y = expression(paste("Swimming speed (cm ", s^{-1} ,")")), 
         fill = "Species", color = "Species") +
    # scale_x_continuous(trans = "log10") +
    # scale_y_continuous(trans = "log10") +
    theme_ppt() +
    theme(legend.position = "none", strip.text = element_blank()) +
    facet_wrap(~species, nrow = 1, scales = "free_x") 
  plot
  
  p <- 
  plot + a + plot_layout(nrow = 2) &
    theme(panel.background = element_rect(fill = cols()[1], color = cols()[1]), 
          plot.background = element_rect(fill = cols()[1], color = cols()[1]))
  save_plot(p, "c7_speed", height = 12)

####### map #############"
loadd(moorea)
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
  mutate(ratio = FMR/SMR, l = SMR * 2) 

moorea_s <- moorea %>%
  filter(!Site_name %in% c("Moorea Entre 2 Baies", "Moorea Haapiti")) %>%
  group_by(Site_name, Method_area, species) %>%
  dplyr::summarise(ab = sum(Abundance)) %>%
  mutate(ab_m2 = ab/Method_area) %>%
  filter(!species %in% c("Chlorurus spilurus", "Chaetodon pelewensis")) %>%
  mutate(species = as.factor(species))


siteplot <- function(site, com, moorea_s){
  
  col <- cols()[c(8,6,3,7,5,4)]
  
  plot1 <-
    ggplot(com[com$Site_name == site, ]) +
    geom_bar(aes(x = " ", y = ratio), fill = "grey" ,
             stat = "identity", alpha = 0.9) +
    theme_bw() +
    theme(axis.title.x=element_blank(), 
          text = element_text(size = 12),
          axis.text = element_text(color = "black"), 
          panel.grid.minor = element_blank())  +
    scale_y_continuous(limits = c(0,2.5)) +
    labs(y ="Ratio FMR:SMR", x = "",
         title = site) +
    theme_ppt()
  
  plot2 <-
    ggplot(moorea_s[moorea_s$Site_name == site,])+ 
    geom_bar(aes(x = "1", y = ab_m2, fill = species, color = species), stat = "identity")+
    xlab("") +
    labs(y = expression(paste("Fish abundance ", (m^{-2}))))+
    scale_color_manual(values = col, drop = FALSE) +
    scale_fill_manual(values = col, drop = FALSE) +
    theme_bw() + labs(x = " ", title = " ") +
    theme(legend.position = "none",  
          axis.title.x = element_blank(), axis.title.y = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.border = element_blank(), plot.margin = unit(c(2,0,4,0), "mm")) +
    theme_ppt() +
    theme(legend.position = "none",  
          axis.title.x = element_blank(), axis.title.y = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.border = element_blank(), plot.margin = unit(c(2,0,4,0), "mm"), 
          axis.line = element_blank()) 
    
  
  plot_grid(plot1, plot2, nrow = 1,ncol = 2, rel_widths = c(3,1))
  
}

a <- siteplot("Gendron", com, moorea_s)
b <- siteplot("Tiahura", com, moorea_s)
c <- siteplot("Entre 2 baies", com, moorea_s)
d <- siteplot("Pihaena", com, moorea_s)
e <- siteplot("Aroa", com, moorea_s)
f <- siteplot("Nuarei", com, moorea_s)
g <- siteplot("Temae", com, moorea_s)
h <- siteplot("Motu Ahi", com, moorea_s)
i <- siteplot("Afareaitu", com, moorea_s)
j <- siteplot("Maatea", com, moorea_s)
k <- siteplot("Haapiti", com, moorea_s)
m <- siteplot("Taotaha", com, moorea_s)
n <- siteplot("Tetaiuo", com, moorea_s)


coast <- sf::st_read("data/coastline.shp") %>%
  st_transform(st_crs("+proj=longlat +datum=WGS84"))

sites <-  moorea %>%
  filter(!Site_name %in% c("Moorea Entre 2 Baies", "Moorea Haapiti")) %>%
  select(Lat, Long, Site_name) %>%
  group_by(Site_name) %>%
  summarize_all(mean)

map <-
  ggplot() +
  geom_sf(data = coast, fill = "grey80", color = "grey80") +
  geom_point(aes(y = Lat, x = Long), data = sites, color = "white") +
  geom_label_repel(aes(y = Lat, x = Long, label = Site_name), 
                   data = sites, color = cols()[1]) +
  coord_sf() +
  theme_ppt() +
  theme(text = element_text(color = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.line = element_blank())

layout <- "
  ABCDE
  NOOOF
  MOOOG
  LKJIH
  "

l <-
  ggplot() +
  geom_text(aes(x = 2, y = 1:6, label = c("Z. scopas", "O. niger", 
                                          "N. lituratus", "C. striatus", 
                                          "C. ornatissimus", "C. argus")),
            size = 3.5, fontface = 3, hjust = 0, color = cols()[3]) +
  add_fishape(family = "Serranidae", option = "Cephalopholis_argus",
              fill = col[1],
              xmin = 1.5, xmax = 2, ymin = 5.5, ymax = 6.5) +
  add_fishape(family = "Chaetodontidae", option = "Chaetodon_ornatissimus",
              fill = col[2],
              xmin = 1.5, xmax = 2, ymin = 4.5, ymax = 5.5) +
  add_fishape(family = "Acanthuridae", option = "Ctenochaetus_striatus",
              fill = col[3],
              xmin = 1.5, xmax = 2, ymin = 3.5, ymax = 4.5) +
  add_fishape(family = "Acanthuridae", option = "Naso_lituratus",
              fill = col[4],
              xmin = 1.5, xmax = 2, ymin = 2.5, ymax = 3.5) +
  add_fishape(family = "Balistidae", option = "Odonus_niger",
              fill = col[5],
              xmin = 1.5, xmax = 2, ymin = 1.5, ymax = 2.5) +
  add_fishape(family = "Acanthuridae", option = "Zebrasoma_scopas",
              fill = col[6],
              xmin = 1.5, xmax = 2, ymin = 0.5, ymax = 1.5) +
  xlim(c(1.5, 3)) + ylim(c(0, 7)) +
  theme_ppt() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
l

combined <-
  a+b+c+d+e+f+g+h+i+j+k+l+m+n+map+ plot_layout(design = layout) & 
  theme(panel.background = element_rect(fill = cols()[1], color = cols()[1]), 
        plot.background = element_rect(fill = cols()[1], color = cols()[1]))

save_plot(combined, "c7_map", width = 2.2*12, height = 2.2*10)



l <-
  ggplot() +
  add_fishape(family = "Serranidae", option = "Cephalopholis_argus",
              fill ="#d7d7d7",
              ymin = 1.5, ymax = 2, xmin = 0.6, xmax = 1.4) +
  add_fishape(family = "Chaetodontidae", option = "Chaetodon_ornatissimus",
              fill = "#d7d7d7",
              ymin = 1.5, ymax = 2, xmin = 1.6, xmax = 2.4) +
  add_fishape(family = "Acanthuridae", option = "Ctenochaetus_striatus",
              fill = "#d7d7d7",
              ymin = 1.5, ymax = 2, xmin = 2.6, xmax = 3.4) +
  add_fishape(family = "Acanthuridae", option = "Naso_lituratus",
              fill = "#d7d7d7",
              ymin = 1.5, ymax = 2, xmin = 3.6, xmax = 4.4) +
  add_fishape(family = "Balistidae", option = "Odonus_niger",
              fill = "#d7d7d7",
              ymin = 1.5, ymax = 2, xmin = 4.6, xmax = 5.4) +
  add_fishape(family = "Acanthuridae", option = "Zebrasoma_scopas",
              fill = "#d7d7d7",
              ymin = 1.5, ymax = 2, xmin = 5.6, xmax = 6.4) +
  ylim(c(1.5, 3)) + xlim(c(0.5, 6.5)) +
  theme_ppt() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
l

save_plot(l, "c7_fishes", height = 4)
