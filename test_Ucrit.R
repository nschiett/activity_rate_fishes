library(tidyverse)

data <- read_csv("data/data_fulton_2007.csv") %>%
  mutate(Ucrit_tls = Ucrit_cms/Length_cm) 

speedmax <- data.frame(
  species = c("Zebrasoma scopas", "Naso lituratus", 
              "Ctenochaetus striatus", "Chaetodon ornatissimus",
              "Cephalopholis argus", "Odonus niger"),
  Ucrit_tls = c(
    mean(data[data$Species == "Zebrasoma scopas", "Ucrit_tls"]),
    mean(data[data$Genus == "Naso", "Ucrit_tls"]),
    mean(data[data$Species == "Ctenochaetus striatus", "Ucrit_tls"]),
    mean(data[data$Species == "Chaetodon lunulatus", "Ucrit_tls"]),
    mean(data[data$Genus == "Cephalopholis", "Ucrit_tls"]),
    6.41))
    
ggplot(data) +
  geom_point(aes(x = Length_cm, y = Ucrit_tls, color = Species))


fit <- brm(Ucrit_cms ~ Length_cm + (1|Species), data = data, backend = "cmdstanr", threads = threading(4))
summary(fit)

nd <- speedmax <- data.frame(
  Species = rep(c("Zebrasoma scopas", "Naso brevirostris", 
              "Ctenochaetus striatus", "Chaetodon lunulatus",
              "Cephalopholis boenak"), each = 20),
  Length_cm = rep(6:25, 5))

pred <- fitted(fit, newdata = nd)  

pred <- cbind(nd, pred)

ggplot(pred) +
  geom_point(aes(x = Length_cm, y = Estimate, color = Species)) +
  geom_abline(slope = 1.77, intercept = 39.5)
  

loadd(mod_speed)

max <- cbind(mod_speed$data, residuals(mod_speed)) %>%
  arrange(desc(Estimate)) %>% 
  dplyr::group_by(Species) %>%
  slice(1:3) %>%
  dplyr::summarize(max = mean(Estimate))



                 