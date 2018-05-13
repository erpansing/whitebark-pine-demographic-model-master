library(dplyr)
library(ggplot2)
library(googlesheets)

load("/Users/elizabethpansing/Box Sync/Yellowstone/88-Fires-Analysis/2017 YNP Data.Rda")

file <- "1990-2017 YNP Non-whitebark Conifers"

others <- gs_title(file) %>% 
  gs_read(., ws = 1) %>% 
  dplyr::rename(., DBH = `Height DBH`)

levels(as.factor(others$DBH))
levels(as.factor(others$Height))

others <- others %>% 
  mutate(., DBH = as.numeric(DBH),
            Height = as.numeric(Height),
            Age = as.numeric(Age)) %>% 
  filter(., !is.na(DBH))%>% 
  filter(., Age < 30 | is.na(Age))

others_1617 <- gs_title(file) %>% 
  gs_read(., ws = 8) 

levels(as.factor(others_1617$DBH))
sum(!is.na(others_1617$CBH))
levels(as.factor(others_1617$Height))


others_1617 <- others_1617 %>% 
  mutate(., DBH = as.numeric(DBH),
            Height = as.numeric(Height)) %>% 
  filter(., Notes != "PREFIRE" | is.na(Notes)) %>% 
  filter(., Notes != "DEAD" | is.na(Notes)) %>% 
  filter(., Notes != "PRE-FIRE"| is.na(Notes))

nonwb <- bind_rows(others, others_1617) %>% 
  filter(., Height >= 137 | !is.na(DBH)) 

sum(!is.na(nonwb$CBH))
sum(!is.na(nonwb$DBH)) == sum(nonwb$Height >= 137, na.rm = T)
sum(!is.na(nonwb$DBH) & is.na(nonwb$Height))
sum(nonwb$Height >= 137 & is.na(nonwb$Height), na.rm = T)


nonwb <- nonwb %>% 
  group_by(., Year, PlotNo) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(., Density = n/20)

stripchart( nonwb$Density~as.factor(nonwb$Year), vertical = T, method = "jitter", pch = 16)

Others <- data.frame(PICO = c( rep(0, 7), 54,  29, 33,   1, 59),
                     PIEN = c( 6,  5, 20,  9, 10, 43, 10,  3, 12, 43, 112,  2),
                     ABLA = c(78, 81, 46, 34, 55, 33, 21, 12,  9,  4,   6, 13),
                     UNK  = c( 6,  4, 14,  4,  0,  1,  1, 25, 32, 27,   6, 43))




other_density <- rowSums(Others)/(300*300) # per plot overall density for all non-whitebark conifers
hist(other_density)

nonwb <- bind_rows(nonwb, data.frame(Year    = rep(1988 + 230, length(other_density)),
                                     PlotNo  = 1:length(other_density),
                                     Density = other_density))

stripchart(nonwb$Density ~ nonwb$Year, vertical = T, method = "jitter", pch = 16)



############################################################################################################
## Calculate the mean, median, and variance of plot density for each year. Don't care about study area or treatment designations

pial_dens_postfire <- pial %>%  ## Calcluate per plot density by year
  filter(., Height > 137) %>% 
  group_by(., PlotNo, Year) %>%
  tally(.) %>% 
  ungroup(.) %>% 
  mutate(., Density = n/20) %>% 
  dplyr::select(., PlotNo, Year, Density)
  
# Plot the histogram of density by year... looks lognormally distributed
ggplot(data = pial_dens_postfire, aes(x = Density))+  
  geom_histogram(binwidth = 0.05)+
  facet_wrap(~ Year)

ggplot(data = pial_dens_postfire, aes(x = Density))+  
  geom_density()+
  facet_wrap(~ Year)

ggplot(data = pial_dens_postfire, aes(x = log(Density)))+  ## But log transforming  it doesn't look normal, so gamma?
  geom_density()+
  facet_wrap(~ Year)

ggplot(data = pial_dens_postfire, aes(x = log(Density)))+  ## Does here, so either?
  geom_histogram(binwidth = 0.3)+
  facet_wrap(~ Year)


pial_dens_postfire_summary <- pial_dens_postfire %>% 
  group_by(., Year) %>% 
  summarise_at(., vars(Density), funs(median, mean, var, min, max)) %>% 
  ungroup(.)
  

ggplot(data = pial_dens_postfire_summary, aes(x = Year, y = median))+
  geom_line() +
  geom_ribbon(data = pial_dens_postfire_summary, aes(x = Year, ymin = min, ymax = max), alpha = 0.2)



pial_dens_prefire <- data.frame(Year    = rep(1988+230, 12),
                                PlotNo  = 1:12,
                                Density = c(0.001, 0.003, 0,0, 0.009, 0.003, 0.001, 0.016, 0, 0.001, 0.001, 0.010))

pial_prefire_summary <- data.frame(Year   = 1988+230, 
                                   median = median(pial_dens_prefire), 
                                   mean   = mean(pial_dens_prefire),
                                   var    = var(pial_dens_prefire),
                                   min    = min(pial_dens_prefire),
                                   max    = max(pial_dens_prefire))


pial_density <- rbind(pial_dens_postfire, pial_dens_prefire)

fit <- glm(Density~ Year, data = pial_density, family = Gamma)
summary(fit)


#################################################################################################################  
  



