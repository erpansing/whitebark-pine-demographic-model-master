library(dplyr)
library(ggplot2)
library(googlesheets)

load("/Users/elizabethpansing/Box Sync/Yellowstone/88-Fires-Analysis/2017 YNP Data.Rda")

file <- "1990-2017 YNP Non-whitebark Conifers"
file2 <- "2017 Non whitebark conifers"

others <- gs_title(file) %>% 
  gs_read(., sheet_title = "Sheet1") %>% 
  dplyr::rename(., DBH       = `Height DBH`,
                   StudyArea = Location) %>% 
  filter(., !is.na(PlotNo)) %>% 
  filter(., Species != "PIAL")

levels(as.factor(others$DBH))
levels(as.factor(others$Height))

others <- others %>% 
  mutate(., DBH = as.numeric(DBH),
            Height = as.numeric(Height),
            Age = as.numeric(Age)) %>% 
  # filter(., !is.na(DBH))%>% 
  filter(., Age < 30 | is.na(Age))


others_1617 <- gs_title(file2) %>% 
  gs_read(.)

levels(as.factor(others_1617$DBH))
sum(!is.na(others_1617$CBH))
levels(as.factor(others_1617$Height))
levels(as.factor(others_1617$Notes))


others_1617 <- others_1617 %>% 
  dplyr::rename(., StudyArea = Location) %>% 
  filter(., Notes != "PREFIRE" | is.na(Notes)) %>% 
  filter(., Notes != "DEAD" | is.na(Notes)) %>% 
  filter(., Notes != "PRE-FIRE"| is.na(Notes)) %>% 
  mutate(., DBH = as.numeric(DBH),
         Height = as.numeric(Height))

nonwb <- bind_rows(others, others_1617) 

sum(!is.na(nonwb$CBH))
sum(!is.na(nonwb$DBH)) == sum(nonwb$Height >= 137, na.rm = T)
sum(!is.na(nonwb$DBH) & is.na(nonwb$Height))
sum(nonwb$Height >= 137 & is.na(nonwb$Height), na.rm = T)

## Create dataframe that contains the number and density of non whitbark pine conifers per plot by year
density <- nonwb %>% 
  filter(., Year != 1996) %>% 
  # filter(., Year != 2005) %>% 
  mutate(., Year = ifelse(Year == 2005, 2001, Year)) %>%    ##MISSING 2001 DATA FOR HENDERSON. GUESSING 2005 DATA IS THAT, BUT NOT SURE
  group_by(., StudyArea, Year, PlotNo) %>%  # want to count by study area, plot no, and year
  tally(.) %>%                              # count the number of observations (i.e., non whitebark conifers by plot and year)
  ungroup(.) %>% 
  do({
    HM <- filter(., StudyArea == "Henderson") %>%   # Create df that has only HM plots to complete cases (i.e., add
      dplyr::select(., -StudyArea) %>%              # plots that didn't have any trees to the data). Zeros matter here
      group_by(., Year) %>% 
      complete(., PlotNo = c(1:150), fill = list(n = 0)) %>% 
      mutate(., StudyArea = "Henderson") %>% 
      ungroup(.)
    
    MW <- filter(., StudyArea == "Washburn") %>%    # Same as above for HM
      dplyr::select(., -StudyArea) %>% 
      group_by(., Year) %>% 
      complete(., PlotNo = c(200:299), fill = list(n = 0)) %>% 
      mutate(., StudyArea = "Washburn") %>% 
      ungroup(.)
    
    bind_rows(HM, MW)                              # Recreate main df
    
  }) %>% 
  mutate(., Density = n/20) %>%    # Create density column
  na.omit() %>% 
  filter(., !(Year == 2017 & PlotNo > 250)) %>%  # No data entry for 2017 plots 251- 299. Will include as these values are entered.
  mutate(., Year = ifelse(Year == 2016, 2017, Year))


stripchart( density$Density~as.factor(density$Year), vertical = T, method = "jitter", pch = 16)

# density %>% 
#   group_by(., StudyArea, Year) %>% 
#   summarise_at(., vars(Density), funs(median, min, max)) %>% 
#   ggplot(., aes(x = Year, y = median, col = StudyArea)) +
#   geom_point()+
#   geom_line() +
#   geom_ribbon(aes(x = Year, ymin = min, ymax = max), alpha = 0.2)+
#   scale_y_continuous(breaks = seq(0, 15, by = 0.5)) +
#   facet_wrap(~StudyArea)


# density %>% 
#   group_by(.,Year) %>% 
#   summarise_at(., vars(Density), funs(median, min, max)) %>% 
#   ggplot(., aes(x = Year, y = median)) +
#   geom_point()+
#   geom_line() +
#   geom_ribbon(aes(x = Year, ymin = min, ymax = max), alpha = 0.2)+
#   scale_y_continuous(breaks = seq(0, 15, by = 0.5)) 




Others <- data.frame(PICO = c( rep(0, 7),                 54, 29, 33,   1, 59),
                     PIEN = c( 6,  5, 20,  9, 10, 43, 10,  3, 12, 43, 112,  2),
                     ABLA = c(78, 81, 46, 34, 55, 33, 21, 12,  9,  4,   6, 13),
                     UNK  = c( 6,  4, 14,  4,  0,  1,  1, 25, 32, 27,   6, 43),
                     StudyArea = c(rep("Henderson", 7), rep("Washburn", 5)))




Others <- Others %>% 
  mutate(n = rowSums(dplyr::select(., -StudyArea))) %>% 
  mutate(Density = n/(30*30)) %>%  # per plot overall density for all non-whitebark conifers
  mutate(., PlotNo = 1:nrow(.)) %>%
  mutate(., Year = 1988 + 230) %>% 
  dplyr::select(., Year, StudyArea, PlotNo, n, Density) %>% 
  mutate(., n = round(Density * 20, digits = 0))
    

hist(Others$Density)

dens <- bind_rows(density,Others) %>% 
  mutate(., firetime = Year - 1989)

# Year <- c(1:3, 5, 6, 12, 28, 229)
# out <- (0.16 * Year)/(1 + Year)

boxplot(dens$Density ~ dens$firetime, vertical = T, method = "jitter", pch = 16)
# points(as.factor(Year), out, col = "red", pch = 19)
# 
# med <- dens %>% 
#   group_by(., Year)%>% 
#   summarise_at(., vars(Density), funs(median)) %>% 
#   ungroup()
# 
# plot(med$Year, med$Density)

library(pscl)

m1 <- zeroinfl(n ~ firetime,
               data = dens, dist = "negbin", EM = TRUE)
summary(m1)

newdata1 <- data.frame(firetime = dens$firetime,
                       phat = predict(m1)) %>% 
  mutate(., density = phat/20)


m2 <- glm(n ~ firetime, data = dens, family = poisson)
summary(m2)

newdata2 <- data.frame(firetime = dens$firetime,
                       phat = predict(m2)) %>%
  mutate(., density = phat/20)


ggplot(newdata1, aes(x = firetime, y = phat)) +
  geom_point() +
  geom_line() +
  labs(x = "Time since fire (Years)", y = "Predicted number of non-wb conifers/20m2")


ggplot2::ggplot(data = dens, aes(x = as.factor(firetime), y = Density)) +
  ggplot2::geom_boxplot() +
  geom_point(data = newdata1, aes(x = as.factor(firetime), y = density), colour = "blue") +
  geom_point(data = newdata2, aes(x = as.factor(firetime), y = density), colour = "red") +
  labs(x = "Time since fire (Years)", y = "Predicted density of non-wb conifers/m2")

# lambda_dens <- m2$coefficients[1] + m2$coefficients[2] * firetime
# LAIb <- rpois(n = 1, lambda = lambda_dens) * area


rm(list=setdiff(ls(), "m2"))

############################################################################################################
## Calculate the mean, median, and variance of plot density for each year. Don't care about study area or treatment designations

# pial_dens_postfire <- pial %>%  ## Calcluate per plot density by year
#   filter(., Height > 137) %>% 
#   group_by(., PlotNo, Year) %>%
#   tally(.) %>% 
#   ungroup(.) %>% 
#   mutate(., Density = n/20) %>% 
#   dplyr::select(., PlotNo, Year, Density)
#   
# # Plot the histogram of density by year... looks lognormally distributed
# ggplot(data = pial_dens_postfire, aes(x = Density))+  
#   geom_histogram(binwidth = 0.05)+
#   facet_wrap(~ Year)
# 
# ggplot(data = pial_dens_postfire, aes(x = Density))+  
#   geom_density()+
#   facet_wrap(~ Year)
# 
# ggplot(data = pial_dens_postfire, aes(x = log(Density)))+  ## But log transforming  it doesn't look normal, so gamma?
#   geom_density()+
#   facet_wrap(~ Year)
# 
# ggplot(data = pial_dens_postfire, aes(x = log(Density)))+  ## Does here, so either?
#   geom_histogram(binwidth = 0.3)+
#   facet_wrap(~ Year)
# 
# 
# pial_dens_postfire_summary <- pial_dens_postfire %>% 
#   group_by(., Year) %>% 
#   summarise_at(., vars(Density), funs(median, mean, var, min, max)) %>% 
#   ungroup(.)
#   
# 
# ggplot(data = pial_dens_postfire_summary, aes(x = Year, y = median))+
#   geom_line() +
#   geom_ribbon(data = pial_dens_postfire_summary, aes(x = Year, ymin = min, ymax = max), alpha = 0.2)
# 
# 
# 
# pial_dens_prefire <- data.frame(Year    = rep(1988+230, 12),
#                                 PlotNo  = 1:12,
#                                 Density = c(0.001, 0.003, 0,0, 0.009, 0.003, 0.001, 0.016, 0, 0.001, 0.001, 0.010))
# 
# pial_prefire_summary <- data.frame(Year   = 1988+230, 
#                                    median = median(pial_dens_prefire), 
#                                    mean   = mean(pial_dens_prefire),
#                                    var    = var(pial_dens_prefire),
#                                    min    = min(pial_dens_prefire),
#                                    max    = max(pial_dens_prefire))
# 
# 
# pial_density <- rbind(pial_dens_postfire, pial_dens_prefire)
# 
# fit <- glm(Density~ Year, data = pial_density, family = Gamma)
# summary(fit)
# 
# 
# data <- density %>% 
#   group_by(., Year) %>% 
#   summarise_at(., vars(Density), funs(median)) %>% 
#   ungroup() %>% 
#   as.data.frame()
# 
# data <- data[-10,]
# 
# library(drm)
# 
# test <- drm(n ~ Year + cluster(id), data = density, fct = MM.2())
# 
# model.nls <- nls(v ~ Vm * S/(K+S), data = mm, 
#                  start = list(K = max(mm$v)/2, Vm = max(mm$v)))

# #################################################################################################################  
#   
# load("/Users/elizabethpansing/Box Sync/Yellowstone/88-Fires-Analysis/2017 YNP Data.Rda")
# file <- "1990-2017 YNP Non-whitebark Conifers"
# others <- gs_title(file) %>%
#   gs_read(., ws = 1) %>%
#   dplyr::rename(., DBH = `Height DBH`)
# others <- others %>%
#   mutate(., DBH = as.numeric(DBH),
#          Height = as.numeric(Height),
#          Age = as.numeric(Age))



