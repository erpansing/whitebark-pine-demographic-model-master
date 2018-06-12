## Convert to nest survival data frame
# rm(list = ls())  # only use when playing with this script, and be sure to comment out before saving as it can ruin dependent scripts

library(dplyr)
library(RMark)
library(ggplot2)
library(reshape2)

load("/Users/elizabethpansing/Box Sync/Yellowstone/88-Fires-Analysis/2017 YNP Data.Rda")

#Omit IDs that only have one entry that is "D"

StatusSummarybyID <- pial %>%        # create table of number of statuses 
  group_by(., IDNumber, Status) %>%  # for each ID
  tally() %>%
  ungroup(.) 

StatusSummarybyID <- dcast(StatusSummarybyID, IDNumber ~ Status, value.var = "n", fill = 0) %>%
  filter(., D > 0 & L == 0)  # ID trees with only one status (D)

IDs <- StatusSummarybyID$IDNumber  # Collect IDNumbers as a vector

pial_SD <- pial %>%
  filter(., !IDNumber %in% IDs) # Remove those trees from dataframe

rm(IDs, pial, StatusSummarybyID) # clean up environment

## Remove IDs with duplicate ID_Year combinations
ID <- unique(pial_SD$IDNumber[duplicated(pial_SD$ID_Year)]) # ID trees
  # with identical ID_Year combinations

pial_SD <- pial_SD %>%
  filter(., !IDNumber %in% ID) # Omit these trees. Once QC is complete, this
  # step will be unnecessary

rm(ID)  # clean up environment

## Remove IDs of pre-fire regen

ID <- unique(pial_SD$IDNumber[pial_SD$YearGerminated < 1990]) #ID prefire trees

pial_SD <- pial_SD %>%
  filter(., !IDNumber %in% ID) # Omit these trees. 

rm(ID)  # clean up environment

pial_SD$Status[is.na(pial_SD$Status)] <-  "D"

# 
# FirstFound <- pial_SD %>%
#   group_by(., IDNumber) %>%
#   filter(., Year == min(Year)) %>% 
#   ungroup(.) %>%
#   dplyr::rename(FirstFound = Year) %>%
#   select(., IDNumber, FirstFound) %>%
#   mutate(., FirstFound = FirstFound - 1989)
#   
# LastChecked <- pial_SD %>%  
#   group_by(., IDNumber) %>%
#   filter(., Year == max(Year)) %>%
#   ungroup(.) %>%
#   dplyr::rename(LastChecked = Year) %>%
#   select(., IDNumber, LastChecked) %>%
#   mutate(., LastChecked = LastChecked - 1989)
# 
# LastPresent <- pial_SD %>%  
#   group_by(., IDNumber) %>%
#   filter(., Status == "L") %>%
#   filter(., Year == max(Year)) %>%
#   ungroup(.) %>%
#   dplyr::rename(LastPresent = Year) %>%
#   select(., IDNumber, LastPresent) %>%
#   mutate(LastPresent = LastPresent - 1989) 
# 
# 
# AgeFound <- pial_SD %>%
#   group_by(., IDNumber) %>%
#   filter(., Year == min(Year)) %>%
#   ungroup(.) %>%
#   mutate(., AgeFound = Age) %>%
#   select(., IDNumber, AgeFound)
# 
# 
# Fate <- pial_SD %>%
#   group_by(., IDNumber) %>%
#   filter(., Year == max(Year)) %>%
#   mutate(., Fate = ifelse(Status == "L", 0, 1)) %>%
#   select(., IDNumber, Fate) 
# 
# pial_SD_ns <- merge(FirstFound, LastChecked, by = "IDNumber") %>%
#   merge(., LastPresent) %>% 
#   merge(., Fate) %>% 
#   merge(., AgeFound) %>%
#   # mutate(., Freq = 1) %>%
#   select(., FirstFound, LastPresent, LastChecked, 
#             Fate, AgeFound, IDNumber) %>%
#   filter(., AgeFound >= 0 | is.na(AgeFound)) %>%
#   filter(., !(FirstFound == LastPresent))
#   
# rm(AgeFound, Fate, FirstFound, LastChecked, LastPresent, pial_SD)
# 
# 
# run.wbp.models=function()
# {
#   # 1. A model of constant survival rate (SR)
#   
#   Dot = mark(pial_SD_ns, nocc = 28, model="Nest",
#            model.parameters = list(S = list(formula = ~1)))
#   
#   # 2. SR varies as a function of time  
#   
#   Time = mark(pial_SD_ns, nocc = 28, model = "Nest",
#             model.parameters = list(S = list(formula = ~ Time)))
#   
#   # # 3. SR varies with discrete time
# 
#   # time = mark(pial_SD_ns, nocc = 28, model = "Nest",
#   #           model.parameters = list(S = list(formula = ~ time)))
# 
#   # Return model table and list of models
#   #
#   return(collect.models() )
# }
# 
# wbp.results=run.wbp.models()
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# # Examine table of model-selection results #
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# wbp.results                        # print model-selection table to screen
# 
# SD_survival_mean <- wbp.results$Dot$results$real$estimate
# SD_survival_var  <- nrow(pial_SD_ns) * (wbp.results$Dot$results$real$se)^2
# 
# # 
# rm(pial_SD_ns,wbp.results, run.wbp.models)

# predictions <- data.frame(time  = seq(from = 1991, to = 2017, by = 1),
#                           prob = wbp.results$Dot$results$real$estimate,
#                           lcl  = wbp.results$Dot$results$real$lcl,
#                           ucl  = wbp.results$Dot$results$real$ucl)
# #
# #
# ggplot(predictions, aes(x = time, y = prob)) +
#   # geom_point() +
#   geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.5, fill = "#ef8a62") +
#   geom_line(size = 0.5) +
#   xlab("Year") +
#   ylab("Annual survival rate") +
#   scale_y_continuous(limits = c(0.75,1)) +
#   scale_x_continuous(limits = c(1990, 2017),
#                      breaks = seq(1990, 2017, 2)) +
#   theme(axis.title = element_text(size = 18, face = "bold")) +
#   theme(axis.text = element_text(size = 15))
  
##--------------------------------------------------------------------##
##                                                                    ##
##                      Age > 10 considered success                   ##
##                                                                    ##
##--------------------------------------------------------------------##

olduns <- pial_SD %>% filter(., YearGerminated < 2006)


FirstFound_10 <- olduns %>%
  dplyr::group_by(., IDNumber) %>%
  dplyr::filter(., Year == min(Year)) %>% 
  dplyr::ungroup(.) %>%
  dplyr::rename(., FirstFound = Year) %>%
  dplyr::select(., IDNumber, FirstFound) %>%
  dplyr::mutate(., FirstFound = FirstFound - 1989)

LastChecked_10 <- olduns %>%  
  dplyr::group_by(., IDNumber) %>%
  dplyr::filter(., Year == max(Year)) %>%
  dplyr::ungroup(.) %>%
  dplyr::rename(LastChecked = Year) %>%
  dplyr::select(., IDNumber, LastChecked) %>%
  dplyr::mutate(., LastChecked = LastChecked - 1989)

LastPresent_10 <- olduns %>%  
  dplyr::group_by(., IDNumber) %>%
  dplyr::filter(., Status == "L") %>%
  dplyr::filter(., Year == max(Year)) %>%
  dplyr::ungroup(.) %>%
  dplyr::rename(LastPresent = Year) %>%
  dplyr::select(., IDNumber, LastPresent) %>%
  dplyr::mutate(LastPresent = LastPresent - 1989) 


AgeFound_10 <- olduns %>%
  dplyr::group_by(., IDNumber) %>%
  dplyr::filter(., Year == min(Year)) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(., AgeFound = Age) %>%
  dplyr::select(., IDNumber, AgeFound)

Fate_10 <- olduns %>%
  dplyr::group_by(., IDNumber) %>%
  dplyr::filter(., Age >= 10) %>%
  dplyr::mutate(., Fate = ifelse(Status == "L", 0, 1)) %>%
  dplyr::select(., IDNumber, Fate) 


pial_SD_ns_10 <- merge(FirstFound_10, LastChecked_10, by = "IDNumber") %>%
  merge(., LastPresent_10) %>% 
  merge(., Fate_10) %>% 
  merge(., AgeFound_10) %>%
  # mutate(., Freq = 1) %>%
  dplyr::select(., FirstFound, LastPresent, LastChecked, 
         Fate, AgeFound, IDNumber) # %>%
  # filter(., AgeFound >= 0 | is.na(AgeFound)) %>%
  # filter(., !(FirstFound == LastPresent))


run.wbp.models_10=function()
{
  # 1. A model of constant survival rate (SR)
  
  Dot = mark(pial_SD_ns_10, nocc = 28, model="Nest",
             model.parameters = list(S = list(formula = ~1)))
  
  # 2. SR varies as a function of time  
  
  Time = mark(pial_SD_ns_10, nocc = 28, model = "Nest",
              model.parameters = list(S = list(formula = ~ Time)))
  
  # 3. SR varies with discrete time

   time = mark(pial_SD_ns_10, nocc = 28, model = "Nest",
             model.parameters = list(S = list(formula = ~ time)))
   
   # 4. SR varies with Age
   
   age = mark(pial_SD_ns_10, nocc = 28, model = "Nest",
              model.parameters = list(S = list(formula =  ~AgeFound)))

  # Return model table and list of models

  return(collect.models() )
}

wbp.results_10=run.wbp.models_10()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Examine table of model-selection results #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

wbp.results_10                        # print model-selection table to screen

SD_survival_mean <- wbp.results_10$Dot$results$real$estimate
SD_survival_var  <- nrow(pial_SD_ns_10) * (wbp.results_10$Dot$results$real$se)^2


# predictions <- data.frame(time  = seq(from = 1991, to = 2017, by = 1),
#                           prob = wbp.results_10$Dot$results$real$estimate,
#                           lcl  = wbp.results_10$Dot$results$real$lcl,
#                           ucl  = wbp.results_10$Dot$results$real$ucl)
# #
# #
# ggplot(predictions, aes(x = time, y = prob)) + 
#   # geom_point() + 
#   geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.5, fill = "#ef8a62") +
#   geom_line(size = 0.5) +
#   xlab("Year") +
#   ylab("Annual survival rate") +
#   scale_y_continuous(limits = c(0.75,1)) +
#   scale_x_continuous(limits = c(1990, 2017),
#                      breaks = seq(1990, 2017, 2)) +
#   theme(axis.title = element_text(size = 18, face = "bold")) +
#   theme(axis.text = element_text(size = 15)) 
# 
# 
# 

# rm(list = ls()[!ls() %in% c("SD_survival_mean", "SD_survival_var")])

del <- paste0("10|old|pial|^del$")

rm(list = ls(pattern = del))

   
      