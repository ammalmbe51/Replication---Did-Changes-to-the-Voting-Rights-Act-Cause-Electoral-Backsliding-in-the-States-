devtools::install_github("synth-inference/synthdid")
library(estimatr)
library(stargazer)
library(tidyverse)
library(devtools)
library(foreign)
library(dplyr)
library(sandwich)
library(lmtest)
library(WVPlots)
library(ggplot2)
library(fixest)
library(DescTools)
library(margins)
library(ggpubr)
library(synthdid)

did_data <- read.csv("did.csv")
View(did_data)

did_data$fullpreclearance <- ifelse(did_data$st == "AL" | did_data$st == "AK" | did_data$st == "AZ" | did_data$st == "GA" | did_data$st == "LA" | did_data$st == "MS" | did_data$st == "SC" | did_data$st == "TX" | did_data$st == "VA", 1, 0)
table(did_data$fullpreclearance) #171 treated, 665 nontreated

#Recode time variable - time since Shelby County
did_data$postperiod <- ifelse(did_data$year >= "2013", 1, 0)

#Pre and post period groups
did_data$treatedpost <- ifelse(did_data$postperiod == 1 & did_data$preclearance == 1, 1, 0)
did_data$treatedpre <- ifelse(did_data$postperiod == 0 & did_data$preclearance == 1, 1, 0)
did_data$untreatedpost <- ifelse(did_data$postperiod == 1 & did_data$preclearance == 0, 1, 0)
did_data$untreatedpre <- ifelse(did_data$postperiod == 0 & did_data$preclearance == 0, 1, 0)

#only fully covered and noncovered states - exclude partially covered
nopartial <- did_data %>% 
  dplyr::filter(!st %in% c("CA", "FL", "MI","NY","NC", "SD"))
table(nopartial$offyear)
dim(nopartial) #leaves 836 obs

#Recode fully covered states w DV 1 if entire state is subject to preclerance
nopartial$preclearance <- ifelse(nopartial$st == "AL" | nopartial$st == "AK" | nopartial$st == "AZ" | nopartial$st == "GA" | nopartial$st == "LA" | nopartial$st == "MS" | nopartial$st == "SC" | nopartial$st == "TX" | nopartial$st == "VA", 1, 0)
table(nopartial$preclearance) #171 treated, 665 nontreated

##Partially covered on noncovered (exclude fully covered)
nofull <- did_data %>% 
  dplyr::filter(!st %in% c("AL", "AK", "AZ","GA","LA", "MS", "SC", "TX", "VA"))
table(nofull$offyear)
dim(nofull) #779 obs

#Partially treated states recoded as 1
nofull$partialpre <- ifelse(nofull$st == "CA" | nofull$st == "FL" | nofull$st == "MI" | nofull$st == "NY" | nofull$st == "NC" | nofull$st == "SD", 1, 0)
table(nofull$partialpre) #665 control, 114 treated

#Pre and post period groups
nofull$treatedpost <- ifelse(nofull$postperiod == 1 & nofull$partialpre == 1, 1, 0)
nofull$treatedpre <- ifelse(nofull$postperiod == 0 & nofull$partialpre == 1, 1, 0)
nofull$untreatedpost <- ifelse(nofull$postperiod == 1 & nofull$partialpre == 0, 1, 0)
nofull$untreatedpre <- ifelse(nofull$postperiod == 0 & nofull$partialpre == 0, 1, 0)

##SynthDiD - Models 1 and 3
#new dataframe for panel matrix and weighted DiD
#fully covered only - Model 1 ATT
weightedshelby2 <- data.frame(nopartial$state, nopartial$year, nopartial$democracy_mcmc_lib, nopartial$treatedpost)
head(weightedshelby2)

#compute weights
shelby.setup <- panel.matrices(weightedshelby2)
summary(shelby.setup$T0)

shelby.setup
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.

#Model 1 ATT
tau.hat2 <- synthdid_estimate(shelby.setup$Y, shelby.setup$N0, shelby.setup$T0)
tau.hat2
se <- 0.173
#confidence intervals
sprintf('95%% CI (%1.2f, %1.2f)', tau.hat2 - 1.96 * se, tau.hat2 + 1.96 * se)

#parallel trends plot for fully covered states - Figure 1 
plot(tau.hat2, overlay=0.5, se.method = 'bootstrap', control.name='Uncovered States', treated.name='Fully Covered States') +
  labs(x = "Year", y = "Average SDI Estimates") 

##partially covered states synthetic DiD - Model 2
weightedshelby3 <- data.frame(nofull$state, nofull$year, nofull$democracy_mcmc_lib, nofull$treatedpost)
shelby.setup.part <- panel.matrices(weightedshelby3)

#ATT estimate for Model 2
tau.hatpart <- synthdid_estimate(shelby.setup.part$Y, shelby.setup.part$N0, shelby.setup.part$T0)
tau.hatpart
se.part <- 0.349
sprintf('95%% CI (%1.2f, %1.2f)', tau.hatpart - 1.96 * se.part, tau.hatpart + 1.96 * se.part)

#Partially covered states parallel trends plot - Figure 2
plot(tau.hatpart, overlay=0.5, se.method = 'bootstrap', control.name='Uncovered States', treated.name='Partially Covered States') +
  labs(x = "Year", y = "Average SDI Estimates") 

#TWFEs - fully covered on noncovered (no partial coverage) - Model 3 ATT
TWFEs <- feols(democracy_mcmc_lib ~ i(treatedpost, preclearance, ref = 0) | year + state,
               cluster = "state", data = nopartial)
summary(TWFEs) #ATT estimate for Model 3

#TWFEs - partial period - Model 4 ATT
TWFEspartial <- feols(democracy_mcmc_lib ~ i(treatedpost, partialpre, ref = 0) | year + state,
                      cluster = "state", data = nofull)
summary(TWFEspartial) #Model 4 ATT

#event study for fully covered states only
esfull <- feols(democracy_mcmc_lib ~ i(year, preclearance, ref = 2013) | year + state,
                cluster = "state", data = nopartial)
eventstudyfull <- iplot(esfull, 
                        xlab = 'Year',
                        main = 'Event study: State-Level Treatment Estimates (with TWFEs) \n Fully Covered States')
eventstudyfull

#event study for partially covered states only
espartial <- feols(democracy_mcmc_lib ~ i(year, partialpre, ref = 2013) | year + state,
                   cluster = "state", data = nofull)
eventstudypart <- iplot(espartial, 
                        xlab = 'Year',
                        main = 'Event study: State-Level Treatment Estimates (with TWFEs) \n Partially Covered States')
eventstudypart

#plot distribution of control unit weights - Fully covered states
#All - Figure B1
fullunitplotall <-  synthdid_units_plot(tau.hat2) + labs(x = "State", y = "Change in Weighted SDI Estimates") +
  theme(legend.background=element_blank(), strip.background=element_blank())
fullunitplotall

#top 10 - Figure B2
fullunitplot <- synthdid_units_plot(tau.hat2, units = rownames(top.controls)) + labs(x = "State", y = "Change in Weighted SDI Estimates") +
  theme(legend.background=element_blank(), strip.background=element_blank(), axis.text = element_text(size = 12, face= "bold"))
fullunitplot

#plot distribution of control unit weights - partially covered states
#All - Figure B3
partialunitplotall <-  synthdid_units_plot(tau.hatpart) + labs(x = "State", y = "Change in Weighted SDI Estimates") +
  theme(legend.background=element_blank(), strip.background=element_blank())
partialunitplotall

#top 10 - Figure B4
partialunitplot <- synthdid_units_plot(tau.hatpart, units = rownames(top.controls)) + labs(x = "State", y = "Change in Weighted SDI Estimates") +
  theme(legend.background=element_blank(), strip.background=element_blank(), axis.text = element_text(size = 12, face= "bold"))
partialunitplot

