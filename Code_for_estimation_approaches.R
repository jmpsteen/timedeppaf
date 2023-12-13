###############################################################################################################################################
##                                                                                                    
##    Supplementary R code for                                                              
##    "Dealing with time-dependent exposures and confounding when defining and estimating attributable fractions       
##     -- Revisiting estimands and estimators"
##    doi: 10.1002/sim.9988
##                                                                                                    
##    Author: Johan Steen                                                                             
##    Date: December 13, 2023           
##
##    ------
##
##    The datasets.RData file contains three data frames
##
##     events: wide format data frame which includes the following event-related variables
##         - id: unique patient identifier
##         - T: time from ICU admission to either ICU death or discharge, whichever occurs first (expressed in days)
##         - epsilon: event type indicator (1 = ICU death, 2 = ICU discharge)
##         - C: time from ICU admission to infection onset (expressed in days)
##         - delta: infection status at ICU death or discharge (1 = infected, 0 = uninfected)
##
##     baseline: wide format data frame which includes the following baseline patient characteristics
##         - id: unique patient identifier
##         - sex: F = female, M = male
##         - age: age (in years) at ICU admission
##         - admissionyr: year of ICU admission (coded as factor)
##         - admissioncat: admission category ("Medicine", "Emergency surgery" or "Scheduled surgery")
##         - apache2: APACHE II score at ICU admission
##         - SOFA_total_0: total SOFA score at ICU admission
##         - CCI: Charlson Comorbidity Index
##
##     timedep: long format data frame which includes the following time-dependent patient information
##         - id: unique patient identifier
##         - tstop: end time of the corresponding 24-h time interval (with follow-up starting at the time of ICU admission)
##         - SOFA_total: daily total SOFA score
##                                                                                                   
###############################################################################################################################################

rm(list=ls())

### load required packages ----
library(dplyr) # for convenient data wrangling tools
library(tidyr) # for LOCF functionality
library(magrittr) # for piping
library(survival) # for survival analysis functions
library(ipw) # for inverse probability of censoring weighting
library(splines) # for fitting splines

# optional packages (for plots)
library(prodlim)
library(ggplot2)
library(ggfortify)


### load data ----
load("./datasets.RData")


### transform the event time data to counting process format ----

## calculate tilde_T and tilde_epsilon as defined in the main text
events %<>% mutate(tilde_T = pmin(T,C), 
                   tilde_epsilon = delta*epsilon)

## calculate discrete time variables (under both competing risk models in Fig 1A and 1B in the main text) 
## to enable transformation to counting process format later on
## and recode event type indicators as labelled factors for convenience and clarity
events %<>% mutate(time1 = ceiling(T),
                   time2 = ceiling(tilde_T),
                   timeC = ifelse(delta == 1, ceiling(tilde_T), floor(tilde_T)),
                   event1 = factor(epsilon, levels = 0:2, labels = c("censor", "death", "discharge")),
                   event2 = factor(tilde_epsilon, levels = -1:2, labels = c("censor", "HAI onset", "HAI-free death", "HAI-free discharge")),
                   eventC = factor(tilde_epsilon, levels = 0:2, labels = c("censor", "death", "discharge")),
                   A_E = 1-delta)
#   Note 1: these new time variables encode the end of the considered 24-h time intervals (with follow-up starting at the time of ICU admission)
# except for timeC, which can be considered a copy of time1 and time2, but with (artificial) censoring due to HAI onset at the start rather than the end of the 24-h time interval
# to enforce the temporal ordering assumption that censoring precedes ICU death or discharge (as described in section 2 of the main text)    
#   Note 2: it is important that the first factor level in the event type variables always corresponds to 'censoring'!
# Even if no censoring occurs, a factor level for censoring should be added and explicitly coded as first factor level.
# If not, the Surv() function may not provide the intended result when dealing with competing events.

# check
# with(events, table(event1, epsilon, useNA = "ifany"))
# with(events, table(event2, tilde_epsilon, useNA = "ifany"))
# with(events, table(eventC, tilde_epsilon, useNA = "ifany"))

## use the survival::survSplit function to transform to counting process format, setting the end of follow-up intervals (day level) as pre-specified cut times  
CPdata1 <- survSplit(Surv(time1, event1) ~ ., data = events %>% select(id, time1, event1, delta, A_E), cut = 0:60, end = "time1") %>% rename(tstop = time1, event = event1)
CPdata2 <- survSplit(Surv(time2, event2) ~ ., data = events %>% select(id, time2, event2), cut = 0:60, end = "time2") %>% rename(tstop = time2, event = event2)
CPdata <- full_join(CPdata1 %>% rename(event1 = event), CPdata2 %>% rename(event2 = event) %>% select(-tstart), by = c("id", "tstop"))
rm(CPdata1, CPdata2)
#   Note: to enable calculation of all the discrete time counting process variables defined in section 2 of the main text, 
# we create a counting process format dataframe for the event times and types of both competing risk models in Figure 1A and 1B,
# which are then joined

## calculate the discrete time counting process variables defined in section 2 of the main text
CPdata %<>% mutate(k = tstop,
                   A_k = as.numeric(event2 == "HAI onset"),
                   D_k = as.numeric(event1 == "discharge"),
                   Y_k = as.numeric(event1 == "death"))
CPdata %<>% group_by(id) %>% fill(A_k) 
#   Note: ensure that A_k is coded 1 at each follow-up interval following the interval at which HAI onset occurs (LOCF for missing A_k)
# CPdata %<>% mutate(tilde_D_k = (1-A_k)*D_k,
#                    tilde_Y_k = (1-A_k)*Y_k)
# calculate lag variables
CPdata %<>% group_by(id) %>% mutate(A_k_lag = lag(A_k, default = 0),
                                    D_k_lag = lag(D_k, default = 0),
                                    Y_k_lag = lag(Y_k, default = 0))


### merge dataframes ----
CPdata %<>% full_join(baseline, by = "id")
CPdata %<>% left_join(timedep, by = c("id", "tstop"))

## create lagged SOFA_total variable
CPdata %<>% group_by(id) %>% mutate(SOFA_total_lag2 = lag(SOFA_total, 2, default = 0))


### estimate observable/factual cumulative incidences ----

# calculate total sample size
n <- length(unique(CPdata$id))

## cumulative incidence of ICU death (competing risk model in Fig 1A)

# Aalen-Johansen estimator via survfit function applied on original (wide) data format
mod1 <- survfit(Surv(time1, event1) ~ 1, data = events)
plot(mod1["death"], xlim = c(0, 30), fun = "event", xlab = "Days from ICU admission", ylab = "Cumulative risk of ICU death")
#ggfortify:::autoplot.survfit(mod1["death"], xlim = c(0, 30), xlab = "Days from ICU admission", ylab = "Cumulative risk of ICU death")
#prodlim::prodlim(prodlim::Hist(time1, event1) ~ 1, data = events) %>% plot(cause = "death", xlim = c(0, 30), ylim = c(0, 0.2))

# Aalen-Johansen estimator via survfit function applied on counting process data format
mod1CP <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata, id = id)
lines(mod1CP["death"], fun = "event", col = "red")

# store cumulative risk estimates based on AJ estimator
E_Y_AJ <- with(mod1CP["death"], data.frame(k = time, E_Y_AJ = pstate))

# compare with empirical cumulative distribution function
E_Y_ecdf <- CPdata %>% group_by(k) %>% summarise(tmp = sum(Y_k*(1-D_k)*(1-Y_k_lag))/n) %>% 
  mutate(E_Y_ecdf = cumsum(tmp)) %>% select(k, E_Y_ecdf)
E_Y <- full_join(E_Y_AJ, E_Y_ecdf, by = "k") %>% arrange(k)
E_Y


## cumulative incidence of HAI-free ICU death (competing risk model in Fig 1B)

# Aalen-Johansen estimator via survfit function applied on original (wide) data format
mod2 <- survfit(Surv(time2, event2) ~ 1, data = events)
plot(mod2["HAI-free death"], fun = "event", xlim = c(0, 30), xlab = "Days from ICU admission", ylab = "Cumulative risk of HAI-free ICU death")

# Aalen-Johansen estimator via survfit function applied on counting process data format
mod2CP <- survfit(Surv(tstart, tstop, event2) ~ 1, data = CPdata, id = id)
lines(mod2CP["HAI-free death"], fun = "event", col = "red")

# store cumulative risk estimates based on AJ estimator
E_tilde_Y_AJ <- with(mod2CP["HAI-free death"], data.frame(k = time, E_tilde_Y_AJ = pstate))

# compare with empirical cumulative distribution function
E_tilde_Y_ecdf <- CPdata %>% group_by(k) %>% summarise(tmp = sum(Y_k*(1-D_k)*(1-A_k)*(1-Y_k_lag))/n) %>% 
  mutate(E_tilde_Y_ecdf = cumsum(tmp)) %>% select(k, E_tilde_Y_ecdf)
E_tilde_Y <- full_join(E_tilde_Y_AJ, E_tilde_Y_ecdf, by = "k") %>% arrange(k)
E_tilde_Y



### Treating (absence of) exposure as a conditioning event or baseline exclusion criterion (see Appendix G1 in Supplementary Data S2) ----

## Aalen-Johansen estimator via survfit function applied on original (wide) data format
mod_star1 <- survfit(Surv(time1, event1) ~ 1, data = events %>% filter(A_E == 0))
mod_star2 <- survfit(Surv(time2, event2) ~ 1, data = events %>% filter(A_E == 0))
mod_star3 <- survfit(Surv(time1, event1) ~ A_E, data = events)
mod_star4 <- survfit(Surv(time2, event2) ~ A_E, data = events)
mod_star5 <- survfit(Surv(time1, event1) ~ 1, data = events, weights = 1-A_E)
mod_star6 <- survfit(Surv(time2, event2) ~ 1, data = events, weights = 1-A_E)
plot(mod_star1["death"], fun = "event", xlim = c(0, 30), main = "Treating HAI as\n a conditioning event", xlab = "Days from ICU admission", ylab = "Counterfactual cumulative risk of ICU death")
lines(mod_star2["HAI-free death"], fun = "event", col = "red")
lines(mod_star3["A_E=0", "death"], fun = "event", col = "green")
lines(mod_star4["A_E=0", "HAI-free death"], fun = "event", col = "purple")
lines(mod_star5["death"], fun = "event", col = "red")
lines(mod_star6["HAI-free death"], fun = "event", col = "green")
#   Note that all the above survfit commands provide the same result because
# 1/ when restricting the data to patients who don't develop an infection during hospitalization,
# for each k, we have Y_k = tilde_Y_k;
# 2/ restricting the data to patients who don't develop an infection is equivalent to stratification 
# and to setting weights of infected patients to zero.

## Aalen-Johansen estimator via survfit function applied on counting process data format
mod_starCP1 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata %>% filter(A_E == 0), id = id)
mod_starCP2 <- survfit(Surv(tstart, tstop, event1) ~ A_E, data = CPdata, id = id)
mod_starCP3 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata, id = id, weights = 1-A_E) 
lines(mod_starCP1["death"], fun = "event", col = "purple")
lines(mod_starCP2["A_E=0", "death"], fun = "event", col = "red")
lines(mod_starCP3["death"], fun = "event", col = "green")

## store cumulative risk estimates based on AJ estimator
E_Y0_star_AJ <- with(mod_starCP1["death"], data.frame(k = time, E_Y0_star_AJ = pstate))

## compare with weighted empirical cumulative distribution function
class(CPdata) <- "data.frame"
#   Note: recode CPdata as a data.frame object to avoid error in ipwpoint function below
ipw_star <- ipwpoint(exposure = A_E, family = "binomial", link = "logit", denominator = ~ 1, data = events)
events$W_star <- ipw_star$ipw.weights
CPdata %<>% left_join(events %>% select(id, W_star))
#   Note: first calculate naive weights (see lines above). 
# As this corresponds to (inappropriately) considering HAI onset as a point exposure at baseline, 
# we can simply use the ipwpoint function from the ipw package for this purpose, specifying an empty covariate set.
# However, in order to obtain the correct weights, this function needs to applied to the original (wide format) dataset, 
# and the resulting weights then need to be added to the counting process format dataset.
E_Y0_star_ecdf <- CPdata %>% group_by(k) %>% summarise(tmp = sum(Y_k*(1-D_k)*(1-A_k)*(1-Y_k_lag)*W_star)/n) %>% 
  mutate(E_Y0_star_ecdf = cumsum(tmp)) %>% select(k, E_Y0_star_ecdf)
E_Y0_star <- full_join(E_Y0_star_AJ, E_Y0_star_ecdf, by = "k") %>% arrange(k)
E_Y0_star

## estimate population attributable fraction
PAF_star <- full_join(E_Y_AJ, E_Y0_star_AJ, by = "k") %>% arrange(k)
PAF_star %<>% mutate(PAF_star = (E_Y_AJ - E_Y0_star_AJ)/E_Y_AJ)
PAF_star


### Treating (absence of) exposure onset as a conditioning event or time-dependent exclusion criterion (see Appendix G2 in Supplementary Data S2) ----

#   Note: because the time-dependent exclusion criterion (time-dependent weights) is (are) indexed by subsequent landmark intervals K at which the cumulative incidence is evaluated
# rather than at follow-up intervals k over which the (weighted) proportion of deaths are summed, 
# this time-dependent exclusion criterion cannot be directly applied in the survfit function, neither by exclusion, stratification nor weighting (as before)
# Instead, we apply the original estimator as proposed by Schumacher et al 2007 as a functional of the estimated cumulative incidence of ICU death and HAI onset, respectively,
# as detailed in Appendix G2 in the Supplementary Material.

Pprime01 <- E_tilde_Y_AJ %>% rename(Pprime01 = E_tilde_Y_AJ)
Pprime03 <- with(mod2CP["HAI onset"], data.frame(k = time, Pprime03 = pstate))

## store cumulative risk estimates based on AJ estimators
E_Y0_dagger_AJ <- full_join(Pprime01, Pprime03, by = "k")
E_Y0_dagger_AJ %<>% mutate(E_Y0_dagger_AJ = Pprime01/(1-Pprime03)) %>% select(k, E_Y0_dagger_AJ)

plot(E_Y0_dagger_AJ ~ k, data = E_Y0_dagger_AJ, type = "s", xlim = c(0, 30), main = "Treating HAI onset as\n a conditioning event", xlab = "Days from ICU admission", ylab = "Counterfactual cumulative risk of ICU death")

## compare with weighted empirical cumulative distribution function
#   Note: calculating the weights defined in section 5.2 in the main text (and in Appendix G2 in the Supplementary Material)
# requires to extend the counting process format data to the maximal follow-up interval E (i.e. beyond the observed event time) for each patient 
E <- 60
CPdata %<>% filter(k <= E)
CPdata_ext <- full_join(expand.grid(k = 1:E, id = unique(CPdata$id)), CPdata, by = c("id", "k"))
CPdata_ext %<>% fill(A_k, D_k, Y_k)
# recalculate lag variables
CPdata_ext %<>% group_by(id) %>% mutate(A_k_lag = lag(A_k, default = 0),
                                        D_k_lag = lag(D_k, default = 0),
                                        Y_k_lag = lag(Y_k, default = 0))

class(CPdata_ext) <- "data.frame"
#   Note: recode CPdata as a data.frame object to avoid error in ipwpoint function below
W_dagger <- data.frame(k = 1:60, W_dagger = NA)
for(K in 1:60) {
  ipw_dagger <- ipwpoint(exposure = A_k, family = "binomial", link = "logit", denominator = ~ 1, data = CPdata_ext %>% filter(k == K))
  W_dagger[K, "W_dagger"] <- ipw_dagger$ipw.weights[1] 
  # Note: make sure to pick the id/index of a subject that is uninfected during the entire hospitalization to extract the weight!
}
CPdata %<>% left_join(W_dagger, by = "k")
CPdata_ext %<>% left_join(W_dagger, by = "k")

#   Note: to calculate the weighted ecdf first sum over intervals k per subject, then multiply by weight indexed by K, and then average over subjects for each k 
# (rather than first averaging over subjects at each k and then taking cumulative sum as before)
E_Y0_dagger_ecdf <- CPdata_ext %>% arrange(id, k) %>% group_by(id) %>% mutate(tmp = cumsum(Y_k*(1-D_k)*(1-A_k)*(1-Y_k_lag))*W_dagger) %>% 
  group_by(k) %>% summarise(E_Y0_dagger_ecdf = sum(tmp)/n) %>% select(k, E_Y0_dagger_ecdf)
E_Y0_dagger <- full_join(E_Y0_dagger_AJ, E_Y0_dagger_ecdf, by = "k") %>% arrange(k)
E_Y0_dagger

## estimate population attributable fraction
PAF_dagger <- full_join(E_Y_AJ, E_Y0_dagger_AJ, by = "k") %>% arrange(k)
PAF_dagger %<>% mutate(PAF_dagger = (E_Y_AJ - E_Y0_dagger_AJ)/E_Y_AJ)
PAF_dagger


### Treating exposure onset as (marginally) independent censoring (see Appendix G3 in Supplementary Data S2) ----

## Aalen-Johansen estimator via survfit function applied on original (wide) data format
mod_circ <- survfit(Surv(timeC, eventC) ~ 1, data = events)
plot(mod_circ["death"], fun = "event", xlim = c(0, 30), main = "Treating HAI onset as\n (marginally) independent censoring", xlab = "Days from ICU admission", ylab = "Counterfactual cumulative risk of ICU death")

## Aalen-Johansen estimator via survfit function applied on counting process data format
mod_circCP1 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata %>% filter(A_k == 0), id = id)
mod_circCP2 <- survfit(Surv(tstart, tstop, event1) ~ A_k, data = CPdata, id = 1:nrow(CPdata))
mod_circCP3 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata, id = 1:nrow(CPdata), weights = 1-A_k) 
lines(mod_circCP1["death"], fun = "event", col = "red")
lines(mod_circCP2["A_k=0", "death"], fun = "event", col = "green")
lines(mod_circCP3["death"], fun = "event", col = "purple")
#   Note that, to guarantee identical results, when stratifying on the time-dependent exposure status A_k or applying time-dependent weights 1-A_k,
# the id argument needs to be specified differently, treating each row in the dataset as a 'pseudo-patient'.

## store cumulative risk estimates based on AJ estimator
E_Y0_circ_AJ <- with(mod_circCP1["death"], data.frame(k = time, E_Y0_circ_AJ = pstate))

## compare with weighted empirical cumulative distribution function
class(CPdata) <- "data.frame"
#   Note: recode CPdata as a data.frame object to avoid error in ipwtm function below
ipw_circ <- ipwtm(exposure = A_k, family = "survival",
                     denominator = ~ 1, id = id,
                     tstart = tstart, timevar = tstop,
                     type = "first", data = CPdata)
CPdata$W_circ <- ipw_circ$ipw.weights
#   Note: first calculate IPC weights based on time-dependent propensity score models conditional on empty covariate sets (see lines above)
E_Y0_circ_ecdf <- CPdata %>% group_by(k) %>% summarise(tmp = sum(Y_k*(1-D_k)*(1-A_k)*(1-Y_k_lag)*W_circ)/n) %>% 
  mutate(E_Y0_circ_ecdf = cumsum(tmp)) %>% select(k, E_Y0_circ_ecdf)
E_Y0_circ <- full_join(E_Y0_circ_AJ, E_Y0_circ_ecdf, by = "k") %>% arrange(k)
E_Y0_circ
#   Note: although estimators are algebraically equivalent, estimates are not identical because of the way weights are calculated using the ipwtm function.
# Identical estimates can be obtained by self-calculated IPC weights (e.g. based on a pooled logit model or the KM estimator). 

## estimate population attributable fraction
PAF_circ <- full_join(E_Y_AJ, E_Y0_circ_AJ, by = "k") %>% arrange(k)
PAF_circ %<>% mutate(PAF_circ = (E_Y_AJ - E_Y0_circ_AJ)/E_Y_AJ)
PAF_circ


### Treating exposure onset as conditionally independent censoring (see Appendix G4 in Supplementary Data S2) ----

## first calculate IPC weights that account for confounder history up to each time
class(CPdata) <- "data.frame"
#   Note: recode CPdata as a data.frame object to avoid error in ipwtm function below
ipw <- ipwtm(exposure = A_k, family = "survival",
             denominator = ~ sex + as.numeric(admissioncat != "Medicine") + admissionyr + ns(age, df = 4) + ns(CCI, df = 2) + ns(SOFA_total_0, df = 4) + ns(SOFA_total_lag2, df = 4), id = id,
             tstart = tstart, timevar = tstop,
             type = "first", data = CPdata)
CPdata$W <- ipw$ipw.weights

## IPC weighted Aalen-Johansen estimator cannot be obtained via survfit function applied on original (wide) data format, because this does not allow accounting for time-dependent confounding

## IPC weighted Aalen-Johansen estimator via survfit function applied on counting process data format
modCP1 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata %>% filter(A_k == 0), id = 1:nrow(CPdata %>% filter(A_k == 0)), weights = W)
modCP2 <- survfit(Surv(tstart, tstop, event1) ~ A_k, data = CPdata, id = 1:nrow(CPdata), weights = W)
modCP3 <- survfit(Surv(tstart, tstop, event1) ~ 1, data = CPdata, id = 1:nrow(CPdata), weights = (1-A_k)*W) 
plot(modCP1["death"], fun = "event", xlim = c(0, 30), main = "Treating HAI onset as\n informative censoring", xlab = "Days from ICU admission", ylab = "Counterfactual cumulative risk of ICU death")
lines(modCP2["A_k=0", "death"], fun = "event", col = "red")
lines(modCP3["death"], fun = "event", col = "green")
#   Note: to guarantee appropriate estimation as intended,
# the id argument needs to be specified such that each row in the dataset is treated as an independent 'pseudo-patient'.
# As a result of this (and the uncertainty in estimation of the weights), confidence intervals can no longer be expected to have nominal coverage.
# Instead, we recommend using the non-parametric bootstrap.

## store cumulative risk estimates based on AJ estimator
E_Y0_AJ <- with(modCP1["death"], data.frame(k = time, E_Y0_AJ = pstate))

## compare with weighted empirical cumulative distribution function
E_Y0_ecdf <- CPdata %>% group_by(k) %>% summarise(tmp = sum(Y_k*(1-D_k)*(1-A_k)*(1-Y_k_lag)*W)/n) %>% 
  mutate(E_Y0_ecdf = cumsum(tmp)) %>% select(k, E_Y0_ecdf)
E_Y0 <- full_join(E_Y0_AJ, E_Y0_ecdf, by = "k") %>% arrange(k)
E_Y0

## estimate population attributable fraction
PAF <- full_join(E_Y_AJ, E_Y0_AJ, by = "k") %>% arrange(k)
PAF %<>% mutate(PAF = (E_Y_AJ - E_Y0_AJ)/E_Y_AJ)


### compare estimated PAFs ----
PAF_star %>% select(k, PAF_star) %>%
  full_join(PAF_dagger %>% select(k, PAF_dagger), by = "k") %>%
  full_join(PAF_circ %>% select(k, PAF_circ), by = "k") %>%
  full_join(PAF %>% select(k, PAF), by = "k") %>% arrange(k)
