#' @title
#' blindness probability function
#' @description
#' calculate blindness probabilities for mean mf counts
#'
#' @param dat this is the data from Little et al. 2004 with blindness probabilities calculated for binned mean mf loads
#'
#' @returns vector of blindness probabilities for mf counts from 0 to 1000
blind_mfcount_prob_func <- function(dat){

  # ep_mfcount_prob_df <- data.frame(prob = c(0.0061, 0.0439, 0.0720, 0.0849, 0.1341, 0.1538, 0.20),
  #                                  mf = c(0, 3, 13, 36, 76, 151, 200)) # data from Chesnais et al. 2018


  fit <- glm(prob ~ log(mf + 1), family = gaussian(link = "log"),
             data = dat) # fitted logarithmic relationship between prob OAE ~ mean mf load;
                         # generalized linear model with log-link function

  newdat <- data.frame(mf = seq(0, 10000, 1)) # generate new dataframe with mf counts from 0 to 1000 by 1

  out <- predict(fit, newdata = newdat, se.fit = T) # use predict function to calculate log(OAE prob) for each mean mf count in newdat

  logpred <- data.frame(fit = out$fit, se = out$se.fit, mf = newdat$mf) # create dataframe with predicted log(blind probs), standard error and mf counts

  logpred <- within(logpred, {
    lwr <- fit - 1.96 * se
    upr <- fit + 1.96 * se
  }) # calculate lwr and upper confidence intervals for log(OAE probs)

  pred <- data.frame(fit = exp(logpred$fit), lwr = exp(logpred$lwr), upr = exp(logpred$upr), mf = logpred$mf) # exponentiate probabilities

  blind_probs <- pred$fit # extract probabilities into a vector

  return(blind_probs)

}




#' @title
#' find individuals which develop new blindness cases
#'
#' @description
#' find individuals with new worm infection and can then develop blindness in the next function
#' (determined by probability associated with individual mf load)
#'
#' @param dat this is the master matrix (all.mats.temp) tracking mf and worm compartments (by age) for each individual
#' @param worms.start starting column for adult worms in master matrix
#' @param tot.worms final column for adult worms in master matrix
#' @param infected_at_all vector of all individuals to identify new infections for blindness
#' @param blindness vector of individuals defined as able to develop blindness to determine whether blindness onset occurs
#' @param tested_blind vector to specify which individuals have been tested (and so not recounted as new incident cases)
#'
find_indiv_blindness_func <- function(dat, worms.start, tot.worms, infected_at_all, blindness, tested_blind, check_ind_blind) {

  worms_all <- rowSums(dat[, worms.start : tot.worms]) # sum all worm categories per host (across 21 age classes)

  ind_new_inf <- which(worms_all > 0 & infected_at_all == 0) # find individuals with at least 1 worm & no previous infection (therefore new)

  infected_at_all[ind_new_inf] <-  1 # where new worm infection, give value 1 at individual position in vector

  # ind_age <- which(round(dat[, 2]) == round(age_to_samp)) # not required for blindness as age of onset not defined (discrete like for OAE 3-18 yrs)
                                                          #  and Little et al. 2004 has binned mf counts for wide age range
                                                          # to obtain blindness probabilities?


  ind_ibf <- which(infected_at_all == 1) # extract position of currently with worm infection to inform blindness (not new infection to inform OAE)

  ind_no_blind <- which(blindness != 1) # extract position of individuals without OAE (so no current OAE?)

  ind_no_samp <- which(tested_blind == 0) # extract position of individuals not (previously?) tested -
                                              # to ensure dont double count individuals when generating new cases


  tot_ind_ep_samp_blind <- Reduce(intersect, list(ind_age, ind_ibf, ind_no_blind, ind_no_samp)) # find common (intersect) where same individual (position number in vector)
                                                                                                        # present in all vectors (any (prior/current) worm infection,
                                                                                                        # no current blindness, not previously tested
                                                                                                        # to generate new incident cases of blindness (current: worm +, blindness -)

  check_ind_blind <- c(check_ind_blind, length(tot_ind_ep_samp_blind)) # ?

  return(list(tot_ind_ep_samp_blind, infected_at_all, check_ind_blind))

}

#' @title
#' calculate blindness incidence and prevalence
#'
#' @description
#' each individual undergoes a Bernoulli trial (using blindness probabilities based on mf counts; new worm infection)
#' to ascertain new incident cases of blindness (indexed by those individuals with new infection?)
#'
#' @param temp.mf vector of all mf per skin snip for each individual
#' @param tot_ind_ep_samp_blind vector with positions of new cases in the population
#' @param blind_probs vector of blindness probabilities for mf counts from 0 to 1000
#' @param dat main matrix tracking mf age compartments and worm age compartments for each individual
#' @param prev_blind vector containing prevalence of blindness to update
#' @param blindness vector of individuals defined as able to develop to blindness to determine whether blindness onset occurs
#' @param tested_blind vector to specify which individuals have been tested
#' @param blind_incidence_DT blindness incidence vector (whole population) to update
#' @param blind_incidence_DT_M blindness incidence vector in males to update
#' @param blind_incidence_DT_F blindness incidence vector in females to update
#'
#' @returns list containing i) updated blindness prevalence, ii) updated blindness incidence in whole population, iii - iv) updated blindness
#' incidence in males and females, v) blindness probabilities for each individual based on mf count, #' vi) those tested/sampled
new_blind_cases_func <- function(temp.mf, tot_ind_ep_samp_blind, blind_probs, dat,
                               prev_blind, blindness, tested_blind,
                               blind_incidence_DT, blind_incidence_DT_M, blind_incidence_DT_F){


  mf_round <- round(temp.mf[[2]][tot_ind_ep_samp_blind]) + 1 # mf count (in those tested ?)
                                                                 # note: if zero the first probability comes from from ep_probs
                                                                 # note: +1 to account for (log) 0


  blind_rates <- blind_probs[mf_round] # get probabilities (based on individual mf count)

  blindness[tot_ind_ep_samp_blind] <- rbinom(length(tot_ind_ep_samp_blind), 1, blind_rates) # for individuals (those positions of individuals as new cases)
                                                                                                        # with an probability of blindness, they undergo a
                                                                                                        # Bernoulli trial to ascertain whether blindness onset realized

  tested_blind[tot_ind_ep_samp_blind] <- 1 # records that they've been tested (those identified in tot_ind_ep_samp_blindness are indexed)

  prev_blind <- c(prev_blind, mean(blind)) # calculate & update prevalence of blindness

  new_inc_blind <- length(which(blind[tot_ind_ep_samp_blind] == 1)) # how many infections (finds total new blind cases in all blind)

  blind_incidence_DT <- c(blind_incidence_DT, new_inc_blind) # record + update number of new OAE cases

  new_inc_blind_M <- length(which(OAE[tot_ind_ep_samp_blind] == 1 & dat[tot_ind_ep_samp_blind ,3] == 1)) # new cases in males
  new_inc_blind_F <- length(which(OAE[tot_ind_ep_samp_blind] == 1 & dat[tot_ind_ep_samp_blind ,3] == 0)) # new cases in females

  blind_incidence_DT_M <- c(blind_incidence_DT_M, new_inc_blind_M) # record & update incidence in males
  blind_incidence_DT_F <- c(blind_incidence_DT_F, new_inc_blind_F) # record & update incidence in females

  return(list(prev_blind, blind_incidence_DT, blind_incidence_DT_F, blind_incidence_DT_M,
              blind_rates, blindness, tested_blind))

}

