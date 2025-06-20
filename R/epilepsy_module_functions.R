#' @title
#' OAE probability function
#' @description
#' calculate onchocerciasis-associated epilepsy (OAE) probabilities for mean mf counts
#'
#' @param dat this is the data from Chesnais et al. 2018 with OAE probabilities for 7 mean mf loads in Cameroon
#'
#' @returns vector of OAE probabilities for mf counts from 0 to 1000
OAE_mfcount_prob_func <- function(dat){

  # ep_mfcount_prob_df <- data.frame(prob = c(0.0061, 0.0439, 0.0720, 0.0849, 0.1341, 0.1538, 0.20),
  #                                  mf = c(0, 3, 13, 36, 76, 151, 200)) # data from Chesnais et al. 2018


  fit <- glm(prob ~ log(mf + 1), family = gaussian(link = "log"),
             data = dat) # fitted logarithmic relationship between prob OAE ~ mean mf load;
                                         # generalized linear model with log-link function

  newdat <- data.frame(mf = seq(0, 10000, 1)) # generate new dataframe with mf counts from 0 to 1000 by 1

  out <- predict(fit, newdata = newdat, se.fit = T) # use predict function to calculate log(OAE prob) for each mean mf count in newdat

  logpred <- data.frame(fit = out$fit, se = out$se.fit, mf = newdat$mf) # create dataframe with predicted log(OAE probs), standard error and mf counts

  logpred <- within(logpred, {
    lwr <- fit - 1.96 * se
    upr <- fit + 1.96 * se
  }) # calculate lwr and upper confidence intervals for log(OAE probs)

  pred <- data.frame(fit = exp(logpred$fit), lwr = exp(logpred$lwr), upr = exp(logpred$upr), mf = logpred$mf) # exponentiate probabilities

  OAE_probs <- pred$fit # extract probabilities into a vector

  return(OAE_probs)

}


#' @title
#' find individuals which develop OAE
#'
#' @description
#' find individuals sampled between specific ages (default is 3 - 10 yrs) with new worm infection
#' and can then develop onchocerciasis-associated epilepsy (OAE) in the next function
#' (determined by probability associated with individual mf load)
#'
#' @param dat this is the master matrix (all.mats.temp) tracking mf and worm compartments (by age) for each individual
#' @param mf.start starting column for mf in master matrix
#' @param mf.end final column for mf in master matrix
#' @param worms.start starting column for adult worms in master matrix
#' @param nfw.start # start of infertile worms
#' @param fw.end # end of fertile worms
#' @param infected_at_all vector of all individuals to identify new infections for OAE disease
#' @param age_to_samp sample of ages for each individual (default between 3 - 10 yrs of age) to match to ages in main matrix
#' @param OAE vector of individuals defined as able to develop to OAE/ to determine whether OAE onset occurs
#' @param tested_OAE vector to specify which individuals have been tested
#'
#' @returns vector with positions of new cases in the population
find_indiv_OAE_func <- function(dat, mf.start, mf.end, worms.start, nfw.start, fw.end,
                                infected_at_all, age_to_samp, OAE, tested_OAE, check_ind){

  mf_all <- rowSums(dat[, mf.start : mf.end]) # sum mf per individual across all 21 age classes

  # worms_all <- rowSums(dat[, worms.start : tot.worms]) # OLD: sum all worm categories per host (across 21 age classes)

  male_worms_all <- rowSums(dat[, worms.start : nfw.start - 1]) # sum all male worms per host (across 21 age classes)
  female_worms_all <- rowSums(dat[, nfw.start : fw.end]) # sum all female worms per host (across 21 age classes)

  # ind_new_inf <- which(worms_all > 0 & infected_at_all == 0) # find individuals with at least 1 worm & no previous infection (therefore new)

  ind_new_inf <- which(male_worms_all > 0 & female_worms_all > 0 & infected_at_all == 0) # find individuals where both male and female worms present
                                                                                         # (i.e. pair of mating worms present in host)

  # print(ind_new_inf)

  infected_at_all[ind_new_inf] <-  1 # where new worm infection, give value 1 at individual position in vector

  ind_age <- which(round(dat[, 2]) == round(age_to_samp)) # find individuals (position in individual vector) where (rounded)
                                                          # age at sampling (e.g. between 3 - 10) matches (rounded) age of individual

  # print(ind_age)

  ind_ibf <- which(infected_at_all == 1) # extract position (vector) of all those infected with a pair of mating worms

  ind_no_OAE <- which(OAE != 1) # extract position of individuals not with OAE (so the same people aren't counted twice)

  ind_no_samp <- which(tested_OAE == 0) # extract position of individuals not previously tested

  # print(age_to_samp)

  tot_ind_ep_samp <- Reduce(intersect, list(ind_age, ind_ibf, ind_no_OAE, ind_no_samp)) # goes through each of the individual lists
                                                                                        # (individual criteria for whether a person)
                                                                                        # If individual has satisfied all those (intersect finds commonalities in columns)
                                                                                        # form pool of individuals which are then sampled for new OAE onset
  check_ind <- c(check_ind, length(tot_ind_ep_samp)) #

  return(list(tot_ind_ep_samp, infected_at_all, check_ind))

}

#' @title
#' calculate OAE incidence and prevalence
#'
#' @description
#' each individual undergoes a Bernoulli trial (using OAE probabilities based on mf counts; new worm infection)
#' to ascertain new incident cases of OAE (indexed by those individuals with new infection?)
#'
#' @param temp.mf vector of all mf per skin snip for each individual
#' @param tot_ind_ep_samp vector with positions of new cases in the population
#' @param OAE_probs vector of OAE probabilities for mf counts from 0 to 1000
#' @param dat main matrix tracking age, sex, mf age compartments and worm age compartments for each individual
#' @param prev_OAE vector containing prevalence of OAE to update
#' @param OAE vector of individuals defined as able to develop to OAE/ to determine whether OAE onset occurs
#' @param tested_OAE vector to specify which individuals have been tested
#' @param OAE_incidence_DT OAE incidence vector (whole population) to update
#' @param OAE_incidence_DT_under_5 OAE incidence vector in under 5 years (3 - 5 years) to update
#' @param OAE_incidence_DT_5_10 OAE incidence vector in 5 - 10 years to update
#' @param OAE_incidence_DT_11_15 OAE incidence vector in 10 - 15 years to update
#' @param OAE_incidence_DT_M OAE incidence vector in males to update
#' @param OAE_incidence_DT_F OAE incidence vector in females to update
#'
#' @returns list containing i) updated OAE prevalence, ii) updated OAE incidence in whole population, iii - vi) updated OAE
#' incidence in under 5 years, 5 - 10 years, 10 - 15 years, and males and females, vii) OAE probabilities for each individual based on mf count,
#' viii) those tested/sampled
new_OAE_cases_func <- function(temp.mf, tot_ind_ep_samp, OAE_probs, dat,
                               OAE, tested_OAE, age_groups){

  mf_round <- round(temp.mf[[2]][tot_ind_ep_samp]) + 1 # mf count (in those tested)
                                                       # note: if zero the first probability comes from from ep_probs
                                                       # note: +1 to account for (log) 0


  OAE_rates <- OAE_probs[mf_round] # get probabilities (based on individual mf count)

  OAE[tot_ind_ep_samp] <- rbinom(length(tot_ind_ep_samp), 1, OAE_rates) # for individuals (those positions of individuals as new cases)
                                                                        # with an probability of OAE, they undergo a
                                                                        # Bernoulli trial to ascertain whether OAE onset realized

  tested_OAE[tot_ind_ep_samp] <- 1 # records that they've been tested (those identified in tot_ind_ep_samp are indexed)

  new_inc <- length(which(OAE[tot_ind_ep_samp] == 1)) # how many infections (finds total new infected/OAE in all OAE)

  new_inc_under_5 <- length(which(OAE[tot_ind_ep_samp] == 1 & dat[tot_ind_ep_samp ,2] >= 3 & dat[tot_ind_ep_samp ,2]< 5 )) # new cases in under 5 age group
  new_inc_5_10 <- length(which(OAE[tot_ind_ep_samp] == 1 & dat[tot_ind_ep_samp ,2] >= 5 & dat[tot_ind_ep_samp ,2]<= 10 )) # new cases in 5 to 10 age group
  new_inc_11_15 <- length(which(OAE[tot_ind_ep_samp] == 1 & dat[tot_ind_ep_samp ,2] >= 10 & dat[tot_ind_ep_samp ,2]<= 15 )) # new cases in 10 to 15 age group

  new_inc_M <- length(which(OAE[tot_ind_ep_samp] == 1 & dat[tot_ind_ep_samp ,3] == 1)) # new cases in males
  new_inc_F <- length(which(OAE[tot_ind_ep_samp] == 1 & dat[tot_ind_ep_samp ,3] == 0)) # new cases in females

  return(list(
    "output_vals" = c(
      mean(OAE), new_inc, new_inc_under_5, new_inc_5_10, new_inc_11_15, new_inc_M, new_inc_F
    ),
    OAE_rates, OAE, tested_OAE
  ))
}
