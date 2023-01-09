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
             data = ep_mfcount_prob_df) # fitted logarithmic relationship between prob OAE ~ mean mf load;
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
#' @param tot.worms final column for adult worms in master matrix
#' @param infected_at_all vector of all individuals to identify new infections for OAE disease
#' @param age_to_samp sample of ages for each individual (default between 3 - 10 yrs of age) to match to ages in main matrix
#' @param OAE vector of individuals defined as able to develop to OAE/ to determine whether OAE onset occurs
#' @param tested_OAE vector to specify which individuals have been tested
#'
#' @returns vector of OAE probabilities for mf counts from 0 to 1000
find_indiv_OAE_func <- function(dat, mf.start, mf.end, worms.start, tot.worms,
                                infected_at_all, age_to_samp, OAE, tested_OAE){

  mf_all <- rowSums(dat[, mf.start : mf.end]) # sum mf per individual across all 21 age classes

  worms_all <- rowSums(dat[, worms.start : tot.worms]) # sum all worm categories per host (across 21 age classes)

  ind_new_inf <- which(worms_all > 0 & infected_at_all == 0) # find individuals with at least 1 worm & no previous infection (therefore new)

  # print(ind_new_inf)

  infected_at_all[ind_new_inf] <-  1 # where new worm infection, give value 1 at individual position in vector

  ind_age <- which(round(all.mats.temp[, 2]) == round(age_to_samp)) # find individuals (position in individual vector) where (rounded)
                                                                    # age at sampling (e.g. between 3 - 10) matches (rounded) age of individual

  ind_sex <- all.mats.temp[, 3] # extract sex of each individual from main matrix

  # print(ind_age)

  ind_ibf <- which(infected_at_all == 1) # extract position of new worm infections (not ever infected)?

  ind_no_OAE <- which(OAE != 1) # extract position of individuals without OAE

  ind_no_samp <- which(tested_OAE == 0) # extract position of individuals not tested

  # print(age_to_samp)

  tot_ind_ep_samp <- Reduce(intersect, list(ind_age, ind_ibf, ind_no_OAE, ind_no_samp, ind_sex)) # find common (intersect) ??? where same
                                                                                                 # individual (position number in vector)
                                                                                                 # present in all vectors (age, OAE, no OAE, sampling sex?)
  check_ind <- c(check_ind, length(tot_ind_ep_samp)) # ?

  return(tot_ind_ep_samp)

}





##INFECTED NEEDS TO BE ADULT WORMS

#mf_round <- round(mf_all[tot_ind_ep_samp]) + 1

temp.mf <- mf.per.skin.snip(ss.wt = 2, num.ss = 2, slope.kmf = 0.0478, int.kMf = 0.313, data = all.mats.temp, nfw.start, fw.end,
                            mf.start, mf.end, pop.size = N)



new_OAE_cases_func <- function(temp.mf, tot_ind_ep_samp, OAE_probs){


  mf_round <- round(temp.mf[[2]][tot_ind_ep_samp]) + 1 # mf count (in those tested ?)
                                                       # note: if zero the first probability comes from from ep_probs
                                                       # note: +1 to account for (log) 0


  OAE_rates <- OAE_probs[mf_round] # get probabilities (based on individual mf count)

  OAE[tot_ind_ep_samp] <- rbinom(length(tot_ind_ep_samp), 1, OAE_rates) # for individuals with an probability of OAE, they undergo a
                                                                        # Bernoulli trial to ascertain whether OAE onset realised

#records that they've been tested
tested_OAE[tot_ind_ep_samp] <- 1

prev_OAE <- c(prev_OAE, mean(OAE))

new_inc <- length(which(OAE[tot_ind_ep_samp] == 1)) #how many infections

OAE_incidence_DT <- c(OAE_incidence_DT, new_inc)


new_inc_3_5 <- length(which(OAE[tot_ind_ep_samp] == 1 & all.mats.temp[tot_ind_ep_samp ,2] >= 3 & all.mats.temp[tot_ind_ep_samp ,2]< 5 ))
new_inc_5_10 <- length(which(OAE[tot_ind_ep_samp] == 1 & all.mats.temp[tot_ind_ep_samp ,2] >= 5 & all.mats.temp[tot_ind_ep_samp ,2]<= 10 ))

new_inc_M <- length(which(OAE[tot_ind_ep_samp] == 1 & all.mats.temp[tot_ind_ep_samp ,3] == 1))
new_inc_F <- length(which(OAE[tot_ind_ep_samp] == 1 & all.mats.temp[tot_ind_ep_samp ,3] == 0))



OAE_incidence_DT_3_5 <- c(OAE_incidence_DT_3_5, new_inc_3_5)
OAE_incidence_DT_5_10 <- c(OAE_incidence_DT_5_10, new_inc_5_10)

OAE_incidence_DT_M <- c(OAE_incidence_DT_M, new_inc_M)
OAE_incidence_DT_F <- c(OAE_incidence_DT_F, new_inc_F)

}
