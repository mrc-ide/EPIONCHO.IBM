
#' @title
#' find the daily probability from the yearly probability
#'
#' @description
#' Take a yearly probability and calculate the daily probability
#'
#' @param prob the yearly probability
#' @param days_per_year the number of days in a year (355 or 366)
#'
#' @returns the daily probability
calc_daily_prob <- function(prob, days) {
  daily_prob <- 1 - (1 - prob) ^ (1 / days)
  return(daily_prob)
}


#' @title
#' find individuals which develop morbidity
#'
#' @description
#' find individuals sampled between specific morb.mat.tmp$Age with current mf and can then develop morbidity in the next function
#' (determined by probability associated with individual mf load)
#'
#' @param dat this is the master matrix (all.mats.temp) tracking mf (by age) for each individual
#' @param mf.start starting column for mf in master matrix
#' @param mf.end final column for mf in master matrix
#' @param morb.mat.tmp dataframe containing columns to determining conditions for testing morbidity
#' @param age_to_samp_vec_reversible 1 year increments of increasing age for when an individual should be tested (matching age), between 5 - 80 yrs
#' @param age_to_samp_vec_nonreversible 1 year increments of increasing age for when an individual should be tested (matching age), between 20 - 80 yrs
#'
#' @returns updated matrix with identifying individuals to test for morbidity
find_indiv_totest_func <- function(dat, mf.start, mf.end, morb.mat.tmp, temp_mf, age_to_samp_vec_reversible,
                                   age_to_samp_vec_nonreversible){

  # update age (and sex for newborns)
  morb.mat.tmp$Age <- dat[,2]
  morb.mat.tmp$Sex <- dat[,3]

  # true number of mf per individual #
  mf_all <- rowSums(dat[, mf.start : mf.end]) # sum mf per individual across all 21 age classes
  morb.mat.tmp$TrueMFCount <- mf_all # true mf count

  # extract number of mf per skin snip per individual
  morb.mat.tmp$ObservedMFCount <- round(temp_mf[[2]])

  # ======================== #
  #  1 ) age sampling        #

  # ================================================================================================================== #
  # determine if age to samp matches for all conditions (only want to sample once per year of age, not every time-step) #
  # morb.mat.tmp[,14] <- ifelse((round(morb.mat.tmp$Age,6) %in% round(age_to_samp_vec_reversible,6)),1,0) # severe itch (age vec to 6 decimla places more most specific matching - only once per year)
  # morb.mat.tmp[,18] <- ifelse((round(morb.mat.tmp$Age,6) %in% round(age_to_samp_vec_reversible,6)),1,0) # RSD (test again single RSD age_to_sample e.g., col 10)

  # ==================================#
  # WHOM TO UNDERGO BERNOULI TRAIL #
  # # based on > 0 true mf
  # morb.mat.tmp[,30] <- ifelse(morb.mat.tmp$TrueMFCount > 0 & morb.mat.tmp$SevereItchStatus == 0 & morb.mat.tmp$Age >= 2, 1, 0) # for SI (if 0 disease state)
  # morb.mat.tmp[,34] <- ifelse(morb.mat.tmp$TrueMFCount > 0 & morb.mat.tmp$RSDStatus == 0 & morb.mat.tmp$Age >= 2, 1, 0) # for RSD (if 0 disease state)
  # morb.mat.tmp[,35] <- ifelse(morb.mat.tmp$TrueMFCount > 0 & morb.mat.tmp$AtrophySampleAges == 1 & morb.mat.tmp$AtrophyStatus == 0, 1, 0) # atrophy (irreversible; only test once in age range)
  # morb.mat.tmp[,36] <- ifelse(morb.mat.tmp$TrueMFCount > 0 & morb.mat.tmp$HangingGroinSampleAges == 1 & morb.mat.tmp$HGStatus == 0, 1, 0) # HG (irreversible; only test once in age range)
  # morb.mat.tmp[,37] <- ifelse(morb.mat.tmp$TrueMFCount > 0 & morb.mat.tmp$DepigSampleAge == 1 & morb.mat.tmp$DepigStatus == 0, 1, 0) # depigm (irreversible; only test once in age range)

  # based on > 0 observed mf
  morb.mat.tmp$ToTestSevereItch <- ifelse(morb.mat.tmp$ObservedMFCount > 0 & morb.mat.tmp$SevereItchStatus == 0 & morb.mat.tmp$Age >= 2, 1, 0) # for SI (if 0 disease state)
  morb.mat.tmp$ToTestRSD <- ifelse(morb.mat.tmp$ObservedMFCount > 0 & morb.mat.tmp$RSDStatus == 0 & morb.mat.tmp$Age >= 2, 1, 0) # for RSD (if 0 disease state)
  morb.mat.tmp$ToTestAtrophy <- ifelse(morb.mat.tmp$ObservedMFCount > 0 & morb.mat.tmp$AtrophyStatus == 0, 1, 0)
  morb.mat.tmp$ToTestHG <- ifelse(morb.mat.tmp$ObservedMFCount > 0 & morb.mat.tmp$HGStatus == 0, 1, 0)
  morb.mat.tmp$ToTestDepig <- ifelse(morb.mat.tmp$ObservedMFCount > 0 & morb.mat.tmp$DepigStatus == 0, 1, 0)


  return(morb.mat.tmp)

}

#' @title
#' individuals designated to be tested undergo Bernoulli trial for realise skin disease status
#'
#' @description
#' each individual undergoes a Bernoulli trial (using probabilities based on mf counts) to ascertain new cases of OAE
#'
#' @param morb.mat.tmp updated dataframe highlighting individuals to test
#' @param temp_mf vector of all mf per skin snip for each individual
#' @param SI_probs probabilities of severe itch for a given mf count
#' @param RSD_probs probabilities of RSD for a given mf count
#' @param Atrp_probs yearly probabilities of Atrophy for a given mf count
#' @param Hg_probs yearly probabilities of hanging groin for a given mf count
#' @param Depigm_probs yearly probabilities of depigmentation for a given mf count
#'
#' @returns updated matrix with disease status updated
new_cases_morbidity_func <- function(morb.mat.tmp, SI_probs, RSD_probs, Atrp_probs, Hg_probs, Depigm_probs){

  ## extract number of mf per skin snip per individual

  # # extract probabilities (rates) to run Bernoulli trial for each condition
  # morb.mat.tmp[,38] <- ifelse(morb.mat.tmp$ObservedMFCount > 0, SI_probs[morb.mat.tmp$ObservedMFCount], 0) # severe itch rates
  # #morb.mat.tmp[,39] <- ifelse(morb.mat.tmp$ObservedMFCount > 0, APOD_probs[morb.mat.tmp$ObservedMFCount], 0) # APOD rates
  # #morb.mat.tmp[,40] <- ifelse(morb.mat.tmp$ObservedMFCount > 0, CPOD_probs[morb.mat.tmp$ObservedMFCount], 0) # CPOD rates
  # #morb.mat.tmp[,41] <- ifelse(morb.mat.tmp$ObservedMFCount > 0, LOD_probs[morb.mat.tmp$ObservedMFCount], 0) # LOD rates
  # morb.mat.tmp[,42] <- ifelse(morb.mat.tmp$ObservedMFCount > 0, RSD_probs[morb.mat.tmp$ObservedMFCount], 0) # RSD rates
  # morb.mat.tmp[,43] <- ifelse(morb.mat.tmp$ObservedMFCount > 0, Atrp_probs[morb.mat.tmp$ObservedMFCount], 0) # atrophy rates
  # morb.mat.tmp[,44] <- ifelse(morb.mat.tmp$ObservedMFCount > 0, Hg_probs[morb.mat.tmp$ObservedMFCount], 0) # hanging groin rates
  # morb.mat.tmp[,45] <- ifelse(morb.mat.tmp$ObservedMFCount > 0, Depigm_probs[morb.mat.tmp$ObservedMFCount], 0) # depigmentation groin rates


  # ======================= #
  # Undergo Bernouli trial  #

  # based on whether true mf present
  morb.mat.tmp$SevereItchStatus <- ifelse(morb.mat.tmp$ToTestSevereItch == 1, rbinom(sum(morb.mat.tmp$ToTestSevereItch), 1, SI_probs), morb.mat.tmp$SevereItchStatus) # severe itch (stay as prior disease condition status if test does not take place)
  morb.mat.tmp$RSDStatus <- ifelse(morb.mat.tmp$ToTestRSD == 1, rbinom(sum(morb.mat.tmp$ToTestRSD), 1, RSD_probs), morb.mat.tmp$RSDStatus) # RSD (stay as prior disease condition status if test does not take place)
  morb.mat.tmp$AtrophyStatus <- ifelse(morb.mat.tmp$ToTestAtrophy == 1, rbinom(sum(morb.mat.tmp$ToTestAtrophy), 1, calc_daily_prob(Atrp_probs, 365)), morb.mat.tmp$AtrophyStatus) # atrophy (non-reversible: only testing 0's
  morb.mat.tmp$HGStatus <- ifelse(morb.mat.tmp$ToTestHG == 1, rbinom(sum(morb.mat.tmp$ToTestHG), 1,  calc_daily_prob(Hg_probs, 365)), morb.mat.tmp$HGStatus) # HG (non-reversible: only testing 0's
  morb.mat.tmp$DepigStatus <- ifelse(morb.mat.tmp$ToTestDepig == 1, rbinom(sum(morb.mat.tmp$ToTestDepig), 1,  calc_daily_prob(Depigm_probs, 365)), morb.mat.tmp$DepigStatus) # depigmentation (non-reversible: only testing 0's

  # # based on whether observed mf present
  # morb.mat.tmp$SevereItchStatus <- ifelse(morb.mat.tmp$ToTestSevereItch == 1, rbinom(sum(morb.mat.tmp$ToTestSevereItch), 1, morb.mat.tmp[,38]), morb.mat.tmp$SevereItchStatus) # severe itch (stay as prior disease condition status if test does not take place)
  # morb.mat.tmp$RSDStatus <- ifelse(morb.mat.tmp$ToTestRSD == 1, rbinom(sum(morb.mat.tmp$ToTestRSD), 1, morb.mat.tmp[,42]), morb.mat.tmp$RSDStatus) # RSD (stay as prior disease condition status if test does not take place)
  # morb.mat.tmp$AtrophyStatus <- ifelse(morb.mat.tmp$ToTestAtrophy == 1, rbinom(sum(morb.mat.tmp$ToTestAtrophy), 1, morb.mat.tmp[,43]), morb.mat.tmp$AtrophyStatus) # atrophy (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp$HGStatus <- ifelse(morb.mat.tmp$ToTestHG == 1, rbinom(sum(morb.mat.tmp$ToTestHG), 1, morb.mat.tmp[,44]), morb.mat.tmp$HGStatus) # HG (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp$DepigStatus <- ifelse(morb.mat.tmp$ToTestDepig == 1, rbinom(sum(morb.mat.tmp$ToTestDepig), 1, morb.mat.tmp[,45]), morb.mat.tmp$DepigStatus) # depigmentation (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)


  return(morb.mat.tmp)


}


#' @title
#' update reversible conditions
#'
#' @description
#' update 3-day delay matrix for tracking how long someone has been disease positive (3rd day return to 0) & updating morbidity matrix
#'
#' @param sequela.postive.mat1 3rd day delay matrix for tracking SI disease state
#' @param sequela.postive.mat2 3rd day delay matrix for tracking RSD disease state
#' @param inds.sequela.mat vector to move delay matrix by one day
#' @param morb.mat.tmp morbidity matrix to update (reversible conditions)
#'
#' @returns a) updated SI 3-day delay matrix b) updated RSD 3rd day delay matrix c) updated morbidity matrix
update_reversible_sequela_func <- function(sequela.postive.mat1, sequela.postive.mat2, inds.sequela.mat,
                                           morb.mat.tmp){

  # # =============== #
  # # for severe itch #

  # note: morb.mat.tmp[,22] = status of whether this is the 3rd day of SI positivity
  # note: morb.mat.tmp[,46] = realized SI state on current day (current iter) - updated to 0 if 4th day of SI positivity
  #  Extract current sequela state for reversible conditions
  sequela.postive.mat1[,inds.sequela.mat] <- sequela.postive.mat1[,(inds.sequela.mat-1)] # move sequela state along one col
  sequela.postive.mat1[,1] <- morb.mat.tmp$SevereItchStatus # update first col with current sequela state on that day

  # new steps : update current sequela state in morb.mat if day of morbidity to 0 #
  morb.mat.tmp$Day3SevereItchStatus <- sequela.postive.mat1[,4] # assign day 3 from sequela delay matrix (4th col of delay mat)
  morb.mat.tmp$SevereItchStatus <- ifelse(morb.mat.tmp$Day3SevereItchStatus == 1, 0, morb.mat.tmp$SevereItchStatus) # update current disease status

  # #  Extract current sequela state for reversible conditions
  # morb.mat.tmp[,22] <- sequela.postive.mat1[,3] # assign day 3 from sequelae delay matrix
  #
  # sequela.postive.mat1[,inds.sequela.mat] <- sequela.postive.mat1[,(inds.sequela.mat-1)] # move sequela state along one col
  #
  # sequela.postive.mat1[,1] <- morb.mat.tmp$SevereItchStatus # update first col with current sequela state on that day
  #
  # morb.mat.tmp[,46] <- ifelse(morb.mat.tmp$DepigSampleAge == 1, 0, morb.mat.tmp$SevereItchStatus) # update current disease status


  # # =============== #
  # #     for RSD     #

  # note: morb.mat.tmp[,23] = status of whether this is the 3rd day of RSD positivity
  # note: morb.mat.tmp[,50] = realized RSD state on current day (current iter) - updated to 0 if 3rd day of RSD positivity

  #  Extract current sequela state for reversible conditions
  sequela.postive.mat2[,inds.sequela.mat] <- sequela.postive.mat2[,(inds.sequela.mat-1)] # move sequela state along one col
  sequela.postive.mat2[,1] <- morb.mat.tmp$RSDStatus # update first col with current sequela state on that day

  # new steps : update current sequela state in morb.mat if 3th day of morbidity to 0 #
  morb.mat.tmp$Day3RSDStatus <- sequela.postive.mat2[,4] # assign day 3 from sequela delay matrix (4th col of delay mat)
  morb.mat.tmp$RSDStatus <- ifelse(morb.mat.tmp$Day3RSDStatus == 1, 0, morb.mat.tmp$RSDStatus) # update current disease status

  # #  Extract current sequela state for reversible conditions
  #
  # morb.mat.tmp$Day3RSDStatus <- sequela.postive.mat2[,3] # assign day 7 from sequela delay matrix
  #
  # sequela.postive.mat2[,inds.sequela.mat] <- sequela.postive.mat2[,(inds.sequela.mat-1)] # move sequela state along one col
  #
  # sequela.postive.mat2[,1] <- morb.mat.tmp$RSDStatus # update first col with current sequela state on that day
  #
  # morb.mat.tmp$RSDStatus <- ifelse(day3RSDStatus == 1, 0, morb.mat.tmp$RSDStatus) # update current disease status
  return(list(sequela.postive.mat1, sequela.postive.mat2, morb.mat.tmp))

}

#' @title
#' disease prevalence calculator
#'
#' @description
#' this is a helper function that takes in a dataframe, age range and the morbidity name, and outputs the prevalence
#'
#' @param data dataframe containing all individuals for a given morbidity name, along with their age
#' @param morbidity_name morbidity name
#' @param age_min minimum age (inclusive)
#' @param age_max maximum age (exclusive)
#'
#' @returns output values for >5 and age-groups for all morbidities, in the order they are passed into the morbidities variable (i.e SI 5-80 + age groups, RSD 5=80 + age groups, etc.)

morbidity_prev_calculator <- function(data, morbidity_name, age_min, age_max = 81) {
  return(
    length(which(data[[morbidity_name]] == 1 & data$Age >= age_min & data$Age < age_max)) /
      length(which(data$Age >= age_min & data$Age < age_max))
  )
}

#' @title
#' disease prevalences
#'
#' @description
#' calculates disease prevalence for each disease state (including age-stratified disease prevalence's)
#'
#' @param non_blindness_morb_dat updated dataframe highlighting individuals to test for morbidities that are not related to blindness
#' @param blidness_dat updated dataframe highlighting individuals to test for morbidities related to blindness
#' @param N human pop size
#'
#' @returns output values for >5 and age-groups for all morbidities, in the order they are passed into the morbidities variable (i.e SI 5-80 + age groups, RSD 5=80 + age groups, etc.)
morbidity_prev_func <- function(
  non_blindness_morb_dat, blidness_dat, oae_status, N, age_groups,
  morbidities = c("SevereItchStatus", "RSDStatus", "AtrophyStatus", "HGStatus", "DepigStatus", "BlindnessStatus", "VisualImpairment", "OAEStatus")
) {
  if (nrow(non_blindness_morb_dat) == length(oae_status)) {
    non_blindness_morb_dat[, "OAEStatus"] <- oae_status
  } else {
    warning(
      "OAE Status population is more than morbidity population. Unable to calculate prevalence."
    )
  }
  # calculate current time-step skin disease state prevalence #
  # 2nd attempt : > 5 yr olds (to match Murdoch et al. 2017 prevs only sampling > 5 yr olds)
  output_vals <- rep(NA, length(morbidities) * length(age_groups))

  for (morb_index in 1:length(morbidities)) {
    morbidity_to_use <- morbidities[morb_index]
    for (age_group_index in 1:length(age_groups)) {
      if (morbidity_to_use == "VisualImpairment" || morbidity_to_use == "BlindnessStatus") {
        multiplier <- 1
        if (morbidity_to_use == "VisualImpairment") {
          multiplier <- 1.78
        }
        output_vals[
          (morb_index - 1) * length(age_groups) + age_group_index
        ] <- morbidity_prev_calculator(
          blidness_dat, "BlindnessStatus",
          age_min = age_groups[[age_group_index]][1],
          age_max = age_groups[[age_group_index]][2]
        ) * multiplier
      } else {
        output_vals[
          (morb_index - 1) * length(age_groups) + age_group_index
        ] <- morbidity_prev_calculator(
          non_blindness_morb_dat, morbidity_to_use,
          age_min = age_groups[[age_group_index]][1],
          age_max = age_groups[[age_group_index]][2]
        )
      }
    }
  }
  return(output_vals)
}


# ================================================================================================================= #
#                                Eye disease functions                                                              #
# ================================================================================================================= #

#' @title
#' find individuals which develop morbidity
#'
#' @description
#' find individuals sampled between specific morb.mat.tmp$Age with current mf and can then develop morbidity in the next function
#' (determined by probability associated with individual mf load or fixed probability)
#'
#' @param dat this is the master matrix (all.mats.temp) tracking mf (by age) for each individual
#' @param mf.start starting column for mf in master matrix
#' @param mf.end final column for mf in master matrix
#' @param temp_mf vector of all mf per skin snip for each individual
#' @param morb.mat.tmp dataframe containing columns to determining conditions for testing morbidity (eye disease)
#' @param age_to_samp_vec_nonreversible 1 year increments of increasing age for when an individual should be tested (matching age), between 20 - 80 yrs
#'
#' @returns updated matrix with identifying individuals to test for morbidity
find_indiv_totest_func2 <- function(dat, mf.start, mf.end, morb.mat.tmp, age_to_samp_vec_nonreversible){

  # update age (and sex for newborns)
  morb.mat.tmp$Age <- dat[,2]
  morb.mat.tmp$Sex <- dat[,3]

  # Decrement the blindness countdown for those who have blindness pending, and update the blindness status for those who are now blind
  morb.mat.tmp$BlindnessCountdown <- ifelse(morb.mat.tmp$BlindnessPending == 1 & morb.mat.tmp$BlindnessCountdown > 0, morb.mat.tmp$BlindnessCountdown - 1, morb.mat.tmp$BlindnessCountdown)
  morb.mat.tmp$BlindnessStatus <- ifelse(morb.mat.tmp$BlindnessPending == 1 & morb.mat.tmp$BlindnessCountdown == 0, 1, morb.mat.tmp$BlindnessStatus)

  # true number of mf per individual #
  mf_all <- rowSums(dat[, mf.start : mf.end]) # sum mf per individual across all 21 age classes
  morb.mat.tmp$TrueMFCount <- mf_all # true mf count

  # find those individuals to test
  morb.mat.tmp$ToTestBlind <- ifelse(morb.mat.tmp$TrueMFCount > 0 & morb.mat.tmp$BlindnessPending == 0 & morb.mat.tmp$BlindnessStatus == 0, 1, 0)


  # # extract number of mf per skin snip per individual
  # morb.mat.tmp[,5] <- round(temp.mf[[2]]) # mf per skin snip for all individuals

 return(morb.mat.tmp)

}



#' @title
#' individuals designated to be tested undergo Bernoulli trial for realize eye disease status
#'
#' @description
#' each individual undergoes a Bernoulli trial (using probabilities based on mf counts) to ascertain new cases of eye disease
#'
#' @param morb.mat.tmp updated dataframe highlighting individuals to test
#' @param temp_mf vector of all mf per skin snip for each individual
#' @param blind_probs yearly probabilities of blindness for a given mf count (based on equation from Little et al. 2004)
#'
#' @returns updated matrix with disease status updated
new_cases_morbidity_func2 <- function(morb.mat.tmp, temp.mf, blind.probs){

  # extract number of mf per skin snip per individual
  morb.mat.tmp$ObservedMFCount <- round(temp.mf[[2]]) + 1 # mf per skin snip for all individuals (+1 because of indexing so that when indexing probabilities goes from 1)

  # extract probabilities (rates) to run Bernoulli trial for each condition
  morb.mat.tmp$BlindnessProb <- ifelse(
    morb.mat.tmp$ObservedMFCount > 0,
    ifelse(morb.mat.tmp$ObservedMFCount <= length(blind.probs), blind.probs[morb.mat.tmp$ObservedMFCount], blind.probs[length(blind.probs)]),
    0
  ) # if observed mf count > 0 then assign blindness probability mapped to observed mf count (if mf count exceeds mapping assign last/highest prob value)

  # ======================= #
  # Undergo Bernouli trial  #

  # Test all users who are not blind or don't have blindness pending
  blind_probs <- calc_daily_prob(morb.mat.tmp$BlindnessProb, 365) # convert to a vector of daily probs
  distribution <- rbinom(length(blind_probs), 1, blind_probs) # create vector of 1 or 0s length of population with daily blindness probs
  morb.mat.tmp$BlindnessPending <- ifelse(morb.mat.tmp$BlindnessPending == 0 & morb.mat.tmp$ToTestBlind == 1, distribution, morb.mat.tmp$BlindnessPending)
  stopifnot(length(which(is.na(distribution))) == 0)
  return(morb.mat.tmp)



}
