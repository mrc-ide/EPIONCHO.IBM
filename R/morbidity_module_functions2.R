
#' @title
#' find individuals which develop morbidity
#'
#' @description
#' find individuals sampled between specific ages with current mf and can then develop morbidity in the next function
#' (determined by probability associated with individual mf load)
#'
#' @param dat this is the master matrix (all.mats.temp) tracking mf (by age) for each individual
#' @param mf.start starting column for mf in master matrix
#' @param mf.end final column for mf in master matrix
#' @param morb.mat.temp matrix containing columns to determining conditions for testing morbidity
#' @param age_to_samp_vec_reversible 1 year increments of increasing age for when an individual should be tested (matching age), between 5 - 80 yrs
#' @param age_to_samp_vec_nonreversible 1 year increments of increasing age for when an individual should be tested (matching age), between 20 - 80 yrs
#'
#' @returns updated matrix with identifying individuals to test for morbidity
find_indiv_totest_func2 <- function(dat, mf.start, mf.end, morb.mat.tmp, age_to_samp_vec_reversible,
                                   age_to_samp_vec_nonreversible){

  # update age (and sex for newborns)
  morb.mat.tmp[,1] <- dat[,2]
  morb.mat.tmp[,2] <- dat[,3]

  # true number of mf per individual #
  mf_all <- rowSums(dat[, mf.start : mf.end]) # sum mf per individual across all 21 age classes
  morb.mat.tmp[,4] <- mf_all # true mf count

  # ======================== #
  #  1 ) age sampling        #

  # ======================================================#
  # # current (only used for non-reversible to sample once) #
  # # determine whether age_to_samp matches age of individual (rounded) for each condition where only test once in age range
  # # morb.mat.tmp[,14] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,6]), 1, 0) # severe itch
  # # morb.mat.tmp[,15] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,10]), 1, 0) # APOD (test again single RSD age_to_sample e.g., col 10)
  # # morb.mat.tmp[,16] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,10]), 1, 0) # CPOD (test again single RSD age_to_sample e.g., col 10)
  # # morb.mat.tmp[,17] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,10]), 1, 0) # LOD (test again single RSD age_to_sample e.g., col 10)
  # # morb.mat.tmp[,18] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,10]), 1, 0) # RSD (test again single RSD age_to_sample e.g., col 10)
  # morb.mat.tmp[,19] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,11]), 1, 0) # atrophy
  # morb.mat.tmp[,20] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,12]), 1, 0) # hanging groin
  # morb.mat.tmp[,21] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,13]), 1, 0) # depigmentation

  # ======================================================== #
  # # new approach: determine if age to samp matches for all conditions (only want to sample once per year of age, not every time-step)
  # morb.mat.tmp[,14] <- ifelse((round(morb.mat.tmp[,1],6) %in% round(age_to_samp_vec_reversible,6)),1,0) # severe itch (age vec to 6 decimla places more most specific matching - only once per year)
  # morb.mat.tmp[,18] <- ifelse((round(morb.mat.tmp[,1],6) %in% round(age_to_samp_vec_reversible,6)),1,0) # RSD (test again single RSD age_to_sample e.g., col 10)
  morb.mat.tmp[,19] <- ifelse((round(morb.mat.tmp[,1],6) %in% round(age_to_samp_vec_nonreversible,6)),1,0) # atrophy
  morb.mat.tmp[,20] <- ifelse((round(morb.mat.tmp[,1],6) %in% round(age_to_samp_vec_nonreversible,6)),1,0) # hanging groin
  morb.mat.tmp[,21] <- ifelse((round(morb.mat.tmp[,1],6) %in% round(age_to_samp_vec_nonreversible,6)),1,0) # depigmentation

  # dont round?
  # morb.mat.tmp[,19] <- ifelse((round(morb.mat.tmp[,1],6) %in% round(age_to_samp_vec_nonreversible,6)),1,0) # atrophy
  # morb.mat.tmp[,20] <- ifelse((round(morb.mat.tmp[,1],6) %in% round(age_to_samp_vec_nonreversible,6)),1,0) # hanging groin
  # morb.mat.tmp[,21] <- ifelse((round(morb.mat.tmp[,1],6) %in% round(age_to_samp_vec_nonreversible,6)),1,0) # depigmentation

  # ==================== #
  # 2) Previously tested #

  # ======= #
  # current #
  # determine whether individual has previously tested for condition (for selecting those where only want to test once for a condition)
  # morb.mat.tmp[,22] <- ifelse(morb.mat.tmp[,30] == 1, 1, morb.mat.tmp[,22]) # severe itch
  # morb.mat.tmp[,23] <- ifelse(morb.mat.tmp[,31] == 1, 1, morb.mat.tmp[,23]) # APOD
  # morb.mat.tmp[,24] <- ifelse(morb.mat.tmp[,32] == 1, 1, morb.mat.tmp[,24]) # CPOD
  # morb.mat.tmp[,25] <- ifelse(morb.mat.tmp[,33] == 1, 1, morb.mat.tmp[,25]) # LOD
  # morb.mat.tmp[,26] <- ifelse(morb.mat.tmp[,34] == 1, 1, morb.mat.tmp[,26]) # RSD
  # morb.mat.tmp[,27] <- ifelse(morb.mat.tmp[,35] == 1, 1, morb.mat.tmp[,27]) # atrophy
  # morb.mat.tmp[,28] <- ifelse(morb.mat.tmp[,36] == 1, 1, morb.mat.tmp[,28]) # hanging groin
  # morb.mat.tmp[,29] <- ifelse(morb.mat.tmp[,37] == 1, 1, morb.mat.tmp[,29]) # depigmentation

  # ==================================#
  # 3) WHOM TO UNDERGO BERNOULI TRAIL #

  # determine whether individual will undergo Bernoulli trial in this time-step for condition (first test)
  # morb.mat.tmp[,30] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # severe itch (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,31] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # APOD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,32] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # CPOD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,33] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # LOD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,34] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # RSD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,19] == 1 & morb.mat.tmp[,27] == 0, 1, 0) # atrophy (irreversible; only test once in age range)
  # morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,20] == 1 & morb.mat.tmp[,28] == 0, 1, 0) # HG (irreversible; only test once in age range)
  # morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,21] == 1 & morb.mat.tmp[,29] == 0, 1, 0) # depigm (irreversible; only test once in age range)
  # morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,19] == 1, 1, 0) # atrophy (irreversible; only test once in age range)
  # morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,20] == 1, 1, 0) # HG (irreversible; only test once in age range)
  # morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,21] == 1, 1, 0) # depigm (irreversible; only test once in age range)


  # determine whether individual will undergo Bernoulli trial in this time-step for condition
  # update: now non-reversible conditions tested daily (only those that do not have the condition)
  # age-truncated for reversible conditions (from 20 yrs) due to age- disease prev profiles (Murdoch et al. 2017)
  # morb.mat.tmp[,30] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # severe itch (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,31] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # APOD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,32] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # CPOD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,33] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # LOD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,34] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # RSD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,51] == 0, 1, 0) # atrophy (irreversible; only test once in age range)
  # morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,52] == 0, 1, 0) # HG (irreversible; only test once in age range)
  # morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,53] == 0, 1, 0) # depigm (irreversible; only test once in age range)
  # morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,51] == 0, 1, 0) # atrophy (irreversible; only test once in age range)
  # morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,52] == 0, 1, 0) # HG (irreversible; only test once in age range)
  # morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,53] == 0, 1, 0) # depigm (irreversible; only test once in age range)

  # ========================= #
  # Version 1                #

  # # determine whether individual will undergo Bernoulli trial in this time-step for condition (first test)
  # # test a) reversible conditions daily and b) non-reversible only once in life-time (at age_to_sample & only not previously tested)
  # morb.mat.tmp[,30] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # severe itch (reversible = tested each dt regardless of presence of condition)
  # # morb.mat.tmp[,31] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # APOD (reversible = tested each dt regardless of presence of condition)
  # # morb.mat.tmp[,32] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # CPOD (reversible = tested each dt regardless of presence of condition)
  # # morb.mat.tmp[,33] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # LOD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,34] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # RSD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,19] == 1 & morb.mat.tmp[,27] == 0, 1, 0) # atrophy (irreversible; only test once in age range)
  # morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,20] == 1 & morb.mat.tmp[,28] == 0, 1, 0) # HG (irreversible; only test once in age range)
  # morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,21] == 1 & morb.mat.tmp[,29] == 0, 1, 0) # depigm (irreversible; only test once in age range)


  # # now both reversible 'tested' at each time-step and non-reversible (only if non-diseases) tested once per year of age
  # morb.mat.tmp[,30] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # severe itch (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,31] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # APOD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,32] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # CPOD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,33] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # LOD (reversible = tested each dt regardless of presence of condition)
  # morb.mat.tmp[,34] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # RSD (reversible = tested each dt regardless of presence of condition)

  morb.mat.tmp[,30] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,46] == 0 & morb.mat.tmp[,1] >= 2, 1, 0) # for SI (if 0 disease state)
  morb.mat.tmp[,34] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,50] == 0 & morb.mat.tmp[,1] >= 2, 1, 0) # for RSD (if 0 disease state)
  morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,19] == 1 & morb.mat.tmp[,51] == 0, 1, 0) # atrophy (irreversible; only test once in age range)
  morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,20] == 1 & morb.mat.tmp[,52] == 0, 1, 0) # HG (irreversible; only test once in age range)
  morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,21] == 1 & morb.mat.tmp[,53] == 0, 1, 0) # depigm (irreversible; only test once in age range)

  # ============ #
  # new approach #

  # # first determine if true mf count for an individual with a reversible condition is 0, then then become sequelea negative
  # morb.mat.tmp[,46] <- ifelse(morb.mat.tmp[,4] == 0, 0, morb.mat.tmp[,46]) # if true mf count is 0, then sequela is 0, else leave
  # morb.mat.tmp[,50] <- ifelse(morb.mat.tmp[,4] == 0, 0, morb.mat.tmp[,50]) # if true mf count is 0, then sequela is 0, else leave

  # # determine whether individual will undergo Bernoulli trial in this time-step for condition
  # # reversible conditions: only tested if sequela 0 once per year of age (cols 14 & 18 if age_to_samp_vec matches age)
  # # update: now non-reversible conditions (only those that do not have the condition) only tested once per year of age (cols 19 - 21)
  # # age-truncated for reversible conditions (from 20 yrs) due to age- disease prev profiles (Murdoch et al. 2017)
  # morb.mat.tmp[,30] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,46] == 0 & morb.mat.tmp[,14] == 1, 1, 0) # severe itch (reversible)
  # morb.mat.tmp[,34] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,50] == 0 & morb.mat.tmp[,18] == 1, 1, 0) # RSD (reversible)
  # morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,19] == 1 & morb.mat.tmp[,51] == 0, 1, 0) # atrophy (irreversible)
  # morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,20] == 1 & morb.mat.tmp[,52] == 0, 1, 0) # HG (irreversible)
  # morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,21] == 1 & morb.mat.tmp[,53] == 0, 1, 0) # depigm (irreversible)


  return(morb.mat.tmp)

}

#' @title
#' individuals designated to be tested undergo Bernoulli trial for realise skin disease status
#'
#' @description
#' each individual undergoes a Bernoulli trial (using probabilities based on mf counts) to ascertain new cases of OAE
#'
#' @param morb.mat.tmp updated matrix highlighting individuals to test
#' @param temp_mf vector of all mf per skin snip for each individual
#' @param SI_probs probabilities of severe itch for a given mf count
#' @param RSD_probs probabilities of RSD for a given mf count
#' @param Atrp_probs probabilities of Atrophy for a given mf count
#' @param Hg_probs probabilities of hanging groin for a given mf count
#' @param Depigm_probs probabilities of depigmentation for a given mf count
#'
#' @returns updated matrix with disease status updated
new_cases_morbidity_func2 <- function(morb.mat.tmp, temp_mf, SI_probs, RSD_probs, Atrp_probs, Hg_probs, Depigm_probs){

  # if need to modiy the probabilities to turn into daily probability (assuming lag of x days between mf count measurement ~ clinical condition test) #
  # SI_probs <- SI_probs / 30 # divide by 5 (if potential lag of up to 5 days between mf count taken and clinical condition measure, so daily rate / 5)
  # RSD_probs <- RSD_probs / 30
  # Atrp_probs <- Atrp_probs / 30
  # Hg_probs <- Hg_probs / 30
  # Depigm_probs <- Depigm_probs / 30

  # extract number of mf per skin snip per individual
  morb.mat.tmp[,5] <- round(temp_mf[[2]]) + 1 # mf per skin snip for all individuals (+1 because of indexing so that when indexing probabilities goes from 1)

  # extract probabilities (rates) to run Bernoulli trial for each condition
  morb.mat.tmp[,38] <- ifelse(morb.mat.tmp[,5] > 0, SI_probs[morb.mat.tmp[,5]], 0) # severe itch rates
  #morb.mat.tmp[,39] <- ifelse(morb.mat.tmp[,5] > 0, APOD_probs[morb.mat.tmp[,5]], 0) # APOD rates
  #morb.mat.tmp[,40] <- ifelse(morb.mat.tmp[,5] > 0, CPOD_probs[morb.mat.tmp[,5]], 0) # CPOD rates
  #morb.mat.tmp[,41] <- ifelse(morb.mat.tmp[,5] > 0, LOD_probs[morb.mat.tmp[,5]], 0) # LOD rates
  morb.mat.tmp[,42] <- ifelse(morb.mat.tmp[,5] > 0, RSD_probs[morb.mat.tmp[,5]], 0) # RSD rates
  # morb.mat.tmp[,43] <- ifelse(morb.mat.tmp[,5] > 0, Atrp_probs[morb.mat.tmp[,5]], 0) # atrophy rates
  # morb.mat.tmp[,44] <- ifelse(morb.mat.tmp[,5] > 0, Hg_probs[morb.mat.tmp[,5]], 0) # hanging groin rates
  # morb.mat.tmp[,45] <- ifelse(morb.mat.tmp[,5] > 0, Depigm_probs[morb.mat.tmp[,5]], 0) # depigmentation groin rates


  # ======================= #
  # Undergo Bernouli trial  #

  # Bernoulli trial for those to be tested for each condition
  # morb.mat.tmp[,46] <- ifelse(morb.mat.tmp[,30] == 1, rbinom(sum(morb.mat.tmp[,30]), 1, morb.mat.tmp[,38]), 0) # severe itch (revert to 0 every time step if not realised)
  # # morb.mat.tmp[,47] <- ifelse(morb.mat.tmp[,31] == 1, rbinom(sum(morb.mat.tmp[,31]), 1, morb.mat.tmp[,39]), 0) # APOD (revert to 0 every time step if not realised)
  # # morb.mat.tmp[,48] <- ifelse(morb.mat.tmp[,32] == 1, rbinom(sum(morb.mat.tmp[,32]), 1, morb.mat.tmp[,40]), 0) # CPOD (revert to 0 every time step if not realised)
  # # morb.mat.tmp[,49] <- ifelse(morb.mat.tmp[,33] == 1, rbinom(sum(morb.mat.tmp[,33]), 1, morb.mat.tmp[,41]), 0) # LOD (revert to 0 every time step if not realised)
  # morb.mat.tmp[,50] <- ifelse(morb.mat.tmp[,34] == 1, rbinom(sum(morb.mat.tmp[,34]), 1, morb.mat.tmp[,42]), 0) # RSD (revert to 0 every time step if not realised)
  # # morb.mat.tmp[,51] <- ifelse(morb.mat.tmp[,35] == 1, rbinom(sum(morb.mat.tmp[,35]), 1, morb.mat.tmp[,43]), morb.mat.tmp[,51]) # atrophy (irreversible to else based on last evaluation)
  # # morb.mat.tmp[,52] <- ifelse(morb.mat.tmp[,36] == 1, rbinom(sum(morb.mat.tmp[,36]), 1, morb.mat.tmp[,44]), morb.mat.tmp[,52]) # HG (irreversible to else based on last evaluation)
  # # morb.mat.tmp[,53] <- ifelse(morb.mat.tmp[,37] == 1, rbinom(sum(morb.mat.tmp[,37]), 1, morb.mat.tmp[,45]), morb.mat.tmp[,53]) # depigmentation (irreversible to else based on last evaluation)
  # # morb.mat.tmp[,51] <- ifelse(morb.mat.tmp[,35] == 1, rbinom(sum(morb.mat.tmp[,35]), 1, morb.mat.tmp[,43]), morb.mat.tmp[,51]) # atrophy (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # # morb.mat.tmp[,52] <- ifelse(morb.mat.tmp[,36] == 1, rbinom(sum(morb.mat.tmp[,36]), 1, morb.mat.tmp[,44]), morb.mat.tmp[,52]) # HG (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # # morb.mat.tmp[,53] <- ifelse(morb.mat.tmp[,37] == 1, rbinom(sum(morb.mat.tmp[,37]), 1, morb.mat.tmp[,45]), morb.mat.tmp[,53]) # depigmentation (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,51] <- ifelse(morb.mat.tmp[,35] == 1, rbinom(1, 1, morb.mat.tmp[,43]), morb.mat.tmp[,51]) # atrophy (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,52] <- ifelse(morb.mat.tmp[,36] == 1, rbinom(1, 1, morb.mat.tmp[,44]), morb.mat.tmp[,52]) # HG (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,53] <- ifelse(morb.mat.tmp[,37] == 1, rbinom(1, 1, morb.mat.tmp[,45]), morb.mat.tmp[,53]) # depigmentation (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)

  # ========================= #
  # Version 1 current working #

  # #Bernoulli trial for those to be tested for each condition
  # morb.mat.tmp[,46] <- ifelse(morb.mat.tmp[,30] == 1, rbinom(sum(morb.mat.tmp[,30]), 1, morb.mat.tmp[,38]), 0) # severe itch (revert to 0 every time step if not realised)
  # morb.mat.tmp[,50] <- ifelse(morb.mat.tmp[,34] == 1, rbinom(sum(morb.mat.tmp[,34]), 1, morb.mat.tmp[,42]), 0) # RSD (revert to 0 every time step if not realised)
  # morb.mat.tmp[,51] <- ifelse(morb.mat.tmp[,35] == 1, rbinom(sum(morb.mat.tmp[,35]), 1, morb.mat.tmp[,43]), morb.mat.tmp[,51]) # atrophy (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,52] <- ifelse(morb.mat.tmp[,36] == 1, rbinom(sum(morb.mat.tmp[,36]), 1, morb.mat.tmp[,44]), morb.mat.tmp[,52]) # HG (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,53] <- ifelse(morb.mat.tmp[,37] == 1, rbinom(sum(morb.mat.tmp[,37]), 1, morb.mat.tmp[,45]), morb.mat.tmp[,53]) # depigmentation (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)

  # ========================= #
  # new approach              #

  # # Bernoulli trial for those to be tested for each condition
  # # for reversible conditions relate disease probability to mf load
  # # for non-reversible conditions relate disease probability to currently infected

  # WORKING #
  # morb.mat.tmp[,46] <- ifelse(morb.mat.tmp[,30] == 1, rbinom(sum(morb.mat.tmp[,30]), 1, morb.mat.tmp[,38]), 0) # severe itch (stay as prior disease condition status if test does not take place)
  # morb.mat.tmp[,50] <- ifelse(morb.mat.tmp[,34] == 1, rbinom(sum(morb.mat.tmp[,34]), 1, morb.mat.tmp[,42]), 0) # RSD (stay as prior disease condition status if test does not take place)
  # morb.mat.tmp[,51] <- ifelse(morb.mat.tmp[,35] == 1, rbinom(sum(morb.mat.tmp[,35]), 1, Atrp_probs), morb.mat.tmp[,51]) # atrophy (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,52] <- ifelse(morb.mat.tmp[,36] == 1, rbinom(sum(morb.mat.tmp[,36]), 1, Hg_probs), morb.mat.tmp[,52]) # HG (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,53] <- ifelse(morb.mat.tmp[,37] == 1, rbinom(sum(morb.mat.tmp[,37]), 1, Depigm_probs), morb.mat.tmp[,53]) # depigmentation (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)

  # # for non-reversible conditions relate disease probability to mf count (but prob divided by average age)
  # morb.mat.tmp[,46] <- ifelse(morb.mat.tmp[,30] == 1, rbinom(sum(morb.mat.tmp[,30]), 1, morb.mat.tmp[,38]), 0) # severe itch (stay as prior disease condition status if test does not take place)
  # morb.mat.tmp[,50] <- ifelse(morb.mat.tmp[,34] == 1, rbinom(sum(morb.mat.tmp[,34]), 1, morb.mat.tmp[,42]), 0) # RSD (stay as prior disease condition status if test does not take place)
  # morb.mat.tmp[,51] <- ifelse(morb.mat.tmp[,35] == 1, rbinom(sum(morb.mat.tmp[,35]), 1, morb.mat.tmp[,43]), morb.mat.tmp[,51]) # atrophy (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,52] <- ifelse(morb.mat.tmp[,36] == 1, rbinom(sum(morb.mat.tmp[,36]), 1, morb.mat.tmp[,44]), morb.mat.tmp[,52]) # HG (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,53] <- ifelse(morb.mat.tmp[,37] == 1, rbinom(sum(morb.mat.tmp[,37]), 1, morb.mat.tmp[,44]), morb.mat.tmp[,53]) # depigmentation (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)

  # CURRENTLY WORKING #
  # morb.mat.tmp[,46] <- ifelse(morb.mat.tmp[,30] == 1, rbinom(sum(morb.mat.tmp[,30]), 1, SI_probs), 0) # severe itch (stay as prior disease condition status if test does not take place)
  # morb.mat.tmp[,50] <- ifelse(morb.mat.tmp[,34] == 1, rbinom(sum(morb.mat.tmp[,34]), 1, RSD_probs), 0) # RSD (stay as prior disease condition status if test does not take place)
  # morb.mat.tmp[,51] <- ifelse(morb.mat.tmp[,35] == 1, rbinom(sum(morb.mat.tmp[,35]), 1, Atrp_probs), morb.mat.tmp[,51]) # atrophy (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,52] <- ifelse(morb.mat.tmp[,36] == 1, rbinom(sum(morb.mat.tmp[,36]), 1, Hg_probs), morb.mat.tmp[,52]) # HG (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,53] <- ifelse(morb.mat.tmp[,37] == 1, rbinom(sum(morb.mat.tmp[,37]), 1, Depigm_probs), morb.mat.tmp[,53]) # depigmentation (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)

  # NOW TESTING #
  morb.mat.tmp[,46] <- ifelse(morb.mat.tmp[,30] == 1, rbinom(sum(morb.mat.tmp[,30]), 1, SI_probs), morb.mat.tmp[,46]) # severe itch (stay as prior disease condition status if test does not take place)
  morb.mat.tmp[,50] <- ifelse(morb.mat.tmp[,34] == 1, rbinom(sum(morb.mat.tmp[,34]), 1, RSD_probs), morb.mat.tmp[,50]) # RSD (stay as prior disease condition status if test does not take place)
  morb.mat.tmp[,51] <- ifelse(morb.mat.tmp[,35] == 1, rbinom(sum(morb.mat.tmp[,35]), 1, Atrp_probs), morb.mat.tmp[,51]) # atrophy (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  morb.mat.tmp[,52] <- ifelse(morb.mat.tmp[,36] == 1, rbinom(sum(morb.mat.tmp[,36]), 1, Hg_probs), morb.mat.tmp[,52]) # HG (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  morb.mat.tmp[,53] <- ifelse(morb.mat.tmp[,37] == 1, rbinom(sum(morb.mat.tmp[,37]), 1, Depigm_probs), morb.mat.tmp[,53]) # depigmentation (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)


  return(morb.mat.tmp)


}


#' @title
#' update reversible conditions
#'
#' @description
#' update 3-day delay matrix for tracking how long someone has been disease positive (3th day return to 0) & updating morbidity matrix
#'
#' @param sequela.postive.mat1 3-day delay matrix for tracking SI disease state
#' @param sequela.postive.mat2 3-day delay matrix for tracking RSD disease state
#' @param inds.sequela.mat vector to move delay matrix by one day
#' @param morb.mat.tmp morbidity matrix to update (reversible conditions)
#'
#' @returns a) updated SI 3-day delay matrix b) updated RSD 3-day delay matrix c) updated morbidity matrix
update_reversible_sequela_func <- function(sequela.postive.mat1, sequela.postive.mat2, inds.sequela.mat,
                                           morb.mat.tmp){

  # =============== #
  # for severe itch #
  #  Extract current sequela state for reversible conditions
  sequela.postive.mat1[,inds.sequela.mat] <- sequela.postive.mat1[,(inds.sequela.mat-1)] # move sequela state along one col
  sequela.postive.mat1[,1] <- morb.mat.tmp[,46] # update first col with current sequela state on that day

  # new steps : update current sequela state in morb.mat if 3th day of morbidity to 0 #
  morb.mat.tmp[,22] <- sequela.postive.mat1[,3] # assign day 3 from sequela delay matrix
  morb.mat.tmp[,46] <- ifelse(morb.mat.tmp[,22] == 1, 0, morb.mat.tmp[,46]) # update current disease status

  # =============== #
  #     for RSD     #
  #  Extract current sequela state for reversible conditions
  sequela.postive.mat2[,inds.sequela.mat] <- sequela.postive.mat2[,(inds.sequela.mat-1)] # move sequela state along one col
  sequela.postive.mat2[,1] <- morb.mat.tmp[,50] # update first col with current sequela state on that day

  # new steps : update current sequela state in morb.mat if 3th day of morbidity to 0 #
  morb.mat.tmp[,23] <- sequela.postive.mat2[,3] # assign day 7 from sequela delay matrix
  morb.mat.tmp[,50] <- ifelse(morb.mat.tmp[,23] == 1, 0, morb.mat.tmp[,50]) # update current disease status

  return(list(sequela.postive.mat1, sequela.postive.mat2, morb.mat.tmp))

}



#' @title
#' disease prevalences
#'
#' @description
#' calculates disease prevalence for each disease state (including age-stratified disease prevalence's)
#'
#' @param morb.mat.tmp updated matrix highlighting individuals to test
#' @param N human pop size
#' @param SI_prev prevalence vector for severe itch to update
#' @param RSD_prev prevalence vector for RSD to update
#' @param Atrp_prev prevalence vector for atrophy to update
#' @param Hg_prev prevalence vector for hanging groin to update
#' @param depigm_prev prevalence vector for depigmentation to update
#' @param SI_prev0_1 prevalence vector for SI (0 - 1 age group) to update
#' @param SI_prev2_4 prevalence vector for SI (2 - 4 age group) to update
#' @param SI_prev5_9 prevalence vector for SI (5 - 9 age group) to update
#' @param SI_prev10_19 prevalence vector for SI (10 - 19 age group) to update
#' @param SI_prev20_29 prevalence vector for SI (20 - 29 age group) to update
#' @param SI_prev30_49 prevalence vector for SI (30 - 49 age group) to update
#' @param SI_prev50_80 prevalence vector for SI (50 - 80 age group) to update
#' @param RSD_prev10_19 prevalence vector for RSD (0 - 1 age group) to update
#'
#' @returns updated matrix with disease status updated
morbidity_prev_func2 <- function(morb.mat.tmp, N, SI_prev, RSD_prev, Atrp_prev, HG_prev, depigm_prev,
                                SI_prev0_1, SI_prev2_4, SI_prev5_9, SI_prev10_19, SI_prev20_29, SI_prev30_49, SI_prev50_80,
                                RSD_prev0_1, RSD_prev2_4, RSD_prev5_9, RSD_prev10_19,RSD_prev20_29, RSD_prev30_49, RSD_prev50_80,
                                Atrp_prev0_1, Atrp_prev2_4, Atrp_prev5_9, Atrp_prev10_19, Atrp_prev20_29, Atrp_prev30_49, Atrp_prev50_80,
                                HG_prev0_1, HG_prev2_4, HG_prev5_9, HG_prev10_19, HG_prev20_29, HG_prev30_49, HG_prev50_80,
                                depigm_prev0_1, depigm_prev2_4, depigm_prev5_9, depigm_prev10_19, depigm_prev20_29, depigm_prev30_49, depigm_prev50_80)
{

  # calculate current time-step skin disease state prevalence #
  # 2nd attempt : > 5 yr olds (to match Murdoch et al. 2017 prevs only sampling > 5 yr olds)

  #SI_prev_temp <- sum(morb.mat.tmp[,46])/N  # severe itch
  SI_prev_temp <- length(which(morb.mat.tmp[,46] == 1 & morb.mat.tmp[,1] >= 5)) /  length(which(morb.mat.tmp[,1] >= 5)) # prev in > 5yrs
  #RSD_prev_temp <- sum(morb.mat.tmp[,50])/N # reactive skin disease
  RSD_prev_temp <- length(which(morb.mat.tmp[,50] == 1 & morb.mat.tmp[,1] >= 5)) /  length(which(morb.mat.tmp[,1] >= 5)) # prev in > 5yrs
  #Atrp_prev_temp <- sum(morb.mat.tmp[,51])/N # skin atrophy
  Atrp_prev_temp <- length(which(morb.mat.tmp[,51] == 1 & morb.mat.tmp[,1] >= 5)) /  length(which(morb.mat.tmp[,1] >= 5)) # prev in > 5yrs
  #HG_prev_temp <- sum(morb.mat.tmp[,52])/N # hanging groin
  HG_prev_temp <- length(which(morb.mat.tmp[,52] == 1 & morb.mat.tmp[,1] >= 5)) /  length(which(morb.mat.tmp[,1] >= 5)) # prev in > 5yrs
  #depigm_prev_temp <- sum(morb.mat.tmp[,53])/N # depigmentation
  depigm_prev_temp <- length(which(morb.mat.tmp[,53] == 1 & morb.mat.tmp[,1] >= 5)) /  length(which(morb.mat.tmp[,1] >= 5)) # prev in > 5yrs

  # update prevalence vectors
  SI_prev <- c(SI_prev, SI_prev_temp)
  RSD_prev <- c(RSD_prev, RSD_prev_temp)
  Atrp_prev <- c(Atrp_prev, Atrp_prev_temp)
  HG_prev <- c(HG_prev, HG_prev_temp)
  depigm_prev <- c(depigm_prev, depigm_prev_temp)

  # update age-prevalence vectors

  # severe itch age age-prev #
  SI_prev0_1_temp <- length(which(morb.mat.tmp[,46] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
  SI_prev0_1 <- c(SI_prev0_1, SI_prev0_1_temp)

  SI_prev2_4_temp <- length(which(morb.mat.tmp[,46] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
  SI_prev2_4 <- c(SI_prev2_4, SI_prev2_4_temp)

  SI_prev5_9_temp <- length(which(morb.mat.tmp[,46] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
  SI_prev5_9 <- c(SI_prev5_9, SI_prev5_9_temp)

  SI_prev10_19_temp <- length(which(morb.mat.tmp[,46] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
  SI_prev10_19 <- c(SI_prev10_19, SI_prev10_19_temp)

  SI_prev20_29_temp <- length(which(morb.mat.tmp[,46] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
  SI_prev20_29 <- c(SI_prev20_29, SI_prev20_29_temp)

  SI_prev30_49_temp <- length(which(morb.mat.tmp[,46] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
  SI_prev30_49 <- c(SI_prev30_49, SI_prev30_49_temp)

  SI_prev50_80_temp <- length(which(morb.mat.tmp[,46] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
  SI_prev50_80 <- c(SI_prev50_80, SI_prev50_80_temp)

  # RSD age age-prev #
  RSD_prev0_1_temp <- length(which(morb.mat.tmp[,50] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
  RSD_prev0_1 <- c(RSD_prev0_1, RSD_prev0_1_temp)

  RSD_prev2_4_temp <- length(which(morb.mat.tmp[,50] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
  RSD_prev2_4 <- c(RSD_prev2_4, RSD_prev2_4_temp)

  RSD_prev5_9_temp <- length(which(morb.mat.tmp[,50] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
  RSD_prev5_9 <- c(RSD_prev5_9, RSD_prev5_9_temp)

  RSD_prev10_19_temp <- length(which(morb.mat.tmp[,50] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
  RSD_prev10_19 <- c(RSD_prev10_19, RSD_prev10_19_temp)

  RSD_prev20_29_temp <- length(which(morb.mat.tmp[,50] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
  RSD_prev20_29 <- c(RSD_prev20_29, RSD_prev20_29_temp)

  RSD_prev30_49_temp <- length(which(morb.mat.tmp[,50] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
  RSD_prev30_49 <- c(RSD_prev30_49, RSD_prev30_49_temp)

  RSD_prev50_80_temp <- length(which(morb.mat.tmp[,50] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
  RSD_prev50_80 <- c(RSD_prev50_80, RSD_prev50_80_temp)

  # Atrophy age age-prev #
  Atrp_prev0_1_temp <- length(which(morb.mat.tmp[,51] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
  Atrp_prev0_1 <- c(Atrp_prev0_1, Atrp_prev0_1_temp)

  Atrp_prev2_4_temp <- length(which(morb.mat.tmp[,51] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
  Atrp_prev2_4 <- c(Atrp_prev2_4, Atrp_prev2_4_temp)

  Atrp_prev5_9_temp <- length(which(morb.mat.tmp[,51] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
  Atrp_prev5_9 <- c(Atrp_prev5_9, Atrp_prev5_9_temp)

  Atrp_prev10_19_temp <- length(which(morb.mat.tmp[,51] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
  Atrp_prev10_19 <- c(Atrp_prev10_19, Atrp_prev10_19_temp)

  Atrp_prev20_29_temp <- length(which(morb.mat.tmp[,51] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
  Atrp_prev20_29 <- c(Atrp_prev20_29, Atrp_prev20_29_temp)

  Atrp_prev30_49_temp <- length(which(morb.mat.tmp[,51] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
  Atrp_prev30_49 <- c(Atrp_prev30_49, Atrp_prev30_49_temp)

  Atrp_prev50_80_temp <- length(which(morb.mat.tmp[,51] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
  Atrp_prev50_80 <- c(Atrp_prev50_80, Atrp_prev50_80_temp)

  # Hanging groin age age-prev #
  HG_prev0_1_temp <- length(which(morb.mat.tmp[,52] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
  HG_prev0_1 <- c(HG_prev0_1, HG_prev0_1_temp)

  HG_prev2_4_temp <- length(which(morb.mat.tmp[,52] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
  HG_prev2_4 <- c(HG_prev2_4, HG_prev2_4_temp)

  HG_prev5_9_temp <- length(which(morb.mat.tmp[,52] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
  HG_prev5_9 <- c(HG_prev5_9, HG_prev5_9_temp)

  HG_prev10_19_temp <- length(which(morb.mat.tmp[,52] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
  HG_prev10_19 <- c(HG_prev10_19, HG_prev10_19_temp)

  HG_prev20_29_temp <- length(which(morb.mat.tmp[,52] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
  HG_prev20_29 <- c(HG_prev20_29, HG_prev20_29_temp)

  HG_prev30_49_temp <- length(which(morb.mat.tmp[,52] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
  HG_prev30_49 <- c(HG_prev30_49, HG_prev30_49_temp)

  HG_prev50_80_temp <- length(which(morb.mat.tmp[,52] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
  HG_prev50_80 <- c(HG_prev50_80, HG_prev50_80_temp)

  # depigmentation age age-prev #
  depigm_prev0_1_temp <- length(which(morb.mat.tmp[,53] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
  depigm_prev0_1 <- c(depigm_prev0_1, depigm_prev0_1_temp)

  depigm_prev2_4_temp <- length(which(morb.mat.tmp[,53] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
  depigm_prev2_4 <- c(depigm_prev2_4, depigm_prev2_4_temp)

  depigm_prev5_9_temp <- length(which(morb.mat.tmp[,53] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
  depigm_prev5_9 <- c(depigm_prev5_9, depigm_prev5_9_temp)

  depigm_prev10_19_temp <- length(which(morb.mat.tmp[,53] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
  depigm_prev10_19 <- c(depigm_prev10_19, depigm_prev10_19_temp)

  depigm_prev20_29_temp <- length(which(morb.mat.tmp[,53] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
  depigm_prev20_29 <- c(depigm_prev20_29, depigm_prev20_29_temp)

  depigm_prev30_49_temp <- length(which(morb.mat.tmp[,53] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
  depigm_prev30_49 <- c(depigm_prev30_49, depigm_prev30_49_temp)

  depigm_prev50_80_temp <- length(which(morb.mat.tmp[,53] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
  depigm_prev50_80 <- c(depigm_prev50_80, depigm_prev50_80_temp)

  return(list(SI_prev, RSD_prev, Atrp_prev, HG_prev, depigm_prev,
              SI_prev0_1, SI_prev2_4, SI_prev5_9, SI_prev10_19, SI_prev20_29, SI_prev30_49, SI_prev50_80,
              RSD_prev0_1, RSD_prev2_4, RSD_prev5_9, RSD_prev10_19,RSD_prev20_29, RSD_prev30_49, RSD_prev50_80,
              Atrp_prev0_1, Atrp_prev2_4, Atrp_prev5_9, Atrp_prev10_19, Atrp_prev20_29, Atrp_prev30_49, Atrp_prev50_80,
              HG_prev0_1, HG_prev2_4, HG_prev5_9, HG_prev10_19, HG_prev20_29, HG_prev30_49, HG_prev50_80,
              depigm_prev0_1, depigm_prev2_4, depigm_prev5_9, depigm_prev10_19, depigm_prev20_29, depigm_prev30_49, depigm_prev50_80))

}


# ================================================================================================================= #
#                                Eye disease functions                                                              #
# ================================================================================================================= #

#' @title
#' find individuals which develop morbidity
#'
#' @description
#' find individuals sampled between specific ages with current mf and can then develop morbidity in the next function
#' (determined by probability associated with individual mf load or fixed probability)
#'
#' @param dat this is the master matrix (all.mats.temp) tracking mf (by age) for each individual
#' @param mf.start starting column for mf in master matrix
#' @param mf.end final column for mf in master matrix
#' @param morb.mat.temp matrix containing columns to determining conditions for testing morbidity (eye disease)
#' @param age_to_samp_vec_nonreversible 1 year increments of increasing age for when an individual should be tested (matching age), between 20 - 80 yrs
#'
#' @returns updated matrix with identifying individuals to test for morbidity
find_indiv_totest_func3 <- function(dat, mf.start, mf.end, morb.mat.tmp, age_to_samp_vec_nonreversible){

  # update age (and sex for newborns)
  morb.mat.tmp[,1] <- dat[,2]
  morb.mat.tmp[,2] <- dat[,3]

  # lagged ages (age 2 year in future) #
  morb.mat.tmp[,10] <- morb.mat.tmp[,1] + 2 # current age + 2 years (in 2 yrs time)
  morb.mat.tmp[,11] <- ifelse(morb.mat.tmp[,10] > 79.99999999, morb.mat.tmp[,10] - 80, morb.mat.tmp[,10]) # if between 78 - 80, will be > 80 yrs in 2 yrs time
                                                                                                          # therefore not alive, so newborn will be future age - 80 yrs (e.g., 81 - 80 = 1 yr old)
                                                                                                          # need to ensure individuals between 0 - 2 yrs with 0 blindness due to lag (below)
  # true number of mf per individual #
  mf_all <- rowSums(dat[, mf.start : mf.end]) # sum mf per individual across all 21 age classes
  morb.mat.tmp[,4] <- mf_all # true mf count

  # ======================== #
  #  1 ) age sampling        #

  # ======================================================== #
  # # new approach: determine if age to samp matches for all conditions (only want to sample once per year of age, not every time-step)
  morb.mat.tmp[,6] <- ifelse((round(morb.mat.tmp[,1],6) %in% round(age_to_samp_vec_nonreversible,6)),1,0) # blindness

  # ==================================#
  # 2) WHOM TO UNDERGO BERNOULI TRAIL #

  morb.mat.tmp[,7] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,6] == 1 & morb.mat.tmp[,9] == 0, 1, 0) # blindness (irreversible; only test once in age range)

 return(morb.mat.tmp)

}



#' @title
#' individuals designated to be tested undergo Bernoulli trial for realize eye disease status
#'
#' @description
#' each individual undergoes a Bernoulli trial (using probabilities based on mf counts) to ascertain new cases of eye disease
#'
#' @param morb.mat.tmp updated matrix highlighting individuals to test
#' @param temp_mf vector of all mf per skin snip for each individual
#' @param SI_probs probabilities of severe itch for a given mf count
#' @param RSD_probs probabilities of RSD for a given mf count
#' @param Atrp_probs probabilities of Atrophy for a given mf count
#' @param Hg_probs probabilities of hanging groin for a given mf count
#' @param Depigm_probs probabilities of depigmentation for a given mf count
#'
#' @returns updated matrix with disease status updated
new_cases_morbidity_func3 <- function(morb.mat.tmp, temp.mf, blind.probs, dat, i, lagged.mat.tmp){

  # extract number of mf per skin snip per individual
  morb.mat.tmp[,5] <- round(temp.mf[[2]]) + 1 # mf per skin snip for all individuals (+1 because of indexing so that when indexing probabilities goes from 1)

  # extract probabilities (rates) to run Bernoulli trial for each condition
  morb.mat.tmp[,8] <- ifelse(morb.mat.tmp[,5] > 0, blind.probs[morb.mat.tmp[,5]], 0) # blindness rate/ prob ~ mf count

  # ======================= #
  # Undergo Bernouli trial  #

  morb.mat.tmp[,9] <- ifelse(morb.mat.tmp[,7] == 1, rbinom(sum(morb.mat.tmp[,7]), 1, morb.mat.tmp[,8]), morb.mat.tmp[,9]) # blindness (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)

  # update blindness status based on age in 2 years (i.e., those 0-2 yrs in 2 yrs time will be 0 blindness) #

  morb.mat.tmp[,12] <- ifelse(morb.mat.tmp[,11] < 2, 0, morb.mat.tmp[,9]) # set to 0 for 0-2 yrs in 2 yrs time, or as in col 9

  # ========================================================================================== #
  #  update lagged matrix (current blindness outcome related to blindness onset 2 years on )   #

  lagged.mat.tmp[,1] <- dat[,2] # update age
  lagged.mat.tmp[,2] <- dat[,3] # update sex

  col.to.select <- 2 + i + 731 # column to update with lagged blindness outcomes :
                               # 2 + because first two columsn are current age and sex, so iter one from thrid col
                               # 731 as want update 2 years - dt on from current iter (e.g., 2 year lag)

  lagged.mat.tmp[,col.to.select] <- morb.mat.tmp[,9] # update lagged iter (column) with current disease state


  return(list(morb.mat.tmp, lagged.mat.tmp))


}




#' @title
#' eye disease prevalence
#'
#' @description
#' calculates disease prevalence for each disease state (including age-stratified disease prevalence's)
#'
#' @param i iteration
#' @param lagged.mat.tmp updated matrix highlighting individuals to test
#' @param N human pop size
#' @param blind_prev prevalence vector for blindness to update
#' @param visual_imp_prev prevalence vector for visual impairment to update
#' @param blind_prev0_1 prevalence vector for blindness (0 - 1 age group) to update
#' @param blind_prev2_4 prevalence vector for blindness (2 - 4 age group) to update
#' @param blind_prev5_9 prevalence vector for blindness (5 - 9 age group) to update
#' @param blind_prev10_19 prevalence vector for blindness (10 - 19 age group) to update
#' @param blind_prev20_29 prevalence vector for blindness (20 - 29 age group) to update
#' @param blind_prev30_49 prevalence vector for blindness (30 - 49 age group) to update
#' @param blind_prev50_80 prevalence vector for blindness (50 - 80 age group) to update
#' @param visual_imp_prev0_1 prevalence vector for visual impairment (0 - 1 age group) to update
#' @param visual_imp_prev2_4 prevalence vector for visual impairment (2 - 4 age group) to update
#' @param visual_imp_prev5_9 prevalence vector for visual impairment (5 - 9 age group) to update
#' @param visual_imp_prev10_19 prevalence vector for visual impairment (10 - 19 age group) to update
#' @param visual_imp_prev20_29 prevalence vector for visual impairment (20 - 29 age group) to update
#' @param visual_imp_prev30_49 prevalence vector for visual impairment (30 - 49 age group) to update
#' @param visual_imp_prev50_80 prevalence vector for visual impairment (50 - 80 age group) to update
#'
#' @returns updated matrix with disease status updated
eye.disease.prev.func <- function(i, lagged.mat.tmp, N, blind_prev, visual_imp_prev,
                                  blind_prev0_1, blind_prev2_4, blind_prev5_9, blind_prev10_19,
                                  blind_prev20_29, blind_prev30_49, blind_prev50_80,
                                  visual_imp_prev0_1, visual_imp_prev2_4, visual_imp_prev5_9,
                                  visual_imp_prev10_19, visual_imp_prev20_29, visual_imp_prev30_49,
                                  visual_imp_prev50_80,
                                  morb.mat.tmp,
                                  blind_prev1, visual_imp_prev1,
                                  blind_prev0_1a, blind_prev2_4a, blind_prev5_9a, blind_prev10_19a, blind_prev20_29a,
                                  blind_prev30_49a, blind_prev50_80a,
                                  visual_imp_prev0_1a, visual_imp_prev2_4a, visual_imp_prev5_9a, visual_imp_prev10_19a,
                                  visual_imp_prev20_29a, visual_imp_prev30_49a, visual_imp_prev50_80a)
{


  # calculate for current time-step, eye disease state prevalence using lagged matrix #
  # overall prevalence calculated in > 5 yr olds

  col_slct <- i + 2 # iter 1 begins from col 3 in lagged matrix

  blind_prev_temp <- length(which(lagged.mat.tmp[,col_slct] == 1 & lagged.mat.tmp[,1] >= 5)) /  length(which(lagged.mat.tmp[,1] >= 5)) # prev in > 5yrs
  visual_imp_prev_temp <- blind_prev_temp * 4/3

  # update prevalence vectors
  blind_prev <- c(blind_prev, blind_prev_temp)
  visual_imp_prev <- c(visual_imp_prev, visual_imp_prev_temp)

  # update age-prevalence vectors

  # blind age age-prev #
  blind_prev0_1_temp <- length(which(lagged.mat.tmp[,col_slct] == 1 & lagged.mat.tmp[,1] >= 0 & lagged.mat.tmp[,1] < 2)) /  length(which(lagged.mat.tmp[,1] >= 0 & lagged.mat.tmp[,1] < 2))# 0 - 1 age
  blind_prev0_1 <- c(blind_prev0_1, blind_prev0_1_temp)

  blind_prev2_4_temp <- length(which(lagged.mat.tmp[,col_slct] == 1 & lagged.mat.tmp[,1] >= 2 & lagged.mat.tmp[,1] < 5)) /  length(which(lagged.mat.tmp[,1] >= 2 & lagged.mat.tmp[,1] < 5))# 2 - 4 age
  blind_prev2_4 <- c(blind_prev2_4, blind_prev2_4_temp)

  blind_prev5_9_temp <- length(which(lagged.mat.tmp[,col_slct] == 1 & lagged.mat.tmp[,1] >= 5 & lagged.mat.tmp[,1] < 10)) /  length(which(lagged.mat.tmp[,1] >= 5 & lagged.mat.tmp[,1] < 10))
  blind_prev5_9 <- c(blind_prev5_9, blind_prev5_9_temp)

  blind_prev10_19_temp <- length(which(lagged.mat.tmp[,col_slct] == 1 & lagged.mat.tmp[,1] >= 10 & lagged.mat.tmp[,1] < 20)) /  length(which(lagged.mat.tmp[,1] >= 10 & lagged.mat.tmp[,1] < 20))
  blind_prev10_19 <- c(blind_prev10_19, blind_prev10_19_temp)

  blind_prev20_29_temp <- length(which(lagged.mat.tmp[,col_slct] == 1 & lagged.mat.tmp[,1] >= 20 & lagged.mat.tmp[,1] < 30)) /  length(which(lagged.mat.tmp[,1] >= 20 & lagged.mat.tmp[,1] < 30))
  blind_prev20_29 <- c(blind_prev20_29, blind_prev20_29_temp)

  blind_prev30_49_temp <- length(which(lagged.mat.tmp[,col_slct] == 1 & lagged.mat.tmp[,1] >= 30 & lagged.mat.tmp[,1] < 50)) /  length(which(lagged.mat.tmp[,1] >= 30 & lagged.mat.tmp[,1] < 50))
  blind_prev30_49 <- c(blind_prev30_49, blind_prev30_49_temp)

  blind_prev50_80_temp <- length(which(lagged.mat.tmp[,col_slct] == 1 & lagged.mat.tmp[,1] >= 50 & lagged.mat.tmp[,1] <= 80)) /  length(which(lagged.mat.tmp[,1] >= 50 & lagged.mat.tmp[,1] < 80))
  blind_prev50_80 <- c(blind_prev50_80, blind_prev50_80_temp)

  # visual impairement age age-prev #
  visual_imp_prev0_1_temp <- blind_prev0_1_temp * 4/3 # 0 - 1 age
  visual_imp_prev0_1 <- c(visual_imp_prev0_1, visual_imp_prev0_1_temp)

  visual_imp_prev2_4_temp <- blind_prev2_4_temp * 4/3 # 2 - 4 age
  visual_imp_prev2_4 <- c(visual_imp_prev2_4, visual_imp_prev2_4_temp)

  visual_imp_prev5_9_temp <- blind_prev5_9_temp * 4/3
  visual_imp_prev5_9 <- c(visual_imp_prev5_9, visual_imp_prev5_9_temp)

  visual_imp_prev10_19_temp <- blind_prev10_19_temp * 4/3
  visual_imp_prev10_19 <- c(visual_imp_prev10_19, visual_imp_prev10_19_temp)

  visual_imp_prev20_29_temp <- blind_prev20_29_temp * 4/3
  visual_imp_prev20_29 <- c(visual_imp_prev20_29, visual_imp_prev20_29_temp)

  visual_imp_prev30_49_temp <- blind_prev30_49_temp * 4/3
  visual_imp_prev30_49 <- c(visual_imp_prev30_49, visual_imp_prev30_49_temp)

  visual_imp_prev50_80_temp <- blind_prev50_80_temp * 4/3
  visual_imp_prev50_80 <- c(visual_imp_prev50_80, visual_imp_prev50_80_temp)


  # ============================================= #
  #    Test by updating a lagged vector approach  #

  # blind_prev_temp1 <- length(which(morb.mat.tmp[,9] == 1 & morb.mat.tmp[,1] >= 5)) /  length(which(morb.mat.tmp[,1] >= 5)) # prev in > 5yrs
  # visual_imp_prev_temp1 <- blind_prev_temp1 * 4/3
  #
  # # update prevalence vectors
  # blind_prev1 <- c(blind_prev1, blind_prev_temp1)
  # visual_imp_prev1 <- c(visual_imp_prev1, visual_imp_prev_temp1)
  #
  # # update age-prevalence vectors
  #
  # # blind age age-prev #
  # blind_prev0_1_temp1 <- length(which(morb.mat.tmp[,9] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
  # blind_prev0_1a <- c(blind_prev0_1a, blind_prev0_1_temp1)
  #
  # blind_prev2_4_temp1 <- length(which(morb.mat.tmp[,9] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
  # blind_prev2_4a <- c(blind_prev2_4a, blind_prev2_4_temp1)
  #
  # blind_prev5_9_temp1 <- length(which(morb.mat.tmp[,9] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
  # blind_prev5_9a <- c(blind_prev5_9a, blind_prev5_9_temp1)
  #
  # blind_prev10_19_temp1 <- length(which(morb.mat.tmp[,9] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
  # blind_prev10_19a <- c(blind_prev10_19a, blind_prev10_19_temp1)
  #
  # blind_prev20_29_temp1 <- length(which(morb.mat.tmp[,9] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
  # blind_prev20_29a <- c(blind_prev20_29a, blind_prev20_29_temp1)
  #
  # blind_prev30_49_temp1 <- length(which(morb.mat.tmp[,9] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
  # blind_prev30_49a <- c(blind_prev30_49a, blind_prev30_49_temp1)
  #
  # blind_prev50_80_temp1 <- length(which(morb.mat.tmp[,9] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
  # blind_prev50_80a <- c(blind_prev50_80a, blind_prev50_80_temp1)
  #
  # # visual impairement age age-prev #
  # visual_imp_prev0_1_temp1 <- blind_prev0_1_temp1 * 4/3 # 0 - 1 age
  # visual_imp_prev0_1a <- c(visual_imp_prev0_1a, visual_imp_prev0_1_temp1)
  #
  # visual_imp_prev2_4_temp1 <- blind_prev2_4_temp1 * 4/3 # 2 - 4 age
  # visual_imp_prev2_4a <- c(visual_imp_prev2_4a, visual_imp_prev2_4_temp1)
  #
  # visual_imp_prev5_9_temp1 <- blind_prev5_9_temp1 * 4/3
  # visual_imp_prev5_9a <- c(visual_imp_prev5_9a, visual_imp_prev5_9_temp1)
  #
  # visual_imp_prev10_19_temp1 <- blind_prev10_19_temp1 * 4/3
  # visual_imp_prev10_19a <- c(visual_imp_prev10_19a, visual_imp_prev10_19_temp1)
  #
  # visual_imp_prev20_29_temp1 <- blind_prev20_29_temp1 * 4/3
  # visual_imp_prev20_29a <- c(visual_imp_prev20_29a, visual_imp_prev20_29_temp1)
  #
  # visual_imp_prev30_49_temp1 <- blind_prev30_49_temp1 * 4/3
  # visual_imp_prev30_49a <- c(visual_imp_prev30_49a, visual_imp_prev30_49_temp1)
  #
  # visual_imp_prev50_80_temp1 <- blind_prev50_80_temp1 * 4/3
  # visual_imp_prev50_80a <- c(visual_imp_prev50_80a, visual_imp_prev50_80_temp1)
  #

  # =============================================================================== #
  # based on age in 2 years (col 11) and updated blindness status in 2 yrs (col 12)

  blind_prev_temp1 <- length(which(morb.mat.tmp[,12] == 1 & morb.mat.tmp[,11] >= 5)) /  length(which(morb.mat.tmp[,11] >= 5)) # prev in > 5yrs
  visual_imp_prev_temp1 <- blind_prev_temp1 * 4/3

  # update prevalence vectors
  blind_prev1 <- c(blind_prev1, blind_prev_temp1)
  visual_imp_prev1 <- c(visual_imp_prev1, visual_imp_prev_temp1)

  # update age-prevalence vectors
  # based on age in 2 years (col 11) and updated blindness status in 2 yrs (col 12)

  # blind age age-prev #
  blind_prev0_1_temp1 <- length(which(morb.mat.tmp[,12] == 1 & morb.mat.tmp[,11] >= 0 & morb.mat.tmp[,11] < 2)) /  length(which(morb.mat.tmp[,11] >= 0 & morb.mat.tmp[,11] < 2))# 0 - 1 age
  blind_prev0_1a <- c(blind_prev0_1a, blind_prev0_1_temp1)

  blind_prev2_4_temp1 <- length(which(morb.mat.tmp[,12] == 1 & morb.mat.tmp[,11] >= 2 & morb.mat.tmp[,11] < 5)) /  length(which(morb.mat.tmp[,11] >= 2 & morb.mat.tmp[,11] < 5))# 2 - 4 age
  blind_prev2_4a <- c(blind_prev2_4a, blind_prev2_4_temp1)

  blind_prev5_9_temp1 <- length(which(morb.mat.tmp[,12] == 1 & morb.mat.tmp[,11] >= 5 & morb.mat.tmp[,11] < 10)) /  length(which(morb.mat.tmp[,11] >= 5 & morb.mat.tmp[,11] < 10))
  blind_prev5_9a <- c(blind_prev5_9a, blind_prev5_9_temp1)

  blind_prev10_19_temp1 <- length(which(morb.mat.tmp[,12] == 1 & morb.mat.tmp[,11] >= 10 & morb.mat.tmp[,11] < 20)) /  length(which(morb.mat.tmp[,11] >= 10 & morb.mat.tmp[,11] < 20))
  blind_prev10_19a <- c(blind_prev10_19a, blind_prev10_19_temp1)

  blind_prev20_29_temp1 <- length(which(morb.mat.tmp[,12] == 1 & morb.mat.tmp[,11] >= 20 & morb.mat.tmp[,11] < 30)) /  length(which(morb.mat.tmp[,11] >= 20 & morb.mat.tmp[,11] < 30))
  blind_prev20_29a <- c(blind_prev20_29a, blind_prev20_29_temp1)

  blind_prev30_49_temp1 <- length(which(morb.mat.tmp[,12] == 1 & morb.mat.tmp[,11] >= 30 & morb.mat.tmp[,11] < 50)) /  length(which(morb.mat.tmp[,11] >= 30 & morb.mat.tmp[,11] < 50))
  blind_prev30_49a <- c(blind_prev30_49a, blind_prev30_49_temp1)

  blind_prev50_80_temp1 <- length(which(morb.mat.tmp[,12] == 1 & morb.mat.tmp[,11] >= 50 & morb.mat.tmp[,11] <= 80)) /  length(which(morb.mat.tmp[,11] >= 50 & morb.mat.tmp[,11] < 80))
  blind_prev50_80a <- c(blind_prev50_80a, blind_prev50_80_temp1)

  # visual impairement age age-prev #
  visual_imp_prev0_1_temp1 <- blind_prev0_1_temp1 * 4/3 # 0 - 1 age
  visual_imp_prev0_1a <- c(visual_imp_prev0_1a, visual_imp_prev0_1_temp1)

  visual_imp_prev2_4_temp1 <- blind_prev2_4_temp1 * 4/3 # 2 - 4 age
  visual_imp_prev2_4a <- c(visual_imp_prev2_4a, visual_imp_prev2_4_temp1)

  visual_imp_prev5_9_temp1 <- blind_prev5_9_temp1 * 4/3
  visual_imp_prev5_9a <- c(visual_imp_prev5_9a, visual_imp_prev5_9_temp1)

  visual_imp_prev10_19_temp1 <- blind_prev10_19_temp1 * 4/3
  visual_imp_prev10_19a <- c(visual_imp_prev10_19a, visual_imp_prev10_19_temp1)

  visual_imp_prev20_29_temp1 <- blind_prev20_29_temp1 * 4/3
  visual_imp_prev20_29a <- c(visual_imp_prev20_29a, visual_imp_prev20_29_temp1)

  visual_imp_prev30_49_temp1 <- blind_prev30_49_temp1 * 4/3
  visual_imp_prev30_49a <- c(visual_imp_prev30_49a, visual_imp_prev30_49_temp1)

  visual_imp_prev50_80_temp1 <- blind_prev50_80_temp1 * 4/3
  visual_imp_prev50_80a <- c(visual_imp_prev50_80a, visual_imp_prev50_80_temp1)


  return(list(blind_prev, visual_imp_prev,
              blind_prev0_1, blind_prev2_4, blind_prev5_9, blind_prev10_19, blind_prev20_29,
              blind_prev30_49, blind_prev50_80,
              visual_imp_prev0_1, visual_imp_prev2_4, visual_imp_prev5_9, visual_imp_prev10_19,
              visual_imp_prev20_29, visual_imp_prev30_49, visual_imp_prev50_80,
              blind_prev1, visual_imp_prev1,
              blind_prev0_1a, blind_prev2_4a, blind_prev5_9a, blind_prev10_19a, blind_prev20_29a,
              blind_prev30_49a, blind_prev50_80a,
              visual_imp_prev0_1a, visual_imp_prev2_4a, visual_imp_prev5_9a, visual_imp_prev10_19a,
              visual_imp_prev20_29a, visual_imp_prev30_49a, visual_imp_prev50_80a))

}


