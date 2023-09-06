
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
#' @param morb_mat dataframe containing columns to determining conditions for testing morbidity
#' @param age_to_samp_vec_reversible 1 year increments of increasing age for when an individual should be tested (matching age), between 5 - 80 yrs
#' @param age_to_samp_vec_nonreversible 1 year increments of increasing age for when an individual should be tested (matching age), between 20 - 80 yrs
#'
#' @returns updated matrix with identifying individuals to test for morbidity
find_indiv_totest_func <- function(dat, mf.start, mf.end, morb_mat, temp_mf, age_to_samp_vec_reversible,
                                   age_to_samp_vec_nonreversible){

  # update age (and sex for newborns)
  morb_mat$Age <- dat[,2]
  morb_mat$Sex <- dat[,3]
  ages <- morb_mat$Age
  sexes <- morb_mat$Sex

  # true number of mf per individual #
  mf_all <- rowSums(dat[, mf.start : mf.end]) # sum mf per individual across all 21 age classes
  morb_mat$TrueMFCount <- mf_all # true mf count
  trueMFCount <- morb_mat$TrueMFCount

  # extract number of mf per skin snip per individual
  morb_mat$ObservedMFCount <- round(temp_mf[[2]])
  observedMFCounts <- morb_mat$ObservedMFCount

  # ======================== #
  #  1 ) age sampling        #

  # ================================================================================================================== #
  # determine if age to samp matches for all conditions (only want to sample once per year of age, not every time-step) #
  # morb.mat.tmp[,14] <- ifelse((round(ages,6) %in% round(age_to_samp_vec_reversible,6)),1,0) # severe itch (age vec to 6 decimla places more most specific matching - only once per year)
  # morb.mat.tmp[,18] <- ifelse((round(ages,6) %in% round(age_to_samp_vec_reversible,6)),1,0) # RSD (test again single RSD age_to_sample e.g., col 10)
  morb_mat$AtrophySampleAges <- ifelse((round(ages,6) %in% round(age_to_samp_vec_nonreversible,6)),1,0) # atrophy
  atrophySampleAges <- morb_mat$AtrophySampleAges
  morb_mat$HangingGroinSampleAges <- ifelse((round(ages,6) %in% round(age_to_samp_vec_nonreversible,6)),1,0) # hanging groin
  hangingGroinSampleAges <- morb_mat$HangingGroinSampleAges
  morb_mat$DepigSampleAge <- ifelse((round(ages,6) %in% round(age_to_samp_vec_nonreversible,6)),1,0) # depigmentation
  depigSampleAges <- morb_mat$DepigSampleAges


  # ==================================#
  # WHOM TO UNDERGO BERNOULI TRAIL #

  severeItchStatuses <- morb_mat$SevereItchStatus
  rsdStatuses <- morb_mat$RSDStatus
  atrophyStatuses <- morb_mat$AtrophyStatus
  hgStatuses <- morb_mat$HGStatus
  depigStatuses <- morb_mat$DepigStatus
  # # based on > 0 true mf
  # morb.mat.tmp[,30] <- ifelse(trueMFCount > 0 & severeItchStatuses == 0 & ages >= 2, 1, 0) # for SI (if 0 disease state)
  # morb.mat.tmp[,34] <- ifelse(trueMFCount > 0 & rsdStatuses == 0 & ages >= 2, 1, 0) # for RSD (if 0 disease state)
  # morb.mat.tmp[,35] <- ifelse(trueMFCount > 0 & atrophySampleAges == 1 & atrophyStatuses == 0, 1, 0) # atrophy (irreversible; only test once in age range)
  # morb.mat.tmp[,36] <- ifelse(trueMFCount > 0 & hangingGroinSampleAges == 1 & hgStatuses == 0, 1, 0) # HG (irreversible; only test once in age range)
  # morb.mat.tmp[,37] <- ifelse(trueMFCount > 0 & depigSampleAges == 1 & depigStatuses == 0, 1, 0) # depigm (irreversible; only test once in age range)

  # based on > 0 observed mf
  morb_mat$ToTestSevereItch <- ifelse(observedMFCounts > 0 & severeItchStatuses == 0 & ages >= 2, 1, 0) # for SI (if 0 disease state)
  morb_mat$ToTestRSD <- ifelse(observedMFCounts > 0 & rsdStatuses == 0 & ages >= 2, 1, 0) # for RSD (if 0 disease state)
  morb_mat$ToTestAtrophy <- ifelse(observedMFCounts > 0 & atrophySampleAges == 1 & atrophyStatuses == 0, 1, 0) # atrophy (irreversible; only test once in age range)
  morb_mat$ToTestHG <- ifelse(observedMFCounts > 0 & hangingGroinSampleAges == 1 & hgStatuses == 0, 1, 0) # HG (irreversible; only test once in age range)
  morb_mat$ToTestDepig <- ifelse(observedMFCounts > 0 & depigSampleAges == 1 & depigStatuses == 0, 1, 0) # depigm (irreversible; only test once in age range)

  return(morb_mat)

}

#' @title
#' individuals designated to be tested undergo Bernoulli trial for realise skin disease status
#'
#' @description
#' each individual undergoes a Bernoulli trial (using probabilities based on mf counts) to ascertain new cases of OAE
#'
#' @param morb_mat updated dataframe highlighting individuals to test
#' @param temp_mf vector of all mf per skin snip for each individual
#' @param SI_probs probabilities of severe itch for a given mf count
#' @param RSD_probs probabilities of RSD for a given mf count
#' @param Atrp_probs probabilities of Atrophy for a given mf count
#' @param Hg_probs probabilities of hanging groin for a given mf count
#' @param Depigm_probs probabilities of depigmentation for a given mf count
#'
#' @returns updated matrix with disease status updated
new_cases_morbidity_func <- function(morb_mat, SI_probs, RSD_probs, Atrp_probs, Hg_probs, Depigm_probs){

  ## extract number of mf per skin snip per individual
  # morb.mat.tmp[,5] <- round(temp_mf[[2]]) # mf per skin snip for all individuals (+1 because of indexing so that when indexing probabilities goes from 1)
  observedMFCounts <- morb_mat$ObservedMFCount

  # # extract probabilities (rates) to run Bernoulli trial for each condition
  # morb.mat.tmp[,38] <- ifelse(observedMFCounts > 0, SI_probs[observedMFCounts], 0) # severe itch rates
  # #morb.mat.tmp[,39] <- ifelse(observedMFCounts > 0, APOD_probs[observedMFCounts], 0) # APOD rates
  # #morb.mat.tmp[,40] <- ifelse(observedMFCounts > 0, CPOD_probs[observedMFCounts], 0) # CPOD rates
  # #morb.mat.tmp[,41] <- ifelse(observedMFCounts > 0, LOD_probs[observedMFCounts], 0) # LOD rates
  # morb.mat.tmp[,42] <- ifelse(observedMFCounts > 0, RSD_probs[observedMFCounts], 0) # RSD rates
  # morb.mat.tmp[,43] <- ifelse(observedMFCounts > 0, Atrp_probs[observedMFCounts], 0) # atrophy rates
  # morb.mat.tmp[,44] <- ifelse(observedMFCounts > 0, Hg_probs[observedMFCounts], 0) # hanging groin rates
  # morb.mat.tmp[,45] <- ifelse(observedMFCounts > 0, Depigm_probs[observedMFCounts], 0) # depigmentation groin rates


  # ======================= #
  # Undergo Bernouli trial  #

  # based on whether true mf present
  toTestSevereItch <- morb_mat$ToTestSevereItch
  severeItchStatuses <- morb_mat$SevereItchStatus
  toTestRSD <- morb_mat$ToTestRSD
  rsdStatuses <- morb_mat$RSDStatus
  toTestAtrophy <- morb_mat$ToTestAtrophy
  atrophyStatuses <- morb_mat$AtrophyStatus
  toTestHG <- morb_mat$ToTestHG
  hgStatuses <- morb_mat$HGStatus
  toTestDepig <- morb_mat$ToTestDepig
  depigStatuses <- morb_mat$DepigStatus
  morb_mat$SevereItchStatus <- ifelse(toTestSevereItch == 1, rbinom(sum(toTestSevereItch), 1, SI_probs), severeItchStatuses) # severe itch (stay as prior disease condition status if test does not take place)
  morb_mat$RSDStatus <- ifelse(toTestRSD == 1, rbinom(sum(toTestRSD), 1, RSD_probs), rsdStatuses) # RSD (stay as prior disease condition status if test does not take place)
  morb_mat$AtrophyStatus <- ifelse(toTestAtrophy == 1, rbinom(sum(toTestAtrophy), 1, Atrp_probs), atrophyStatuses) # atrophy (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  morb_mat$HGStatus <- ifelse(toTestHG == 1, rbinom(sum(toTestHG), 1, Hg_probs), hgStatuses) # HG (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  morb_mat$DepigStatus <- ifelse(toTestDepig == 1, rbinom(sum(toTestDepig), 1, Depigm_probs), depigStatuses) # depigmentation (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)

  # # based on whether observed mf present
  # morb.mat.tmp[,46] <- ifelse(toTestSevereItch == 1, rbinom(sum(toTestSevereItch), 1, morb.mat.tmp[,38]), severeItchStatuses) # severe itch (stay as prior disease condition status if test does not take place)
  # morb.mat.tmp[,50] <- ifelse(toTestRSD == 1, rbinom(sum(toTestRSD), 1, morb.mat.tmp[,42]), rsdStatuses) # RSD (stay as prior disease condition status if test does not take place)
  # morb.mat.tmp[,51] <- ifelse(toTestAtrophy == 1, rbinom(sum(toTestAtrophy), 1, morb.mat.tmp[,43]), atrophyStatuses) # atrophy (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,52] <- ifelse(toTestHG == 1, rbinom(sum(toTestHG), 1, morb.mat.tmp[,44]), hgStatuses) # HG (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
  # morb.mat.tmp[,53] <- ifelse(toTestDepig == 1, rbinom(sum(toTestDepig), 1, morb.mat.tmp[,45]), depigStatuses) # depigmentation (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)


  return(morb_mat)


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
#' @param morb_mat morbidity matrix to update (reversible conditions)
#'
#' @returns a) updated SI 3-day delay matrix b) updated RSD 3rd day delay matrix c) updated morbidity matrix
update_reversible_sequela_func <- function(sequela.postive.mat1, sequela.postive.mat2, inds.sequela.mat,
                                           morb_mat){

  # # =============== #
  # # for severe itch #

  # note: morb.mat.tmp[,22] = status of whether this is the 3rd day of SI positivity
  # note: morb.mat.tmp[,46] = realized SI state on current day (current iter) - updated to 0 if 4th day of SI positivity
  severeItchStatuses <- morb_mat$SevereItchStatus
  #  Extract current sequela state for reversible conditions
  sequela.postive.mat1[,inds.sequela.mat] <- sequela.postive.mat1[,(inds.sequela.mat-1)] # move sequela state along one col
  sequela.postive.mat1[,1] <- severeItchStatuses # update first col with current sequela state on that day

  # new steps : update current sequela state in morb.mat if day of morbidity to 0 #
  morb_mat$Day3SevereItchStatus <- sequela.postive.mat1[,4] # assign day 3 from sequela delay matrix (4th col of delay mat)
  depigSampleAges <- morb_mat$DepigSampleAge
  day3SevereItchStatus <- morb_mat$Day3SevereItchStatus
  morb_mat$SevereItchStatus <- ifelse(depigSampleAges == 1, 0, severeItchStatuses) # update current disease status

  # #  Extract current sequela state for reversible conditions
  # morb.mat.tmp[,22] <- sequela.postive.mat1[,3] # assign day 3 from sequelae delay matrix
  #
  # sequela.postive.mat1[,inds.sequela.mat] <- sequela.postive.mat1[,(inds.sequela.mat-1)] # move sequela state along one col
  #
  # sequela.postive.mat1[,1] <- severeItchStatuses # update first col with current sequela state on that day
  #
  # morb.mat.tmp[,46] <- ifelse(depigSampleAges == 1, 0, severeItchStatuses) # update current disease status


  # # =============== #
  # #     for RSD     #

  # note: morb.mat.tmp[,23] = status of whether this is the 3rd day of RSD positivity
  # note: morb.mat.tmp[,50] = realized RSD state on current day (current iter) - updated to 0 if 3rd day of RSD positivity

  #  Extract current sequela state for reversible conditions
  rsdStatuses <- morb_mat$RSDStatus
  sequela.postive.mat2[,inds.sequela.mat] <- sequela.postive.mat2[,(inds.sequela.mat-1)] # move sequela state along one col
  sequela.postive.mat2[,1] <- rsdStatuses # update first col with current sequela state on that day

  # new steps : update current sequela state in morb.mat if 3th day of morbidity to 0 #
  morb_mat$Day3RSDStatus <- sequela.postive.mat2[,4] # assign day 3 from sequela delay matrix (4th col of delay mat)
  day3RSDStatus <- morb_mat$Day3RSDStatus
  morb_mat$RSDStatus <- ifelse(day3RSDStatus == 1, 0, rsdStatuses) # update current disease status

  # #  Extract current sequela state for reversible conditions
  #
  # morb.mat.tmp[,23] <- sequela.postive.mat2[,3] # assign day 7 from sequela delay matrix
  #
  # sequela.postive.mat2[,inds.sequela.mat] <- sequela.postive.mat2[,(inds.sequela.mat-1)] # move sequela state along one col
  #
  # sequela.postive.mat2[,1] <- rsdStatuses # update first col with current sequela state on that day
  #
  # morb.mat.tmp[,50] <- ifelse(day3RSDStatus == 1, 0, rsdStatuses) # update current disease status

  return(list(sequela.postive.mat1, sequela.postive.mat2, morb_mat))

}



#' @title
#' disease prevalences
#'
#' @description
#' calculates disease prevalence for each disease state (including age-stratified disease prevalence's)
#'
#' @param morb_mat updated dataframe highlighting individuals to test
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
morbidity_prev_func <- function(morb_mat, N, SI_prev, RSD_prev, Atrp_prev, HG_prev, depigm_prev,
                                SI_prev0_1, SI_prev2_4, SI_prev5_9, SI_prev10_19, SI_prev20_29, SI_prev30_49, SI_prev50_80,
                                RSD_prev0_1, RSD_prev2_4, RSD_prev5_9, RSD_prev10_19,RSD_prev20_29, RSD_prev30_49, RSD_prev50_80,
                                Atrp_prev0_1, Atrp_prev2_4, Atrp_prev5_9, Atrp_prev10_19, Atrp_prev20_29, Atrp_prev30_49, Atrp_prev50_80,
                                HG_prev0_1, HG_prev2_4, HG_prev5_9, HG_prev10_19, HG_prev20_29, HG_prev30_49, HG_prev50_80,
                                depigm_prev0_1, depigm_prev2_4, depigm_prev5_9, depigm_prev10_19, depigm_prev20_29, depigm_prev30_49, depigm_prev50_80)
{
  ages <- morb_mat$Ages
  severeItchStatuses <- morb_mat$SevereItchStatus
  rsdStatuses <- morb_mat$RSDStatus
  atrophyStatuses <- morb_mat$AtrophyStatus
  hgStatuses <- morb_mat$HGStatus
  depigStatuses <- morb_mat$DepigStatus
  # calculate current time-step skin disease state prevalence #
  # 2nd attempt : > 5 yr olds (to match Murdoch et al. 2017 prevs only sampling > 5 yr olds)

  #SI_prev_temp <- sum(severeItchStatuses)/N  # severe itch
  SI_prev_temp <- length(which(severeItchStatuses == 1 & ages >= 5)) /  length(which(ages >= 5)) # prev in > 5yrs
  #RSD_prev_temp <- sum(rsdStatuses)/N # reactive skin disease
  RSD_prev_temp <- length(which(rsdStatuses == 1 & ages >= 5)) /  length(which(ages >= 5)) # prev in > 5yrs
  #Atrp_prev_temp <- sum(atrophyStatuses)/N # skin atrophy
  Atrp_prev_temp <- length(which(atrophyStatuses == 1 & ages >= 5)) /  length(which(ages >= 5)) # prev in > 5yrs
  #HG_prev_temp <- sum(hgStatuses)/N # hanging groin
  HG_prev_temp <- length(which(hgStatuses == 1 & ages >= 5)) /  length(which(ages >= 5)) # prev in > 5yrs
  #depigm_prev_temp <- sum(depigStatuses)/N # depigmentation
  depigm_prev_temp <- length(which(depigStatuses == 1 & ages >= 5)) /  length(which(ages >= 5)) # prev in > 5yrs

  # update prevalence vectors
  SI_prev <- c(SI_prev, SI_prev_temp)
  RSD_prev <- c(RSD_prev, RSD_prev_temp)
  Atrp_prev <- c(Atrp_prev, Atrp_prev_temp)
  HG_prev <- c(HG_prev, HG_prev_temp)
  depigm_prev <- c(depigm_prev, depigm_prev_temp)

  # update age-prevalence vectors

  # severe itch age age-prev #
  SI_prev0_1_temp <- length(which(severeItchStatuses == 1 & ages >= 0 & ages < 2)) /  length(which(ages >= 0 & ages < 2))# 0 - 1 age
  SI_prev0_1 <- c(SI_prev0_1, SI_prev0_1_temp)

  SI_prev2_4_temp <- length(which(severeItchStatuses == 1 & ages >= 2 & ages < 5)) /  length(which(ages >= 2 & ages < 5))# 2 - 4 age
  SI_prev2_4 <- c(SI_prev2_4, SI_prev2_4_temp)

  SI_prev5_9_temp <- length(which(severeItchStatuses == 1 & ages >= 5 & ages < 10)) /  length(which(ages >= 5 & ages < 10))
  SI_prev5_9 <- c(SI_prev5_9, SI_prev5_9_temp)

  SI_prev10_19_temp <- length(which(severeItchStatuses == 1 & ages >= 10 & ages < 20)) /  length(which(ages >= 10 & ages < 20))
  SI_prev10_19 <- c(SI_prev10_19, SI_prev10_19_temp)

  SI_prev20_29_temp <- length(which(severeItchStatuses == 1 & ages >= 20 & ages < 30)) /  length(which(ages >= 20 & ages < 30))
  SI_prev20_29 <- c(SI_prev20_29, SI_prev20_29_temp)

  SI_prev30_49_temp <- length(which(severeItchStatuses == 1 & ages >= 30 & ages < 50)) /  length(which(ages >= 30 & ages < 50))
  SI_prev30_49 <- c(SI_prev30_49, SI_prev30_49_temp)

  SI_prev50_80_temp <- length(which(severeItchStatuses == 1 & ages >= 50 & ages <= 80)) /  length(which(ages >= 50 & ages < 80))
  SI_prev50_80 <- c(SI_prev50_80, SI_prev50_80_temp)

  # RSD age age-prev #
  RSD_prev0_1_temp <- length(which(rsdStatuses == 1 & ages >= 0 & ages < 2)) /  length(which(ages >= 0 & ages < 2))# 0 - 1 age
  RSD_prev0_1 <- c(RSD_prev0_1, RSD_prev0_1_temp)

  RSD_prev2_4_temp <- length(which(rsdStatuses == 1 & ages >= 2 & ages < 5)) /  length(which(ages >= 2 & ages < 5))# 2 - 4 age
  RSD_prev2_4 <- c(RSD_prev2_4, RSD_prev2_4_temp)

  RSD_prev5_9_temp <- length(which(rsdStatuses == 1 & ages >= 5 & ages < 10)) /  length(which(ages >= 5 & ages < 10))
  RSD_prev5_9 <- c(RSD_prev5_9, RSD_prev5_9_temp)

  RSD_prev10_19_temp <- length(which(rsdStatuses == 1 & ages >= 10 & ages < 20)) /  length(which(ages >= 10 & ages < 20))
  RSD_prev10_19 <- c(RSD_prev10_19, RSD_prev10_19_temp)

  RSD_prev20_29_temp <- length(which(rsdStatuses == 1 & ages >= 20 & ages < 30)) /  length(which(ages >= 20 & ages < 30))
  RSD_prev20_29 <- c(RSD_prev20_29, RSD_prev20_29_temp)

  RSD_prev30_49_temp <- length(which(rsdStatuses == 1 & ages >= 30 & ages < 50)) /  length(which(ages >= 30 & ages < 50))
  RSD_prev30_49 <- c(RSD_prev30_49, RSD_prev30_49_temp)

  RSD_prev50_80_temp <- length(which(rsdStatuses == 1 & ages >= 50 & ages <= 80)) /  length(which(ages >= 50 & ages < 80))
  RSD_prev50_80 <- c(RSD_prev50_80, RSD_prev50_80_temp)

  # Atrophy age age-prev #
  Atrp_prev0_1_temp <- length(which(atrophyStatuses == 1 & ages >= 0 & ages < 2)) /  length(which(ages >= 0 & ages < 2))# 0 - 1 age
  Atrp_prev0_1 <- c(Atrp_prev0_1, Atrp_prev0_1_temp)

  Atrp_prev2_4_temp <- length(which(atrophyStatuses == 1 & ages >= 2 & ages < 5)) /  length(which(ages >= 2 & ages < 5))# 2 - 4 age
  Atrp_prev2_4 <- c(Atrp_prev2_4, Atrp_prev2_4_temp)

  Atrp_prev5_9_temp <- length(which(atrophyStatuses == 1 & ages >= 5 & ages < 10)) /  length(which(ages >= 5 & ages < 10))
  Atrp_prev5_9 <- c(Atrp_prev5_9, Atrp_prev5_9_temp)

  Atrp_prev10_19_temp <- length(which(atrophyStatuses == 1 & ages >= 10 & ages < 20)) /  length(which(ages >= 10 & ages < 20))
  Atrp_prev10_19 <- c(Atrp_prev10_19, Atrp_prev10_19_temp)

  Atrp_prev20_29_temp <- length(which(atrophyStatuses == 1 & ages >= 20 & ages < 30)) /  length(which(ages >= 20 & ages < 30))
  Atrp_prev20_29 <- c(Atrp_prev20_29, Atrp_prev20_29_temp)

  Atrp_prev30_49_temp <- length(which(atrophyStatuses == 1 & ages >= 30 & ages < 50)) /  length(which(ages >= 30 & ages < 50))
  Atrp_prev30_49 <- c(Atrp_prev30_49, Atrp_prev30_49_temp)

  Atrp_prev50_80_temp <- length(which(atrophyStatuses == 1 & ages >= 50 & ages <= 80)) /  length(which(ages >= 50 & ages < 80))
  Atrp_prev50_80 <- c(Atrp_prev50_80, Atrp_prev50_80_temp)

  # Hanging groin age age-prev #
  HG_prev0_1_temp <- length(which(hgStatuses == 1 & ages >= 0 & ages < 2)) /  length(which(ages >= 0 & ages < 2))# 0 - 1 age
  HG_prev0_1 <- c(HG_prev0_1, HG_prev0_1_temp)

  HG_prev2_4_temp <- length(which(hgStatuses == 1 & ages >= 2 & ages < 5)) /  length(which(ages >= 2 & ages < 5))# 2 - 4 age
  HG_prev2_4 <- c(HG_prev2_4, HG_prev2_4_temp)

  HG_prev5_9_temp <- length(which(hgStatuses == 1 & ages >= 5 & ages < 10)) /  length(which(ages >= 5 & ages < 10))
  HG_prev5_9 <- c(HG_prev5_9, HG_prev5_9_temp)

  HG_prev10_19_temp <- length(which(hgStatuses == 1 & ages >= 10 & ages < 20)) /  length(which(ages >= 10 & ages < 20))
  HG_prev10_19 <- c(HG_prev10_19, HG_prev10_19_temp)

  HG_prev20_29_temp <- length(which(hgStatuses == 1 & ages >= 20 & ages < 30)) /  length(which(ages >= 20 & ages < 30))
  HG_prev20_29 <- c(HG_prev20_29, HG_prev20_29_temp)

  HG_prev30_49_temp <- length(which(hgStatuses == 1 & ages >= 30 & ages < 50)) /  length(which(ages >= 30 & ages < 50))
  HG_prev30_49 <- c(HG_prev30_49, HG_prev30_49_temp)

  HG_prev50_80_temp <- length(which(hgStatuses == 1 & ages >= 50 & ages <= 80)) /  length(which(ages >= 50 & ages < 80))
  HG_prev50_80 <- c(HG_prev50_80, HG_prev50_80_temp)

  # depigmentation age age-prev #
  depigm_prev0_1_temp <- length(which(depigStatuses == 1 & ages >= 0 & ages < 2)) /  length(which(ages >= 0 & ages < 2))# 0 - 1 age
  depigm_prev0_1 <- c(depigm_prev0_1, depigm_prev0_1_temp)

  depigm_prev2_4_temp <- length(which(depigStatuses == 1 & ages >= 2 & ages < 5)) /  length(which(ages >= 2 & ages < 5))# 2 - 4 age
  depigm_prev2_4 <- c(depigm_prev2_4, depigm_prev2_4_temp)

  depigm_prev5_9_temp <- length(which(depigStatuses == 1 & ages >= 5 & ages < 10)) /  length(which(ages >= 5 & ages < 10))
  depigm_prev5_9 <- c(depigm_prev5_9, depigm_prev5_9_temp)

  depigm_prev10_19_temp <- length(which(depigStatuses == 1 & ages >= 10 & ages < 20)) /  length(which(ages >= 10 & ages < 20))
  depigm_prev10_19 <- c(depigm_prev10_19, depigm_prev10_19_temp)

  depigm_prev20_29_temp <- length(which(depigStatuses == 1 & ages >= 20 & ages < 30)) /  length(which(ages >= 20 & ages < 30))
  depigm_prev20_29 <- c(depigm_prev20_29, depigm_prev20_29_temp)

  depigm_prev30_49_temp <- length(which(depigStatuses == 1 & ages >= 30 & ages < 50)) /  length(which(ages >= 30 & ages < 50))
  depigm_prev30_49 <- c(depigm_prev30_49, depigm_prev30_49_temp)

  depigm_prev50_80_temp <- length(which(depigStatuses == 1 & ages >= 50 & ages <= 80)) /  length(which(ages >= 50 & ages < 80))
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
#' @param temp_mf vector of all mf per skin snip for each individual
#' @param morb_mat dataframe containing columns to determining conditions for testing morbidity (eye disease)
#' @param age_to_samp_vec_nonreversible 1 year increments of increasing age for when an individual should be tested (matching age), between 20 - 80 yrs
#'
#' @returns updated matrix with identifying individuals to test for morbidity
find_indiv_totest_func2 <- function(dat, mf.start, mf.end, morb_mat, age_to_samp_vec_nonreversible){

  # update age (and sex for newborns)
  morb_mat$Age <- dat[,2]
  morb_mat$Sex <- dat[,3]
  ages <- morb_mat$Age
  sexes <- morb_mat$Sex

  # lagged ages (age 2 year in future) #

  morb_mat$LaggedAges <- ages + 2 # current age + 2 years (in 2 yrs time)
  laggedAges <- morb_mat$LaggedAges
  morb_mat$LaggedAgeOver80 <- ifelse(laggedAges > 79.99999999, laggedAges - 80, laggedAges) # if between 78 - 80, will be > 80 yrs in 2 yrs time
                                                                                                          # therefore not alive, so newborn will be future age - 80 yrs (e.g., 81 - 80 = 1 yr old)
                                                                                                          # need to ensure individuals between 0 - 2 yrs with 0 blindness due to lag (below)
  # true number of mf per individual #
  mf_all <- rowSums(dat[, mf.start : mf.end]) # sum mf per individual across all 21 age classes
  morb_mat$TrueMFCount <- mf_all # true mf count
  trueMFCount <- morb_mat$TrueMFCount

  # # extract number of mf per skin snip per individual
  # morb.mat.tmp[,5] <- round(temp.mf[[2]]) # mf per skin snip for all individuals
  observedMFCounts <- morb_mat$ObservedMFCount

  # ======================== #
  #  1 ) age sampling        #

  # ======================================================== #
  # # new approach: determine if age to samp matches for all conditions (only want to sample once per year of age, not every time-step)
  morb_mat$AgeToSampleEyeDist <- ifelse((round(ages,6) %in% round(age_to_samp_vec_nonreversible,6)),1,0) # blindness
  ageToSampleEyeDist <- morb_mat$AgeToSampleEyeDist

  # ==================================#
  # 2) WHOM TO UNDERGO BERNOULI TRAIL #

  # selection based on true mf count
  blindnessStatuses <- morb_mat$BlindnessStatus
  morb_mat$ToTestBlindness <- ifelse(trueMFCount > 0 & ageToSampleEyeDist == 1 & blindnessStatuses == 0, 1, 0) # blindness (irreversible; only test once in age range)
  toTestBlindness <- morb_mat$ToTestBlindness

  # # selection based on observed mf count
  # toTestBlindness <- ifelse(observedMFCounts > 0 & ageToSampleEyeDist == 1 & blindnessStatuses == 0, 1, 0) # blindness (irreversible; only test once in age range)


 return(morb_mat)

}



#' @title
#' individuals designated to be tested undergo Bernoulli trial for realize eye disease status
#'
#' @description
#' each individual undergoes a Bernoulli trial (using probabilities based on mf counts) to ascertain new cases of eye disease
#'
#' @param morb_mat updated dataframe highlighting individuals to test
#' @param temp_mf vector of all mf per skin snip for each individual
#' @param blind_probs probabilities of blindness for a given mf count (based on equation from Little et al. 2004)
#'
#' @returns updated matrix with disease status updated
new_cases_morbidity_func2 <- function(morb_mat, temp.mf, blind.probs){

  # extract number of mf per skin snip per individual
  morb_mat$ObservedMFCount <- round(temp.mf[[2]]) + 1 # mf per skin snip for all individuals (+1 because of indexing so that when indexing probabilities goes from 1)
  observedMFCounts <- morb_mat$ObservedMFCount

  # extract probabilities (rates) to run Bernoulli trial for each condition
  morb_mat$BlindnessProb <- ifelse(observedMFCounts > 0, blind.probs[observedMFCounts], 0) # blindness rate/ prob ~ mf count
  blindnessProbs <- morb_mat$BlindnessProb
  # ======================= #
  # Undergo Bernouli trial  #

  toTestBlindness <- morb_mat$ToTestBlindness #Supposed to be morb.mat.tmp[,7]
  blindnessStatuses <- morb_mat$BlindnessStatus
  morb_mat$BlindnessStatus <- ifelse(toTestBlindness == 1, rbinom(sum(toTestBlindness), 1, blindnessProbs), blindnessStatuses) # blindness (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)

  # update blindness status based on age in 2 years (i.e., those 0-2 yrs in 2 yrs time will be 0 blindness) #

  morb_mat$BlindnessStatus2Yrs <- ifelse(morb_mat$LaggedAgeOver80 < 2, 0, blindnessStatuses) # set to 0 for 0-2 yrs in 2 yrs time, or as in col 9

  return(morb_mat)



}

#' @title
#' eye disease prevalence
#'
#' @description
#' calculates disease prevalence for each disease state (including age-stratified disease prevalence's)
#'
#' @param N human pop size
#' @param morb_mat updated dataframe highlighting individuals to test
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
eye.disease.prev.func <- function(N, morb_mat,
                                  blind_prev, visual_imp_prev,
                                  blind_prev0_1, blind_prev2_4, blind_prev5_9, blind_prev10_19,
                                  blind_prev20_29, blind_prev30_49, blind_prev50_80,
                                  visual_imp_prev0_1, visual_imp_prev2_4, visual_imp_prev5_9,
                                  visual_imp_prev10_19, visual_imp_prev20_29, visual_imp_prev30_49,
                                  visual_imp_prev50_80)
{

  # ========================================================================================================= #
  # Approach: based on age in 2 years (col 11 of morb.mat.tmp) and updated blindness status in 2 yrs (col 12) #

  blind_prev_temp <- length(which(morb_mat$BlindnessStatus2Yrs == 1 & morb_mat$LaggedAgeOver80 >= 5)) /  length(which(morb_mat$LaggedAgeOver80 >= 5)) # prev in > 5yrs
  visual_imp_prev_temp <- blind_prev_temp * 1.78

  # update prevalence vectors
  blind_prev <- c(blind_prev, blind_prev_temp)
  visual_imp_prev <- c(visual_imp_prev, visual_imp_prev_temp)

  # blind age age-prev #
  blind_prev0_1_temp <- length(which(morb_mat$BlindnessStatus2Yrs == 1 & morb_mat$LaggedAgeOver80 >= 0 & morb_mat$LaggedAgeOver80 < 2)) /  length(which(morb_mat$LaggedAgeOver80 >= 0 & morb_mat$LaggedAgeOver80 < 2))# 0 - 1 age
  blind_prev0_1 <- c(blind_prev0_1, blind_prev0_1_temp)

  blind_prev2_4_temp <- length(which(morb_mat$BlindnessStatus2Yrs == 1 & morb_mat$LaggedAgeOver80 >= 2 & morb_mat$LaggedAgeOver80 < 5)) /  length(which(morb_mat$LaggedAgeOver80 >= 2 & morb_mat$LaggedAgeOver80 < 5))# 2 - 4 age
  blind_prev2_4 <- c(blind_prev2_4, blind_prev2_4_temp)

  blind_prev5_9_temp <- length(which(morb_mat$BlindnessStatus2Yrs == 1 & morb_mat$LaggedAgeOver80 >= 5 & morb_mat$LaggedAgeOver80 < 10)) /  length(which(morb_mat$LaggedAgeOver80 >= 5 & morb_mat$LaggedAgeOver80 < 10))
  blind_prev5_9 <- c(blind_prev5_9, blind_prev5_9_temp)

  blind_prev10_19_temp <- length(which(morb_mat$BlindnessStatus2Yrs == 1 & morb_mat$LaggedAgeOver80 >= 10 & morb_mat$LaggedAgeOver80 < 20)) /  length(which(morb_mat$LaggedAgeOver80 >= 10 & morb_mat$LaggedAgeOver80 < 20))
  blind_prev10_19 <- c(blind_prev10_19, blind_prev10_19_temp)

  blind_prev20_29_temp <- length(which(morb_mat$BlindnessStatus2Yrs == 1 & morb_mat$LaggedAgeOver80 >= 20 & morb_mat$LaggedAgeOver80 < 30)) /  length(which(morb_mat$LaggedAgeOver80 >= 20 & morb_mat$LaggedAgeOver80 < 30))
  blind_prev20_29 <- c(blind_prev20_29, blind_prev20_29_temp)

  blind_prev30_49_temp <- length(which(morb_mat$BlindnessStatus2Yrs == 1 & morb_mat$LaggedAgeOver80 >= 30 & morb_mat$LaggedAgeOver80 < 50)) /  length(which(morb_mat$LaggedAgeOver80 >= 30 & morb_mat$LaggedAgeOver80 < 50))
  blind_prev30_49 <- c(blind_prev30_49, blind_prev30_49_temp)

  blind_prev50_80_temp <- length(which(morb_mat$BlindnessStatus2Yrs == 1 & morb_mat$LaggedAgeOver80 >= 50 & morb_mat$LaggedAgeOver80 <= 80)) /  length(which(morb_mat$LaggedAgeOver80 >= 50 & morb_mat$LaggedAgeOver80 < 80))
  blind_prev50_80 <- c(blind_prev50_80, blind_prev50_80_temp)

  # visual impairement age age-prev #
  visual_imp_prev0_1_temp <- blind_prev0_1_temp * 1.78 # 0 - 1 age
  visual_imp_prev0_1 <- c(visual_imp_prev0_1, visual_imp_prev0_1_temp)

  visual_imp_prev2_4_temp <- blind_prev2_4_temp * 1.78 # 2 - 4 age
  visual_imp_prev2_4 <- c(visual_imp_prev2_4, visual_imp_prev2_4_temp)

  visual_imp_prev5_9_temp <- blind_prev5_9_temp * 1.78
  visual_imp_prev5_9 <- c(visual_imp_prev5_9, visual_imp_prev5_9_temp)

  visual_imp_prev10_19_temp <- blind_prev10_19_temp * 1.78
  visual_imp_prev10_19 <- c(visual_imp_prev10_19, visual_imp_prev10_19_temp)

  visual_imp_prev20_29_temp <- blind_prev20_29_temp * 1.78
  visual_imp_prev20_29 <- c(visual_imp_prev20_29, visual_imp_prev20_29_temp)

  visual_imp_prev30_49_temp <- blind_prev30_49_temp * 1.78
  visual_imp_prev30_49 <- c(visual_imp_prev30_49, visual_imp_prev30_49_temp)

  visual_imp_prev50_80_temp <- blind_prev50_80_temp * 1.78
  visual_imp_prev50_80 <- c(visual_imp_prev50_80, visual_imp_prev50_80_temp)


  return(list(blind_prev, visual_imp_prev,
              blind_prev0_1, blind_prev2_4, blind_prev5_9, blind_prev10_19, blind_prev20_29,
              blind_prev30_49, blind_prev50_80,
              visual_imp_prev0_1, visual_imp_prev2_4, visual_imp_prev5_9, visual_imp_prev10_19,
              visual_imp_prev20_29, visual_imp_prev30_49, visual_imp_prev50_80))

}
