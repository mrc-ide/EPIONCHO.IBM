# ========================================================================================================================== #
#           External testing of morbidity model (mf count ~ onset approach)  - 01 / 06/ 2023                                 #
# ========================================================================================================================== #
DT <- 1/366

morbidity_matrix_out_eq <- output_equilibrium_morbidity$morbidity.matrix # after running 130 yrs

all.mats.temp_out <- output_equilibrium_morbidity$all_equilibrium_outputs[[1]] # to get mf out

# specify how matrices called in functions ~
dat <- all.mats.temp_out
morb.mat.tmp <- morbidity_matrix_out_eq

mf.start <- 7
mf.end <- 6 + 21

# update age
dat[,2] <- dat[,2] + 1/366
dat[,2] <- dat[,2] + 3


# getting a vector of years
# all.morb.temp[,6] <- sample(seq(5, 80, DT), size = N, replace = TRUE)
#
# # age_to_samp_vec <- sample(seq(5, 80, DT), size = N, replace = TRUE)
#
# age_seq <- round(seq(5, 80, DT), digits = 2)
# morb.mat.tmp[,1] <- round(morb.mat.tmp[,1], digits = 2)
#
# match(morb.mat.tmp[,1],age_seq)
#
# age_to_samp_vec <- seq(5.01,79.01,1)
# match(morb.mat.tmp[,1],age_to_samp_vec)
#
# test_vec <- ifelse(match(morb.mat.tmp[,1],age_to_samp_vec),1,0)

# ============================================================================================================== #
# 1)                                           find_indiv_totest_func                                            #
# ============================================================================================================== #

# update age (and sex for newborns)
morb.mat.tmp[,1] <- dat[,2] + 3
morb.mat.tmp[,2] <- dat[,3]

# true number of mf per individual #
mf_all <- rowSums(dat[, mf.start : mf.end]) # sum mf per individual across all 21 age classes
mf_all_dummy <- mf_all + 1
morb.mat.tmp[,4] <- mf_all_dummy # true mf count

# testing cols to add
mat.to.add = matrix(NA, nrow = nrow(morb.mat.tmp), ncol = 2)
morb.mat.tmp = cbind(morb.mat.tmp, mat.to.add)

morb.mat.tmp[,54] <- dat[,2]
morb.mat.tmp[,55] <- dat[,2] + 3


# update age (and sex for newborns)
morb.mat.tmp[,1] <- dat[,2]
morb.mat.tmp[,2] <- dat[,3]

# true number of mf per individual #
mf_all <- rowSums(dat[, mf.start : mf.end]) # sum mf per individual across all 21 age classes
morb.mat.tmp[,4] <- mf_all # true mf count

age_to_samp_vec_reversible <- seq(5+1/366, 79+1/366, 1) # between 5 and 80, sample once year year of age
age_to_samp_vec_nonreversible <- seq(20+1/366, 79+1/366, 1) # between 5 and 80, sample once year year of age

# age_to_samp_vec_reversible <- seq(5.01,79.01,1) # between 5 and 80, sample once year year of age
# age_to_samp_vec_nonreversible <- seq(20.01,79.01,1) # between 5 and 80, sample once year year of age

all_ages_vec <- seq(5+1/366, 79+1/366, 1/366) # between 5 and 80, sample once year year of age

match(all_ages_vec, age_to_samp_vec_reversible)
match(round(all_ages_vec,2), age_to_samp_vec_reversible)
match(round(all_ages_vec,6), round(age_to_samp_vec_reversible,6))

morb.mat.tmp[,1] <- morb.mat.tmp[,1] + DT

morb.mat.tmp[,1] <- round(morb.mat.tmp[,1],6) # needs to be rounded to 5 decimal places

any(morb.mat.tmp[,1] %in% round(age_to_samp_vec_reversible,6))
morb.mat.tmp[,1] %in% round(age_to_samp_vec_reversible,6)

# ======================== #
#  1 ) age sampling        #

# ======================================================#
# current (only used for non-reversible to sample once) #
# determine whether age_to_samp matches age of individual (rounded) for each condition where only test once in age range
# # morb.mat.tmp[,14] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,6]), 1, 0) # severe itch
# # morb.mat.tmp[,15] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,10]), 1, 0) # APOD (test again single RSD age_to_sample e.g., col 10)
# # morb.mat.tmp[,16] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,10]), 1, 0) # CPOD (test again single RSD age_to_sample e.g., col 10)
# # morb.mat.tmp[,17] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,10]), 1, 0) # LOD (test again single RSD age_to_sample e.g., col 10)
# # morb.mat.tmp[,18] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,10]), 1, 0) # RSD (test again single RSD age_to_sample e.g., col 10)
# morb.mat.tmp[,19] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,11]), 1, 0) # atrophy
# morb.mat.tmp[,20] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,12]), 1, 0) # hanging groin
# morb.mat.tmp[,21] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,13]), 1, 0) # depigmentation

round(morb.mat.tmp[,1], 3)
# match(round(morb.mat.tmp[,1],3), age_to_samp_vec_reversible)

# ======================================================== #
# # new approach: determine if age to samp matches for all conditions (only want to sample once per year of age, not every time-step)
morb.mat.tmp[,14] <- ifelse((round(morb.mat.tmp[,1],2) %in% age_to_samp_vec_reversible),1,0) # severe itch
# morb.mat.tmp[,15] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,10]), 1, 0) # APOD (test again single RSD age_to_sample e.g., col 10)
# morb.mat.tmp[,16] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,10]), 1, 0) # CPOD (test again single RSD age_to_sample e.g., col 10)
# morb.mat.tmp[,17] <- ifelse(round(morb.mat.tmp[,1]) == round(morb.mat.tmp[,10]), 1, 0) # LOD (test again single RSD age_to_sample e.g., col 10)
morb.mat.tmp[,18] <- ifelse((round(morb.mat.tmp[,1],2) %in% age_to_samp_vec_reversible),1,0) # RSD (test again single RSD age_to_sample e.g., col 10)
morb.mat.tmp[,19] <- ifelse((round(morb.mat.tmp[,1],2) %in% age_to_samp_vec_nonreversible),1,0) # atrophy
morb.mat.tmp[,20] <- ifelse((round(morb.mat.tmp[,1],2) %in% age_to_samp_vec_nonreversible),1,0) # hanging groin
morb.mat.tmp[,21] <- ifelse((round(morb.mat.tmp[,1],2) %in% age_to_samp_vec_nonreversible),1,0) # depigmentation

morb.mat.tmp[,14] <- ifelse((morb.mat.tmp[,1] %in% age_to_samp_vec_reversible),1,0) # severe itch
morb.mat.tmp[,18] <- ifelse((morb.mat.tmp[,1] %in% age_to_samp_vec_reversible),1,0) # RSD (test again single RSD age_to_sample e.g., col 10)
morb.mat.tmp[,19] <- ifelse((morb.mat.tmp[,1] %in% age_to_samp_vec_nonreversible),1,0) # atrophy
morb.mat.tmp[,20] <- ifelse((morb.mat.tmp[,1] %in% age_to_samp_vec_nonreversible),1,0) # hanging groin
morb.mat.tmp[,21] <- ifelse((morb.mat.tmp[,1] %in% age_to_samp_vec_nonreversible),1,0) # depigmentation


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
morb.mat.tmp[,27] <- ifelse(morb.mat.tmp[,35] == 1, 1, morb.mat.tmp[,27]) # atrophy
morb.mat.tmp[,28] <- ifelse(morb.mat.tmp[,36] == 1, 1, morb.mat.tmp[,28]) # hanging groin
morb.mat.tmp[,29] <- ifelse(morb.mat.tmp[,37] == 1, 1, morb.mat.tmp[,29]) # depigmentation

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
# Version 1 current working #

# # determine whether individual will undergo Bernoulli trial in this time-step for condition (first test)
# # test a) reversible conditions daily and b) non-reversible only once in life-time (at age_to_sample & only not previously tested)
# morb.mat.tmp[,30] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # severe itch (reversible = tested each dt regardless of presence of condition)
# morb.mat.tmp[,31] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # APOD (reversible = tested each dt regardless of presence of condition)
# morb.mat.tmp[,32] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # CPOD (reversible = tested each dt regardless of presence of condition)
# morb.mat.tmp[,33] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # LOD (reversible = tested each dt regardless of presence of condition)
# morb.mat.tmp[,34] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # RSD (reversible = tested each dt regardless of presence of condition)
# morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,19] == 1 & morb.mat.tmp[,27] == 0, 1, 0) # atrophy (irreversible; only test once in age range)
# morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,20] == 1 & morb.mat.tmp[,28] == 0, 1, 0) # HG (irreversible; only test once in age range)
# morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,21] == 1 & morb.mat.tmp[,29] == 0, 1, 0) # depigm (irreversible; only test once in age range)

# ============ #
# new approach #

# test #
# morb.mat.tmp[4,4] = 0
# morb.mat.tmp[7,4] = 0

# # first determine if true mf count for an individual with a reversible condition is 0, then then become sequelea negative
morb.mat.tmp[,46] <- ifelse(morb.mat.tmp[,4] == 0, 0, morb.mat.tmp[,46]) # if true mf count is 0, then sequela is 0, else leave
morb.mat.tmp[,47] <- ifelse(morb.mat.tmp[,4] == 0, 0, morb.mat.tmp[,47]) # if true mf count is 0, then sequela is 0, else leave
morb.mat.tmp[,48] <- ifelse(morb.mat.tmp[,4] == 0, 0, morb.mat.tmp[,48]) # if true mf count is 0, then sequela is 0, else leave
morb.mat.tmp[,49] <- ifelse(morb.mat.tmp[,4] == 0, 0, morb.mat.tmp[,49]) # if true mf count is 0, then sequela is 0, else leave
morb.mat.tmp[,50] <- ifelse(morb.mat.tmp[,4] == 0, 0, morb.mat.tmp[,50]) # if true mf count is 0, then sequela is 0, else leave

# determine whether individual will undergo Bernoulli trial in this time-step for condition
# reversible conditions: only tested if sequela 0 once per year of age (cols 14 & 18 if age_to_samp_vec matches age)
# update: now non-reversible conditions (only those that do not have the condition) only tested once per year of age (cols 19 - 21)
# age-truncated for reversible conditions (from 20 yrs) due to age- disease prev profiles (Murdoch et al. 2017)
morb.mat.tmp[,30] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,46] == 0 & morb.mat.tmp[,14] == 1, 1, 0) # severe itch (reversible)
morb.mat.tmp[,34] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,50] == 0 & morb.mat.tmp[,18] == 1, 1, 0) # RSD (reversible)
morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,19] == 1 & morb.mat.tmp[,51] == 0, 1, 0) # atrophy (irreversible)
morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,20] == 1 & morb.mat.tmp[,52] == 0, 1, 0) # HG (irreversible)
morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,21] == 1 & morb.mat.tmp[,53] == 0, 1, 0) # depigm (irreversible)



# ======================================================================================================================== #
#  2)                                               new_cases_morbidity_func                                               #

# extract number of mf per skin snip per individual
morb.mat.tmp[,5] <- round(temp_mf[[2]]) + 1 # mf per skin snip for all individuals (+1 because of indexing so that when indexing probabilities goes from 1)

# extract probabilities (rates) to run Bernoulli trial for each condition
morb.mat.tmp[,38] <- ifelse(morb.mat.tmp[,5] > 0, SI_probs[morb.mat.tmp[,5]], 0) # severe itch rates
#morb.mat.tmp[,39] <- ifelse(morb.mat.tmp[,5] > 0, APOD_probs[morb.mat.tmp[,5]], 0) # APOD rates
#morb.mat.tmp[,40] <- ifelse(morb.mat.tmp[,5] > 0, CPOD_probs[morb.mat.tmp[,5]], 0) # CPOD rates
#morb.mat.tmp[,41] <- ifelse(morb.mat.tmp[,5] > 0, LOD_probs[morb.mat.tmp[,5]], 0) # LOD rates
morb.mat.tmp[,42] <- ifelse(morb.mat.tmp[,5] > 0, RSD_probs[morb.mat.tmp[,5]], 0) # RSD rates
morb.mat.tmp[,43] <- ifelse(morb.mat.tmp[,5] > 0, Atrp_probs[morb.mat.tmp[,5]], 0) # atrophy rates
morb.mat.tmp[,44] <- ifelse(morb.mat.tmp[,5] > 0, Hg_probs[morb.mat.tmp[,5]], 0) # hanging groin rates
morb.mat.tmp[,45] <- ifelse(morb.mat.tmp[,5] > 0, Depigm_probs[morb.mat.tmp[,5]], 0) # depigmentation groin rates

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

#Bernoulli trial for those to be tested for each condition
morb.mat.tmp[,46] <- ifelse(morb.mat.tmp[,30] == 1, rbinom(sum(morb.mat.tmp[,30]), 1, morb.mat.tmp[,38]), 0) # severe itch (revert to 0 every time step if not realised)
morb.mat.tmp[,50] <- ifelse(morb.mat.tmp[,34] == 1, rbinom(sum(morb.mat.tmp[,34]), 1, morb.mat.tmp[,42]), 0) # RSD (revert to 0 every time step if not realised)
morb.mat.tmp[,51] <- ifelse(morb.mat.tmp[,35] == 1, rbinom(sum(morb.mat.tmp[,35]), 1, morb.mat.tmp[,43]), morb.mat.tmp[,51]) # atrophy (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
morb.mat.tmp[,52] <- ifelse(morb.mat.tmp[,36] == 1, rbinom(sum(morb.mat.tmp[,36]), 1, morb.mat.tmp[,44]), morb.mat.tmp[,52]) # HG (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
morb.mat.tmp[,53] <- ifelse(morb.mat.tmp[,37] == 1, rbinom(sum(morb.mat.tmp[,37]), 1, morb.mat.tmp[,45]), morb.mat.tmp[,53]) # depigmentation (non-reversible: only testing 0's - stay as previous if tested in this time-step i.e, currently diseases (1) stay as 1)
