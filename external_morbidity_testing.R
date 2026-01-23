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

# ======================================== #
# 1)    find_indiv_totest_func             #

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

# determine whether age_to_samp matches age of individual (rounded) for each condition where only test once in age range
morb.mat.tmp[,14] <- ifelse(round(morb.mat.tmp[,2]) == round(morb.mat.tmp[,6]), 1, 0) # severe itch
morb.mat.tmp[,15] <- ifelse(round(morb.mat.tmp[,2]) == round(morb.mat.tmp[,10]), 1, 0) # APOD (test again single RSD age_to_sample e.g., col 10)
morb.mat.tmp[,16] <- ifelse(round(morb.mat.tmp[,2]) == round(morb.mat.tmp[,10]), 1, 0) # CPOD (test again single RSD age_to_sample e.g., col 10)
morb.mat.tmp[,17] <- ifelse(round(morb.mat.tmp[,2]) == round(morb.mat.tmp[,10]), 1, 0) # LOD (test again single RSD age_to_sample e.g., col 10)
morb.mat.tmp[,18] <- ifelse(round(morb.mat.tmp[,2]) == round(morb.mat.tmp[,10]), 1, 0) # RSD (test again single RSD age_to_sample e.g., col 10)
morb.mat.tmp[,19] <- ifelse(round(morb.mat.tmp[,2]) == round(morb.mat.tmp[,11]), 1, 0) # atrophy
morb.mat.tmp[,20] <- ifelse(round(morb.mat.tmp[,2]) == round(morb.mat.tmp[,12]), 1, 0) # hanging groin
morb.mat.tmp[,21] <- ifelse(round(morb.mat.tmp[,2]) == round(morb.mat.tmp[,13]), 1, 0) # depigmentation

# determine whether individual has previously tested for condition (for selecting those where only want to test once for a condition)
morb.mat.tmp[,22] <- ifelse(morb.mat.tmp[,30] == 1, 1, morb.mat.tmp[,22]) # severe itch
morb.mat.tmp[,23] <- ifelse(morb.mat.tmp[,31] == 1, 1, morb.mat.tmp[,23]) # APOD
morb.mat.tmp[,24] <- ifelse(morb.mat.tmp[,32] == 1, 1, morb.mat.tmp[,24]) # CPOD
morb.mat.tmp[,25] <- ifelse(morb.mat.tmp[,33] == 1, 1, morb.mat.tmp[,25]) # LOD
morb.mat.tmp[,26] <- ifelse(morb.mat.tmp[,34] == 1, 1, morb.mat.tmp[,26]) # RSD
morb.mat.tmp[,27] <- ifelse(morb.mat.tmp[,35] == 1, 1, morb.mat.tmp[,27]) # atrophy
morb.mat.tmp[,28] <- ifelse(morb.mat.tmp[,36] == 1, 1, morb.mat.tmp[,28]) # hanging groin
morb.mat.tmp[,29] <- ifelse(morb.mat.tmp[,37] == 1, 1, morb.mat.tmp[,29]) # depigmentation

# determine whether individual will undergo Bernoulli trial in this time-step for condition (first test)
# morb.mat.tmp[,30] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # severe itch (reversible = tested each dt regardless of presence of condition)
# morb.mat.tmp[,31] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # APOD (reversible = tested each dt regardless of presence of condition)
# morb.mat.tmp[,32] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # CPOD (reversible = tested each dt regardless of presence of condition)
# morb.mat.tmp[,33] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # LOD (reversible = tested each dt regardless of presence of condition)
# morb.mat.tmp[,34] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # RSD (reversible = tested each dt regardless of presence of condition)
# morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,19] == 1 & morb.mat.tmp[,27] == 0, 1, 0) # atrophy (irreversible; only test once in age range)
# morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,20] == 1 & morb.mat.tmp[,28] == 0, 1, 0) # HG (irreversible; only test once in age range)
# morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,21] == 1 & morb.mat.tmp[,29] == 0, 1, 0) # depigm (irreversible; only test once in age range)


# determine whether individual will undergo Bernoulli trial in this time-step for condition
# update: now non-reversible conditions tested daily (only those that do not have the condition)
# age-truncated for reversible conditions (from 20 yrs) due to age- disease prev profiles (Murdoch et al. 2017)
morb.mat.tmp[,30] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # severe itch (reversible = tested each dt regardless of presence of condition)
morb.mat.tmp[,31] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # APOD (reversible = tested each dt regardless of presence of condition)
morb.mat.tmp[,32] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # CPOD (reversible = tested each dt regardless of presence of condition)
morb.mat.tmp[,33] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # LOD (reversible = tested each dt regardless of presence of condition)
morb.mat.tmp[,34] <- ifelse(morb.mat.tmp[,4] > 0, 1, 0) # RSD (reversible = tested each dt regardless of presence of condition)
# morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,51] == 0, 1, 0) # atrophy (irreversible; only test once in age range)
# morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,52] == 0, 1, 0) # HG (irreversible; only test once in age range)
# morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,53] == 0, 1, 0) # depigm (irreversible; only test once in age range)
morb.mat.tmp[,35] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,51] == 0, 1, 0) # atrophy (irreversible; only test once in age range)
morb.mat.tmp[,36] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,52] == 0, 1, 0) # HG (irreversible; only test once in age range)
morb.mat.tmp[,37] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,53] == 0, 1, 0) # depigm (irreversible; only test once in age range)




rbinom(sum(morb.mat.tmp[,35]), 1, morb.mat.tmp[,43])

0.04034684

rbinom(sum(morb.mat.tmp[,35]), 1, 0.04034684)
rbinom(1, 1, 0.04034684)

rbinom(1, 1, 0.14456026)


ggplot(data = morb.mat.tmp)+



# ======================================== #
# 2)    find_indiv_totest_func             #




test_list <-
  list("SI_prev", "RSD_prev", "Atrp_prev", "HG_prev", "depigm_prev",
     "SI_prev0_1", "SI_prev2_4", "SI_prev5_9", "SI_prev10_19", "SI_prev20_29", "SI_prev30_49", "SI_prev50_80",
     "RSD_prev0_1", "RSD_prev2_4", "RSD_prev5_9", "RSD_prev10_19", "RSD_prev20_29", "RSD_prev30_49", "RSD_prev50_80",
     "Atrp_prev0_1", "Atrp_prev2_4", "Atrp_prev5_9", "Atrp_prev10_19", "Atrp_prev20_29", "Atrp_prev30_49", "Atrp_prev50_80",
     "HG_prev0_1", "HG_prev2_4", "HG_prev5_9", "HG_prev10_19", "HG_prev20_29", "HG_prev30_49", "HG_prev50_80",
     "depigm_prev0_1", "depigm_prev2_4", "depigm_prev5_9", "depigm_prev10_19", "depigm_prev20_29", "depigm_prev30_49", "depigm_prev50_80")

list( SI_prev0_1 = morbidity_prev_out[[6]], SI_prev2_4 = morbidity_prev_out[[7]], SI_prev5_9 = morbidity_prev_out[[8]],
      SI_prev10_19 = morbidity_prev_out[[9]], SI_prev20_29 = morbidity_prev_out[[10]], SI_prev30_49 = morbidity_prev_out[[11]],
      SI_prev50_80 =  morbidity_prev_out[[12]],
      RSD_prev0_1 =  morbidity_prev_out[[13]], RSD_prev2_4 =  morbidity_prev_out[[14]], RSD_prev5_9 =  morbidity_prev_out[[15]],
      RSD_prev10_19 =  morbidity_prev_out[[16]], RSD_prev20_29 =  morbidity_prev_out[[17]], RSD_prev30_49 =  morbidity_prev_out[[18]],
      RSD_prev50_80 =  morbidity_prev_out[[17]],
      Atrp_prev0_1 =  morbidity_prev_out[[18]], Atrp_prev2_4 =  morbidity_prev_out[[19]], Atrp_prev5_9 =  morbidity_prev_out[[20]],
      Atrp_prev10_19 =  morbidity_prev_out[[21]], Atrp_prev20_29 =  morbidity_prev_out[[22]], Atrp_prev30_49 =  morbidity_prev_out[[23]],
      Atrp_prev50_80 =  morbidity_prev_out[[24]],
      HG_prev0_1 =  morbidity_prev_out[[25]], HG_prev2_4 =  morbidity_prev_out[[26]], HG_prev5_9 =  morbidity_prev_out[[27]],
      HG_prev10_19 =  morbidity_prev_out[[28]], HG_prev20_29 =  morbidity_prev_out[[29]], HG_prev30_49 =  morbidity_prev_out[[30]],
      HG_prev50_80 =  morbidity_prev_out[[31]],
      depigm_prev0_1 =  morbidity_prev_out[[32]], depigm_prev2_4 =  morbidity_prev_out[[33]], depigm_prev5_9 =  morbidity_prev_out[[34]],
      depigm_prev10_19 =  morbidity_prev_out[[35]], depigm_prev20_29 =  morbidity_prev_out[[36]], depigm_prev30_49 =  morbidity_prev_out[[37]],
      depigm_prev50_80 =  morbidity_prev_out[[38]]),
