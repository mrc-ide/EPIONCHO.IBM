#' @title
#' mf dying calculation
#' @description
#' calculates number of mf that die per individual (summed over mf age classes)
#' @param mat.in mf.dying matrix in (columns for current number of mf per age class in each individual) & columns to calculate number of mf that die in each age class
#' @param num.mf.comps number of mf columns (e.g., 21 age classes)
#' @param N human population size
#' @param mort.rates.mf mortality rates for mf by age class according to weibull distribution (note DT = per month e.g., 30/366 rather than 1/366)
#' @param mf.start columns where mf age classes start
#' @param mf.end column where mf age classes end
#' @param all.mats.temp main matrix tracking mf numbers
#'
#' @returns list including 1) number of mf dying per individual (summed over mf age classes), 2) number of mf alive, 3) final matrix
number_mfdie_perindiv_func <- function(mat.in, num.mf.comps, N, mort.rates.mf,
                                       mf.start, mf.end, all.mats.temp){

  mat.in[,c(1:21)] <- all.mats.temp[,c(mf.start:mf.end)] # make first 21 current equilibrium number of mf per individual x age class

  for (k in 1 : num.mf.comps) {

    cl <- k # define col to select for mf age class

    cur.mf <- mat.in[, cl] # select mf age class column from matrix

    cl_to_upd <- cl + num.mf.comps # column to update number of mf dying in matrix

    mat.in[,cl_to_upd] <- mat.in[,cl] * mort.rates.mf[cl] # calculate number of mf dying in that age class
                                                                                    # based on mf.mort.rate for that age class

  }

  mf.dying.mat <- mat.in[,c(22:42)]

  indiv.mf.dying <- rowSums(mf.dying.mat)


  return(list(indiv.mf.dying, mat.in))
}


#' @title
#' mf dying per month calculation
#' @description
#' calculate accrued number of mf dying per individual every 30 days
#' @param mf.die.mat delay matrix tracking mf dying in each individual per day over a rolling 30 day period (updated each time step)
#' @param inds.mfdie.mat for moving columns along with time (each time-step)
#' @param mf.die.day.in calculated number of mf that have died in current day per individual
#'
#' @returns list including 1) vector of mf that have died per individual over past 30 days, 2) updated delay matrix
update_mfdying_month_mat_func <- function(mf.die.mat, inds.mfdie.mat, mf.die.day.in){

  mf.die.mat[,inds.mfdie.mat] <- mf.die.mat[,(inds.mfdie.mat-1)] # move mf dying along one col

  mf.die.mat[,1] <- mf.die.day.in # update first col with number of mf dying on that day

  mf.die.month <- c()

  mf.die.month <- rowSums(mf.die.mat) # calculate number of mf dying in last 30 days

  return(list(mf.die.month, mf.die.mat))

  }


#' @title
#' disease status update matrix
#' @description
#' calculates accrued tissue damage per month due to mf dying in each individual and whether disease threshold(s) passed
#' @param morb.mat matrix with cols representing components to calculate presence of disease for each individual
#' @param all.mats.temp main matrix tracking mf numbers
#' @param mf.die.lastmnth number of mf dying per individual in last 30 days
#' @param SI_regression_rate regression rate per month for severe itch (default value is 0.015)
#' @param RSD_regression_rate regression rate per month for RSD (default value is 0.030)
#' @param SI_threshold disease threshold required to develop severe itch (default value is 0.255 * 1000 ?)
#' @param RSD_threshold disease threshold required to develop reactive skin disease (default value is 0.210 * 1000 ?)
#' @param atrophy_threshold disease threshold required to develop atrophy (default value is 11.3 * 1000 ?)
#' @param HG_threshold disease threshold required to develop hanging groin (default value is 21.4 * 1000 ?)
#' @param mdpigm_threshold disease threshold required to develop mild depigmentation (default value is 2.35 * 1000 ?)
#' @param sdpigm_threshold disease threshold required to develop severe depigmentation (default value is 4.3 * 1000 ?)
#'
#' @returns matrix
disease_update_function <- function(morb.mat, all.mats.temp, mf.die.lastmnth,
                                    SI_regression_rate, RSD_regression_rate, SI_threshold, RSD_threshold,
                                    atrophy_threshold, HG_threshold, mdpigm_threshold, sdepigm_threshold){

  # update age / sex cols (if need to cal prev by age and sex)
  morb.mat[,1] <- all.mats.temp[,2]
  morb.mat[,2] <- all.mats.temp[,3]

  # update column 9 of morbidity matrix with total mf dying per individual in time-step
  morb.mat[,8] <- mf.die.lastmnth

  # calculate (new) change in amount of tissue activation/damage (change in Dix_T) for each condition #
  # this is calculated by: individual.susceptibility.condition * number.mf.dying - regression.rate * current tissue damage
  morb.mat[,10] <- morb.mat[,3] * morb.mat[,8] - (SI_regression_rate) * morb.mat[,9] # severe itch
  morb.mat[,10] <- ifelse(morb.mat[,10] < 0, 0, morb.mat[,10]) # if any values negative, change to 0 (as cant have negative tissue dmg)
  morb.mat[,12] <- morb.mat[,4] * morb.mat[,8] - (RSD_regression_rate) * morb.mat[,11] # reactive skin disease
  morb.mat[,12] <- ifelse(morb.mat[,12] < 0, 0, morb.mat[,12]) # if any values negative, change to 0 (as cant have negative tissue dmg)

  morb.mat[,14] <- morb.mat[,5] * morb.mat[,8] - 0 * morb.mat[,13] # atrophy (non-reversible)
  morb.mat[,16] <- morb.mat[,6] * morb.mat[,8] - 0 * morb.mat[,15] # hanging groin (non-reversible)
  morb.mat[,18] <- morb.mat[,7] * morb.mat[,8] - 0 * morb.mat[,17] # All depigmentation (non-reversible)

  # want the current amount of tissue damage to be updated with change in tissue damage col
  morb.mat[,9] <- morb.mat[,10] # severe itch
  morb.mat[,11] <- morb.mat[,12] # reactive skin disease
  morb.mat[,13] <- morb.mat[,14] # atrophy (non-reversible)
  morb.mat[,15] <- morb.mat[,16] # hanging groin (non-reversible)
  morb.mat[,17] <- morb.mat[,18] # All depigmentation (non-reversible)

  # check if threshold exceeded and disease condition present (thresholds are on original month scale)
  morb.mat[,19] <- ifelse(morb.mat[,10] > SI_threshold, 1, 0) # severe itch (threshold = 0.255)
  morb.mat[,20] <- ifelse(morb.mat[,12] > RSD_threshold, 1, 0) # reactive skin disease (threshold = 0.255)
  morb.mat[,21] <- ifelse(morb.mat[,21] == 0 & morb.mat[,14] > atrophy_threshold, 1, morb.mat[,21]) # atrophy (threshold = 11.3) - need to only test 0's as 1's permanent
  morb.mat[,22] <- ifelse(morb.mat[,22] == 0 & morb.mat[,16] > HG_threshold, 1, morb.mat[,22]) # hanging groin (threshold = 21.4) - need to only test 0's as 1's permanent

  morb.mat[,23] <- ifelse(morb.mat[,23] == 0 & morb.mat[,18] > mdpigm_threshold, 1, morb.mat[,23]) # temporary mild depigmentation (threshold = 2.35) - need to only test 0's as 1's permanent
  morb.mat[,25] <- ifelse(morb.mat[,25] == 0 & morb.mat[,18] > sdepigm_threshold, 1, morb.mat[,25]) # severe depigmentation (threshold = 4.3) - need to only test 0's as 1's permanent
  morb.mat[,24] <- ifelse(morb.mat[,25] == 1, 0, morb.mat[,23]) # update mild depigm i.e., with 0 if severe depigm is 1 or leave (this used for prev calculation for mild depigmentation)

return(morb.mat)

}

#' @title
#' prevalence of OSD sequelae
#' @description
#' calculates OSD sequelae prevalence both overall across all age groups and within specific age groups (to test vs. data)
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
#' @returns list including 1)

OSD_prevalence_function <- function(morb.mat.tmp, N, SI_prev, RSD_prev, Atrp_prev, HG_prev, MDp_prev, SDp_prev,
                                    SI_prev0_1, SI_prev2_4, SI_prev5_9, SI_prev10_19, SI_prev20_29, SI_prev30_49, SI_prev50_80,
                                    RSD_prev0_1, RSD_prev2_4, RSD_prev5_9, RSD_prev10_19,RSD_prev20_29, RSD_prev30_49, RSD_prev50_80,
                                    Atrp_prev0_1, Atrp_prev2_4, Atrp_prev5_9, Atrp_prev10_19, Atrp_prev20_29, Atrp_prev30_49, Atrp_prev50_80,
                                    HG_prev0_1, HG_prev2_4, HG_prev5_9, HG_prev10_19, HG_prev20_29, HG_prev30_49, HG_prev50_80,
                                    MDp_prev0_1, MDp_prev2_4, MDp_prev5_9, MDp_prev10_19, MDp_prev20_29, MDp_prev30_49, MDp_prev50_80,
                                    SDp_prev0_1, SDp_prev2_4, SDp_prev5_9, SDp_prev10_19, SDp_prev20_29, SDp_prev30_49, SDp_prev50_80){

SI_prev_temp <- sum(morb.mat.tmp[,19])/N  # severe itch
RSD_prev_temp <- sum(morb.mat.tmp[,20])/N # reactive skin disease
Atrp_prev_temp <- sum(morb.mat.tmp[,21])/N # skin atrophy
HG_prev_temp <- sum(morb.mat.tmp[,22])/N # hanging groin
MDp_prev_temp <- sum(morb.mat.tmp[,24])/N # mild depigmentation
SDp_prev_temp <- sum(morb.mat.tmp[,25])/N # severe depigmentation

# update prevalence vectors
SI_prev <- c(SI_prev, SI_prev_temp)
RSD_prev <- c(RSD_prev, RSD_prev_temp)
Atrp_prev <- c(Atrp_prev, Atrp_prev_temp)
HG_prev <- c(HG_prev, HG_prev_temp)
MDp_prev <- c(MDp_prev, MDp_prev_temp)
SDp_prev <- c(SDp_prev, SDp_prev_temp)

# update age-prevalence vectors

# severe itch age age-prev #
SI_prev0_1_temp <- length(which(morb.mat.tmp[,19] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
SI_prev0_1 <- c(SI_prev0_1, SI_prev0_1_temp)

SI_prev2_4_temp <- length(which(morb.mat.tmp[,19] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
SI_prev2_4 <- c(SI_prev2_4, SI_prev2_4_temp)

SI_prev5_9_temp <- length(which(morb.mat.tmp[,19] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
SI_prev5_9 <- c(SI_prev5_9, SI_prev5_9_temp)

SI_prev10_19_temp <- length(which(morb.mat.tmp[,19] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
SI_prev10_19 <- c(SI_prev10_19, SI_prev10_19_temp)

SI_prev20_29_temp <- length(which(morb.mat.tmp[,19] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
SI_prev20_29 <- c(SI_prev20_29, SI_prev20_29_temp)

SI_prev30_49_temp <- length(which(morb.mat.tmp[,19] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
SI_prev30_49 <- c(SI_prev30_49, SI_prev30_49_temp)

SI_prev50_80_temp <- length(which(morb.mat.tmp[,19] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
SI_prev50_80 <- c(SI_prev50_80, SI_prev50_80_temp)

# RSD age age-prev #
RSD_prev0_1_temp <- length(which(morb.mat.tmp[,20] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
RSD_prev0_1 <- c(RSD_prev0_1, RSD_prev0_1_temp)

RSD_prev2_4_temp <- length(which(morb.mat.tmp[,20] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
RSD_prev2_4 <- c(RSD_prev2_4, RSD_prev2_4_temp)

RSD_prev5_9_temp <- length(which(morb.mat.tmp[,20] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
RSD_prev5_9 <- c(RSD_prev5_9, RSD_prev5_9_temp)

RSD_prev10_19_temp <- length(which(morb.mat.tmp[,20] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
RSD_prev10_19 <- c(RSD_prev10_19, RSD_prev10_19_temp)

RSD_prev20_29_temp <- length(which(morb.mat.tmp[,20] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
RSD_prev20_29 <- c(RSD_prev20_29, RSD_prev20_29_temp)

RSD_prev30_49_temp <- length(which(morb.mat.tmp[,20] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
RSD_prev30_49 <- c(RSD_prev30_49, RSD_prev30_49_temp)

RSD_prev50_80_temp <- length(which(morb.mat.tmp[,20] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
RSD_prev50_80 <- c(RSD_prev50_80, RSD_prev50_80_temp)

# Atrophy age age-prev #
Atrp_prev0_1_temp <- length(which(morb.mat.tmp[,21] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
Atrp_prev0_1 <- c(Atrp_prev0_1, Atrp_prev0_1_temp)

Atrp_prev2_4_temp <- length(which(morb.mat.tmp[,21] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
Atrp_prev2_4 <- c(Atrp_prev2_4, Atrp_prev2_4_temp)

Atrp_prev5_9_temp <- length(which(morb.mat.tmp[,21] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
Atrp_prev5_9 <- c(Atrp_prev5_9, Atrp_prev5_9_temp)

Atrp_prev10_19_temp <- length(which(morb.mat.tmp[,21] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
Atrp_prev10_19 <- c(Atrp_prev10_19, Atrp_prev10_19_temp)

Atrp_prev20_29_temp <- length(which(morb.mat.tmp[,21] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
Atrp_prev20_29 <- c(Atrp_prev20_29, Atrp_prev20_29_temp)

Atrp_prev30_49_temp <- length(which(morb.mat.tmp[,21] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
Atrp_prev30_49 <- c(Atrp_prev30_49, Atrp_prev30_49_temp)

Atrp_prev50_80_temp <- length(which(morb.mat.tmp[,21] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
Atrp_prev50_80 <- c(Atrp_prev50_80, Atrp_prev50_80_temp)

# Hanging groin age age-prev #
HG_prev0_1_temp <- length(which(morb.mat.tmp[,22] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
HG_prev0_1 <- c(HG_prev0_1, HG_prev0_1_temp)

HG_prev2_4_temp <- length(which(morb.mat.tmp[,22] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
HG_prev2_4 <- c(HG_prev2_4, HG_prev2_4_temp)

HG_prev5_9_temp <- length(which(morb.mat.tmp[,22] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
HG_prev5_9 <- c(HG_prev5_9, HG_prev5_9_temp)

HG_prev10_19_temp <- length(which(morb.mat.tmp[,22] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
HG_prev10_19 <- c(HG_prev10_19, HG_prev10_19_temp)

HG_prev20_29_temp <- length(which(morb.mat.tmp[,22] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
HG_prev20_29 <- c(HG_prev20_29, HG_prev20_29_temp)

HG_prev30_49_temp <- length(which(morb.mat.tmp[,22] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
HG_prev30_49 <- c(HG_prev30_49, HG_prev30_49_temp)

HG_prev50_80_temp <- length(which(morb.mat.tmp[,22] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
HG_prev50_80 <- c(HG_prev50_80, HG_prev50_80_temp)

# Mild depigmentation age age-prev #
MDp_prev0_1_temp <- length(which(morb.mat.tmp[,24] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
MDp_prev0_1 <- c(MDp_prev0_1, MDp_prev0_1_temp)

MDp_prev2_4_temp <- length(which(morb.mat.tmp[,24] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
MDp_prev2_4 <- c(MDp_prev2_4, MDp_prev2_4_temp)

MDp_prev5_9_temp <- length(which(morb.mat.tmp[,24] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
MDp_prev5_9 <- c(MDp_prev5_9, MDp_prev5_9_temp)

MDp_prev10_19_temp <- length(which(morb.mat.tmp[,24] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
MDp_prev10_19 <- c(MDp_prev10_19, MDp_prev10_19_temp)

MDp_prev20_29_temp <- length(which(morb.mat.tmp[,24] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
MDp_prev20_29 <- c(MDp_prev20_29, MDp_prev20_29_temp)

MDp_prev30_49_temp <- length(which(morb.mat.tmp[,24] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
MDp_prev30_49 <- c(MDp_prev30_49, MDp_prev30_49_temp)

MDp_prev50_80_temp <- length(which(morb.mat.tmp[,24] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
MDp_prev50_80 <- c(MDp_prev50_80, MDp_prev50_80_temp)

# Severe depigmentation age age-prev #
SDp_prev0_1_temp <- length(which(morb.mat.tmp[,25] == 1 & morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2)) /  length(which(morb.mat.tmp[,1] >= 0 & morb.mat.tmp[,1] < 2))# 0 - 1 age
SDp_prev0_1 <- c(SDp_prev0_1, SDp_prev0_1_temp)

SDp_prev2_4_temp <- length(which(morb.mat.tmp[,25] == 1 & morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5)) /  length(which(morb.mat.tmp[,1] >= 2 & morb.mat.tmp[,1] < 5))# 2 - 4 age
SDp_prev2_4 <- c(SDp_prev2_4, SDp_prev2_4_temp)

SDp_prev5_9_temp <- length(which(morb.mat.tmp[,25] == 1 & morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10)) /  length(which(morb.mat.tmp[,1] >= 5 & morb.mat.tmp[,1] < 10))
SDp_prev5_9 <- c(SDp_prev5_9, SDp_prev5_9_temp)

SDp_prev10_19_temp <- length(which(morb.mat.tmp[,25] == 1 & morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20)) /  length(which(morb.mat.tmp[,1] >= 10 & morb.mat.tmp[,1] < 20))
SDp_prev10_19 <- c(SDp_prev10_19, SDp_prev10_19_temp)

SDp_prev20_29_temp <- length(which(morb.mat.tmp[,25] == 1 & morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30)) /  length(which(morb.mat.tmp[,1] >= 20 & morb.mat.tmp[,1] < 30))
SDp_prev20_29 <- c(SDp_prev20_29, SDp_prev20_29_temp)

SDp_prev30_49_temp <- length(which(morb.mat.tmp[,25] == 1 & morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50)) /  length(which(morb.mat.tmp[,1] >= 30 & morb.mat.tmp[,1] < 50))
SDp_prev30_49 <- c(SDp_prev30_49, SDp_prev30_49_temp)

SDp_prev50_80_temp <- length(which(morb.mat.tmp[,25] == 1 & morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] <= 80)) /  length(which(morb.mat.tmp[,1] >= 50 & morb.mat.tmp[,1] < 80))
SDp_prev50_80 <- c(SDp_prev50_80, SDp_prev50_80_temp)

return(list(SI_prev, RSD_prev, Atrp_prev, HG_prev, MDp_prev, SDp_prev,
            SI_prev0_1, SI_prev2_4, SI_prev5_9, SI_prev10_19, SI_prev20_29, SI_prev30_49, SI_prev50_80,
            RSD_prev0_1, RSD_prev2_4, RSD_prev5_9, RSD_prev10_19,RSD_prev20_29, RSD_prev30_49, RSD_prev50_80,
            Atrp_prev0_1, Atrp_prev2_4, Atrp_prev5_9, Atrp_prev10_19, Atrp_prev20_29, Atrp_prev30_49, Atrp_prev50_80,
            HG_prev0_1, HG_prev2_4, HG_prev5_9, HG_prev10_19, HG_prev20_29, HG_prev30_49, HG_prev50_80,
            MDp_prev0_1, MDp_prev2_4, MDp_prev5_9, MDp_prev10_19, MDp_prev20_29, MDp_prev30_49, MDp_prev50_80,
            SDp_prev0_1, SDp_prev2_4, SDp_prev5_9, SDp_prev10_19, SDp_prev20_29, SDp_prev30_49, SDp_prev50_80))

}



