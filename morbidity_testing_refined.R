# ================================================================================================================================== #
#                             Refined morbidity code (external to dynamic model) - 19/05/23                                          #
# ================================================================================================================================== #

# ===========================#
#   Onchosim approach        #

# ========================== #
#    TO DO AT ITER = 1       #
# ========================== #

# ========================== #
# For main morbidity matrix #
# things already present and needed for matrix in model_wrappers:
N <- 440
cur.age <- rep(0, N)
DT <- 1/366
mean.age <- 50
real.max.age <- 80

# get number of worms from iter = 1 and and endemic eq run

eq_outputs <- output_equilibrium$all_equilibrium_outputs[[1]]

# =================== #
# for mf dying matrix #
num.mf.comps = 21
all.mats.mfdie.temp <- matrix(nrow=N, ncol=num.mf.comps*2)

# if extracting from eq #
mf.start <- 7
mf.end <- 6 + num.mf.comps

# 1) calculate age-dep mort rates for mf age classes (1 - 21)

age.cats.mf <- seq(0, 2.5, length = num.mf.comps) #up to 2.5 years old (assume mf die after age 2.5 years)

weibull.mortality <- function(DT, par1, par2, age.cats)

{
  out <- DT  * (par1 ^ par2) * par2 * (age.cats ^ (par2-1))

  return(out)
}

mu.mf1 = 1.089; mu.mf2 = 1.428
#DT not relevent here because RK4 is used to calculate change in mf

#mort.rates.mf <- weibull.mortality(DT = 1, par1 = mu.mf1, par2 = mu.mf2, age.cats = age.cats.mf) # DT = 1 for RK4 - switch to DT = DT

DT <- 1/366
DT <- 30/366
mort.rates.mf <- weibull.mortality(DT = DT, par1 = mu.mf1, par2 = mu.mf2, age.cats = age.cats.mf) # DT = 1 for RK4 - switch to DT = DT


#================================================#
#              TO DO AT EACH ITER                #
#================================================#

# ===================== #
#     Mf dying matrix   #

# needed in each iter #
all.mats.mfdie.temp <- matrix(nrow=N, ncol=num.mf.comps*2)
all.mats.mfdie.temp[,c(1:21)] <- eq_outputs[,c(mf.start:mf.end)] # make first 21 current eq number of mf per individual x age class
all.mats.mfdie.temp[,c(1:21)] <- floor(all.mats.mfdie.temp[,c(1:21)]) # make mf integers: round down decimals (so mf 0 - 1 counted as 0)

# loop through to identify number of mf to die for each individual #

test_func_mfdie <- function(all.mats.mfdie.temp, num.mf.comps, N, mort.rates.mf){
  for (k in 1 : num.mf.comps) {

    # cl <- (ws-1) + compartment #calculate which column to use depending on sex, type (fertile or infertile) and compartment

    cl <- k # define col to select for mf age class

    cur.mf <- all.mats.mfdie.temp[, cl] # select mf age class column from matrix

    cl_to_upd <- cl + num.mf.comps # column to update number of mf dying in matrix

    all.mats.mfdie.temp[,cl_to_upd] <- rbinom(N, all.mats.mfdie.temp[,cl], rep(mort.rates.mf[cl], N))


  }

  mf.dying.mat <- all.mats.mfdie.temp[,c(22:42)]
  mf.alive.mat <- all.mats.mfdie.temp[,c(1:21)]

  indiv.mf.dying <- rowSums(mf.dying.mat)
  indiv.mf.alive <- rowSums(mf.alive.mat)

  return(list(indiv.mf.dying, indiv.mf.alive, all.mats.mfdie.temp))
}

mf.die.out <- test_func_mfdie(all.mats.mfdie.temp = all.mats.mfdie.temp, num.mf.comps = num.mf.comps, N = N,
                              mort.rates.mf = mort.rates.mf)

all.mats.mfdie.temp <- mf.die.out[[3]] # updated matrix
# =========================== #
#   Morbidity matrix          #

# morbidity model dims to define morbidity matrix (for N individuals)
num_indv_suscept_cols <- 6
tissue_dmg_cols <- 2 + (6 * 2) # need * 2 as need a last step tissue damage (dTx) and change in tissue damage (DTx)
disease_cond_cols <- 6

num.cols.morb <- 2 + num_indv_suscept_cols + tissue_dmg_cols + disease_cond_cols

# define morbidity matrix (skin disease states and blindness/VI)
all.morb.temp2 <- matrix(nrow=N, ncol=num.cols.morb)

#all.morb.temp2[,1] <- cur.age
#all.morb.temp2[,2] <- sex
all.morb.temp2[,c(10,12,14,16,18,22,23,24,25,26,28)] <- 0


# all.morb.temp2[,10] <- 0 # assign Dix_T as 0 initially (this will then be the old) - severe itch
# all.morb.temp2[,12] <- 0 # assign Dix_T as 0 initially (this will then be the old) - reactive skin
# all.morb.temp2[,14] <- 0 # assign Dix_T as 0 initially (this will then be the old) - atrophy
# all.morb.temp2[,16] <- 0 # assign Dix_T as 0 initially (this will then be the old) - hanging groin
# all.morb.temp2[,18] <- 0 # assign Dix_T as 0 initially (this will then be the old) - mild dipigmentation
# #all.morb.temp[,20] <- 0 # assign Dix_T as 0 initially (this will then be the old) - severe dipigmentation
#
# all.morb.temp2[,22] <- 0 # assign disease condition initially - severe itch
# all.morb.temp2[,23] <- 0 # assign disease condition initially - reactive skin
# all.morb.temp2[,24] <- 0 # assign disease condition initially - atrophy
# all.morb.temp2[,25] <- 0 # assign disease condition initially - hanging groin
# all.morb.temp2[,26] <- 0 # assign disease condition initially - mild dipigmentation
# all.morb.temp2[,28] <- 0 # assign disease condition initially - severe dipigmentation

# define individual variation in susceptibility
sevr_itch_suscpt <- rgamma(n=N, shape = 0.316, rate = 0.316) # matt: randomly assigned (binomal dist) with equal probability (e.g., psi_F = 0.5, psi_m = 0.5)
all.morb.temp2[,3] <- sevr_itch_suscpt

react_skin_suscpt <- rgamma(n=N, shape = 0.425, rate = 0.425) # matt: randomly assigned (binomal dist) with equal probability (e.g., psi_F = 0.5, psi_m = 0.5)
all.morb.temp2[,4] <- react_skin_suscpt

atrophy_suscpt <- rgamma(n=N, shape = 0.279, rate = 0.279) # matt: randomly assigned (binomal dist) with equal probability (e.g., psi_F = 0.5, psi_m = 0.5)
all.morb.temp2[,5] <- atrophy_suscpt

hanging_grn_suscpt <- rgamma(n=N, shape = 0.857, rate = 0.857) # matt: randomly assigned (binomal dist) with equal probability (e.g., psi_F = 0.5, psi_m = 0.5)
all.morb.temp2[,6] <- hanging_grn_suscpt

# mild_depgm_suscpt <- rgamma(n=N, shape = 0.246, rate = 0.246) # mild depigmentation (gamma dist: shape = 0.246)
# all.morb.temp[,7] <- mild_depgm_suscpt
#
# sevr_depgm_suscpt <- rgamma(n=N, shape = 0.246, rate = 0.246) # severe depigmentation (gamma dist: shape = 0.246)
# all.morb.temp[,8] <- sevr_depgm_suscpt

# only one individual susceptibility for pigmentation (so same for mild and severe)?
depgm_suscpt <- rgamma(n=N, shape = 0.246, rate = 0.246) # mild depigmentation (gamma dist: shape = 0.246)
all.morb.temp2[,7] <- depgm_suscpt

# ====================================== #
#  this should all be done AFTER first i #
# ====================================== #

# 1) add number of female worms (same as accured tissue damage di_t)
#all.morb.temp2[,9] <- indiv.mf.dying # assigning from equilibirum ouputs after 30 years (this is di_T)
all.morb.temp2[,9] <- mf.die.out[[1]] # assigning from equilibirum ouputs after 30 years (this is di_T)

hist(all.morb.temp[,9], breaks = 30)
mean(all.morb.temp[,9])
median(all.morb.temp[,9])

hist(all.morb.temp2[,9], breaks = 30)
mean(all.morb.temp2[,9])
median(all.morb.temp2[,9])

# =====================================================================================================#
# 2) calculate (new) change in amount of tissue activation/damage (change in Dix_T) for each condition #
# these are used to update the Dix_T col (matrix col 10) at the end, so then used in this equation as all.morb.temp[,10] in next time step

all.morb.temp2[,11] <- all.morb.temp2[,3] * all.morb.temp2[,9] - (0.015/30) * all.morb.temp2[,10] # severe itch

all.morb.temp2[,13] <- all.morb.temp2[,4] * all.morb.temp2[,9] - (0.035/30) * all.morb.temp2[,12] # reactive skin disease

all.morb.temp2[,15] <- all.morb.temp2[,5] * all.morb.temp2[,9] - 0 * all.morb.temp2[,14] # atrophy

all.morb.temp2[,17] <- all.morb.temp2[,6] * all.morb.temp2[,9] - 0 * all.morb.temp2[,16] # hanging groin

# all.morb.temp[,19] <- all.morb.temp[,7] * all.morb.temp[,9] - 0 * all.morb.temp[,18] # mild dipigmentation

# all.morb.temp[,21] <- all.morb.temp[,8] * all.morb.temp[,9] - 0 * all.morb.temp[,20] # severe dipigmentation

all.morb.temp2[,19] <- all.morb.temp2[,7] * all.morb.temp2[,9] - 0 * all.morb.temp2[,18] # ALL dipigmentation

# =============================================================#
# 3) check if threshold exceeded and disease condition present

all.morb.temp2[,22] <- ifelse(all.morb.temp2[,11] > 0.255, 1, 0) # severe itch (threshold = 0.255)

all.morb.temp2[,23] <- ifelse(all.morb.temp2[,13] > 0.210, 1, 0) # reactive skin disease (threshold = 0.255)

# need to select only 0's, if 1 then permanent in next disease conditions (so not chnaged, will remain as 1)

# all.morb.temp[,24] <- test_condition_col # to test if only 0's updated due to permanence of 1
all.morb.temp2[,24] <- ifelse(all.morb.temp2[,24] == 0 & all.morb.temp2[,15] > 11.3, 1, all.morb.temp2[,24]) # atrophy (threshold = 11.3) - need to only test 0's as 1's permanent?
#all.morb.temp[,24] <- ifelse(all.morb.temp[,15] > 11.3, 1, 0) # atrophy (threshold = 11.3) - need to only test 0's as 1's permanent?

all.morb.temp2[,25] <- ifelse(all.morb.temp2[,25] == 0 & all.morb.temp2[,17] > 21.4, 1, all.morb.temp2[,25]) # hanging groin (threshold = 21.4) - need to only test 0's as 1's permanent?
#all.morb.temp[,25] <- ifelse(all.morb.temp[,17] > 21.4, 1, 0)

# dipigmentation (two-stage: need to specify if mild or severe)

all.morb.temp2[,26] <- ifelse(all.morb.temp2[,26] == 0 & all.morb.temp2[,19] > 2.35, 1, all.morb.temp2[,26]) # temporary mild dipigmentation (threshold = 2.35) - need to only test 0's as 1's permanent?

# all.morb.temp[,27] <- ifelse(all.morb.temp[,21] > 4.3, 1, 0) # severe dipigmentation (threshold = 4.3) - need to only test 0's as 1's permanent?

all.morb.temp2[,28] <- ifelse(all.morb.temp2[,28] == 0 & all.morb.temp2[,19] > 4.3, 1, all.morb.temp2[,28]) # severe dipigmentation (threshold = 4.3) - need to only test 0's as 1's permanent?

all.morb.temp2[,27] <- ifelse(all.morb.temp2[,28] == 1, 0, all.morb.temp2[,26]) # updated)mild dipigm 0 if severe dipigm is 1 or leave (this used for prev calc)

# prev #
SI_prev2 <- sum(all.morb.temp2[,22])/N
RSD_prev2 <- sum(all.morb.temp2[,23])/N
Atrp_prev2 <- sum(all.morb.temp2[,24])/N
HG_prev2 <- sum(all.morb.temp2[,25])/N
MDp_prev2 <- sum(all.morb.temp2[,27])/N
SDp_prev2 <- sum(all.morb.temp2[,28])/N

prev2_out <- c(SI_prev2, RSD_prev2, Atrp_prev2, HG_prev2, MDp_prev2, SDp_prev2)

prev2_out*100


# =========================================================== #
#   Updating  the all.morb.temp when people die testing       #

num.comps.worm <- 21

to.die.test <- c(4,7)

#cols.to.zero <- seq(from = 1, to = (6 + num.mf.comps + 3*num.comps.worm))
#cols.to.zero <- cols.to.zero[-c(1,5, 6)] #compliance, L2 and L3 do not become zero when an individual dies

cols.to.zero.morb <- seq(from = 1, to = ncol(all.morb.temp)) # want all cols included in selecting whether becomes 0 if individual dies

all.morb.temp1[to.die.test, cols.to.zero.morb] <- 0 #sset all o 0 (check if individual susceptibilities should be reset?)

# redraw individual suscpetibilities in newly entering individuals

all.morb.temp1[to.die.test,7] <- rgamma(n=length(to.die.test), shape = 0.246, rate = 0.246) # mild depigmentation (gamma dist: shape = 0.246)

#  mf dying matrix #
to.die.test.mf <- c(2,3)
cols.to.zero.mfdie <- seq(from = 1, to = ncol(all.mats.mfdie.temp))
all.mats.mfdie.temp[to.die.test.mf, cols.to.zero.mfdie] <- 0


#=============================================================================================================#
#                               New Development 14-06-23                                                      #
# ============================================================================================================#

# =============================================================== #
#    Implementing a delay matrix for number of mf dying (30 days) #

N <- 440

# at i = 1 #
num.days.cols <- 30
mf.die.mat <- matrix(0, ncol = num.days.cols, nrow = N)
inds.mfdie.mat <- seq(2,(length(mf.die.mat[1,]))) # for moving columns along with time

# i > 1 #

mf.die.mat[,inds.mfdie.mat] <- mf.die.mat[,(inds.mfdie.mat-1)] # move mf dying along one col

mf_die_per_day <- c() # new vec for this day to fill

dummy_mf_die_day_vec <- ifelse(runif(440) < 0.2, 0, (runif(440) ^ 5)) # generate list of numbers (dummy mf dying per day numbers)
                                                                      # bias towards small numbers and 0 values to simulate mf dying numbers?
mf_die_per_day <- c(dummy_mf_die_day_vec) # extract sum of mf dying per individual for that day

mf.die.mat[,1] <- mf_die_per_day # update first col with number of mf dying on that day

mf.die.month <- c()
mf.die.month <- rowSums(mf.die.mat) # calculate number of mf dying in last 30 days

# all.morb.temp2[,9] <- mf.die.out[[1]] # update col in main morb.mat to perform accumulated tissue damage calcs

# update any individuals that die (set row to 0)
to.die <- c(2,50,400)
mf.die.mat[to.die, ] <- 0


# ================================================================== #
#     Setting up a new morb.matrix structure                         #

# at i = 1 #

# morbidity model dims to define morbidity matrix (for N individuals)
demographic_cols <- 2 # age and sex cols
num_indv_suscept_cols <- 5
accured_tissue_dmg_mnth_col <- 1
tissue_dmg_amount_cols <- 5 * 2 # need * 2 as need a last step tissue damage (dTx) and change in tissue damage (DTx)
disease_cond_cols <- 6 + 1 # 6 for all OSD conditions, + 1 for two-stage depigm process
num.cols.morb <- demographic_cols + num_indv_suscept_cols + accured_tissue_dmg_mnth_col + tissue_dmg_amount_cols + disease_cond_cols # total cols required

# define Morbidity matrix #
all.morb.temp <- matrix(nrow = N, ncol = num.cols.morb)

all.morb.temp[,1] <- cur.age
all.morb.temp[,2] <- sex

# define individual variation in susceptibility
all.morb.temp[,3] <- rgamma(n = N, shape = 0.316, rate = 0.316) # severe itch (gamma dist: shape = 0.316)
all.morb.temp[,4] <- rgamma(n = N, shape = 0.425, rate = 0.425) # reactive skin disease (gamma dist: shape = 0.425)
all.morb.temp[,5] <- rgamma(n = N, shape = 0.279, rate = 0.279) # atrophy (gamma dist: shape = 0.279)
all.morb.temp[,6] <- rgamma(n = N, shape = 0.857, rate = 0.857) # hanging groin (gamma dist: shape = 0.857)
all.morb.temp[,7] <- rgamma(n = N, shape = 0.246, rate = 0.246) # mild depigmentation (gamma dist: shape = 0.246)

#assign Dix_T as 0 initially (this will then be the old) & assign disease condition initially as 0
all.morb.temp[,c(8,9,11,13,15,17,19,20,21,22,23,24,25)] <- 0


# ==================================================================== #
#            Updated mf mortality matrix by age class                  #



# ===================== #
#   at i = 1            #
# =================== #
# for mf dying matrix #
num.mf.comps = 21
mf.start <- 7
mf.end <- 6 + num.mf.comps

all.mats.mfdie.temp <- matrix(nrow = N, ncol = num.mf.comps * 2)


# ===================== #
#       i > 1           #
# ===================== #
#     Mf dying matrix   #

all.mats.temp <- output_equilibrium$all_equilibrium_outputs[[1]]

# needed in each iter #

all.mats.mfdie.temp[,c(1:21)] <- all.mats.temp[,c(mf.start:mf.end)] # make first 21 current eq number of mf per individual x age class
#all.mats.mfdie.temp[,c(1:21)] <- floor(all.mats.mfdie.temp[,c(1:21)]) # make mf integers: round down decimals (so mf 0 - 1 counted as 0)

# loop through to identify number of mf to die for each individual #

test_func_mfdie <- function(all.mats.mfdie.temp, num.mf.comps, N, mort.rates.mf){
  for (k in 1 : num.mf.comps) {

    # cl <- (ws-1) + compartment #calculate which column to use depending on sex, type (fertile or infertile) and compartment

    cl <- k # define col to select for mf age class

    cur.mf <- all.mats.mfdie.temp[, cl] # select mf age class column from matrix

    cl_to_upd <- cl + num.mf.comps # column to update number of mf dying in matrix

    # all.mats.mfdie.temp[,cl_to_upd] <- rbinom(N, all.mats.mfdie.temp[,cl], rep(mort.rates.mf[cl], N))
    all.mats.mfdie.temp[,cl_to_upd] <- all.mats.mfdie.temp[,cl] * mort.rates.mf[cl]


  }

  mf.dying.mat <- all.mats.mfdie.temp[,c(22:42)]
  #mf.alive.mat <- all.mats.mfdie.temp[,c(1:21)]

  indiv.mf.dying <- rowSums(mf.dying.mat)
  #indiv.mf.alive <- rowSums(mf.alive.mat)

  return(list(indiv.mf.dying, all.mats.mfdie.temp))
}

mf.die.out <- test_func_mfdie(all.mats.mfdie.temp = all.mats.mfdie.temp, num.mf.comps = num.mf.comps, N = N,
                              mort.rates.mf = mort.rates.mf)

range(mf.die.out[[1]])

# =================================================== #
#       Updating main.morb.mat (i > 1)                #
# =================================================== #

#mat_test <- matrix(nrow=440, ncol=25)
# test_vec <- numeric_vector <- seq(1, -1, length.out = 440)
# mat_test[,10] <- test_vec
mat_test[,10] <- ifelse(mat_test[,10] < 0, 0, mat_test[,10])

morb.mat <- all.morb.temp

# update age / sex cols (if need to cal prev by age and sex)
morb.mat[,1] <- all.mats.temp[,2]
morb.mat[,2] <- all.mats.temp[,3]

# update column 9 of morbidity matrix with total mf dying per individual in time-step
morb.mat[,8] <- mf.die.out[[1]]

# calculate (new) change in amount of tissue activation/damage (change in Dix_T) for each condition #
# this is calculated by: individual.susceptibility.condition * number.mf.dying - regression.rate * current tissue damage
SI_regression_rate <- 0.015
RSD_regression_rate <- 0.030

morb.mat[,10] <- morb.mat[,3] * morb.mat[,8] - (SI_regression_rate) * morb.mat[,9] # severe itch

morb.mat[,12] <- morb.mat[,4] * morb.mat[,8] - (RSD_regression_rate) * morb.mat[,11] # reactive skin disease
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
SI_threshold <- 0.255
RSD_threshold <- 0.210
atrophy_threshold <- 11.3
HG_threshold <- 21.4
mdpigm_threshold <- 2.35
sdepigm_threshold <- 4.3

morb.mat[,19] <- ifelse(morb.mat[,10] > SI_threshold, 1, 0) # severe itch (threshold = 0.255)
morb.mat[,20] <- ifelse(morb.mat[,12] > RSD_threshold, 1, 0) # reactive skin disease (threshold = 0.255)
morb.mat[,21] <- ifelse(morb.mat[,21] == 0 & morb.mat[,14] > atrophy_threshold, 1, morb.mat[,21]) # atrophy (threshold = 11.3) - need to only test 0's as 1's permanent
morb.mat[,22] <- ifelse(morb.mat[,22] == 0 & morb.mat[,16] > HG_threshold, 1, morb.mat[,22]) # hanging groin (threshold = 21.4) - need to only test 0's as 1's permanent

morb.mat[,23] <- ifelse(morb.mat[,23] == 0 & morb.mat[,18] > mdpigm_threshold, 1, morb.mat[,23]) # temporary mild depigmentation (threshold = 2.35) - need to only test 0's as 1's permanent
morb.mat[,25] <- ifelse(morb.mat[,25] == 0 & morb.mat[,18] > sdepigm_threshold, 1, morb.mat[,25]) # severe depigmentation (threshold = 4.3) - need to only test 0's as 1's permanent
morb.mat[,24] <- ifelse(morb.mat[,25] == 1, 0, morb.mat[,23]) # update mild depigm i.e., with 0 if severe depigm is 1 or leave (this used for prev calculation for mild depigmentation)

# check if threshold exceeded and disease condition present (thresholds are on original month scale)
# morb.mat[,22] <- ifelse(morb.mat[,11] > SI_threshold, 1, 0) # severe itch (threshold = 0.255)
# morb.mat[,23] <- ifelse(morb.mat[,13] > RSD_threshold, 1, 0) # reactive skin disease (threshold = 0.255)
# morb.mat[,24] <- ifelse(morb.mat[,24] == 0 & morb.mat[,15] > atrophy_threshold, 1, morb.mat[,24]) # atrophy (threshold = 11.3) - need to only test 0's as 1's permanent
# morb.mat[,25] <- ifelse(morb.mat[,25] == 0 & morb.mat[,17] > HG_threshold, 1, morb.mat[,25]) # hanging groin (threshold = 21.4) - need to only test 0's as 1's permanent
#
# morb.mat[,26] <- ifelse(morb.mat[,26] == 0 & morb.mat[,19] > mdpigm_threshold, 1, morb.mat[,26]) # temporary mild depigmentation (threshold = 2.35) - need to only test 0's as 1's permanent
# morb.mat[,28] <- ifelse(morb.mat[,28] == 0 & morb.mat[,19] > sdepigm_threshold, 1, morb.mat[,28]) # severe depigmentation (threshold = 4.3) - need to only test 0's as 1's permanent
# morb.mat[,27] <- ifelse(morb.mat[,28] == 1, 0, morb.mat[,26]) # update mild depigm i.e., with 0 if severe depigm is 1 or leave (this used for prev calculation for mild depigmentation)

# ========================================================== #
#   making a time vector to run OSD model once per month     #

DT <- 1/366
time.its <- 200 # 100 years
time.its <- round(time.its / (DT)) # how many total iterations to run (on a per-day scale, 1 iter = 1 day)

# DT.month <- 1/30
# time.its2 <- 200 # 100 years
# time.its.month <- round(time.its2 / (DT.month))

# how many its househoulds there be (individual iterations) e.g., for 200 years
366*200 # how many iters where each iter = 1 day
12*200 # how many iters where each iter = 1 month
1*200 # how many iters where each iter = 1 year

day_its <- seq(0, time.its, 1)
month_its <- seq(0, time.its, 30)
year_its <- seq(0, time.its, 366)

# final.yr.tmept <- 200
# day.iter <- seq(0, (final.yr.tmept / (1/366)), by = 1/366)
#
#
# day.iter2 <- seq(0, 50, 1/366)

# ================== #
# to set up vector   #

# day <- 1
# yr <- 366 * day
n.yrs <- 200
# tot.days <- yr * n.yrs
iters.day <- seq(1, 366*n.yrs, 1)

mnth_iter_select <- iters.day[seq(30, length(iters.day), by = 30)] # select every 30th element (1 month)

