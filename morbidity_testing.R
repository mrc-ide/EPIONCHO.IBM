# ==================================================================================================================== #

# things already present and needed for matrix in model_wrappers:

N <- 440
cur.age <- rep(0, N)
DT <- 1/366
mean.age <- 50
real.max.age <- 80


#(the approach below must be used, drawing human lifespans from an exponential distribution eventually leads to a non-exponential distribution)
for(i in 1 : 75000) #if at equilibrium you saved the age at which inds die and simulated further, you should get an exponential distribution
{
  cur.age <- cur.age + DT

  death.vec <- rbinom(N, 1, (1/mean.age) * DT) # Matt: human mortality (constant with age) - no. of deaths at time step t is a random variable drawn from binomial distribution (N = human pop size?)

  cur.age[which(death.vec == 1)] <- 0 #set individuals which die to age 0
  cur.age[which(cur.age >= real.max.age)] <- 0 #all individuals >= maximum imposed age die (matt: distribution truncated to prevent excessively long life spans - a_max)
}

sex.rat <- 0.5
sex <- rbinom(N, 1, sex.rat)

# get number of worms from iter = 1 and and endemic eq run

#
eq_outputs <- output_equilibrium$all_equilibrium_outputs[[1]]
#nfw.start <- 7 + 21 + 21 # start of infertile worms
#fw.end <- 90 # end of fertile worms
# female_worms_all <- rowSums(eq_outputs[, nfw.start : fw.end])

female_worms_all <- c(33, 31, 1, 1, 1, 20, 43, 2, 0, 11, 38, 25, 0, 0, 0, 0, 12, 2, 2, 11, 16, 1, 0, 31, 4, 0, 45,
                      0, 0, 7, 4, 8, 2, 14, 4, 0, 0, 0, 17, 6, 7, 0, 8, 11, 97, 17, 0, 2, 1, 0, 0, 2, 2, 8, 2, 0,
                      0, 47, 7, 0, 0, 40, 1, 89, 0, 0, 0, 25, 1, 0, 0, 18, 17, 4, 6, 2, 77, 0, 30, 28, 0, 1, 8, 0,
                      1, 0, 1, 110, 0, 18, 0, 2, 0, 22, 19, 0, 0, 0, 0, 4, 0, 70, 23, 45, 0, 5, 2, 5, 2, 0, 0, 1,
                      0, 0, 0, 0, 26, 48, 4, 1, 36, 9, 37, 12, 0, 19, 0, 17, 1, 20, 0, 0, 0, 4, 7, 6, 10, 0, 5, 47,
                      7, 30, 16, 0, 27, 0, 7, 3, 1, 42, 7, 0, 39, 0, 26, 0, 1, 0, 4, 0, 0, 1, 0, 14, 1, 0, 21, 19,
                      1, 29, 73, 1, 0, 2, 14, 8, 0, 9, 0, 25, 1, 43, 0, 0, 81, 2, 0, 22, 2, 0, 6, 1, 15, 4, 6, 0,
                      1, 91, 8, 0, 46, 3, 0, 8, 0, 5, 12, 13, 7, 1, 0, 0, 0, 1, 0, 9, 0, 8, 11, 0, 4, 2, 0, 0, 49,
                      21, 24, 29, 6, 6, 7, 0, 0, 2, 9, 13, 1, 3, 0, 4, 11, 0, 4, 55, 9, 0, 28, 0, 3, 20, 0, 11, 21,
                      19, 5, 1, 0, 11, 6, 27, 2, 0, 15, 18, 7, 3, 0, 0, 35, 2, 0, 0, 12, 14, 0, 0, 4, 0, 67, 1, 0,
                      13, 15, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 62, 64, 0, 3, 7, 0, 7, 0, 0, 51, 4, 91, 0, 4, 52,
                      5, 27, 0, 0, 20, 23, 3, 72, 0, 36, 0, 0, 5, 7, 8, 15, 1, 0, 53, 15, 61, 8, 2, 74, 6, 2, 0, 1,
                      9, 47, 1, 4, 5, 14, 0, 0, 60, 0, 10, 14, 8, 0, 0, 0, 30, 0, 4, 70, 11, 6, 0, 0, 0, 0, 30, 3,
                      0, 9, 1, 7, 19, 0, 16, 3, 0, 27, 22, 77, 1, 0, 1, 26, 0, 3, 1, 9, 0, 4, 2, 0, 5, 4, 0, 0, 5,
                      5, 12, 7, 39, 0, 12, 19, 2, 81, 0, 13, 23, 34, 41, 0, 40, 0, 0, 2, 19, 0, 15, 88, 47, 14, 9,
                      0, 3, 3, 0, 0, 10, 2, 10, 7, 9, 5, 60, 54, 0, 10, 9, 0, 0, 1, 1) # taken from EQ at 30 yrs

#write.csv(female_worms_all,file="female_worms_all_EQ.csv",row.names=F)

# ============================================================================= #
#    need to track number of mf killed per individual (instead of female worms) #


# ========================== #
#   New Morbidity            #

# =================================== #
#  this should all be done at first i #

# morbidity model dims to define morbidity matrix (for N individuals)
num_indv_suscept_cols <- 6
tissue_dmg_cols <- 2 + (6 * 2) # need * 2 as need a last step tissue damage (dTx) and change in tissue damage (DTx)
disease_cond_cols <- 6

num.cols.morb <- 2 + num_indv_suscept_cols + tissue_dmg_cols + disease_cond_cols

# define morbidity matrix (skin disease states and blindness/VI)
all.morb.temp <- matrix(nrow=N, ncol=num.cols.morb)

all.morb.temp[,1] <- cur.age
all.morb.temp[,2] <- sex

all.morb.temp[,10] <- 0 # assign Dix_T as 0 initially (this will then be the old) - severe itch
all.morb.temp[,12] <- 0 # assign Dix_T as 0 initially (this will then be the old) - reactive skin
all.morb.temp[,14] <- 0 # assign Dix_T as 0 initially (this will then be the old) - atrophy
all.morb.temp[,16] <- 0 # assign Dix_T as 0 initially (this will then be the old) - hanging groin
all.morb.temp[,18] <- 0 # assign Dix_T as 0 initially (this will then be the old) - mild dipigmentation
#all.morb.temp[,20] <- 0 # assign Dix_T as 0 initially (this will then be the old) - severe dipigmentation

all.morb.temp[,22] <- 0 # assign disease condition initially - severe itch
all.morb.temp[,23] <- 0 # assign disease condition initially - reactive skin
all.morb.temp[,24] <- 0 # assign disease condition initially - atrophy
all.morb.temp[,25] <- 0 # assign disease condition initially - hanging groin
all.morb.temp[,26] <- 0 # assign disease condition initially - mild dipigmentation
all.morb.temp[,28] <- 0 # assign disease condition initially - severe dipigmentation



# define individual variation in susceptibility
sevr_itch_suscpt <- rgamma(n=N, shape = 0.316, rate = 0.316) # matt: randomly assigned (binomal dist) with equal probability (e.g., psi_F = 0.5, psi_m = 0.5)
all.morb.temp[,3] <- sevr_itch_suscpt

react_skin_suscpt <- rgamma(n=N, shape = 0.425, rate = 0.425) # matt: randomly assigned (binomal dist) with equal probability (e.g., psi_F = 0.5, psi_m = 0.5)
all.morb.temp[,4] <- react_skin_suscpt

atrophy_suscpt <- rgamma(n=N, shape = 0.279, rate = 0.279) # matt: randomly assigned (binomal dist) with equal probability (e.g., psi_F = 0.5, psi_m = 0.5)
all.morb.temp[,5] <- atrophy_suscpt

hanging_grn_suscpt <- rgamma(n=N, shape = 0.857, rate = 0.857) # matt: randomly assigned (binomal dist) with equal probability (e.g., psi_F = 0.5, psi_m = 0.5)
all.morb.temp[,6] <- hanging_grn_suscpt

# mild_depgm_suscpt <- rgamma(n=N, shape = 0.246, rate = 0.246) # mild depigmentation (gamma dist: shape = 0.246)
# all.morb.temp[,7] <- mild_depgm_suscpt
#
# sevr_depgm_suscpt <- rgamma(n=N, shape = 0.246, rate = 0.246) # severe depigmentation (gamma dist: shape = 0.246)
# all.morb.temp[,8] <- sevr_depgm_suscpt

# only one individual susceptibility for pigmentation (so same for mild and severe)?
depgm_suscpt <- rgamma(n=N, shape = 0.246, rate = 0.246) # mild depigmentation (gamma dist: shape = 0.246)
all.morb.temp[,7] <- depgm_suscpt

# ====================================== #
#  this should all be done AFTER first i #
# ====================================== #

# 1) add number of female worms (same as accured tissue damage di_t)
all.morb.temp[,9] <- female_worms_all # assigning from equilibirum ouputs after 30 years (this is di_T)

# =====================================================================================================#
# 2) calculate (new) change in amount of tissue activation/damage (change in Dix_T) for each condition #
# these are used to update the Dix_T col (matrix col 10) at the end, so then used in this equation as all.morb.temp[,10] in next time step

all.morb.temp[,11] <- all.morb.temp[,3] * all.morb.temp[,9] - (0.015/30) * all.morb.temp[,10] # severe itch

all.morb.temp[,13] <- all.morb.temp[,4] * all.morb.temp[,9] - (0.035/30) * all.morb.temp[,12] # reactive skin disease

all.morb.temp[,15] <- all.morb.temp[,5] * all.morb.temp[,9] - 0 * all.morb.temp[,14] # atrophy

all.morb.temp[,17] <- all.morb.temp[,6] * all.morb.temp[,9] - 0 * all.morb.temp[,16] # hanging groin

# all.morb.temp[,19] <- all.morb.temp[,7] * all.morb.temp[,9] - 0 * all.morb.temp[,18] # mild dipigmentation

# all.morb.temp[,21] <- all.morb.temp[,8] * all.morb.temp[,9] - 0 * all.morb.temp[,20] # severe dipigmentation

all.morb.temp[,19] <- all.morb.temp[,7] * all.morb.temp[,9] - 0 * all.morb.temp[,18] # ALL dipigmentation

# =============================================================#
# 3) check if threshold exceeded and disease condition present

all.morb.temp[,22] <- ifelse(all.morb.temp[,11] > 0.255, 1, 0) # severe itch (threshold = 0.255)

all.morb.temp[,23] <- ifelse(all.morb.temp[,13] > 0.210, 1, 0) # reactive skin disease (threshold = 0.255)


# need to select only 0's, if 1 then permanent in next disease conditions (so not chnaged, will remain as 1)

# all.morb.temp[,24] <- test_condition_col # to test if only 0's updated due to permanence of 1
all.morb.temp[,24] <- ifelse(all.morb.temp[,24] == 0 & all.morb.temp[,15] > 11.3, 1, all.morb.temp[,24]) # atrophy (threshold = 11.3) - need to only test 0's as 1's permanent?
#all.morb.temp[,24] <- ifelse(all.morb.temp[,15] > 11.3, 1, 0) # atrophy (threshold = 11.3) - need to only test 0's as 1's permanent?

all.morb.temp[,25] <- ifelse(all.morb.temp[,25] == 0 & all.morb.temp[,17] > 21.4, 1, all.morb.temp[,25]) # hanging groin (threshold = 21.4) - need to only test 0's as 1's permanent?
#all.morb.temp[,25] <- ifelse(all.morb.temp[,17] > 21.4, 1, 0)

# dipigmentation (two-stage: need to specify if mild or severe)

all.morb.temp[,26] <- ifelse(all.morb.temp[,26] == 0 & all.morb.temp[,19] > 2.35, 1, all.morb.temp[,26]) # temporary mild dipigmentation (threshold = 2.35) - need to only test 0's as 1's permanent?

# all.morb.temp[,27] <- ifelse(all.morb.temp[,21] > 4.3, 1, 0) # severe dipigmentation (threshold = 4.3) - need to only test 0's as 1's permanent?

all.morb.temp[,28] <- ifelse(all.morb.temp[,28] == 0 & all.morb.temp[,19] > 4.3, 1, all.morb.temp[,28]) # severe dipigmentation (threshold = 4.3) - need to only test 0's as 1's permanent?

all.morb.temp[,27] <- ifelse(all.morb.temp[,28] == 1, 0, all.morb.temp[,26]) # updated)mild dipigm 0 if severe dipigm is 1 or leave (this used for prev calc)



# random testing #
sum(all.morb.temp[,22])
210/440
sum(all.morb.temp[,23])
sum(all.morb.temp[,24])
68/440
sum(all.morb.temp[,25])
sum(all.morb.temp[,26])

test_condition_col <- sample(c(0,1), replace=TRUE, size=440)

all.morb.temp1 <- all.morb.temp # store


# ================================================================================================================= #
#                       Tracking mf dying                                                                           #

num.mf.comps = 21
all.mats.mfdie.temp <- matrix(nrow=N, ncol=num.mf.comps*2)

# if extracting from eq #
 mf.start <- 7
 mf.end <- 6 + num.mf.comps
# all.mats.mfdie.temp <- eq_outputs[,c(mf.start:mf.end)]
# all.mf.t1 <- format(rowSums(all.mats.mfdie.temp), scientific = FALSE)

 all.mats.mfdie.temp[,c(1:21)] <- eq_outputs[,c(mf.start:mf.end)] # make first 21 current eq number of mf per individual x age class

 # =========================================================== #
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

# age.cats.worms <- seq(0, 20, length = 21)
# mort.rates.worms <- weibull.mortality(DT = DT, par1 = 0.09953, par2 = 6.00569, age.cats = age.cats.worms) # DT = 1 for RK4 - switch to DT = DT


# ============================================================= #
# 1b) diff way of calculating no. of mf dead based on worm approach #

# cl <- (ws-1) + compartment #calculate which column to use depending on sex, type (fertile or infertile) and compartment
#
# cur.Wm <- total.dat[, cl] #take current number of worms from matrix
#
# worm.dead.males <- rbinom(N, cur.Wm, rep(mort.rates[compartment], N))
# worm.loss.males <- rbinom(N, (cur.Wm - worm.dead.males), rep((DT / time.each.comp), N))

# ============================================================================================================= #
# 2) apply age-dep mortality rates to each respective col to calculate no. mf expected to die in each age class #

range(all.mats.mfdie.temp[,4])

# need to loop through each colum (age class) and apply specific age-dep mortality rates
# all.mats.mfdie.temp[,22] <- mort.rates.mf[1] * all.mats.mfdie.temp[,1] # age group 1 mf dying
# all.mats.mfdie.temp[,23] <- mort.rates.mf[2] * all.mats.mfdie.temp[,2] # age group 2 mf dying
# all.mats.mfdie.temp[,24] <- mort.rates.mf[3] * all.mats.mfdie.temp[,3] # age group 3 mf dying
# all.mats.mfdie.temp[,25] <- mort.rates.mf[4] * all.mats.mfdie.temp[,4] # age group 4 mf dying
# all.mats.mfdie.temp[,26] <- mort.rates.mf[5] * all.mats.mfdie.temp[,5] # age group 5 mf dying
# all.mats.mfdie.temp[,27] <- mort.rates.mf[6] * all.mats.mfdie.temp[,6] # age group 6 mf dying
# all.mats.mfdie.temp[,28] <- mort.rates.mf[7] * all.mats.mfdie.temp[,7] # age group 7 mf dying
# all.mats.mfdie.temp[,29] <- mort.rates.mf[8] * all.mats.mfdie.temp[,8] # age group 4 mf dying
# all.mats.mfdie.temp[,30] <- mort.rates.mf[9] * all.mats.mfdie.temp[,9] # age group 4 mf dying
# all.mats.mfdie.temp[,31] <- mort.rates.mf[10] * all.mats.mfdie.temp[,10]
# all.mats.mfdie.temp[,32] <- mort.rates.mf[11] * all.mats.mfdie.temp[,11]
# all.mats.mfdie.temp[,33] <- mort.rates.mf[12] * all.mats.mfdie.temp[,12]
# all.mats.mfdie.temp[,34] <- mort.rates.mf[13] * all.mats.mfdie.temp[,13]
# all.mats.mfdie.temp[,35] <- mort.rates.mf[14] * all.mats.mfdie.temp[,14]
# all.mats.mfdie.temp[,36] <- mort.rates.mf[15] * all.mats.mfdie.temp[,15]
# all.mats.mfdie.temp[,37] <- mort.rates.mf[16] * all.mats.mfdie.temp[,16]
# all.mats.mfdie.temp[,38] <- mort.rates.mf[17] * all.mats.mfdie.temp[,17]
# all.mats.mfdie.temp[,39] <- mort.rates.mf[18] * all.mats.mfdie.temp[,18]
# all.mats.mfdie.temp[,40] <- mort.rates.mf[19] * all.mats.mfdie.temp[,19]
# all.mats.mfdie.temp[,41] <- mort.rates.mf[20] * all.mats.mfdie.temp[,20]
# all.mats.mfdie.temp[,42] <- mort.rates.mf[21] * all.mats.mfdie.temp[,21] # age group 21 mf dying

# 1) need to extract mf age-class cols from matrix

eq_outputs <- output_equilibrium$all_equilibrium_outputs[[1]]
mf.start <- 7
mf.end <- 6 + num.mf.comps
num.mf.comps = 21

# needed in each iter #
all.mats.mfdie.temp <- matrix(nrow=N, ncol=num.mf.comps*2)
all.mats.mfdie.temp[,c(1:21)] <- eq_outputs[,c(mf.start:mf.end)] # make first 21 current eq number of mf per individual x age class
all.mats.mfdie.temp[,c(1:21)] <- floor(all.mats.mfdie.temp[,c(1:21)]) # make mf integers: round down decimals (so mf 0 - 1 counted as 0)

# make loop #

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


# all.mats.mfdie.temp[,22] <- rbinom(N, all.mats.mfdie.temp[,1], rep(mort.rates.mf[1], N))
# all.mats.mfdie.temp[,23] <- rbinom(N, all.mats.mfdie.temp[,2], rep(mort.rates.mf[2], N))
# # rbinom(N, all.mats.mfdie.temp[,2], rep(0.05, N)) # playing around with probability (changing 0.001809674 to 0.5)
# all.mats.mfdie.temp[,24] <- rbinom(N, all.mats.mfdie.temp[,3], rep(mort.rates.mf[3], N))
# all.mats.mfdie.temp[,25] <- rbinom(N, all.mats.mfdie.temp[,4], rep(mort.rates.mf[4], N))
# all.mats.mfdie.temp[,26] <- rbinom(N, all.mats.mfdie.temp[,5], rep(mort.rates.mf[5], N))
# all.mats.mfdie.temp[,27] <- rbinom(N, all.mats.mfdie.temp[,6], rep(mort.rates.mf[6], N))
# all.mats.mfdie.temp[,28] <- rbinom(N, all.mats.mfdie.temp[,7], rep(mort.rates.mf[7], N))
# all.mats.mfdie.temp[,29] <- rbinom(N, all.mats.mfdie.temp[,8], rep(mort.rates.mf[8], N))
# all.mats.mfdie.temp[,30] <- rbinom(N, all.mats.mfdie.temp[,9], rep(mort.rates.mf[9], N))
# all.mats.mfdie.temp[,31] <- rbinom(N, all.mats.mfdie.temp[,10], rep(mort.rates.mf[10], N))
# all.mats.mfdie.temp[,32] <- rbinom(N, all.mats.mfdie.temp[,11], rep(mort.rates.mf[11], N))
# all.mats.mfdie.temp[,33] <- rbinom(N, all.mats.mfdie.temp[,12], rep(mort.rates.mf[12], N))
# all.mats.mfdie.temp[,34] <- rbinom(N, all.mats.mfdie.temp[,13], rep(mort.rates.mf[13], N))
# all.mats.mfdie.temp[,35] <- rbinom(N, all.mats.mfdie.temp[,14], rep(mort.rates.mf[14], N))
# all.mats.mfdie.temp[,36] <- rbinom(N, all.mats.mfdie.temp[,15], rep(mort.rates.mf[15], N))
# all.mats.mfdie.temp[,37] <- rbinom(N, all.mats.mfdie.temp[,16], rep(mort.rates.mf[16], N))
# all.mats.mfdie.temp[,38] <- rbinom(N, all.mats.mfdie.temp[,17], rep(mort.rates.mf[17], N))
# all.mats.mfdie.temp[,39] <- rbinom(N, all.mats.mfdie.temp[,18], rep(mort.rates.mf[18], N))
# all.mats.mfdie.temp[,40] <- rbinom(N, all.mats.mfdie.temp[,19], rep(mort.rates.mf[19], N))
# all.mats.mfdie.temp[,41] <- rbinom(N, all.mats.mfdie.temp[,20], rep(mort.rates.mf[20], N))
# all.mats.mfdie.temp[,42] <- rbinom(N, all.mats.mfdie.temp[,21], rep(mort.rates.mf[21], N))

# range(all.mats.mfdie.temp[,2]) # current mf in age group 2
# sum(all.mats.mfdie.temp[,2])
# range(all.mats.mfdie.temp[,23]) # mf to die in age group 2
# sum(all.mats.mfdie.temp[,23])
#
# range(all.mats.mfdie.temp[,3]) # current mf in age group 3
# sum(all.mats.mfdie.temp[,3])
# range(all.mats.mfdie.temp[,24]) # mf to die in age group 3
# sum(all.mats.mfdie.temp[,24])
#
# range(all.mats.mfdie.temp[,4]) # current mf in age group 4
# sum(all.mats.mfdie.temp[,4])
# range(all.mats.mfdie.temp[,25]) # mf to die in age group 4
# sum(all.mats.mfdie.temp[,25])
#
# range(all.mats.mfdie.temp[,21]) # current mf in age group 21
# sum(all.mats.mfdie.temp[,21])
# range(all.mats.mfdie.temp[,42]) # mf to die in age group 21
# sum(all.mats.mfdie.temp[,42])


mf.dying.mat <- all.mats.mfdie.temp[,c(22:42)]
mf.alive.mat <- all.mats.mfdie.temp[,c(1:21)]

#indiv.mf.dying <- format(rowSums(mf.dying.mat), scientific = FALSE)
#indiv.mf.alive <- format(rowSums(mf.alive.mat), scientific = FALSE)

indiv.mf.dying <- rowSums(mf.dying.mat)
indiv.mf.alive <- rowSums(mf.alive.mat)

# indiv.mf.dying2 <- indiv.mf.alive * (0.8/365)


# =========================================================================== #
# try in the morbidity matrix (indiv.mf.dying instead of num.female_worms)    #

# =================================== #
#  this should all be done at first i #

# morbidity model dims to define morbidity matrix (for N individuals)
num_indv_suscept_cols <- 6
tissue_dmg_cols <- 2 + (6 * 2) # need * 2 as need a last step tissue damage (dTx) and change in tissue damage (DTx)
disease_cond_cols <- 6

num.cols.morb <- 2 + num_indv_suscept_cols + tissue_dmg_cols + disease_cond_cols

# define morbidity matrix (skin disease states and blindness/VI)
all.morb.temp2 <- matrix(nrow=N, ncol=num.cols.morb)

all.morb.temp2[,1] <- cur.age
all.morb.temp2[,2] <- sex

all.morb.temp2[,10] <- 0 # assign Dix_T as 0 initially (this will then be the old) - severe itch
all.morb.temp2[,12] <- 0 # assign Dix_T as 0 initially (this will then be the old) - reactive skin
all.morb.temp2[,14] <- 0 # assign Dix_T as 0 initially (this will then be the old) - atrophy
all.morb.temp2[,16] <- 0 # assign Dix_T as 0 initially (this will then be the old) - hanging groin
all.morb.temp2[,18] <- 0 # assign Dix_T as 0 initially (this will then be the old) - mild dipigmentation
#all.morb.temp[,20] <- 0 # assign Dix_T as 0 initially (this will then be the old) - severe dipigmentation

all.morb.temp2[,22] <- 0 # assign disease condition initially - severe itch
all.morb.temp2[,23] <- 0 # assign disease condition initially - reactive skin
all.morb.temp2[,24] <- 0 # assign disease condition initially - atrophy
all.morb.temp2[,25] <- 0 # assign disease condition initially - hanging groin
all.morb.temp2[,26] <- 0 # assign disease condition initially - mild dipigmentation
all.morb.temp2[,28] <- 0 # assign disease condition initially - severe dipigmentation

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
SI_prev <- sum(all.morb.temp[,22])/N
RSD_prev <- sum(all.morb.temp[,23])/N
Atrp_prev <- sum(all.morb.temp[,24])/N
HG_prev <- sum(all.morb.temp[,25])/N
MDp_prev <- sum(all.morb.temp[,27])/N
SDp_prev <- sum(all.morb.temp[,28])/N

prev1_out <- c(SI_prev, RSD_prev, Atrp_prev, HG_prev, MDp_prev, SDp_prev)

SI_prev2 <- sum(all.morb.temp2[,22])/N
RSD_prev2 <- sum(all.morb.temp2[,23])/N
Atrp_prev2 <- sum(all.morb.temp2[,24])/N
HG_prev2 <- sum(all.morb.temp2[,25])/N
MDp_prev2 <- sum(all.morb.temp2[,27])/N
SDp_prev2 <- sum(all.morb.temp2[,28])/N

prev2_out <- c(SI_prev2, RSD_prev2, Atrp_prev2, HG_prev2, MDp_prev2, SDp_prev2)

prev1_out*100
prev2_out*100

test_condition_col <- sample(c(0,1), replace=TRUE, size=440)

all.morb.temp2 <- all.morb.temp2 # store


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

# repeat for all individ_suscpt cols




# # ================= #
# # function 1 needed #
# rotate <- function(x) {
#   x_rotated <- t(apply(x, 2, rev))
#
#   return(x_rotated)
# }
#
# # ================= #
# # function 2 needed #
# change.micro <- function(dat, num.comps, mf.cpt, num.mf.comps, ws, DT, time.each.comp, mu.rates.mf, fec.rates, mf.move.rate,
#                          up, kap, iteration, treat.vec, give.treat, treat.start)
#
# {
#   N <- length(dat[,1])
#
#   # indexes for fertile worms (to use in production of mf)
#   fert.worms.start <-  ws + num.comps * 2
#   fert.worms.end <-  (ws - 1) + num.comps * 3
#
#   # indexes to check if there are males (males start is just 'ws')
#   # there must be >= 1 male worm for females to produce microfilariae
#   mal.worms.end <- (ws - 1) + num.comps
#   mf.mu <- rep(mu.rates.mf[mf.cpt], N)
#   fert.worms <- dat[, fert.worms.start:fert.worms.end] #number of fertile females worms
#
#   # increases microfilarial mortality if treatment has started
#   if(give.treat == 1 & iteration >= treat.start)
#   {
#     tao <- ((iteration - 1) * DT) - treat.vec # tao is zero if treatment has been given at this timestep
#
#     mu.mf.prime <- ((tao + up) ^ (- kap)) # additional mortality due to ivermectin treatment
#
#     mu.mf.prime[which(is.na(mu.mf.prime) == TRUE)] <- 0
#
#     mf.mu <- mf.mu + mu.mf.prime
#
#   }
#
#   # if the first age class of microfilariae
#   if(mf.cpt == 1)
#   {
#     mp <- rep(0, N)
#
#     inds.fec <- which(rowSums(dat[, ws : mal.worms.end]) > 0); mp[inds.fec] <- 1 # need to check there is more than one male
#
#     k1 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = 0)  #fert worms and epin are vectors
#     k2 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k1/2)
#     k3 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k2/2)
#     k4 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k3)
#
#     out <- dat[, 6 + mf.cpt] + DT/6 * (k1 + 2 * k2 + 2* k3 + k4)
#
#   }
#
#   # if age class of microfilariae is >1
#   if(mf.cpt > 1)
#   {
#     k1 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = 0)
#     k2 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k1/2)
#     k3 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k2/2)
#     k4 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k3)
#
#     out <- dat[, 6 + mf.cpt] + DT/6 * (k1 + 2 * k2 + 2 * k3 + k4)
#
#   }
#
#   return(out)
# }
#
# # ================= #
# # function 3 needed #
#
# derivmf.one <- function(fert.worms, mf.in, ep.in, mf.mort, mf.move, mp, k.in)  #fert worms and epin are vectors
# {
#   new.in <- (rotate(fert.worms) * ep.in) #need to rotate matrix to each column is multiplied by respective fecundity rate, not each row
#   new.in <- rotate(rotate(rotate(new.in)))
#   new.in <- rowSums(new.in)
#
#   mort.temp <- mf.mort * (mf.in + k.in)
#   move.temp <- mf.move * (mf.in + k.in)
#
#   mort.temp[which(mort.temp < 0)] <- 0
#   move.temp [which(move.temp < 0)] <- 0
#
#   if(length(which(mf.mort * (mf.in + k.in) < 0)) > 0) {print('MF NEGATIVE1')}
#   if(length(which(mf.move * (mf.in + k.in) < 0)) > 0) {print('MF NEGATIVE2')}
#
#
#   out <- mp * new.in - mort.temp - move.temp
#
#   return(out)
# }
#
# # ================= #
# # function 4 needed #
#
# derivmf.rest <- function(mf.in, mf.mort, mf.move, mf.comp.minus.one, k.in)
# {
#   move.last <- mf.comp.minus.one * mf.move
#   mort.temp <- mf.mort * (mf.in + k.in)
#   move.temp <- mf.move * (mf.in + k.in)
#
#   move.last[which(move.last < 0)] <- 0
#   mort.temp[which(mort.temp < 0)] <- 0
#   move.temp [which(move.temp < 0)] <- 0
#
#   if(length(which(mf.mort * (mf.in + k.in) < 0)) > 0) {print('WARNING MF NEGATIVE3')}
#   if(length(which(mf.move * (mf.in + k.in) < 0)) > 0) {print('WARNING MF NEGATIVE4')}
#
#
#   out <- move.last - mort.temp - move.temp
#
#   return(out)
# }
#
# # ================================== #
# #         loop to run                #
#
# for(mf.c in 1 : num.mf.comps)
#
# {
#
#   res.mf <- change.micro(dat = all.mats.cur, num.comps =num.comps.worm, mf.cpt = mf.c,
#                          num.mf.comps = num.mf.comps, ws=worms.start, DT=DT, time.each.comp = time.each.comp.mf,
#                          mu.rates.mf = mort.rates.mf, fec.rates = fec.rates.worms, mf.move.rate = mf.move.rate, up = up, kap = kap, iteration = i,
#                          treat.vec = treat.vec.in, give.treat = give.treat, treat.start = treat.start)
#
#   all.mats.temp[, 6 + mf.c] <- res.mf
# }
