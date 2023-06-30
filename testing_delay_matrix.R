# =============================================================================================================== #
#                      testing delay matrix for remaining sequelae positive for 5 or 7 days                       #
#                                    28 - 06 - 23                                                                 #
# =============================================================================================================== #

# ========================== #
# For main morbidity matrix #
# things already present and needed for matrix in model_wrappers:
N <- 440
cur.age <- rep(0, N)
DT <- 1/366
mean.age <- 50
real.max.age <- 80

# get number of worms from iter = 1 and and endemic eq run

eq_outputs <- output_equilibrium_morbidity$all_equilibrium_outputs[[1]]

morb.eq.outputs <- output_equilibrium_morbidity$morbidity.matrix

# morb.eq.outputs[,46] # SI state
# morb.eq.outputs[,50] # RSD state

morb.mat.tmp <- morb.eq.outputs

# ==================== #
# set- up delay matrix #
N <- 440

# ======== #
# at i = 1 #
num.days.cols <- 7
sequela.postive.mat <- matrix(0, ncol = num.days.cols, nrow = N)
inds.sequela.mat <- seq(2,(length(sequela.postive.mat[1,]))) # for moving columns along with time

# =============== #
# i > 1 each loop #

# for SI #
# steps to go through to update states #


# 2) whether should be tested (old step updated)
morb.mat.tmp[,30] <- ifelse(morb.mat.tmp[,4] > 0 & morb.mat.tmp[,46] == 0, 1, 0) # whether to be tested: severe itch

# 3) do the Bernoulli trial

# after updating the states in the morbidity matrix #

#  4) extract current sequela state for reversible conditions
sequela.postive.mat[,inds.sequela.mat] <- sequela.postive.mat[,(inds.sequela.mat-1)] # move sequela state along one col
sequela.postive.mat[,1] <- morb.eq.outputs[,46] # update first col with current sequela state on that day

# new steps : update current sequela state in morb.mat if 7th day of morbidity to 0 #
morb.mat.tmp[,22] <- sequela.postive.mat[,7] # NEW STEP: assign day 7 from sequela delay matrix
morb.mat.tmp[,46] <- ifelse(morb.mat.tmp[,22] == 1, 0, morb.mat.tmp[,46]) # NEW STEP: update current disease status

# need to update the sequela matrix = reset any individuals that have reached sequela on day 7 as 0 across row
current_day7_sequela <- which(sequela.postive.mat[,7] == 1)

if(length(current_day7_sequela) > 0) {

  sequela.postive.mat[current_day7_sequela, ] <- 0

}

to_die <- c(12)

sequela.postive.mat[to_die, ] <- 0

# morb.mat.tmp[1,46] = 1
# sequela.postive.mat[1,1] = 1
