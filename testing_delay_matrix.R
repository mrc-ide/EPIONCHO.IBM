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
num.days.cols <- 3
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

# ================================================================================================================= #
#                         Testing matrix for lagged blindness (07 - 07 -23)                                         #
# ================================================================================================================= #

N <- 440
cur.age <- rep(0, N)
DT <- 1/366
mean.age <- 50
real.max.age <- 80

time.its <- 80 # same as timesteps e.g., 80 yrs
time.its <- round(time.its / (DT))
year_its <- seq(0, time.its, 366)

# get number of worms from iter = 1 and and endemic eq run

eq_outputs <- output_equilibrium_morbidity$all_equilibrium_outputs[[1]]

morb.eq.outputs <- output_equilibrium_morbidity$morbidity.outputs[[1]]

# eye disease matrix (specified at i = 1)

num.cols.blindness <- 9

# define Morbidity matrix #
all.blind.temp <- matrix(nrow = N, ncol = num.cols.blindness)

all.blind.temp[,1] <- eq_outputs[,2] # current age
all.blind.temp[,2] <- eq_outputs[,3] # sex
all.blind.temp[,3] <- eq_outputs[,1] # compliance
all.blind.temp[,4] <- 0 # current TRUE mf count
all.blind.temp[,5] <- 0 # current mf skin snip count

# set other cols to 0 at i = 1
all.blind.temp[,c(6:9)] <- 0

# ================================================================================= #
# set up 'lagged' matrix for recording number of blind positive 2 yrs after 'trial' #

lag.iters <- (366 * 2) - 1 # 2 yrs - timestep

n.iter <- time.its

ncol.lag.mat <- 2 + lag.iters + n.iter # col 1 = curr.age, col 2 = sex

blind.lag.mat <- matrix(nrow = N, ncol = ncol.lag.mat) # large size?

blind.lag.mat[,c(3:733)] <- 0 # these should all be 0 up to 2 years (iter = 732)

# at iter = 1

blind.lag.mat[,1] <- eq_outputs[,2] # current age
blind.lag.mat[,2] <- eq_outputs[,3] # sex

i <- 1
col.to.select <- 2 + i + 731 # 2 + because first two columsn are current age and sex, so iter one from thrid col

blind.lag.mat[,col.to.select] <- morb.eq.outputs[,51] # update lagged iter with current disease state

# prev cals #
lagged.mat.tmp <- blind.lag.mat

col_slct <- i + 2 # iter 1 begins from col 3 in lagged matrix
blind_prev_temp <- length(which(lagged.mat.tmp[,col_slct] == 1 & lagged.mat.tmp[,1] >= 5)) /  length(which(lagged.mat.tmp[,1] >= 5)) # prev in > 5yrs

#blind_prev_temp <- 0.07
visual_imp_prev_temp <- blind_prev_temp * 4/3

# blind.lag.mat[,c(732:ncol(blind.lag.mat))] <- morb.eq.outputs[,51]

# for iter 1, do prevalence calcs etc based on column 3, iter 2 do prev calcs based on column 4 etc etc
# for iter 1, disease state updated in column 734 (2 + 732 iters), as disease state at iter 1 relevant for 2 yrs later

# at iter = 2

blind.lag.mat[,1] <- eq_outputs[,2] + 1/366 # current age
blind.lag.mat[,2] <- eq_outputs[,3] # sex

i <- 2
col.to.select <- 2 + i + 731 # 2 + because first two columsn are current age and sex, so iter one from thrid col

blind.lag.mat[,col.to.select] <- morb.eq.outputs[,51] # update lagged iter with current disease state
