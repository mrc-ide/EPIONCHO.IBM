library(devtools)
devtools::load_all()
library(dplyr)

iter <- 1#as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
set.seed(iter + (iter*3758))

DT.in <- 1/366

# kE = 0.2
# delta.hz.in =  0.385,
# delta.hinf.in = 0.003,
# c.h.in = 0.008,
# gam.dis.in = 0.2,

# kE = 0.3
# delta.hz.in = 0.186,
# delta.hinf.in = 0.003,
# c.h.in = 0.005,
# gam.dis.in = 0.2,


# Gabon Elisa Dataset

# kE = 0.2
#ABR.in <- 73
#mda.val <- 0
# kE = 0.3
# ABR.in <- 176
# mda.val <- 0
#
# # Gabon RDT Dataset
# # kE = 0.2
# ABR.in <- 76
# mda.val <- 0
# # kE = 0.3
# ABR.in <- 180
# mda.val <- 0
ABR.in <- round(rgamma(1, 18.9, .0072)) # 70%

# For Keran
mda.val <- 26

vctr.control.strt <- 80
vctr.control.duration <- 31
vector.control.efficacies <- rep(c(.60, .75, .95), 3300)
vctr.control.efficacy <- vector.control.efficacies[iter]

# treat.strt.yrs = 1989
treat.len = mda.val; treat.strt.yrs = 93; yrs.post.treat = 10

treat.strt = treat.strt.yrs; treat.stp = treat.strt + treat.len
timesteps = treat.stp + yrs.post.treat #final duration
cstm_treat_params <- list(start_biannual=treat.strt.yrs+14, coverage_changes=c(treat.strt.yrs+7, treat.strt.yrs+14), coverage_change_values=c(0.55, 0.75, 0.85))


give.treat.in = 1; trt.int = 1

output <- ep.equi.sim(time.its = timesteps,
                      ABR = ABR.in,
                      treat.int = trt.int,
                      treat.prob = 0.80,
                      give.treat = give.treat.in,
                      treat.start = treat.strt,
                      treat.stop = treat.stp,
                      treat.timing = NA,
                      pnc = 0.01,
                      min.mont.age = 5,
                      vector.control.strt = vctr.control.strt,
                      vector.control.duration = vctr.control.duration,
                      vector.control.efficacy = vctr.control.efficacy,
                      #delta.hz.in =  0.385,
                      #delta.hinf.in = 0.003,
                      #c.h.in = 0.008,
                      #gam.dis.in = 0.2,
                      delta.hz.in = 0.186,
                      delta.hinf.in = 0.003,
                      c.h.in = 0.005,
                      gam.dis.in = 0.3,
                      #delta.hz.in = 0.118,
                      #delta.hinf.in = 0.002,
                      #c.h.in = 0.004,
                      #gam.dis.in = 0.4,
                      run_equilibrium = FALSE,
                      #equilibrium = output_equilibrium$all_equilibrium_outputs,
                      print_progress=TRUE,
                      calc_ov16 = TRUE,
                      #ov16_equilibrium = output_equilibrium$ov16_equilibrium,
                      no_prev_run=TRUE,
                      custom_treat_params=cstm_treat_params,
                      seroreversion=FALSE)

params <- list(mda.val, ABR.in)
names(params) <- c('MDA', 'ABR')
output <- append(output, params)

saveRDS(output, paste("/rds/general/user/ar722/home/ov16_test/ov16_output/ov16_any_worm_output",iter,".rds", sep=""))
