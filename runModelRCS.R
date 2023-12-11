library(dplyr)

iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
set.seed(iter + (iter*3758))

DT.in <- 1/366

kE = 0.2

#### Current file: runModelRCS.R
#iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
set.seed(parallel_iterator + (parallel_iterator*3758))

DT.in <- 1/366

if(kE == 0.2) {
  delta.hz.in.val =  0.385
  delta.hinf.in.val = 0.003
  c.h.in.val = 0.008
  gam.dis.in.val = 0.2
} else {
  delta.hz.in.val =  0.186
  delta.hinf.in.val = 0.003
  c.h.in.val = 0.005
  gam.dis.in.val = 0.3
}


# Gabon Elisa Dataset
if (kE == 0.2) {
  ABR.in <- 73
  mda.val <- 0
} else {
  ABR.in <- 176
  mda.val <- 0
}

# try 13/14 as well
treat.len = mda.val; treat.strt.yrs = 100; yrs.post.treat = 5

treat.strt = treat.strt.yrs; treat.stp = treat.strt + treat.len
timesteps = treat.stp + yrs.post.treat #final duration

give.treat.in = 0; trt.int = 1

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
                      vector.control.strt = NA,
                      delta.hz.in =  0.385,
                      delta.hinf.in = 0.003,
                      c.h.in = 0.008,
                      gam.dis.in = 0.2,
                      run_equilibrium = FALSE,
                      #equilibrium = output_equilibrium$all_equilibrium_outputs,
                      print_progress=TRUE,
                      calc_ov16 = TRUE,
                      #ov16_equilibrium = output_equilibrium$ov16_equilibrium,
                      no_prev_run=TRUE)

params <- list(mda.val, ABR.in)
names(params) <- c('MDA', 'ABR')
output <- append(output, params)

saveRDS(output, paste("/rds/general/user/ar722/home/ov16_test/ov16_output/ov16_any_worm_output",iter,".rds", sep=""))

