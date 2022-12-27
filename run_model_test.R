# ====================================================================================================================#
#                       Compiling/loading package, testing etc                                                        #

devtools::has_devel()

library(devtools)

# use_mit_license("EPIONCHO.IBM")

devtools::load_all()

devtools::check()

devtools::install()

library("EPIONCHO.IBM") # when updated/ to install locally


# ==================================================== #
DT.in <- 1/366 #timestep must be one day

treat.len <- 1 #treatment duration in years
treat.len <- 25 #treatment duration in years

treat.strt  = round(1 / (DT.in )); treat.stp = treat.strt + round(treat.len / (DT.in )) #treatment start and stop
treat.strt  = round(80 / (DT.in )); treat.stp = treat.strt + round(treat.len / (DT.in )) #treatment start and stop

timesteps = treat.stp + round(1 / (DT.in )) #final duration (3 is number of years after treatment stops to continue running model)
timesteps = treat.stp + round(10 / (DT.in )) #final duration (3 is number of years after treatment stops to continue running model)

gv.trt = 1
trt.int = 1 #treatment interval (years, 0.5 gives biannual)


ABR.in <- 647 #annual biting rate for mf 50% prevalence with k_E 0.2

output <-  ep.equi.sim(time.its = timesteps,
                       ABR = ABR.in,
                       DT = DT.in,
                       treat.int = trt.int,
                       treat.prob = 80,
                       give.treat = gv.trt,
                       treat.start = treat.strt,
                       treat.stop = treat.stp,
                       pnc = 0.01, min.mont.age = 5,
                       delta.hz.in = 0.385,
                       delta.hinf.in = 0.003,
                       c.h.in = 0.008,
                       gam.dis.in = 0.2)
