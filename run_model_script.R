# ====================================================================================================================#
#                       Compiling/loading package, testing etc                                                        #

# ==================== #
#    install first     #

#install.packages("remotes")
#remotes::install_github("Mad206/EPIONCHO.IBM") # need to change from Mad206 to mrc-ide

# load locally #

library(EPIONCHO.IBM) # when updated/ to install locally

# ==================================================== #
#     baseline model run                               #
# ==================================================== #

DT.in <- 1/366 #timestep must be one day

treat.len <- 0 #treatment duration in years

treat.strt  = round(0 / (DT.in )); treat.stp = treat.strt + round(treat.len / (DT.in )) #treatment start and stop

timesteps = treat.stp + round(80 / (DT.in )) #final duration (3 is number of years after treatment stops to continue running model)

gv.trt = 0
trt.int = 0 #treatment interval (years, 0.5 gives biannual)


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

# original parameters
year_seq <- seq(from= 1/366, to= 79.99727, by = 1/366)
prev_out <- output[[2]]

plot(x = year_seq, y = prev_out)
lines(year_seq, prev_out, col="red", lwd=2)


# ==================================================== #
#     model run with treatment                         #
# ==================================================== #

DT.in <- 1/366 #timestep must be one day

treat.len <- 25 #treatment duration in years

treat.strt  = round(80 / (DT.in )); treat.stp = treat.strt + round(treat.len / (DT.in )) #treatment start and stop

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

# original parameters
year_seq <- seq(from= 1/366, to= 114.99727, by = 1/366)
prev_out <- output[[2]]

plot(x = year_seq, y = prev_out)
lines(year_seq, prev_out, col="red", lwd=2)

# ==================================================== #
#     epilepsy module run to equilibrium               #
# ==================================================== #

DT.in <- 1/366 #timestep must be one day


# ============== #
# no treatment #

treat.len <- 0 #treatment duration in years
#treat.strt  = round(200 / (DT.in )); treat.stp = treat.strt + round(treat.len / (DT.in )) #treatment start and stop
treat.strt  = round(0 / (DT.in )); treat.stp = treat.strt + round(treat.len / (DT.in )) #treatment start and stop
timesteps = treat.stp + round(200 / (DT.in )) #final duration
gv.trt = 0
trt.int = 0

# ============== #
# with treatment #

treat.len <- 25 #treatment duration in years
#treat.strt  = round(200 / (DT.in )); treat.stp = treat.strt + round(treat.len / (DT.in )) #treatment start and stop
treat.strt  = round(200 / (DT.in )); treat.stp = treat.strt + round(treat.len / (DT.in )) #treatment start and stop
timesteps = treat.stp + round(1 / (DT.in )) #final duration
gv.trt = 1
trt.int = 1 #treatment interval


ABR.in <- 41922 #annual biting rate (from Jacob's report)


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
                       gam.dis.in = 0.2,
                       epilepsy_module = "YES") # need to specify epilepsy module run here



prev_out <- output[[2]]
OAE_prev_out <- output[[6]]
OAE_incid_out <- output[[8]]

year_seq <- seq(from= 1/366, to= 200 - 1/366, by = 1/366) # change from 200 to 226 if treatment run
plot(x = year_seq, y = prev_out)
lines(year_seq, prev_out, col="red", lwd=2)

year_seq2 <- seq(from= 1/366, to= 200, by = 1/366) # change from 200 to 226 if treatment run
plot(x = year_seq2, y = OAE_prev_out, ylim=c(0,0.1))
plot(x = year_seq2, y = OAE_prev_out, ylim=c(0,1))
lines(year_seq2, OAE_prev_out, col="red", lwd=2)

year_seq_incidence <- seq(from= 1/366, to= 199.9973, by = 1/366) # change from 200 to 225.9973 if treatment run
plot(x = year_seq_incidence, y = OAE_incid_out)
lines(year_seq_incidence, OAE_incid_out, col="red", lwd=2)

