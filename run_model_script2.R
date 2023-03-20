# ========================================================================================================================== #
#               Testing new approach to modelling blindness (same structure as Hugo's e.g., post-hoc analysis)               #
# ========================================================================================================================== #

devtools::load_all()


# ==================================================== #
DT.in <- 1/366 #timestep must be one day

treat.len <- 1 #treatment duration in years
#treat.len <- 25 #treatment duration in years

treat.strt  = round(1 / (DT.in )); treat.stp = treat.strt + round(treat.len / (DT.in )) #treatment start and stop
treat.strt  = round(40 / (DT.in )); treat.stp = treat.strt + round(treat.len / (DT.in )) #treatment start and stop

timesteps = treat.stp + round(1 / (DT.in )) #final duration (3 is number of years after treatment stops to continue running model)
#timesteps = treat.stp + round(10 / (DT.in )) #final duration (3 is number of years after treatment stops to continue running model)

gv.trt = 1
gv.trt = 0
trt.int = 1 #treatment interval (years, 0.5 gives biannual)


ABR.in <- 1000 #annual biting rate for mf 50% prevalence with k_E 0.2

output <-  ep.equi.sim(time.its = timesteps,
                       ABR = ABR.in,
                       DT = DT.in,
                       treat.int = trt.int,
                       treat.prob = 80,
                       give.treat = gv.trt,
                       treat.start = treat.strt,
                       treat.stop = treat.stp,
                       pnc = 0.01, min.mont.age = 5,
                       delta.hz.in = 0.186,
                       delta.hinf.in = 0.003,
                       c.h.in = 0.005,
                       gam.dis.in = 0.3)

# original parameters
year_seq <- seq(from= 1/366, to= 114.99727, by = 1/366)
prev_out <- output[[2]]

plot(x = year_seq, y = prev_out)
lines(year_seq, prev_out, col="red", lwd=2)


# extract strata prev outputs (dataframes)

male_prev_df <- output[[5]][[1]]

female_prev_df <- output[[5]][[2]]

plot(x = year_seq, y = male_prev_df$X0.5)
lines(year_seq, male_prev_df$X0.5, col="red", lwd=2)
lines(year_seq, male_prev_df$X5.10, col="blue", lwd=2)
lines(year_seq, male_prev_df$X10.15, col="green", lwd=2)
lines(year_seq, male_prev_df$X15.20, col="purple", lwd=2)
lines(year_seq, male_prev_df$X20.25, col="brown", lwd=2)
lines(year_seq, male_prev_df$X25.30, col="pink", lwd=2)
lines(year_seq, male_prev_df$X30.40, col="yellow", lwd=2)
lines(year_seq, male_prev_df$X60.80, col="grey", lwd=2)

male_df <- output[[5]]

male_df$time <- time_vec_test
