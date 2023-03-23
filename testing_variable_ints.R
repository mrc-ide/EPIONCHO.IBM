# ========================================================================================================= #
# current code for simple intervention: treat start ( in years), treat stop and interval (e.g., 1 = 1 year) #

DT.t <- 1/366
#time.its <- round(time.its / (DT))
#year_its <- seq(0, time.its, 366)

#length of simulation in years
timesteps.t = 80

#should treatment be given, when and how often
give.treat.t = 0
treat.start.t = 1; treat.stop.t = 16
treat.int.t = 1


#if(give.treat == 1)
#{
treat.stop.t <- round(treat.stop.t / (DT.t))
treat.start.t <-  round( (treat.start.t) / (DT.t)) + 1
#}



times.of.treat.in.t <- seq(treat.start.t, treat.stop.t - (treat.int.t / DT.t), treat.int.t / DT.t)
times.of.treat.in.t <- seq(treat.start.t, treat.stop.t, treat.int.t / DT.t)

# ============================================================== #
# additional code to specify more complex intervention histories #

treat.start.t = 31; treat.stop.t = 46

treat.timing.in.t <- c(1,2,3,4,5,6,7,8,9,10,10.5,11,11.5,12, 13, 14, 15)

treat.timing.in.t2 <- treat.timing.in.t + (treat.start.t - 1)

times.of.treat.in.t2 <- round((treat.timing.in.t2) / (DT.t)) + 1


# =============================================================== #
#             Test in model  (variable timing/frequency)          #

devtools::load_all()

# ============================ #
#length of simulation in years
timesteps = 47

#should treatment be given, when and how often
give.treat.in = 1
treat.strt = 31; treat.stp = 46
#trt.int = 1
treat.timing.in <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10.5, 11, 11.5, 12, 13, 14, 15)

# =========================== #
# specify large gap in middle #
timesteps = 47

#should treatment be given, when and how often
give.treat.in = 1
treat.strt = 31; treat.stp = 46
#trt.int = 1
treat.timing.in <- c(1, 2, 3, 4, 10, 11, 11.5, 12, 13, 14, 15)


#annual biting rate, which determines infection prevalence
ABR.in = 10000

output_MDAvariablefreq2 <- ep.equi.sim(time.its = timesteps,
                                  ABR = ABR.in,
                                  N.in = 440,
                                  treat.timing = treat.timing.in,
                                  treat.prob = 0.65,
                                  give.treat = give.treat.in,
                                  treat.start = treat.strt,
                                  treat.stop = treat.stp,
                                  pnc = 0.05,
                                  min.mont.age = 5,
                                  vector.control.strt = NA,
                                  delta.hz.in = 0.186,
                                  delta.hinf.in = 0.003,
                                  c.h.in = 0.005,
                                  gam.dis.in = 0.3,
                                  run_equilibrium = TRUE,
                                  equilibrium = NA,
                                  print_progress = TRUE)

names(output_MDAvariablefreq2)


tme <- seq(1, 47*366-1)/366

plot(tme, output_MDAvariablefreq2$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))
plot(tme, output_MDAvariablefreq2$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1), xlim = c(30,47), xaxt = "n")
axis(1, at = seq(1, 47, by = 1), las=2)


# ========================================================================= #
#             Test in model  (both VC & variable timing/frequency)          #

devtools::load_all()

# ============================ #
#length of simulation in years
timesteps = 81

#should treatment be given, when and how often
give.treat.in = 1
treat.strt = 61; treat.stp = 80
#trt.int = 1
treat.timing.in <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20)

#annual biting rate, which determines infection prevalence
ABR.in = 1000

output_MDAvariablefreq_VC <- ep.equi.sim(time.its = timesteps,
                                       ABR = ABR.in,
                                       N.in = 440,
                                       treat.timing = treat.timing.in,
                                       treat.prob = 0.65,
                                       give.treat = give.treat.in,
                                       treat.start = treat.strt,
                                       treat.stop = treat.stp,
                                       pnc = 0.05,
                                       min.mont.age = 5,
                                       #vector.control.strt = NA,
                                       vector.control.strt = 40,
                                       vector.control.duration = 30,
                                       vector.control.efficacy = 0.7,
                                       delta.hz.in = 0.186,
                                       delta.hinf.in = 0.003,
                                       c.h.in = 0.005,
                                       gam.dis.in = 0.3,
                                       run_equilibrium = TRUE,
                                       equilibrium = NA,
                                       print_progress = TRUE)

names(output_MDAvariablefreq_VC)


tme <- seq(1, 81*366-1)/366

plot(tme, output_MDAvariablefreq_VC$ABR_recorded, type = 'l', xlab = 'time (years)', ylab = 'ABR', ylim = c(0, 1000))
plot(tme, output_MDAvariablefreq_VC$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))
plot(tme, output_MDAvariablefreq_VC$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1), xlim = c(40,81), xaxt = "n")
axis(1, at = seq(1, 81, by = 1), las=2)

