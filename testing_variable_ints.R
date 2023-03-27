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
treat.start.t = 1/366


#if(give.treat == 1)
#{
treat.stop.t <- round(treat.stop.t / (DT.t))
treat.start.t <-  round( (treat.start.t) / (DT.t)) + 1

treat.start.t <-  round( (treat.start) / (DT)) + 1
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
timesteps = 45

#should treatment be given, when and how often
give.treat.in = 1
treat.strt = 31; treat.stp = 45
treat.strt = 30; treat.stp = 44
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


tme <- seq(1, 45*366-1)/366

plot(tme, output_MDAvariablefreq2$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))
plot(tme, output_MDAvariablefreq2$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1), xlim = c(30,47), xaxt = "n")
axis(1, at = seq(1, 45, by = 1), las=2)


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
ABR.in = 10000

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
                                       vector.control.duration = 20,
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

plot(tme, output_MDAvariablefreq_VC$ABR_recorded, type = 'l', xlab = 'time (years)', ylab = 'ABR', ylim = c(0, 10000))
plot(tme, output_MDAvariablefreq_VC$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))
plot(tme, output_MDAvariablefreq_VC$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1), xlim = c(40,81), xaxt = "n")
axis(1, at = seq(1, 81, by = 1), las=2)



# ================================================================================================================================== #
#                             Variable coverage (vector of coverages)                                                                #

test.function.extracting.cov <- function(i, treat.prob.variable.t, times.of.treat.in.t) {

if(all(!is.na(treat.prob.variable.t))){

  if(any(i == times.of.treat.in.t)) {
    index.iter.treat <- match(i, times.of.treat.in.t)
    treat.prob.out <- treat.prob.variable.t[index.iter.treat]} # find element where iteration number matches a time in times.of.treat vector

  #if(any(i == times.of.treat.in.t)) {treat.prob.out <- treat.prob.variable.t[index.iter.treat]}

  return(treat.prob.out)

  }

}

treat.prob.variable.t <- c(0.65, 0.7, 0.8, 0.65)

times.of.treat.in.t <- c(2, 7, 20, 100)

test.function.extracting.cov(i = 5, treat.prob.variable.t = treat.prob.variable.t, times.of.treat.in.t = times.of.treat.in.t)

# ================== #
# test in main model #

devtools::load_all()

timesteps = 45

#should treatment be given, when and how often
# give.treat.in = 1
# treat.strt = 3; treat.stp = 15
# #trt.int = 1
# treat.timing.in <- c(3, 6, 9) # not years, this will instead directly infer the
#                               # (normally this would be years when MDA occur then converetd)
#                               # so this is directly the iters rather than the years (ONLY FOR TESTING!)
#
# treat.prob.variable.in <- c(0.65, 0.75, 0.8)


#should treatment be given, when and how often
give.treat.in = 1
treat.strt = 31; treat.stp = 45

treat.timing.in <- c(1, 2, 3, 4, 10, 11, 11.5, 12, 13, 14, 15)
treat.prob.variable.in <- c(0.65, 0.75, 0.8, 0.85, 0.9, 0.85, 0.5, 0.65, 0.9, 0.8, 0.6)


timesteps = 52
treat.strt = 31
treat.stp = 51 # note this must be year after last treatment year in vector of treatment round timings,
               # so if last one is at 20 years, last year of treatment is 31 + (20 - 1) = year 50.
               # therefore treat.stop needs to be at year 51.
treat.timing.in <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20) # last round at 20 year will be 50 years in above
treat.prob.variable.in <- c(0.65, 0.75, 0.8, 0.85, 0.9, 0.85, 0.5, 0.65, 0.9, 0.8, 0.6, 0.95, 0.95, 0.9, 0.8, 0.5, 0.55, 0.85, 0.9)

#annual biting rate, which determines infection prevalence
ABR.in = 10000

output_MDAvarCov <- ep.equi.sim(time.its = timesteps,
                                       ABR = ABR.in,
                                       N.in = 440,
                                       treat.timing = treat.timing.in,
                                       treat.prob.variable = treat.prob.variable.in,
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

names(output_MDAvarCov)


tme <- seq(1, 52*366-1)/366

plot(tme, output_MDAvarCov$coverage.recorded, type = 'l', xlab = 'time (years)', ylab = 'coverage (0-1)', ylim = c(0, 1), xaxt = "n")
axis(1, at = seq(1, 52, by = 1), las=2)

plot(tme, output_MDAvarCov$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))
plot(tme, output_MDAvarCov$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1), xlim = c(30,51), xaxt = "n")
axis(1, at = seq(1, 52, by = 1), las=2)


# =============================== #
#    testing fre, cov & VC        #

timesteps = 52
treat.strt = 31
treat.stp = 51 # note this must be year after last treatment year in vector of treatment round timings,
# so if last one is at 20 years, last year of treatment is 31 + (20 - 1) = year 50.
# therefore treat.stop needs to be at year 51.
treat.timing.in <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20) # last round at 20 year will be 50 years in above
treat.prob.variable.in <- c(0.65, 0.75, 0.8, 0.85, 0.9, 0.85, 0.5, 0.65, 0.9, 0.8, 0.6, 0.95, 0.95, 0.9, 0.8, 0.5, 0.55, 0.85, 0.9)

#annual biting rate, which determines infection prevalence
ABR.in = 10000

output_MDA_freqcovVC <- ep.equi.sim(time.its = timesteps,
                                ABR = ABR.in,
                                N.in = 440,
                                treat.timing = treat.timing.in,
                                treat.prob.variable = treat.prob.variable.in,
                                give.treat = give.treat.in,
                                treat.start = treat.strt,
                                treat.stop = treat.stp,
                                pnc = 0.05,
                                min.mont.age = 5,
                                #vector.control.strt = NA,
                                vector.control.strt = 30,
                                vector.control.duration = 10,
                                vector.control.efficacy = 0.7,
                                delta.hz.in = 0.186,
                                delta.hinf.in = 0.003,
                                c.h.in = 0.005,
                                gam.dis.in = 0.3,
                                run_equilibrium = TRUE,
                                equilibrium = NA,
                                print_progress = TRUE)

names(output_MDA_freqcovVC)


tme <- seq(1, 52*366-1)/366

plot(tme, output_MDA_freqcovVC$coverage.recorded, type = 'l', xlab = 'time (years)', ylab = 'coverage (0-1)', ylim = c(0, 1), xaxt = "n")
axis(1, at = seq(1, 52, by = 1), las=2)
