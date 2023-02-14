# ====================================================================================================================#
#                       Compiling/loading package, testing etc                                                        #

# ==================== #
#    install first     #

#install.packages("remotes")
#remotes::install_github("Mad206/EPIONCHO.IBM") # need to change from Mad206 to mrc-ide

devtools::load_all()

# load locally #

library(EPIONCHO.IBM) # when updated/ to install locally


devtools::build_rmd("vignettes/Running_EPIONCHO_IBM.Rmd")

# ==================================================== #
#     baseline model run                               #
# ==================================================== #

#length of simulation in years
timesteps = 30

#should treatment be given, when and how often
give.treat.in = 0
treat.strt = 1; treat.stp = 16
trt.int = 1

#annual biting rate, which determines infection prevalence
ABR.in = 1082

output_equilibrium <- ep.equi.sim(time.its = timesteps,
                                  ABR = ABR.in,
                                  N.in = 440,
                                  treat.int = trt.int,
                                  treat.prob = 0.65,
                                  give.treat = give.treat.in,
                                  treat.start = treat.strt,
                                  treat.stop = treat.stp,
                                  pnc = 0.05,
                                  min.mont.age = 5,
                                  delta.hz.in = 0.186,
                                  delta.hinf.in = 0.003,
                                  c.h.in = 0.005,
                                  gam.dis.in = 0.3,
                                  run_equilibrium = TRUE,
                                  equilibrium,
                                  print_progress = TRUE)

names(output_equilibrium)

tme <- seq(1, 30*366-1)/366

plot(tme, output_equilibrium$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))

plot(tme, output_equilibrium$mf_intens, type = 'l', xlab = 'time (years)', ylab = 'mean microfilarial intensity')

plot(tme, output_equilibrium$L3, type = 'l', xlab = 'time (years)', ylab = 'mean L3 in fly population')

# ======================================== #
# simulating MDA from endemic equilibirium #

timesteps = 20
treat.strt = 1; treat.stp = 16 #if treatment stops in year 26, the last round is at the beginning of year 15
gv.trt = 1
trt.int = 1

timesteps = 30
treat.strt = 1; treat.stp = 26 #if treatment stops in year 26, the last round is at the beginning of year 15
gv.trt = 1
trt.int = 1

#treat.strt = 1/366

output_treat_annual <- ep.equi.sim(time.its = timesteps,
                                   ABR = ABR.in,
                                   N.in = 440,
                                   treat.int = trt.int,
                                   treat.prob = 0.8,
                                   give.treat = gv.trt,
                                   treat.start = treat.strt,
                                   treat.stop = treat.stp,
                                   pnc = 0.01,
                                   min.mont.age = 5,
                                   delta.hz.in = 0.186,
                                   delta.hinf.in = 0.003,
                                   c.h.in = 0.005,
                                   gam.dis.in = 0.3,
                                   run_equilibrium = FALSE,
                                   equilibrium = output_equilibrium[[4]],
                                   print_progress = TRUE)

tme <- seq(0, 20*366-1)/366
tme <- seq(0, 30*366-1)/366

plot(tme, output_treat_annual$mf_prev, type = 'l', xlab = 'time', ylab = 'microfilarial prevalence', ylim = c(0, 1))

abline(v = seq(1, 25), col = 'grey', lwd = 0.1)

abline(v = seq(1, 25, 0.5), col = 'grey', lwd = 0.1)

plot(tme, output_treat_annual$mf_intens, type = 'l', xlab = 'time', ylab = 'mean microfilarial intensity')

abline(v = seq(0, 20), col = 'black', lwd = 0.1)

# ================================================= #
# simulating BIANNUAL MDA from endemic equilibirium #

timesteps = 20
treat.strt = 1; treat.stp = 16 #if treatment stops in year 16, the last round is at the beginning of year 15
gv.trt = 1
trt.int = 0.5


output_treat_biannual <- ep.equi.sim(time.its = timesteps,
                                     ABR = ABR.in,
                                     N.in = 440,
                                     treat.int = trt.int,
                                     treat.prob = 0.65,
                                     give.treat = gv.trt,
                                     treat.start = treat.strt,
                                     treat.stop = treat.stp,
                                     pnc = 0.05,
                                     min.mont.age = 5,
                                     delta.hz.in = 0.385,
                                     delta.hinf.in = 0.003,
                                     c.h.in = 0.008,
                                     gam.dis.in = 0.2,
                                     run_equilibrium = FALSE,
                                     equilibrium = output_equilibrium[[4]],
                                     print_progress = FALSE)


tme <- seq(0, 20*366-1)/366

plot(tme, output_treat_biannual$mf_prev, type = 'l', xlab = 'time', ylab = 'microfilarial prevalence', ylim = c(0, 1))

abline(v = seq(0, 20), col = 'black', lwd = 0.1)

plot(tme, output_treat_biannual$mf_intens, type = 'l', xlab = 'time', ylab = 'mean microfilarial intensity')

abline(v = seq(0, 20), col = 'black', lwd = 0.1)
