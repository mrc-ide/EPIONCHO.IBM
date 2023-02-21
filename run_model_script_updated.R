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


#=====================================================================================================================================#
#                           EPILEPSY MODEL WITH NEW STRUCTURE                                                                         #
# =================================================================================================================================== #

# =========================================== #
# with equilbiirum option selected - 20/02/23 #

#length of simulation in years
timesteps = 200

#should treatment be given, when and how often
give.treat.in = 0
treat.strt = 1; treat.stp = 16
trt.int = 1

#annual biting rate, which determines infection prevalence
ABR.in = 41922

output_equilibrium_OAE <- ep.equi.sim(time.its = timesteps,
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
                                      print_progress = TRUE,
                                      epilepsy_module = "YES")

names(output_equilibrium_OAE)

tme <- seq(1, 200*366-1)/366

plot(tme, output_equilibrium_OAE$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))

tme2 <- seq(1, 200*366)/366
tme3 <- seq(1, 200*366)

plot(tme2, output_equilibrium_OAE[[5]], type = 'l', xlab = 'time (years)', ylab = 'OAE prevalence', ylim = c(0, 1))
plot(tme3, output_equilibrium_OAE[[5]], type = 'l', xlab = 'time (days)', ylab = 'OAE prevalence', ylim = c(0, 0.1), xlim = c(25000,80000))
plot(tme3, output_equilibrium_OAE[[5]], type = 'l', xlab = 'time (days)', ylab = 'OAE prevalence', ylim = c(0.02, 0.06), xlim = c(25000,80000))

tail(output_equilibrium_OAE[[5]])


# =========================================== #
# supplying equilbiirum option     - 20/02/23 #

timesteps = 26
treat.strt = 1; treat.stp = 26 #if treatment stops in year 26, the last round is at the beginning of year 15
gv.trt = 1
trt.int = 1

#treat.strt = 1/366

# make object with OAE infection inputs (if wish to take directly)
# OAE_input <- list(tail(output_equilibrium_OAE$OAE_prev,1), tail(output_equilibrium_OAE$OAE_incidence,1),
#                   tail(output_equilibrium_OAE$OAE_incidence_3_5yrs, 1), tail(output_equilibrium_OAE$OAE_incidence_5_10yrs, 1),
#                   tail(output_equilibrium_OAE$OAE_incidence_males, 1), tail(output_equilibrium_OAE$OAE_incidence_females,1))

output_treat_annual_OAE <- ep.equi.sim(time.its = timesteps,
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
                                   equilibrium = output_equilibrium_OAE[[4]],
                                   print_progress = TRUE,
                                   epilepsy_module = "YES",
                                   OAE_equilibrium = output_equilibrium_OAE[[12]])

# plot
tme <- seq(0, 26*366-1)/366

plot(tme, output_treat_annual_OAE$mf_prev, type = 'l', xlab = 'time', ylab = 'microfilarial prevalence', ylim = c(0, 1))

abline(v = seq(1, 25), col = 'grey', lwd = 0.1)

abline(v = seq(1, 25, 0.5), col = 'grey', lwd = 0.1)

plot(tme, output_treat_annual_OAE$mf_intens, type = 'l', xlab = 'time', ylab = 'mean microfilarial intensity')

abline(v = seq(0, 20), col = 'black', lwd = 0.1)

# OAE
tme3 <- seq(1, 26*366)

plot(tme, output_treat_annual_OAE$OAE_prev, type = 'l', xlab = 'time', ylab = 'OAE prevalence', ylim = c(0, 0.1))

abline(v = seq(1, 25), col = 'grey', lwd = 0.1)

head(output_treat_annual_OAE$OAE_prev)
tail(output_treat_annual_OAE$OAE_prev)

(head(output_treat_annual_OAE$OAE_prev, 1) - tail(output_treat_annual_OAE$OAE_prev, 1)) / head(output_treat_annual_OAE$OAE_prev, 1)

# ============================================================================ #
# new structure but without using the equilibirum input option - OLD: 17/02/23 #

output_equilibrium_OAE <- ep.equi.sim(time.its = timesteps,
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
                                  print_progress = TRUE,
                                  epilepsy_module = "YES")

names(output_equilibrium_OAE)

output_equilibrium_OAE$mf_prev

output_equilibrium_OAE$all_equilibrium_outputs[[1]] # this is all_mats_temp in
# OAE_prev_out <- output_equilibrium_OAE[[2]]
# OAE_prev_out <- output_equilibrium_OAE[[6]]
# OAE_incid_out <- output_equilibrium_OAE[[8]]

tme <- seq(1, 200*366-1)/366

plot(tme, output_equilibrium_OAE[[2]], type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))

tme2 <- seq(1, 200*366)/366
tme3 <- seq(1, 200*366)

plot(tme2, output_equilibrium_OAE[[6]], type = 'l', xlab = 'time (years)', ylab = 'OAE prevalence', ylim = c(0, 1))
plot(tme3, output_equilibrium_OAE[[6]], type = 'l', xlab = 'time (days)', ylab = 'OAE prevalence', ylim = c(0, 0.1), xlim = c(25000,80000))



# ============================================ #
#   With MDA (& running equilibirium)          #

#length of simulation in years
timesteps = 235

#should treatment be given, when and how often
give.treat.in = 1
treat.strt = 200; treat.stp = 225
trt.int = 1

#annual biting rate, which determines infection prevalence
ABR.in = 41922

output_equilibrium_OAE_MDA <- ep.equi.sim(time.its = timesteps,
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
                                      print_progress = TRUE,
                                      epilepsy_module = "YES")

names(output_equilibrium_OAE_MDA)


tme <- seq(1, 235*366-1)/366

plot(tme, output_equilibrium_OAE_MDA[[2]], type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))
plot(tme, output_equilibrium_OAE_MDA[[2]], type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1), xlim = c(190, 235))

tme2 <- seq(1, 235*366)/366
tme3 <- seq(1, 235*366)

plot(tme2, output_equilibrium_OAE_MDA[[6]], type = 'l', xlab = 'time (years)', ylab = 'OAE prevalence', ylim = c(0, 1))
plot(tme3, output_equilibrium_OAE_MDA[[6]], type = 'l', xlab = 'time (days)', ylab = 'OAE prevalence', ylim = c(0, 0.1), xlim = c(25000,366*200))
plot(tme3, output_equilibrium_OAE_MDA[[6]], type = 'l', xlab = 'time (days)', ylab = 'OAE prevalence', ylim = c(0, 0.1), xlim = c(366*190,366*235))



