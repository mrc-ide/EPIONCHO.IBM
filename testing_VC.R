

test.vc.func <- function(vector.control.strt, DT.in, vector.control.duration)

if(!is.na(vector.control.strt)){

  vc.iter.strt <- round(vector.control.strt / (DT.in ));
  vector.control.stp = vector.control.strt + round(vector.control.duration / (DT.in ))

  return(c(vc.iter.strt, vector.control.stp))
}


vc.test <- test.vc.func(vector.control.strt = 1, DT.in = 1/366, vector.control.duration = 20)


# treat.strt.t  = round(1 / (1/366 ));
# treat.stp.t  = 1 + round(20  / (1 / 366 )) #treatment start and stop

test.vc.ABR.func <- function(vc.iter.strt, vector.control.stp, vector.control.efficacy, i, ABR){

  if (i >= vc.iter.strt && i < vector.control.stp) {

  ABR_updated <- ABR * vector.control.efficacy # multiply ABR by effiacy to reduce ABR during VC
  ABR <- ABR_updated
}


m = ABR * ((1/104) / 0.63) # update m

return(c(ABR, m))

}


abr.test <- test.vc.ABR.func(vc.iter.strt = vc.test[1], vector.control.stp = vc.test[2], vector.control.efficacy = 0.7,
                             i = 366, ABR = 20000)

abr.test2 <- test.vc.ABR.func(vc.iter.strt = vc.test[1], vector.control.stp = vc.test[2], vector.control.efficacy = 0.7,
                             i = 365, ABR = 20000)

# ================================================================================================================================== #
#                                           Testing vector control in full model                                                     #


# library(EPIONCHO.IBM) # when updated/ to install locally

devtools::load_all()

# ==================================================== #
#     VC  model run                                    #
# ==================================================== #

#length of simulation in years
timesteps = 80

#should treatment be given, when and how often
give.treat.in = 0
treat.strt = 1; treat.stp = 16
trt.int = 1

#annual biting rate, which determines infection prevalence
ABR.in = 1082

output_VC <- ep.equi.sim(time.its = timesteps,
                                  ABR = ABR.in,
                                  N.in = 440,
                                  treat.int = trt.int,
                                  treat.prob = 0.65,
                                  give.treat = give.treat.in,
                                  treat.start = treat.strt,
                                  treat.stop = treat.stp,
                                  pnc = 0.05,
                                  min.mont.age = 5,
                                  #vector.control.strt = 0.005464481,
                                  vector.control.strt = 40,
                                  vector.control.duration = 40,
                                  vector.control.efficacy = 0.7,
                                  delta.hz.in = 0.186,
                                  delta.hinf.in = 0.003,
                                  c.h.in = 0.005,
                                  gam.dis.in = 0.3,
                                  run_equilibrium = TRUE,
                                  equilibrium = NA,
                                  print_progress = TRUE)

names(output_VC)

tme <- seq(1, 80*366-1)/366

plot(tme, output_VC$ABR_recorded, type = 'l', xlab = 'time (years)', ylab = 'ABR', ylim = c(0, 6000))

plot(tme, output_VC$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))

eq_output <- output_VC$all_equilibrium_outputs

plot(tme, output_VC$all_equilibrium_outputs, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))

plot(tme, output_equilibrium$mf_intens, type = 'l', xlab = 'time (years)', ylab = 'mean microfilarial intensity')

plot(tme, output_equilibrium$L3, type = 'l', xlab = 'time (years)', ylab = 'mean L3 in fly population')


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
                                  vector.control.strt = NA,
                                  delta.hz.in = 0.186,
                                  delta.hinf.in = 0.003,
                                  c.h.in = 0.005,
                                  gam.dis.in = 0.3,
                                  run_equilibrium = TRUE,
                                  equilibrium = NA,
                                  print_progress = TRUE)

names(output_equilibrium)

eq_outputs_baseline <- output_equilibrium$all_equilibrium_outputs

tme <- seq(1, 30*366-1)/366

plot(tme, output_equilibrium$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))

plot(tme, output_equilibrium$mf_intens, type = 'l', xlab = 'time (years)', ylab = 'mean microfilarial intensity')

plot(tme, output_equilibrium$L3, type = 'l', xlab = 'time (years)', ylab = 'mean L3 in fly population')



