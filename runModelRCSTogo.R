library(dplyr)

iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
set.seed(iter + (iter*3758))

kEs = c(rep(0.3, 4500), rep(0.4, 4500))
seroreversions = rep("no_infection", 9000)

kE = kEs[iter]
sero_val <- seroreversions[iter]

DT.in <- 1/366

if(kE == 0.4) {
    delta.hz.in.val = 0.118,
    delta.hinf.in.val = 0.002,
    c.h.in.val = 0.004,
    gam.dis.in.val = 0.4,
} else {
    delta.hz.in.val =  0.186
    delta.hinf.in.val = 0.003
    c.h.in.val = 0.005
    gam.dis.in.val = 0.3
}

vctr.control.strt <- 80
vctr.control.duration <- 31
vector.control.efficacies <- rep(rep(c(.60, .75, .95), 4500), 2)
vctr.control.efficacy <- vector.control.efficacies[iter]

prefecture = "oti"
if(kE = 0.3) {
    if(prefecture == "bassar") {
        ABR.in <- round(rgamma(1, 20.12, .0077)) # 70% Bassar
    }
    if(prefecture == "oti") {
        ABR.in <- round(rgamma(1, 14.69, .0032)) # 75% Oti
    }
    if(prefecture == "keran") {
        ABR.in <- round(rgamma(1, 7.09, .00029)) # 85% Keran
    }
} else {
    if(prefecture == "bassar") {
        ABR.in <- round(rgamma(1, 38.81, .020)) # 70% Bassar
    }
    if(prefecture == "oti") {
        ABR.in <- round(rgamma(1, 27.08, .0094)) # 75% Oti
    }
    if(prefecture == "keran") {
        ABR.in <- round(rgamma(1, 12.28, .0014)) # 85% Keran
    }
}

if(prefecture == "bassar") {
    # treat.strt.yrs = 1989
    mda.val <- 26
    treat.len = mda.val; treat.strt.yrs = 93; yrs.post.treat = 10

    treat.strt = treat.strt.yrs; treat.stp = treat.strt + treat.len
    timesteps = treat.stp + yrs.post.treat #final duration
    cstm_treat_params <- list(start_biannual=treat.strt.yrs+14, coverage_changes=c(treat.strt.yrs+7, treat.strt.yrs+14), coverage_change_values=c(0.60, 0.75, 0.85))
}
if(prefecture == "oti") {
    # treat.strt.yrs = 1996
    mda.val <- 19
    treat.len = mda.val; treat.strt.yrs = 100; yrs.post.treat = 10

    treat.strt = treat.strt.yrs; treat.stp = treat.strt + treat.len
    timesteps = treat.stp + yrs.post.treat #final duration
    cstm_treat_params <- list(start_biannual=treat.strt.yrs+7, coverage_changes=c(treat.strt.yrs+7), coverage_change_values=c(0.75, 0.80))
}
if(prefecture == "keran") {
    # treat.strt.yrs = 1989
    mda.val <- 26
    treat.len = mda.val; treat.strt.yrs = 93; yrs.post.treat = 10

    treat.strt = treat.strt.yrs; treat.stp = treat.strt + treat.len
    timesteps = treat.stp + yrs.post.treat #final duration
    cstm_treat_params <- list(start_biannual=treat.strt.yrs+14, coverage_changes=c(treat.strt.yrs+7, treat.strt.yrs+14), coverage_change_values=c(0.55, 0.75, 0.85))
}



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
                      delta.hz.in =  delta.hz.in.val,
                      delta.hinf.in = delta.hinf.in.val,
                      c.h.in = c.h.in.val,
                      gam.dis.in = gam.dis.in.val,
                      N.in = 500,
                      run_equilibrium = FALSE,
                      print_progress=TRUE,
                      calc_ov16 = TRUE,
                      no_prev_run=TRUE,
                      custom_treat_params=cstm_treat_params,
                      seroreversion=sero_val)

params <- list(mda.val, ABR.in, kE)
names(params) <- c('MDA', 'ABR', 'Ke')
output <- append(output, params)

saveRDS(output, paste("/rds/general/user/ar722/home/ov16_test/ov16_output/ov16_any_worm_output", kE, "_", iter,".rds", sep=""))
