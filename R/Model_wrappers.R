#EPIONCHO-IBM
#30/04/2020
#Jonathan Hamley


#' @title
#' Run EPIONCHO-IBM epidemiological model with or without interventions
#' @description
#' Runs the individual-based Onchocerciasis transmission model (EPIONCHO-IBM), based on Hamley et al. 2019
#' ep.equi.sim will run one repeat of the model, typically the mean of 500 repeats is required
#' model must be run to equilibrium (100 years), before treatment can begin
#' treatment is parameterised based on ivermectin
#' if the mf prevalence is zero 50 years after the final treatment, we assume elimination has occured
#' code is available which saves the equilibrium and receives it as an input
#'
#' @param time.its number of iterations (this input is a single value)
#' @param ABR annual biting rate (this input is a single value)
#' @param N.in human population size
#' @param DT timestep (must be one day e.g., 1/366)
#' @param treat.int treatment interval in years e.g., 1 is every 1 year, 0.5 is every 6 months (this input is a single value)
#' @param treat.timing specific timing of treatment rounds (default is NA)
#' @param treat.prob total population coverage (this input is a single value between 0 - 1)
#' @param treat.prob.variable variable total population coverage (value between 0 - 1) specified as a vector of coverages for each treatment round
#' @param give.treat takes 1 (MDA) or 0 (no MDA)
#' @param treat.start iteration where treatment starts
#' @param treat.stop iteration where treatment stops
#' @param pnc proportion of population which never receive treatment (single input value between 0 - 1)
#' @param min.mont.age minimum age for giving a skin snip (single input value, default is 5)
#' @param vector.control.strt start year for vector control
#' @param vector.control.duration duration in years for vector control
#' @param vector.control.efficacy efficacy of vector (proportion of original ABR value)
#' @param delta.hz.in this is a new user-input for density-dependence in humans (proportion of L3 larvae establishing/ developing to adult stage within human host, per bit, when ATP tends to 0)
#' @param delta.hinf.in this is a new user-input for density-dependence in humans (proportion of L3 larvae establishing/ developing to adult stage within human host, per bit, when ATP tends to infinity)
#' @param c.h.in this is a new user-input for density-dependence in humans (severity of transmission intensity - dependent parasite establishment within humans)
#' @param gam.dis.in this is a new user-input for individual-level exposure in humans (overdispersion parameter k_E determining degree of individual exposure heterogeneity; default value is 0.3)
#' @param kM.const.toggle specifies whether kM set to constant (if yes, then overdispersion in mf in skin set to 15)
#' @param run_equilibrium specify whether to run to equilibrium first
#' @param equilibrium equilibrium input given to continue model
#' @param print_progess print the current status of the model run
#' @param epilepsy_module this element determines whether the epilepsy model is turned on ("YES" will activate this)
#' @param OAE_equilibrium OAE equilibrium input given to continue model
#' @param OAE_infection OAE prevalence and incidence inputs at equilibrium
#'
#' @export

ep.equi.sim <- function(time.its,
                        ABR,
                        N.in,
                        treat.int,
                        treat.timing,
                        treat.prob,
                        treat.prob.variable = NA,
                        give.treat,
                        treat.start,
                        treat.stop,
                        pnc,
                        min.mont.age,
                        vector.control.strt,
                        vector.control.duration,
                        vector.control.efficacy,
                        delta.hz.in, # these inputs are new (matt) for testing DD
                        delta.hinf.in,
                        c.h.in,
                        gam.dis.in,
                        kM.const.toggle = FALSE,
                        run_equilibrium,
                        equilibrium,
                        print_progress = TRUE,
                        epilepsy_module = "NO",
                        OAE_equilibrium,
                        calc_ov16=FALSE,
                        ov16_equilibrium=NA,
                        ov16_store_times = c(),
                        no_prev_run=FALSE,
                        custom_treat_params = list(),
                        seroreversion = FALSE)


{
  # ====================== #
  # Set-up time parameters #

  DT <- 1/366
  time.its <- round(time.its / (DT))
  year_its <- seq(0, time.its, 366)
  # if(give.treat == 1) #calculate timesteps at which treatment is given
  # {times.of.treat.in <- seq(treat.start, treat.stop - (treat.int / DT), treat.int / DT)}
  # else {times.of.treat.in <- 0}


  if(give.treat == 1)
  {
    treat.stop <- round(treat.stop / (DT))
    if(treat.start >= 1) {treat.start <-  round( (treat.start) / (DT)) + 1}
    if(treat.start == 0) {treat.start <-  1}
    ov16_store_times <- c(treat.start-1, treat.stop, treat.stop+(1/DT)+1)
  } else {
    ov16_store_times <- c(time.its-1)
  }
  print("Sero Store Times")
  print(ov16_store_times)

  # vector control #
  if(!is.na(vector.control.strt)){

    vc.iter.strt <- round(vector.control.strt / (DT))
    vc.iter.stp <- vc.iter.strt + round(vector.control.duration / (DT))

  } else {
    vc.iter.strt <- NA
    vc.iter.stp <- NA
    vector.control.duration <- NA
  }


  # ================ #
  # hard coded parms #

  # density dep pars (worm establishment in humans)
  delta.hz <- delta.hz.in
  delta.hinf <- delta.hinf.in
  c.h <- c.h.in

  m = ABR * ((1/104) / 0.63) # matt: m = vector to host ratio (V/H) ?; ABR = beta * V/H, where V/H = ABR/beta or V/H = ABR/(h/g), or V/H = ABR * (g/h) which you see here (note, V/H is inferred from the ABR - KEY INPUT for adjusting endemicity level of EPIONCHO-IBM sims)
  beta = 0.63 / (1/104) # matt: beta = per blackfly biting rate on humans (h/g; where h is the human blood index, g is duration of the gonotrophic cycle)
  mu.v = 26
  int.mf = 0
  sigma.L0 = 52
  a.H = 0.8
  g = 0.0096
  a.v = 0.39
  real.max.age = 80 #no humans live longer than 80 years
  N = 2700 #human population size
  mean.age = 50 #mean human age (before truncation)
  int.L3 = 0.03; int.L2 = 0.03; int.L1 = 0.03
  lambda.zero = 0.33 # (matt:) per-capita rate that female worms lose their fertility (W_FF) & return to non-fertile state (W_FN)
  omeg = 0.59 # (matt:) per-capita rate that female worms progress from non-fertile (W_FN) to fertile (W_FF)
  delta.vo = 0.0166 # matt : delta V0, density dependence when microfilarae tend to 0
  c.v = 0.0205 # matt: severity of constraining density-dep larval development (per dermal mf) : Table F Supp
  num.mf.comps = 21; num.comps.worm = 21 #number of age classes for mf and adult worms (matt: c_max ?)
  time.each.comp.worms = 1; time.each.comp.mf = 0.125; mf.move.rate = 8.133333 #for aging in parasites (matt: time.each.comp.worms = q_W & time.each.comp.mf = q_M in supp table E )
  int.worms=1 #initial number of worms in each worm age compartment
  ss.wt = 2; num.ss = 2 #skin snip parameters (matt weight and number)
  slope.kmf = 0.0478 # matt: parameter associated with decreasing degree of aggregation of skin mf with increasing no. of adult female worms (slope in linear model -> k_M = 0.0478 * W_F + 0.313)
  int.kMf = 0.313 # matt: parameter associated with decreasing degree of aggregation of skin mf with increasing no. of adult female worms (inital y value in linear model -> k_M = 0.0478 * W_F + 0.313)

  sex.rat = 0.5 #sex ratio (matt: inidividual assigned a sex randomly - equal probability psi_F = 0.5, psi_m = 0.5)

  nuone = 201.6189; nutwo = 207.7384 #movement of fly parasite life stages

  mu.w1 = 0.09953; mu.w2 = 6.00569 #parameters controlling age-dependent mortality in adult worms (matt: these are y_l = y_w and d_l = d_w in equation S6/S7 & Table E)
  mu.mf1 = 1.089; mu.mf2 = 1.428 #parameters controlling age-dependent mortality in mf (matt: these are y_l = y_m and d_l = d_m in equation S6/S7 & Table E)
  fec.w.1 = 70; fec.w.2 = 0.72 #parameters controlling age-dependent fecundity in adult worms (matt: fec.w.1 = F and fec.w.2 = G in Supp table E)
  l3.delay = 10; dt.days = DT*366 #delay in worms entering humans and joining the first adult worm age class (dt.days = DT.in*366)
  lam.m = 32.4; phi = 19.6 #effects of ivermectin (matt: embryostatic effect - lam.m is the max rate of treatment-induced sterility; phi is the rate of decay of this effect - Table G in Supp)
  cum.infer= 0.345 # permanent infertility in worms due to ivermectin (irreversible sterlising effect- "global") - this could be changed as a macrofilaricidal to 0.9 (90%)
  up = 0.0096; kap = 1.25 #effects of ivermectin (matt: parameters u (up) and k (kap) define the microfilaricidal effect curve, u = finite effect follwoed by decline (rebound) = k - table G in Supp)
  # gam.dis = 0.3 #individual level exposure heterogeneity (matt: shape par in gamma dist, K_E)
  gam.dis <- gam.dis.in # when specifying user input (K_E)
  E0 = 0; q = 0; m.exp = 1.08; f.exp = 0.9; age.exp.m = 0.007; age.exp.f = -0.023 #age-dependent exposure to fly bites age.exp.m or .f = alpha_m or alpha_f)


  # print error messages when incorrect inputs / combination of inputs specified #

  if(give.treat == 1)
  {
    if(all(is.na(treat.prob.variable))){ if(treat.prob > 1) stop('treatment probability must be between 0 & 1') }
    else {if(any(treat.prob.variable > 1)) stop('treatment probability must be between 0 & 1')}

    if(treat.stop > time.its) stop('not enough time for requested MDA duration')

    #times.of.treat.in <- seq(treat.start, treat.stop - (treat.int / DT), treat.int / DT)
    if(length(custom_treat_params) > 0) {
      times.of.treat.in <- c(seq(treat.start, round(custom_treat_params$start_biannual / DT) - treat.int, treat.int / DT),
                             seq(round(custom_treat_params$start_biannual / DT), treat.stop, 0.5 / DT))

      treat.prob.variable <- ifelse(times.of.treat.in < round(custom_treat_params$coverage_changes[1] / DT),
                                    custom_treat_params$coverage_change_values[1],
                                    ifelse(times.of.treat.in < round(custom_treat_params$coverage_changes[2] / DT),
                                           custom_treat_params$coverage_change_values[2],
                                           custom_treat_params$coverage_change_values[3]))
    } else {
      if(all(!is.na(treat.timing))) {treat.timing <- treat.timing + ((treat.start - 367)/ 366)}
      if(all(is.na(treat.timing)))
      {times.of.treat.in <- seq(treat.start, treat.stop, treat.int / DT)}
      else {times.of.treat.in <- round((treat.timing) / (DT)) + 1}
    }

    print(paste(length(times.of.treat.in), 'MDA rounds to be given', sep = ' '))

    print('MDA given at')
    print(paste(round(times.of.treat.in / 366, digits = 2), 'yrs', sep = ''))

    print(times.of.treat.in)

    print('Coverage at each round')
    if(all(!is.na(treat.prob.variable))) {print(paste(treat.prob.variable*100, "%", sep = ''))}
    else{print(paste(treat.prob*100, "%", sep = ''))}

    print('ABR is')
    print(paste(ABR, 'bites person-1 yr-1 at endemic equilibria'))
    if(!is.na(vector.control.strt)){
      print('ABR changes to the following due to vector control')
      print(paste((ABR - (ABR * vector.control.efficacy)), 'bites person-1 yr-1 at', vector.control.strt, 'yrs'))}

  }
  else {times.of.treat.in <- 0; print('no MDA to be simulated')}

  if(is.logical(run_equilibrium) == FALSE) stop('no input indicating if model needs to be run to equilibrium')


  #columns to set to zero when an individual dies
  cols.to.zero <- seq(from = 1, to = (6 + num.mf.comps + 3*num.comps.worm))
  cols.to.zero <- cols.to.zero[-c(1,5, 6)] #compliance, L2 and L3 do not become zero when an individual dies

  #columns, used to perform operations on different worm and mf compartments
  tot.worms <- num.comps.worm*3
  ov16.col <- ifelse(calc_ov16, 1, 0)
  num.cols <- 6 + num.mf.comps + tot.worms
  worms.start <- 7 + num.mf.comps

  nfw.start <- 7 + num.mf.comps + num.comps.worm # start of infertile worms
  fw.start <- nfw.start + num.comps.worm
  fw.end <- num.cols # end of fertile worms
  mf.start <- 7
  mf.end <- 6 + num.mf.comps

  #age-dependent mortality and fecundity rates of parasite life stages

  age.cats <- seq(0, 20, length = num.comps.worm) #up to 20 years old (assume all worms die after age 20 years)

  mort.rates.worms <- weibull.mortality(DT = DT, par1 = mu.w1, par2 = mu.w2, age.cats = age.cats)

  fec.rates.worms <- 1.158305 * fec.w.1 / (fec.w.1 + (fec.w.2 ^ -age.cats) - 1) #no DT - Rk4

  age.cats.mf <- seq(0, 2.5, length = num.mf.comps) #up to 2.5 years old (assume mf die after age 2.5 years)

  #DT not relevent here because RK4 is used to calculate change in mf
  mort.rates.mf <- weibull.mortality(DT = 1, par1 = mu.mf1, par2 = mu.mf2, age.cats = age.cats.mf)

  # ========================================================================================================== #
  # create inital age distribution and simulate stable age distribution (where equilibrium input not provided) #
  if(isTRUE(run_equilibrium) | no_prev_run){
    cur.age <- rep(0, N)

    #(the approach below must be used, drawing human lifespans from an exponential distribution eventually leads to a non-exponential distribution)
    for(i in 1 : 75000) #if at equilibrium you saved the age at which inds die and simulated further, you should get an exponential distribution
    {
      cur.age <- cur.age + DT

      death.vec <- rbinom(N, 1, (1/mean.age) * DT) # Matt: human mortality (constant with age) - no. of deaths at time step t is a random variable drawn from binomial distribution (N = human pop size?)

      cur.age[which(death.vec == 1)] <- 0 #set individuals which die to age 0
      cur.age[which(cur.age >= real.max.age)] <- 0 #all individuals >= maximum imposed age die (matt: distribution truncated to prevent excessively long life spans - a_max)
    }


    ex.vec <- rgamma(N, gam.dis, gam.dis) #individual level exposure to fly bites (matt: individual-specific exposure factor assigned at birth - drawn from gamma dist, with shape par (K_E = gam.dis, and rate par set to this))

    ###############################################
    #matrix for delay in L3 establishment in humans
    num.delay.cols <- l3.delay * (28 / dt.days)
    l.extras <- matrix(0, ncol= num.delay.cols, nrow= N)
    inds.l.mat <- seq(2,(length(l.extras[1,]))) #for moving columns along with time

    ################################################
    #L1 delay in flies
    l1.delay <- rep(int.L1, N)

    ###############################################
    #matrix for tracking mf for L1 delay
    num.mfd.cols <- 4 / dt.days
    mf.delay <- matrix(int.mf, ncol= num.mfd.cols, nrow= N)
    inds.mfd.mats <- seq(2,(length(mf.delay[1,])))

    ###############################################
    #matrix for exposure (to fly bites) for L1 delay
    num.exp.cols <- 4 / dt.days
    exposure.delay <- matrix(ex.vec, ncol= num.exp.cols, nrow= N)
    inds.exp.mats <- seq(2,(length(exposure.delay[1,])))

    #matrix for first timestep, contains all parasite values, human age, sex and compliance

    #all.mats.temp <- matrix(, nrow=N, ncol=num.cols) # error here? (remove the ,)
    all.mats.temp <- matrix(nrow=N, ncol=num.cols+ov16.col) # error here? (remove the ,)

    all.mats.temp[,  (worms.start) : num.cols] <- int.worms

    all.mats.temp[, 4] <- int.L1

    all.mats.temp[, 5] <- int.L2

    all.mats.temp[, 6] <- int.L3

    all.mats.temp[, 7 : (7 + (num.mf.comps-1))] <- int.mf

    all.mats.temp[,1] <- rep(0, N) #column used during treatment
    all.mats.temp[,2] <- cur.age

    #assign sex to humans

    sex <- rbinom(N, 1, sex.rat) # matt: randomly assigned (binomal dist) with equal probability (e.g., psi_F = 0.5, psi_m = 0.5)

    all.mats.temp[,3] <- sex

    #non-compliant people
    non.comp <- ceiling(N * pnc)
    out.comp <- rep(0, N)
    s.comp <- sample(N, non.comp)
    out.comp[s.comp] <- 1
    all.mats.temp[,1] <- out.comp

    treat.vec.in <- rep(NA, N) #for time since treatment calculations

    prev <-  c()
    mean.mf.per.snip <- c()
    L3_vec <- vector()
    ABR_recorded <- c()
    coverage.recorded <- c()

    # i <- 1

    if(epilepsy_module == "YES"){

      # new inputs required #
      infected_at_all <- rep(1, N)
      #age_to_samp <- runif(N, min = 3, max = 10)
      age_to_samp <- sample(seq(3, 15, DT), size = N, replace = TRUE)
      OAE <- rep(1, N)
      prev_OAE <- 1
      tested_OAE <- rep(0, N)
      check_ind <- c()
      OAE_incidence_DT <- c()
      OAE_incidence_DT_under_5 <- c()
      OAE_incidence_DT_5_10 <- c()
      OAE_incidence_DT_11_15 <- c()
      OAE_incidence_DT_M <- c()
      OAE_incidence_DT_F <- c()

      # data and function to obtain OAE probability for a given mf count
      Chesnais_dat <- data.frame(prob = c(0.0061, 0.0439, 0.0720, 0.0849, 0.1341, 0.1538, 0.20),
                                 mf = c(0, 3, 13, 36, 76, 151, 200))

      OAE_probs <- OAE_mfcount_prob_func(dat = Chesnais_dat)
    }
    if(calc_ov16) {
      Ov16_Seropositive <- rep(0, N)
      Ov16_Seropositive_L3 <- rep(0, N)
      Ov16_Seropositive_L4 <- rep(0, N)
      Ov16_Seropositive_mating_no_mf <- rep(0, N)
      Ov16_Seropositive_mating_detectable_mf <- rep(0, N)
      Ov16_Seropositive_mating_any_mf <- rep(0, N)
      time_since_mf_neg <- matrix(-1, nrow=N, ncol=6)


      mf_indv_prev <- rep(0, N)


      Ov16_Seropositive_matrix <- matrix(0, nrow=N, ncol=length(ov16_store_times)*9)
      Ov16_Seropositive_Serorevert_matrix <- matrix(0, nrow=N, ncol=length(ov16_store_times)*9)
      old_matrix_index <- -1
      matrix_index <- 1

      prev_Ov16 <- 0

      # 80% of pop is able to mount Antibody response to Ov16
      all.mats.temp[,num.cols+ov16.col] <- sample(rep(c(1,1,1,1,1,1,1,1,0,0), N*0.1), N)
    }
  }

  # =============================================================================================#
  # set-up main structures for tracking human/parasite stages if equilibrium dataframe given     #
  if(isFALSE(run_equilibrium) & !no_prev_run)

  {
    if(is.list(equilibrium) == FALSE) stop('equilibrium condition not in correct format')
    # if(epilepsy_module == "YES") stop(' cannot run epilepsy model unless model run at/to equilibrium first')

    ex.vec <- equilibrium[[2]] #exposure

    ###############################################
    #matrix for delay in l3 establishment in humans
    num.delay.cols <- l3.delay * (28 / dt.days)
    l.extras <- equilibrium[[4]]
    inds.l.mat <- seq(2,(length(l.extras[1,]))) #for moving columns along with time

    ################################################
    l1.delay <- equilibrium[[6]]

    ###############################################
    #matrix for tracking mf for l1 delay
    num.mfd.cols <- 4 / dt.days
    mf.delay <- equilibrium[[5]]
    inds.mfd.mats <- seq(2,(length(mf.delay[1,])))

    ###############################################
    #matrix for exposure for L1 delay
    num.exp.cols <- 4 / dt.days

    #matrix for first time step
    all.mats.temp <- equilibrium[[1]]


    exposure.delay <- equilibrium[[8]]
    inds.exp.mats <- seq(2,(length(exposure.delay[1,])))

    temp <- mf.per.skin.snip(ss.wt = 2, num.ss = 2, slope.kmf = 0.0478, int.kMf = 0.313, data = all.mats.temp, nfw.start, fw.end,
                             mf.start, mf.end, pop.size = N, kM.const.toggle)


    prev <-  prevalence.for.age(age = 5, ss.in = temp, main.dat = all.mats.temp)

    mean.mf.per.snip <- mean(temp[[2]][which(all.mats.temp[,2] >= 5)])

    L3_vec <- mean(all.mats.temp[, 6])

    treat.vec.in <- equilibrium[[3]]

    ABR_recorded <- c()
    coverage.recorded <- c()

    if(epilepsy_module == "YES"){

      OAE <- OAE_equilibrium[[1]]
      age_to_samp <- OAE_equilibrium[[2]]
      tested_OAE <- OAE_equilibrium[[3]]
      infected_at_all <- OAE_equilibrium[[4]]
      check_ind <- OAE_equilibrium[[5]]
      tot_ind_ep_samp <- OAE_equilibrium[[6]]
      OAE_probs <- OAE_equilibrium[[7]]

      prev_OAE <- mean(OAE) # calculate & update prevalence of OAE

      new_inc <- length(which(OAE[tot_ind_ep_samp] == 1)) # how many infections (finds total new infected/OAE in all OAE)

      OAE_incidence_DT <- new_inc # record + update number of new OAE cases


      new_inc_under_5 <- length(which(OAE[tot_ind_ep_samp] == 1 & all.mats.temp[tot_ind_ep_samp ,2] >= 3 & all.mats.temp[tot_ind_ep_samp ,2]< 5 )) # new cases in 3 to 5 age group
      new_inc_5_10 <- length(which(OAE[tot_ind_ep_samp] == 1 & all.mats.temp[tot_ind_ep_samp ,2] >= 5 & all.mats.temp[tot_ind_ep_samp ,2]<= 10 )) # new cases in 5 to 10 age group
      new_inc_11_15 <- length(which(OAE[tot_ind_ep_samp] == 1 & all.mats.temp[tot_ind_ep_samp ,2] >= 10 & all.mats.temp[tot_ind_ep_samp ,2]<= 15 )) # new cases in 10 to 15 age group

      new_inc_M <- length(which(OAE[tot_ind_ep_samp] == 1 & all.mats.temp[tot_ind_ep_samp ,3] == 1)) # new cases in males
      new_inc_F <- length(which(OAE[tot_ind_ep_samp] == 1 & all.mats.temp[tot_ind_ep_samp ,3] == 0)) # new cases in females

      OAE_incidence_DT_under_5 <- new_inc_under_5 # record & update incidence in 3 to 5 age group
      OAE_incidence_DT_5_10 <- new_inc_5_10 # record & update incidence in 5 to 10 age group
      OAE_incidence_DT_11_15 <- new_inc_11_15 # record & update incidence in 10 to 15 age group

      OAE_incidence_DT_M <- new_inc_M # record & update incidence in males
      OAE_incidence_DT_F <- new_inc_F # record & update incidence in females

      # # if taking straight from eq input
      # prev_OAE <- OAE_infection[[1]] # calculate & update prevalence of OAE
      #
      # OAE_incidence_DT <- OAE_infection[[2]] # record + update number of new OAE cases
      #
      # OAE_incidence_DT_3_5 <- OAE_infection[[3]] # record & update incidence in 3 to 5 age group
      # OAE_incidence_DT_5_10 <- OAE_infection[[4]] # record & update incidence in 5 to 10 age group
      #
      # OAE_incidence_DT_M <- OAE_infection[[5]] # record & update incidence in males
      # OAE_incidence_DT_F <- OAE_infection[[6]] # record & update incidence in females

    }

    if(calc_ov16) {
      Ov16_Seropositive <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive
      Ov16_Seropositive_L3 <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive_l3
      Ov16_Seropositive_L4 <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive_l4
      Ov16_Seropositive_mating_no_mf <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive_mating_no_mf
      Ov16_Seropositive_mating_detectable_mf <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive_mating_detectable_mf
      Ov16_Seropositive_mating_any_mf <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive_mating_any_mf

      mf_indv_prev <- ov16_equilibrium$mf_indv_prev

      Ov16_Seropositive_matrix <- matrix(0, nrow=N, ncol=length(ov16_store_times)*9)
      matrix_index <- 1

      prev_Ov16 <- sum(Ov16_Seropositive)/N
    }

  }

  i <- 1

  # ================================================================================================= #
  #               Main loop for running through processes in model with i > 1                         #
  while(i < time.its) #over time

  {
    #print(paste(round(i * DT, digits = 2), 'yrs', sep = ' '))

    # if(isTRUE(print_progress)) {print(paste(round(i * DT, digits = 2), 'yrs', sep=' '))}

    # if(isTRUE(print_progress) & (any(i == year_its))) {print(paste(round(i * DT, digits = 2), 'yrs', sep=' '))}
    # if(isTRUE(print_progress) & (any(i == year_its))) {print(paste(round(i/time.its * 100, digits = 1), '%', sep=' '))}

    if(isTRUE(print_progress) & (any(i == year_its))) {print(paste(round(i * DT, digits = 2), 'yrs;',
                                                                   (paste(round(i/time.its * 100, digits = 1), '%', sep=' '))))}

    # extract ov16 results before treatment
    if (calc_ov16 & i == treat.start) {

      all.mats.temp_pre_treat <- all.mats.temp
      Ov16_Seropositive_pre_treat <- Ov16_Seropositive
      mf_indv_prev_pre_treat <- mf_indv_prev

    }

    #stores mean L3 and adult worms from previous timesteps

    all.mats.cur <- all.mats.temp

    # extract element from treat.prob.variable depending on specific treatment round (iteration in times.of.treat.in) #
    if(all(!is.na(treat.prob.variable))){

      if(any(i == times.of.treat.in)) {
        index.iter.treat <- match(i, times.of.treat.in) # find element where iteration number matches a time in times.of.treat vector
        treat.prob <- treat.prob.variable[index.iter.treat]} # index prob value from treat.prob.variable vector

    }

    # to track variable coverage
    if(i >= treat.start & i <= treat.stop & give.treat == 1){
      coverage.upd <- treat.prob
    } else
    {coverage.upd <- 0}

    #which individuals will be treated if treatment is given
    if(i >= treat.start & give.treat==1) {cov.in <- os.cov(all.dt = all.mats.cur, pncomp = pnc, covrg = treat.prob, N = N)}

    #sex and age dependent exposure, mean exposure must be 1, so ABR is meaningful

    mls <- which(all.mats.cur[,3] == 1) # matt : ?
    fmls <- which(all.mats.cur[,3] == 0) # matt: ?

    s.a.exp <- rep(0, N)

    s.a.exp[mls] <- m.exp * exp(-age.exp.m * (all.mats.cur[mls, 2]))

    gam.m <- 1 / mean(s.a.exp[mls]) #normalize so mean = 1 (matt: is this equivalent to including the gamma_s term?)
    s.a.exp[mls] <- s.a.exp[mls] * gam.m

    s.a.exp[fmls] <- f.exp * exp(-age.exp.f * (all.mats.cur[fmls, 2]))

    gam.f <- 1 / mean(s.a.exp[fmls]) #normalize so mean = 1
    s.a.exp[fmls] <- s.a.exp[fmls] * gam.f

    ex.vec <- ex.vec * (1 / mean(ex.vec)) #normalize so mean = 1 (matt: normalising the indvidual-specific exposure from line 565)

    tot.ex.ai <- s.a.exp * ex.vec # matt: combine sex/age specific exposure + individual specific exposure (total exposure to blackfly bites)
    tot.ex.ai <- tot.ex.ai * (1 / mean(tot.ex.ai)) #normalize so mean = 1

    #increase age (for next time step)

    all.mats.temp[,2] <- (all.mats.cur[,2]) + DT #increase age for all individuals

    death.vec <- rbinom(N, 1, (1/mean.age) * DT) #select individuals to die

    to.die <- which(death.vec == 1)

    at.ab.max <- which(all.mats.temp[,2] >= real.max.age)

    to.die <- c(to.die, at.ab.max)

    to.die <- unique(to.die) #may have repeated indivudals i.e selected by binom and >80

    ##################
    #delay calculations
    ##################

    #there is a delay in new parasites entering humans (from fly bites) and entering the first adult worm age class

    new.worms.m <- c()
    new.worms.nf <- c()

    new.worms.m <- rbinom(N, size = l.extras[,length(l.extras[1,])], prob = 0.5) #draw males and females from last column of delay matrix
    new.worms.nf <- l.extras[,length(l.extras[1,])] - new.worms.m

    #move individuals along
    l.extras[,inds.l.mat] <- l.extras[,(inds.l.mat-1)]

    #mean number of L3 in fly population
    L3.in <- mean(all.mats.cur[, 6])


    # change m based on ABR change due to vector control if called (during vector control duration iteration period)
    if(!is.na(vector.control.strt)){

      if (i >= vc.iter.strt && i < vc.iter.stp) {

        ABR_updated <- ABR - (ABR * vector.control.efficacy) # proportional reduction in ABR (x efficacy) during VC

        m = ABR_updated * ((1/104) / 0.63) # update m
      }
    }

    # change m back to original after VC finishes (can change when this occurs i.e, one year after VC ends)
    if(!is.na(vector.control.strt)){

      if (i >= vc.iter.stp) {

        m = ABR * ((1/104) / 0.63) # update m
      }
    }


    # to track #
    if (!is.na(vector.control.strt) && i >= vc.iter.strt && i < vc.iter.stp) {

      ABR_upd <- ABR_updated

    } else {

      ABR_upd <- ABR
    }


    #rate of infections in humans
    #delta.hz, delta.hinf, c.h are density dependence parameters, expos is the exposure of each person to bites
    nw.rate <- Wplus1.rate(delta.hz, delta.hinf, c.h, L3 = L3.in, m ,
                           beta, expos = tot.ex.ai, DT)


    new.worms <- rpois(N, nw.rate) #total new establishing L3 for each individual

    l.extras[,1] <- new.worms


    for(k in 1 : num.comps.worm) #go through each adult worm compartment

    {

      if(k == 1) {from.last <- rep(0, N)} #create vector for worms coming from previous compartment (needs to be 0 when k ==1)


      res.w1 <- change.worm.per.ind1(treat.vec = treat.vec.in, lambda.zero = lambda.zero, DT = DT, omeg = omeg,
                                     ws = worms.start, compartment = k, total.dat = all.mats.cur, mort.rates = mort.rates.worms,
                                     time.each.comp = time.each.comp.worms, new.worms.m = new.worms.m, w.f.l.c = from.last,
                                     num.comps = num.comps.worm)

      res.w.treat <- change.worm.per.ind.treat(give.treat = give.treat, iteration = i, treat.start = treat.start, times.of.treat = times.of.treat.in, treat.stop = treat.stop,
                                               onchosim.cov = cov.in, treat.vec = treat.vec.in, DT = DT, cum.infer = cum.infer, lam.m = lam.m, phi = phi, N = res.w1[[3]],
                                               mort.fems = res.w1[[2]], lambda.zero.in = res.w1[[1]])

      res.w2 <- change.worm.per.ind2(DT = DT, time.each.comp = time.each.comp.worms, compartment = k, new.worms.nf.fo = new.worms.nf, w.f.l.c = from.last,
                                     N = res.w1[[3]], cur.Wm.nf = res.w1[[4]], mort.fems = res.w.treat[[3]], cur.Wm.f = res.w1[[5]], omeg = res.w1[[7]],
                                     male.tot.worms = res.w1[[8]], worm.loss.males = res.w1[[9]],
                                     lambda.zero.in = res.w.treat[[1]], treat.vec = res.w.treat[[2]])


      res <- res.w2 # (matt: re-label the final result output to res so do not have to change res below)

      from.last <- res # (matt: re-label the final result output to res so do not have to change res below)

      # from.last <- res #assign output to use at next iteration, indexes 2, 5, 6 (worms moving through compartments)

      # update male worms in matrix for compartment k

      all.mats.temp[, (6 + num.mf.comps + k)] <- res[[1]]

      # update females worms in matrix

      all.mats.temp[, (6 + num.mf.comps + num.comps.worm + k)] <- res[[3]] # infertile, num.comps.worm skips over males
      all.mats.temp[, (6 + num.mf.comps + 2*num.comps.worm + k)] <- res[[4]] # fertile, num.comps.worm skips over males and infertile females

      if(give.treat == 1 & i >= treat.start & k == num.comps.worm) {treat.vec.in <- res[[7]]} #treated individuals
    }

    for(mf.c in 1 : num.mf.comps)

    {

      res.mf <- change.micro(dat = all.mats.cur, num.comps =num.comps.worm, mf.cpt = mf.c,
                             num.mf.comps = num.mf.comps, ws=worms.start, DT=DT, time.each.comp = time.each.comp.mf,
                             mu.rates.mf = mort.rates.mf, fec.rates = fec.rates.worms, mf.move.rate = mf.move.rate, up = up, kap = kap, iteration = i,
                             treat.vec = treat.vec.in, give.treat = give.treat, treat.start = treat.start)

      all.mats.temp[, 6 + mf.c] <- res.mf
    }


    #inputs for delay in L1
    exp.delay.temp <- exposure.delay[, length(exposure.delay[1,])]
    mf.delay.temp <- mf.delay[, length(mf.delay[1,])]
    l1.delay.temp <- l1.delay #L1 from previous timestep

    #move values along
    exposure.delay[, inds.exp.mats] <- exposure.delay[, (inds.exp.mats -1)]
    mf.delay[, inds.mfd.mats] <- mf.delay[, (inds.mfd.mats - 1)]

    #update L1, L2 and L3

    #total number of mf in each person
    mf.temp <- rowSums(all.mats.cur[, 7 : (6 + num.mf.comps)]) #sum mf over compartments, mf start in column 7

    all.mats.temp[, 4] <- calc.L1(beta, mf = mf.temp, mf.delay.in = mf.delay.temp, expos = tot.ex.ai, delta.vo, c.v, nuone, mu.v, a.v, expos.delay = exp.delay.temp)
    all.mats.temp[, 5] <- calc.L2(nuone, L1.in = l1.delay.temp, mu.v, nutwo, mf = mf.delay.temp, a.v, expos = exp.delay.temp)
    all.mats.temp[, 6] <- calc.L3(nutwo, L2.in = all.mats.cur[, 5], a.H, g, mu.v, sigma.L0)

    #new values for delay parts
    l1.delay <- all.mats.temp[, 4]
    mf.delay[, 1] <- rowSums(all.mats.cur[, 7 : (6 + num.mf.comps)])
    exposure.delay[, 1] <- tot.ex.ai + 0

    #===========================#
    #     OAE module funcs      #

    if(epilepsy_module == "YES"){

      if(i == 1){

        OAE_out1 <- find_indiv_OAE_func(dat = all.mats.temp, mf.start = mf.start, mf.end = mf.end, worms.start = worms.start, nfw.start = nfw.start, fw.end = fw.end,
                                        infected_at_all = infected_at_all, age_to_samp = age_to_samp, OAE = OAE, tested_OAE = tested_OAE, check_ind = check_ind) # step 1

        infected_at_all = OAE_out1[[2]] # updated (when i = 1)
        check_ind = OAE_out1[[3]] # updated (when i = 1)

        temp.mf <- mf.per.skin.snip(ss.wt = 2, num.ss = 2, slope.kmf = 0.0478, int.kMf = 0.313, data = all.mats.temp, nfw.start, fw.end,
                                    mf.start, mf.end, pop.size = N, kM.const.toggle)

        OAE_out2 <- new_OAE_cases_func(temp.mf = temp.mf, tot_ind_ep_samp = OAE_out1[[1]], OAE_probs = OAE_probs, dat = all.mats.temp,
                                       OAE = OAE, tested_OAE = tested_OAE,
                                       prev_OAE = prev_OAE, OAE_incidence_DT = OAE_incidence_DT,
                                       OAE_incidence_DT_under_5 = OAE_incidence_DT_under_5, OAE_incidence_DT_5_10 = OAE_incidence_DT_5_10, OAE_incidence_DT_11_15 = OAE_incidence_DT_11_15,
                                       OAE_incidence_DT_M = OAE_incidence_DT_M, OAE_incidence_DT_F = OAE_incidence_DT_F) # step 2

        OAE = OAE_out2[[9]] # updated (when i = 1)
        tested_OAE = OAE_out2[[10]] # updated (when i = 1)
      }

      if(i > 1){

        OAE_out1 <- find_indiv_OAE_func(dat = all.mats.temp, mf.start = mf.start, mf.end = mf.end, worms.start = worms.start, nfw.start = nfw.start, fw.end = fw.end,
                                        infected_at_all = infected_at_all, age_to_samp = age_to_samp, OAE = OAE, tested_OAE = tested_OAE, check_ind = check_ind) # step 1

        infected_at_all = OAE_out1[[2]] # updated (when i > 1)
        check_ind = OAE_out1[[3]] # updated (when i > 1)

        temp.mf <- mf.per.skin.snip(ss.wt = 2, num.ss = 2, slope.kmf = 0.0478, int.kMf = 0.313, data = all.mats.temp, nfw.start, fw.end,
                                    mf.start, mf.end, pop.size = N, kM.const.toggle)

        OAE_out2 <- new_OAE_cases_func(temp.mf = temp.mf, tot_ind_ep_samp = OAE_out1[[1]], OAE_probs = OAE_probs, dat = all.mats.temp,
                                       OAE = OAE, tested_OAE = tested_OAE,
                                       prev_OAE = OAE_out2[[1]], OAE_incidence_DT = OAE_out2[[2]],
                                       OAE_incidence_DT_under_5 = OAE_out2[[3]], OAE_incidence_DT_5_10 = OAE_out2[[4]], OAE_incidence_DT_11_15 = OAE_out2[[5]],
                                       OAE_incidence_DT_M = OAE_out2[[6]], OAE_incidence_DT_F = OAE_out2[[7]]) # step 2

        OAE = OAE_out2[[9]] # updated (when i > 1)
        tested_OAE = OAE_out2[[10]] # updated (when i > 1)
      }
    }

    #save prevalence at current timestep
    temp.mf <- mf.per.skin.snip(ss.wt = 2, num.ss = 2, slope.kmf = 0.0478, int.kMf = 0.313, data = all.mats.temp, nfw.start, fw.end,
                                mf.start, mf.end, pop.size = N, kM.const.toggle)

    prev <-  c(prev, prevalence.for.age(age = min.mont.age, ss.in = temp.mf, main.dat = all.mats.temp))


    mean.mf.per.snip <- c(mean.mf.per.snip, mean(temp.mf[[2]][which(all.mats.temp[,2] >= min.mont.age)]))

    mf.per.skin.snp.out <- temp.mf[[2]] #to extract out mf per skin snip for each individual?

    L3_vec <- c(L3_vec, mean(all.mats.temp[, 6]))

    #ov16_sero
    if(calc_ov16) {
      any_juvy_worms <- (rowSums(all.mats.temp[, worms.start:num.cols]) > 0)
      any_l3_exposure <- (l.extras[,1] > 0)
      any_larvae <- (rowSums(l.extras) > 0)
      l4_development <- (l.extras[,floor(length(l.extras[1,])/2)] > 0)
      any_worms <- rowSums(all.mats.temp[,worms.start:fw.end])
      mating_worm <- ((rowSums(all.mats.temp[,worms.start:nfw.start])) > 0 & (rowSums(all.mats.temp[, fw.start:fw.end]) > 0))
      mating_worm_detectable_mf <- (mating_worm & temp.mf[[2]] > 0)
      mating_worm_any_mf <- (mating_worm == TRUE & (rowSums(all.mats.temp[,mf.start:mf.end]) > 0))

      indv_antibody_response <- all.mats.temp[,91]

      findPositives <- function(exposure_array, curr_array, antibody_resp, mf_neg_index, doSerorevert=FALSE) {
        curr_array[which(exposure_array == TRUE & curr_array == 0 & antibody_resp == 1)] <- 1
        # hard seroreversion
        if(doSerorevert & seroreversion) {
          curr_array[which(curr_array == 1 & any_larvae == FALSE & exposure_array == FALSE & any_worms == FALSE & rowSums(all.mats.temp[,mf.start:mf.end]) == 0)] <- 0
        }

        # From kamlan et al.
        # actual_mf_ov16_data <- c(1826-1825, 1875-1735, 1323-1323, 1987-1938)
        # seroreversion_start <- which(curr_array == 1 & exposure_array == FALSE & any_worms == FALSE & rowSums(temp.mf[[2]]) == 0)
        # time_since_mf_neg[seroreversion_start,mf_neg_index] <- time_since_mf_neg + 1
        # should_serorevert <- rnorm(N, mean(actual_mf_ov16_data), sd(actual_mf_ov16_data))
        # curr_array[which(time_since_mf_neg[,1] >0 & should_serorevert <= time_since_mf_neg[,1])] <- 0
        return(curr_array)
      }

      Ov16_Seropositive <- findPositives(any_juvy_worms, Ov16_Seropositive, indv_antibody_response, 1)
      Ov16_Seropositive_L3 <- findPositives(any_l3_exposure, Ov16_Seropositive_L3, indv_antibody_response, 2)
      Ov16_Seropositive_L4 <- findPositives(l4_development, Ov16_Seropositive_L4, indv_antibody_response, 3)
      Ov16_Seropositive_mating_no_mf <- findPositives(mating_worm, Ov16_Seropositive_mating_no_mf, indv_antibody_response, 4)
      Ov16_Seropositive_mating_detectable_mf <- findPositives(mating_worm_detectable_mf, Ov16_Seropositive_mating_detectable_mf, indv_antibody_response, 5)
      Ov16_Seropositive_mating_any_mf <- findPositives(mating_worm_any_mf, Ov16_Seropositive_mating_any_mf, indv_antibody_response, 6)

      index_to_use <- matrix_index
      if(old_matrix_index != -1) {
        index_to_use <- index_to_use - 1
        old_matrix_index <- -1
      }

      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-5] <- findPositives(Ov16_Seropositive_Serorevert_matrix[,9*index_to_use-5], Ov16_Seropositive, indv_antibody_response, 1, doSerorevert=TRUE)
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-4] <- findPositives(Ov16_Seropositive_Serorevert_matrix[,9*index_to_use-4], Ov16_Seropositive_L3, indv_antibody_response, 2, doSerorevert=TRUE)
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-3] <- findPositives(Ov16_Seropositive_Serorevert_matrix[,9*index_to_use-3], Ov16_Seropositive_L4, indv_antibody_response, 3, doSerorevert=TRUE)
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-2] <- findPositives(Ov16_Seropositive_Serorevert_matrix[,9*index_to_use-2], Ov16_Seropositive_mating_no_mf, indv_antibody_response, 4, doSerorevert=TRUE)
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-1] <- findPositives(Ov16_Seropositive_Serorevert_matrix[,9*index_to_use-1], Ov16_Seropositive_mating_detectable_mf, indv_antibody_response, 5, doSerorevert=TRUE)
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index] <- findPositives(Ov16_Seropositive_Serorevert_matrix[,9*index_to_use], Ov16_Seropositive_mating_any_mf, indv_antibody_response, 6, doSerorevert=TRUE)

      mf_indv_prev <- as.integer(temp.mf[[2]] > 0)
      prev_Ov16 <- c(prev_Ov16, sum(Ov16_Seropositive)/N)
    }

    if(!is.na(vector.control.strt)) {
      ABR_recorded <- c(ABR_recorded, ABR_upd) # tracking changing ABR
    }

    if(isTRUE(run_equilibrium) & i >= treat.start & i <= treat.stop & give.treat == 1) {
      coverage.recorded <- c(coverage.recorded, coverage.upd) # track changing coverage if specified
    }

    # new individual exposure for newborns, clear rows for new borns
    if(length(to.die) > 0)
    {
      ex.vec[to.die] <- rgamma(length(to.die), gam.dis, gam.dis)

      l.extras[to.die, ] <- 0 #establishing adult worms

      mf.delay[to.die, 1] <- 0 #individual dies so no contribution to L1s at this timestep

      l1.delay[to.die] <- 0

      treat.vec.in[to.die] <- NA

      all.mats.temp[to.die, cols.to.zero] <- 0 #set age, sex and parasites to 0 (includes L1, but not L2 L3)
      all.mats.temp[to.die, 3] <- rbinom(length(to.die), 1, 0.5) #draw sex

      if(calc_ov16) {
        Ov16_Seropositive[to.die] <- 0
        Ov16_Seropositive_L3[to.die] <- 0
        Ov16_Seropositive_L4[to.die] <- 0
        Ov16_Seropositive_mating_no_mf[to.die] <- 0
        Ov16_Seropositive_mating_detectable_mf[to.die] <- 0
        Ov16_Seropositive_mating_any_mf[to.die] <- 0
        mf_indv_prev[to.die] <- 0
      }

      if(epilepsy_module == "YES"){

        infected_at_all[to.die] <- 0 # index those individuals to die as no longer ever infected

        age_to_samp[to.die] <- sample(seq(3, 15, DT), size = length(to.die), replace = TRUE) # for those individuals set to die, resample

        OAE[to.die] <-  0 # index those individuals to die as no longer with OAE

        tested_OAE[to.die] <-  0 # index those individuals to die as no longer tested

      }

    }
    if(calc_ov16 & !is.na(match(i, ov16_store_times))) {
      Ov16_Seropositive_matrix[,9*matrix_index-8] <- all.mats.temp[,2]
      Ov16_Seropositive_matrix[,9*matrix_index-7] <- all.mats.temp[,3]
      Ov16_Seropositive_matrix[,9*matrix_index-6] <- as.integer(temp.mf[[2]] > 0)
      Ov16_Seropositive_matrix[,9*matrix_index-5] <- Ov16_Seropositive
      Ov16_Seropositive_matrix[,9*matrix_index-4] <- Ov16_Seropositive_L3
      Ov16_Seropositive_matrix[,9*matrix_index-3] <- Ov16_Seropositive_L4
      Ov16_Seropositive_matrix[,9*matrix_index-2] <- Ov16_Seropositive_mating_no_mf
      Ov16_Seropositive_matrix[,9*matrix_index-1] <- Ov16_Seropositive_mating_detectable_mf
      Ov16_Seropositive_matrix[,9*matrix_index] <- Ov16_Seropositive_mating_any_mf

      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-8] <- all.mats.temp[,2]
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-7] <- all.mats.temp[,3]
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-6] <- as.integer(temp.mf[[2]] > 0)

      old_matrix_index <- matrix_index
      matrix_index <- matrix_index + 1
    }


    i <- i + 1

  }

  if(epilepsy_module == "YES"){

    # return(list(all.mats.temp, prev, mean.mf.per.snip, mf.per.skin.snp.out, OAE, prev_OAE = OAE_out2[[1]], check_ind = OAE_out1[[3]],
    #             OAE_incidence_DT = OAE_out2[[2]], OAE_incidence_DT_3_5 = OAE_out2[[3]], OAE_incidence_DT_5_10 = OAE_out2[[4]],
    #             OAE_incidence_DT_M = OAE_out2[[5]], OAE_incidence_DT_F = OAE_out2[[6]])) #[[2]] is mf prevalence, [[3]] is intensity

    if(isTRUE(run_equilibrium)){
      outp <- (list(prev, mean.mf.per.snip, L3_vec,
                    list(all.mats.temp, ex.vec, treat.vec.in, l.extras, mf.delay, l1.delay, ABR, exposure.delay),
                    prev_OAE = OAE_out2[[1]], OAE_incidence_DT = OAE_out2[[2]],
                    OAE_incidence_DT_under_5 = OAE_out2[[3]], OAE_incidence_DT_5_10 = OAE_out2[[4]], OAE_incidence_DT_11_15 = OAE_out2[[5]],
                    OAE_incidence_DT_M = OAE_out2[[6]], OAE_incidence_DT_F = OAE_out2[[7]],
                    list(OAE = OAE, age_to_samp = age_to_samp, tested_OAE = tested_OAE, infected_at_all = infected_at_all,
                         check_ind = OAE_out1[[3]], tot_ind_ep_samp = OAE_out1[[1]], OAE_probs = OAE_probs),
                    ABR_recorded, coverage.recorded))

      names(outp) <- c('mf_prev', 'mf_intens', 'L3', 'all_equilibrium_outputs', 'OAE_prev','OAE_incidence',
                       'OAE_incidence_under_5yrs','OAE_incidence_5_10yrs','OAE_incidence_10_15yrs',
                       'OAE_incidence_males','OAE_incidence_females','all_OAE_equilibirum_ouputs',
                       'ABR_recorded', 'coverage.recorded')
      return(outp)
    }

    if(isFALSE(run_equilibrium))
    {
      outp <- (list(prev, mean.mf.per.snip, L3_vec, ABR, all.mats.temp,
                    prev_OAE = OAE_out2[[1]], OAE_incidence_DT = OAE_out2[[2]],
                    OAE_incidence_DT_under_5 = OAE_out2[[3]], OAE_incidence_DT_5_10 = OAE_out2[[4]], OAE_incidence_DT_11_15 = OAE_out2[[5]],
                    OAE_incidence_DT_M = OAE_out2[[6]], OAE_incidence_DT_F = OAE_out2[[7]],
                    list(OAE = OAE, age_to_samp = age_to_samp, tested_OAE = tested_OAE, infected_at_all = infected_at_all,
                         check_ind = OAE_out1[[3]], tot_ind_ep_samp = OAE_out1[[1]], OAE_probs = OAE_probs),
                    ABR_recorded, coverage.recorded))

      names(outp) <- c('mf_prev', 'mf_intens', 'L3', 'ABR','all_equilibrium_outputs', 'OAE_prev','OAE_incidence',
                       'OAE_incidence_under_5yrs','OAE_incidence_5_10yrs','OAE_incidence_10_15yrs',
                       'OAE_incidence_males','OAE_incidence_females','all_OAE_equilibirum_ouputs',
                       'ABR_recorded', 'coverage.recorded')
      return(outp)
    }

  } else {

    # ================================#
    #  When epilpsy module not called #


    #enough outputs to restart sims
    if(isTRUE(run_equilibrium))
    {
      outp <- list(prev, mean.mf.per.snip, L3_vec, list(all.mats.temp, ex.vec, treat.vec.in, l.extras, mf.delay, l1.delay, ABR, exposure.delay), ABR_recorded, coverage.recorded)
      names(outp) <- c('mf_prev', 'mf_intens', 'L3', 'all_equilibrium_outputs', 'ABR_recorded', 'coverage.recorded')
      if(calc_ov16) {
        ov16_seropos_outputs <- list(Ov16_Seropositive, Ov16_Seropositive_L3, Ov16_Seropositive_L4, Ov16_Seropositive_mating_no_mf, Ov16_Seropositive_mating_detectable_mf, Ov16_Seropositive_mating_any_mf)
        names(ov16_seropos_outputs) <- c('ov16_seropositive', 'ov16_seropositive_l3', 'ov16_seropositive_l4', 'ov16_seropositive_mating_no_mf', 'ov16_seropositive_mating_detectable_mf', 'ov16_seropositive_mating_any_mf')
        ov16_equilibrium_outputs <- list(ov16_seropos_outputs, mf_indv_prev)
        names(ov16_equilibrium_outputs) <- c('ov16_seropos_outputs', 'mf_indv_prev')
        ov16_output <- list(ov16_equilibrium_outputs)
        names(ov16_output) <- c('ov16_equilibrium')
        outp <- append(outp, ov16_output)
      }
      return(outp)
    }

    #assuming output will not be used for further sims
    if(isFALSE(run_equilibrium))
    {
      outp <- list(prev, mean.mf.per.snip, L3_vec, ABR, all.mats.temp, ABR_recorded, coverage.recorded)
      names(outp) <-  c('mf_prev', 'mf_intens', 'L3', 'ABR', 'all_infection_burdens', 'ABR_recorded', 'coverage.recorded')
      if(calc_ov16) {
        ov16_output <- list(Ov16_Seropositive, mf_indv_prev, Ov16_Seropositive_pre_treat, mf_indv_prev_pre_treat, all.mats.temp_pre_treat, Ov16_Seropositive_matrix, Ov16_Seropositive_Serorevert_matrix)
        names(ov16_output) <- c('ov16_seropositive', 'mf_indv_prevalence', 'ov16_seropositive_pre_treatment', 'mf_indv_prevalence_pre_treatment', 'all_infection_burdens_pre_treatment', 'ov16_seropositive_matrix', "ov16_seropositive_matrix_serorevert")
        outp <- append(outp, ov16_output)
      }
      return(outp)
    }
  }


}
