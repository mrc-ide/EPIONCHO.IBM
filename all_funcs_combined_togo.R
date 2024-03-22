

#### Current file: R/adult_worm_dynamics_functions.R 

#' @title
#' delta.h
#' @description
#' proportion of L3 larvae (final life stage in the fly population) developing into adult worms in humans (Pi_H in Hamley et al. 2019)
#' @param delta.hz proportion of L3 larvae establishing/ developing to adult stage within human host, per bit, when ATP tends to 0
#' @param delta.hinf proportion of L3 larvae establishing/ developing to adult stage within human host, per bit, when ATP tends to infinity
#' @param c.h severity of transmission intensity - dependent parasite establishment within humans
#' @param L3 mean number of L3 in fly population
#' @param m  vector to human host ratio
#' @param beta per blackfly biting rate on humans (product of proportion of blackfly bites take on humans (h) & reciprocal of the duration of the gonotrophic cycle (g))
#' @param expos individuals total (age/sex-specific and individual-specific) exposure to blackfly bites
#'
#' @returns vector of density-dependent (delta_h) values per individual in the population
delta.h <- function(delta.hz, delta.hinf, c.h, L3, m , beta, expos)

{
  out <- (delta.hz + delta.hinf * c.h * m * beta *  L3 * expos) / (1 + c.h * m * beta * L3 * expos)
  return(out)
}


#' @title
#' Wplus1.rate
#' @description
#' Individual rate of acquisition of new infections (male and female adult worms) in humans
#' @param delta.hz proportion of L3 larvae establishing/ developing to adult stage within human host, per bit, when ATP tends to 0
#' @param delta.hinf proportion of L3 larvae establishing/ developing to adult stage within human host, per bit, when ATP tends to infinity
#' @param c.h severity of transmission intensity - dependent parasite establishment within humans
#' @param L3 mean number of L3 in fly population
#' @param m  vector to human host ratio
#' @param beta per blackfly biting rate on humans (product of proportion of blackfly bites take on humans (h) & reciprocal of the duration of the gonotrophic cycle (g))
#' @param expos individuals total (age/sex-specific and individual-specific) exposure to blackfly bites
#' @param DT timestep
#'
#' @returns vector of acquisition rates (probabilities) per individual in the population
Wplus1.rate <- function(delta.hz, delta.hinf, c.h, L3, m , beta, expos, DT)

{
  dh <- delta.h(delta.hz, delta.hinf, c.h, L3, m , beta, expos) # matt: density-dep (L3 establishing in human host)

  out <- DT * m * beta * dh * expos * L3 # matt: individual rate of acquisition of male / female adult worms - distrecte-time stochastic process with this probablility (pg.7 SI) ?

  return(out)
}

#' @title
#' weibull.mortality
#' @description
#' age-dependent mortality for adult worms and microfilariae (mortality rates assumed to increase as a funcction of parasite age, according to a Weibull distribution of survival times)
#' @param DT timestep
#' @param par1 shape parameter 1 in Weibull distribution (y_l in Hamley et al. 2019 supplementary material; equations S6 and S7)
#' @param par2 shape parameter 2 in Weibull distribution (d_l in Hamley et al. 2019 supplementary material; equations S6 and S7)
#' @param age.cats vector of age categories for mf or adults worms (default is 21)
#'
#' @returns vector of mortality rates (mf or adult worms) for each age class/category
weibull.mortality <- function(DT, par1, par2, age.cats)

{
  out <- DT  * (par1 ^ par2) * par2 * (age.cats ^ (par2-1))

  return(out)
}


#' @title
#' change.worm.per.ind1
#' @description
#' Extracting key columns from overall matrix & evaluating different worm categories (for a specific age class of worms in each individual population).
#' These vectors are important for the next functions to calculate change in number of adults in one adult worm age class for all people.
#' @param treat.vec vector of treatment values (1 or 0) for each individual in the population
#' @param lambda.zero per-capita rate that female worms lose their fertility (W_FF) & return to non-fertile state (W_FN); default value of 0.33 year-1
#' @param DT timestep
#' @param omeg per-capita rate that female worms progress from non-fertile (W_FN) to fertile (W_FF); default is 0.59 year-1
#' @param ws starting worm compartment (first age class) - this is column 28 in the main matrix
#' @param compartment which worm age-compartment (iteration k)
#' @param total.dat main matrix containing different worm compartments for age-classes
#' @param mort.rates vector of mortality rates (adult worms) for each age class/category
#' @param time.each.comp Duration of each age class for adult worms (q_W in Hamley et al. 2019 supp); 1 year default
#' @param new.worms.m vector of (binomial) drawn adult male worms from last column of delay matrix
#' @param w.f.l.c vector for worms coming from previous compartment (0 when k ==1)
#' @param num.comps number of age compartments for each worm category
#'
#' @returns list containing vector/values of a) lambda values per individual, b) mortality value for females per individual c) buman pop size
#' d) number of current non-fertile female worms per individual e) number of current fertile female worms per individual
#' f) treatment or no treatment value per individual g) omega value per individual h) total male worms per individual i) male worms lost per individual
change.worm.per.ind1 <- function(treat.vec, lambda.zero, DT, omeg, ws, compartment, total.dat,
                                 mort.rates, time.each.comp, new.worms.m, w.f.l.c, num.comps)
{
  N <- length(treat.vec)

  lambda.zero.in <- rep(lambda.zero * DT, N) #loss of fertility
  omeg <- rep(omeg * DT, N) #becoming fertile

  #male worms

  cl <- (ws-1) + compartment #calculate which column to use depending on sex, type (fertile or infertile) and compartment

  cur.Wm <- total.dat[, cl] #take current number of worms from matrix

  worm.dead.males <- rbinom(N, cur.Wm, rep(mort.rates[compartment], N))
  worm.loss.males <- rbinom(N, (cur.Wm - worm.dead.males), rep((DT / time.each.comp), N))


  if(compartment == 1)

  {
    male.tot.worms <- cur.Wm + new.worms.m - worm.loss.males - worm.dead.males
  }

  if(compartment > 1)

  {
    male.tot.worms <- cur.Wm + w.f.l.c[[2]] - worm.loss.males - worm.dead.males
  }

  #female worms

  clnf <- (ws - 1) + num.comps + compartment #column for infertile females, num.comps skips over males

  clf <- (ws - 1) + 2*num.comps + compartment #column for fertile females, 2*num.comps skips over males and infertile females

  cur.Wm.nf <- total.dat[, clnf] #take current worm number from matrix infertile females

  cur.Wm.f <- total.dat[, clf] #take current worm number from matrix fertile females

  mort.fems <- rep(mort.rates[compartment], N)


  return(list(lambda.zero.in, mort.fems, N, cur.Wm.nf, cur.Wm.f, treat.vec,
              omeg, male.tot.worms, worm.loss.males))
}


#' @title
#' change.worm.per.ind.treat
#' @description
#' updating worm compartments based on treatment (if this occurs) i.e., considering additional mortality and loss of fertility due to ivermectin treatment;
#' approach assumes individuals which are moved from fertile to non and fertile class due to treatment re-enter fertile class at standard rate
#' @param give.treat value (1 or 0) indicating whether treatment occurs
#' @param iteration iteration along vector of times
#' @param treat.start iteration where treatment starts
#' @param times.of.treat sequence of treatment times
#' @param treat.stop iteration where treatment stops
#' @param onchosim.cov vector of which individuals will be treated if treatment is given
#' @param treat.vec vector of times since treatment for each individual
#' @param DT timestep
#' @param cum.infer proportion of adult female worms made permanently infertility due to ivermectin at each round (lambda'_p in Hamley et al. 2019 supp); default is 0.345
#' @param lam.m embryostatic effects of ivermectin; lam.m is the max rate of transient treatment-induced sterility (lambda_max in Hamley et al. 2019 supp): default is 32.4 year-1
#' @param phi rate of decay of ivermectin-induced female worm sterlisation; default is 19.6 year-1
#' @param N total human population size
#' @param mort.fems vector of mortality rates (adult worms) for each age class/category; updated in change.worm.per.ind1
#' @param lambda.zero.in updated vector of per-capita rate that female worms lose their fertility & return to non-fertile state
#'
#' @returns list containing vector of a) updated lambda values per individual, b) updated (or not) vector of times since treatment for each individual
change.worm.per.ind.treat <- function(give.treat, iteration, treat.start, times.of.treat, treat.stop,
                                 onchosim.cov, treat.vec, DT, cum.infer, lam.m, phi, N,
                                 mort.fems, lambda.zero.in)
{

  if(give.treat == 1 & iteration >= treat.start)
  {

    if((sum(times.of.treat == iteration) == 1) & iteration <= treat.stop) #is it a treatment time

    {
      #print('TREATMENT GIVEN')

      inds.to.treat <- which(onchosim.cov == 1) #which individuals will received treatment

      treat.vec[inds.to.treat]  <-  (iteration-1) * DT #alter time since treatment
      #cum.infer is the proportion of female worms made permanently infertile, killed for simplicity
      if(iteration > treat.start) {mort.fems[inds.to.treat] <- mort.fems[inds.to.treat] + (cum.infer)} #alter mortality
    }


    tao <- ((iteration-1)*DT) - treat.vec #vector of toas, some will be NA

    lam.m.temp <- rep(0, N); lam.m.temp[which(is.na(treat.vec) != TRUE)] <- lam.m #individuals which have been treated get additional infertility rate

    f.to.nf.rate <- DT * (lam.m.temp * exp(-phi * tao)) #account for time since treatment (fertility reduction decays with time since treatment)

    f.to.nf.rate[which(is.na(treat.vec) == TRUE)] <- 0 #these entries in f.to.nf.rate will be NA, lambda.zero.in cannot be NA

    lambda.zero.in <- lambda.zero.in + f.to.nf.rate #update 'standard' fertile to non fertile rate to account for treatment

  }

  return(list(lambda.zero.in, treat.vec, mort.fems))

}

#' @title
#' change.worm.per.ind2
#' @description
#' final function to calculate the change in the number of adult worms in one adult worm age class for all people
#' @param DT timestep
#' @param time.each.comp Duration of each age class for adult worms (q_W in Hamley et al. 2019 supp); 1 year default
#' @param compartment which worm age-compartment (iteration k)
#' @param new.worms.nf.fo new non-fertile female worms (binomial draw)
#' @param w.f.l.c vector for worms coming from previous compartment (0 when k ==1)
#' @param N total human population size
#' @param cur.Wm.nf number of current non-fertile female worms per individual
#' @param mort.fems mortality value for females per individual (from change.worm.per.ind1 function)
#' @param cur.Wm.f number of current fertile female worms per individual
#' @param omeg omega value per individual (from change.worm.per.ind1 function)
#' @param male.tot.worms vector of total male worms per individual (from change.worm.per.ind1 function)
#' @param worm.loss.males vector of male worms lost per individual (from change.worm.per.ind1 function)
#' @param lambda.zero.in updated lambda values per individual (updated fertile to non fertile rate to account for treatment or non-updated) from change.worm.per.ind.treat
#' @param treat.vec updated (or not) vector of times since treatment for each individual from change.worm.per.ind.treat
#'
#' @returns list containing vector of a) total male worms per individual (not updated on treatment, so from change.worm.per.ind1),
#' b) male worms lost per individual (not updated on treatment, so from change.worm.per.ind1), c) vector of non-fertile female worms per individual
#' d) vector of fertile female worms per individual, e) vector of number of non-fertile worms lost per individual,
#' f) vector of number fertile worms lost per individual, g) updated (or not) vector of times since treatment for each individual (taken from change.worm.per.ind.treat)
change.worm.per.ind2 <- function(DT, time.each.comp, compartment, new.worms.nf.fo, w.f.l.c,
                                 N, cur.Wm.nf, mort.fems, cur.Wm.f, omeg,
                                 male.tot.worms, worm.loss.males,
                                 lambda.zero.in, treat.vec)

{
  worm.dead.nf <- rbinom(N, cur.Wm.nf, mort.fems) #movement to next compartment

  worm.dead.f <- rbinom(N, cur.Wm.f, mort.fems)

  worm.loss.age.nf <- rbinom(N, (cur.Wm.nf - worm.dead.nf), rep((DT / time.each.comp), N))

  worm.loss.age.f <- rbinom(N, (cur.Wm.f - worm.dead.f), rep((DT / time.each.comp), N))


  #calculate worms moving between fertile and non fertile, deaths and aging

  #females from fertile to infertile

  new.worms.nf.fi <- rep(0, N)

  trans.fc <- which((cur.Wm.f - worm.dead.f - worm.loss.age.f) > 0)

  #individuals which still have fertile worms in an age compartment after death and aging
  if(length(trans.fc) > 0)
  {
    new.worms.nf.fi[trans.fc] <- rbinom(length(trans.fc), (cur.Wm.f[trans.fc] - worm.dead.f[trans.fc] - worm.loss.age.f[trans.fc]), lambda.zero.in[trans.fc])
  }


  #females worms from infertile to fertile, this happens independent of males, but production of mf depends on males

  #individuals which still have non fertile worms in an age compartment after death and aging
  new.worms.f.fi <- rep(0, N)

  trans.fc <-  which((cur.Wm.nf - worm.dead.nf - worm.loss.age.nf) > 0)
  if(length(trans.fc) > 0)
  {
    new.worms.f.fi[trans.fc] <- rbinom(length(trans.fc), (cur.Wm.nf[trans.fc] - worm.dead.nf[trans.fc] - worm.loss.age.nf[trans.fc]), omeg[trans.fc])#females moving from infertile to fertile
  }


  if(compartment == 1) #if it's the first adult worm age compartment

  {
    nf.out <- cur.Wm.nf + new.worms.nf.fo + new.worms.nf.fi - worm.loss.age.nf - new.worms.f.fi - worm.dead.nf #final number of infertile worms

    f.out <- cur.Wm.f + new.worms.f.fi - worm.loss.age.f - new.worms.nf.fi - worm.dead.f #final number of fertile worms
  }

  if(compartment > 1)

  {
    nf.out <- cur.Wm.nf + new.worms.nf.fi - worm.loss.age.nf - new.worms.f.fi + w.f.l.c[[5]] - worm.dead.nf#w.f.l.c = worms from previous compartment

    f.out <- cur.Wm.f + new.worms.f.fi - worm.loss.age.f - new.worms.nf.fi + w.f.l.c[[6]] - worm.dead.f
  }


  return(list(male.tot.worms,
              worm.loss.males,
              nf.out,
              f.out,
              worm.loss.age.nf,
              worm.loss.age.f, treat.vec))

}




#### Current file: R/epilepsy_module_functions.R 

#' @title
#' OAE probability function
#' @description
#' calculate onchocerciasis-associated epilepsy (OAE) probabilities for mean mf counts
#'
#' @param dat this is the data from Chesnais et al. 2018 with OAE probabilities for 7 mean mf loads in Cameroon
#'
#' @returns vector of OAE probabilities for mf counts from 0 to 1000
OAE_mfcount_prob_func <- function(dat){

  # ep_mfcount_prob_df <- data.frame(prob = c(0.0061, 0.0439, 0.0720, 0.0849, 0.1341, 0.1538, 0.20),
  #                                  mf = c(0, 3, 13, 36, 76, 151, 200)) # data from Chesnais et al. 2018


  fit <- glm(prob ~ log(mf + 1), family = gaussian(link = "log"),
             data = dat) # fitted logarithmic relationship between prob OAE ~ mean mf load;
                                         # generalized linear model with log-link function

  newdat <- data.frame(mf = seq(0, 10000, 1)) # generate new dataframe with mf counts from 0 to 1000 by 1

  out <- predict(fit, newdata = newdat, se.fit = T) # use predict function to calculate log(OAE prob) for each mean mf count in newdat

  logpred <- data.frame(fit = out$fit, se = out$se.fit, mf = newdat$mf) # create dataframe with predicted log(OAE probs), standard error and mf counts

  logpred <- within(logpred, {
    lwr <- fit - 1.96 * se
    upr <- fit + 1.96 * se
  }) # calculate lwr and upper confidence intervals for log(OAE probs)

  pred <- data.frame(fit = exp(logpred$fit), lwr = exp(logpred$lwr), upr = exp(logpred$upr), mf = logpred$mf) # exponentiate probabilities

  OAE_probs <- pred$fit # extract probabilities into a vector

  return(OAE_probs)

}


#' @title
#' find individuals which develop OAE
#'
#' @description
#' find individuals sampled between specific ages (default is 3 - 10 yrs) with new worm infection
#' and can then develop onchocerciasis-associated epilepsy (OAE) in the next function
#' (determined by probability associated with individual mf load)
#'
#' @param dat this is the master matrix (all.mats.temp) tracking mf and worm compartments (by age) for each individual
#' @param mf.start starting column for mf in master matrix
#' @param mf.end final column for mf in master matrix
#' @param worms.start starting column for adult worms in master matrix
#' @param nfw.start # start of infertile worms
#' @param fw.end # end of fertile worms
#' @param infected_at_all vector of all individuals to identify new infections for OAE disease
#' @param age_to_samp sample of ages for each individual (default between 3 - 10 yrs of age) to match to ages in main matrix
#' @param OAE vector of individuals defined as able to develop to OAE/ to determine whether OAE onset occurs
#' @param tested_OAE vector to specify which individuals have been tested
#'
#' @returns vector with positions of new cases in the population
find_indiv_OAE_func <- function(dat, mf.start, mf.end, worms.start, nfw.start, fw.end,
                                infected_at_all, age_to_samp, OAE, tested_OAE, check_ind){

  mf_all <- rowSums(dat[, mf.start : mf.end]) # sum mf per individual across all 21 age classes

  # worms_all <- rowSums(dat[, worms.start : tot.worms]) # OLD: sum all worm categories per host (across 21 age classes)

  male_worms_all <- rowSums(dat[, worms.start : nfw.start - 1]) # sum all male worms per host (across 21 age classes)
  female_worms_all <- rowSums(dat[, nfw.start : fw.end]) # sum all female worms per host (across 21 age classes)

  # ind_new_inf <- which(worms_all > 0 & infected_at_all == 0) # find individuals with at least 1 worm & no previous infection (therefore new)

  ind_new_inf <- which(male_worms_all > 0 & female_worms_all > 0 & infected_at_all == 0) # find individuals where both male and female worms present
                                                                                         # (i.e. pair of mating worms present in host)

  # print(ind_new_inf)

  infected_at_all[ind_new_inf] <-  1 # where new worm infection, give value 1 at individual position in vector

  ind_age <- which(round(dat[, 2]) == round(age_to_samp)) # find individuals (position in individual vector) where (rounded)
                                                          # age at sampling (e.g. between 3 - 10) matches (rounded) age of individual

  # print(ind_age)

  ind_ibf <- which(infected_at_all == 1) # extract position (vector) of all those infected with a pair of mating worms

  ind_no_OAE <- which(OAE != 1) # extract position of individuals not with OAE (so the same people aren't counted twice)

  ind_no_samp <- which(tested_OAE == 0) # extract position of individuals not previously tested

  # print(age_to_samp)

  tot_ind_ep_samp <- Reduce(intersect, list(ind_age, ind_ibf, ind_no_OAE, ind_no_samp)) # goes through each of the individual lists
                                                                                        # (individual criteria for whether a person)
                                                                                        # If individual has satisfied all those (intersect finds commonalities in columns)
                                                                                        # form pool of individuals which are then sampled for new OAE onset
  check_ind <- c(check_ind, length(tot_ind_ep_samp)) #

  return(list(tot_ind_ep_samp, infected_at_all, check_ind))

}

#' @title
#' calculate OAE incidence and prevalence
#'
#' @description
#' each individual undergoes a Bernoulli trial (using OAE probabilities based on mf counts; new worm infection)
#' to ascertain new incident cases of OAE (indexed by those individuals with new infection?)
#'
#' @param temp.mf vector of all mf per skin snip for each individual
#' @param tot_ind_ep_samp vector with positions of new cases in the population
#' @param OAE_probs vector of OAE probabilities for mf counts from 0 to 1000
#' @param dat main matrix tracking age, sex, mf age compartments and worm age compartments for each individual
#' @param prev_OAE vector containing prevalence of OAE to update
#' @param OAE vector of individuals defined as able to develop to OAE/ to determine whether OAE onset occurs
#' @param tested_OAE vector to specify which individuals have been tested
#' @param OAE_incidence_DT OAE incidence vector (whole population) to update
#' @param OAE_incidence_DT_under_5 OAE incidence vector in under 5 years (3 - 5 years) to update
#' @param OAE_incidence_DT_5_10 OAE incidence vector in 5 - 10 years to update
#' @param OAE_incidence_DT_11_15 OAE incidence vector in 10 - 15 years to update
#' @param OAE_incidence_DT_M OAE incidence vector in males to update
#' @param OAE_incidence_DT_F OAE incidence vector in females to update
#'
#' @returns list containing i) updated OAE prevalence, ii) updated OAE incidence in whole population, iii - vi) updated OAE
#' incidence in under 5 years, 5 - 10 years, 10 - 15 years, and males and females, vii) OAE probabilities for each individual based on mf count,
#' viii) those tested/sampled
new_OAE_cases_func <- function(temp.mf, tot_ind_ep_samp, OAE_probs, dat,
                               prev_OAE, OAE, tested_OAE,
                               OAE_incidence_DT, OAE_incidence_DT_under_5, OAE_incidence_DT_5_10, OAE_incidence_DT_11_15,
                               OAE_incidence_DT_M, OAE_incidence_DT_F){


  mf_round <- round(temp.mf[[2]][tot_ind_ep_samp]) + 1 # mf count (in those tested)
                                                       # note: if zero the first probability comes from from ep_probs
                                                       # note: +1 to account for (log) 0


  OAE_rates <- OAE_probs[mf_round] # get probabilities (based on individual mf count)

  OAE[tot_ind_ep_samp] <- rbinom(length(tot_ind_ep_samp), 1, OAE_rates) # for individuals (those positions of individuals as new cases)
                                                                        # with an probability of OAE, they undergo a
                                                                        # Bernoulli trial to ascertain whether OAE onset realized

  tested_OAE[tot_ind_ep_samp] <- 1 # records that they've been tested (those identified in tot_ind_ep_samp are indexed)

  prev_OAE <- c(prev_OAE, mean(OAE)) # calculate & update prevalence of OAE

  new_inc <- length(which(OAE[tot_ind_ep_samp] == 1)) # how many infections (finds total new infected/OAE in all OAE)

  OAE_incidence_DT <- c(OAE_incidence_DT, new_inc) # record + update number of new OAE cases

  new_inc_under_5 <- length(which(OAE[tot_ind_ep_samp] == 1 & dat[tot_ind_ep_samp ,2] >= 3 & dat[tot_ind_ep_samp ,2]< 5 )) # new cases in under 5 age group
  new_inc_5_10 <- length(which(OAE[tot_ind_ep_samp] == 1 & dat[tot_ind_ep_samp ,2] >= 5 & dat[tot_ind_ep_samp ,2]<= 10 )) # new cases in 5 to 10 age group
  new_inc_11_15 <- length(which(OAE[tot_ind_ep_samp] == 1 & dat[tot_ind_ep_samp ,2] >= 10 & dat[tot_ind_ep_samp ,2]<= 15 )) # new cases in 10 to 15 age group

  new_inc_M <- length(which(OAE[tot_ind_ep_samp] == 1 & dat[tot_ind_ep_samp ,3] == 1)) # new cases in males
  new_inc_F <- length(which(OAE[tot_ind_ep_samp] == 1 & dat[tot_ind_ep_samp ,3] == 0)) # new cases in females

  OAE_incidence_DT_under_5 <- c(OAE_incidence_DT_under_5, new_inc_under_5) # record & update incidence in under 5 age group
  OAE_incidence_DT_5_10 <- c(OAE_incidence_DT_5_10, new_inc_5_10) # record & update incidence in 5 to 10 age group
  OAE_incidence_DT_11_15 <- c(OAE_incidence_DT_11_15, new_inc_11_15) # record & update incidence in 10 to 15 age group

  OAE_incidence_DT_M <- c(OAE_incidence_DT_M, new_inc_M) # record & update incidence in males
  OAE_incidence_DT_F <- c(OAE_incidence_DT_F, new_inc_F) # record & update incidence in females

  return(list(prev_OAE, OAE_incidence_DT, OAE_incidence_DT_under_5, OAE_incidence_DT_5_10, OAE_incidence_DT_11_15,
              OAE_incidence_DT_F, OAE_incidence_DT_M,
              OAE_rates, OAE, tested_OAE))

}


#### Current file: R/Intervention_functions.R 

#' @title
#' Treat individuals function
#' @description
#' function to calculate which individuals will be treated.
#'
#' @param all.dt matrix containing (among other things) information on compliance to treatment (whether or not an individual can be treated) and the age of each individual.
#' @param pncomp proportion non-compliers.
#' @param covrg total coverage of the population.
#' @param N human population size (default = 400).
#'
#' @return matrix of individuals to be treated (?)
os.cov <- function(all.dt, pncomp, covrg, N)

{
  pop.ages <- all.dt[,2] #age of each individual in population

  iny <- which(pop.ages < 5 | all.dt[,1] == 1)

  nc.age <- length(iny) / length(pop.ages)

  covrg <- covrg / (1 - nc.age) # probability a complying individual will be treated

  out.cov <- rep(covrg, length(pop.ages))

  out.cov[iny] <- 0 # non-compliers get probability 0

  f.cov <- rep(0, N)

  r.nums <- runif(N, 0, 1)

  inds.c <- which(r.nums < out.cov)

  f.cov[inds.c] <- 1

  return(f.cov)

}


#### Current file: R/larval_dynamics_functions.R 

#' @title
#' Density-dependence in vector
#' @description
#' proportion of mf per mg developing into infective larvae within the vector
#' @param delta.vo density dependence when microfilarae tend to 0
#' @param c.v severity of constraining density-dep larval development (per dermal mf)
#' @param mf vector with total number of mf in each person (extracted from main matrix, summed over each mf age cat per individual)
#' @param expos vector with total exposure (age/sex-specific and individual) per individual in population
#'
#' @returns vector for density-dependent (delat_V) value for each individual
delta.v <- function(delta.vo, c.v, mf, expos)

{
  out <- delta.vo / (1 + c.v * mf *expos)

  return(out)
}


#' @title
#' calc.L1
#' @description
#' L1 (parasite life stages) dynamics in the fly population, assumed to be at equilibrium (modelled deterministically).
#' @param beta per blackfly biting rate on humans (product of proportion of blackfly bites take on humans (h) & reciprocal of the duration of the gonotrophic cycle (g))
#' @param mf vector with total number of mf in each person (extracted from main matrix, summed over each mf age cat per individual)
#' @param mf.delay.in vector of ??
#' @param expos vector with total exposure (age/sex-specific and individual) per individual in population
#' @param delta.vo density dependence when microfilarae tend to 0
#' @param c.v severity of constraining density-dep larval development (per dermal mf)
#' @param nuone per-capita development rate from L1 to L2 larvae (default 201.6 year-1)
#' @param mu.v per-capita mortality rate of blackfly vectors (default 26 year-1)
#' @param a.v per-capita microfilaria-induced mortality of blackfly vectors (default 0.39 year-1)
#' @param expos.delay vector of values from final (4th) column of matrix for exposure (to fly bites) (for L1 delay)exposure.delay) as each value for individual in population
#'
#' @returns vector of mean L1 per individual
calc.L1 <- function(beta, mf, mf.delay.in, expos, delta.vo, c.v, nuone, mu.v, a.v, expos.delay)

{
  delta.vv <- delta.v(delta.vo, c.v, mf, expos)#density dependent establishment

  out <- (delta.vv * beta * expos *  mf)  / ((mu.v + a.v * mf*expos) + (nuone * exp (-(4/366) * (mu.v + (a.v * mf.delay.in*expos.delay)))))

  return(out)
}


#' @title
#' calc.L2
#' @description
#' L2 (parasite life stages) dynamics in the fly population, assumed to be at equilibrium (modelled deterministically).
#' delay of 4 days for parasites moving from L1 to L2
#'
#' @param nuone per-capita development rate from L1 to L2 larvae (default 201.6 year-1)
#' @param L1.in vector of in initial L1 or L1 from last timestep
#' @param mu.v per-capita mortality rate of blackfly vectors (default 26 year-1)
#' @param nutwo per-capita development rate from L2 to L3 larvae (default 201.6 year-1)
#' @param mf vector with delayed mf (?)
#' @param a.v per-capita microfilaria-induced mortality of blackfly vectors (default 0.39 year-1)
#' @param expos vector of values from final (4th) column of matrix for exposure (to fly bites) (for L1 delay)exposure.delay) as each value for individual in population
#'
#' @returns vector of mean L2 per individual
calc.L2 <- function(nuone, L1.in, mu.v, nutwo, mf, a.v, expos)

{
  out <- (L1.in * (nuone * exp (-(4/366) * (mu.v + (a.v * mf * expos))))) / (mu.v + nutwo)

  return(out)
}


#' @title
#' calc.L3
#' @description
#' L3 (parasite life stages) dynamics in the fly population, assumed to be at equilibrium (modelled deterministically).
#'
#' @param nutwo per-capita development rate from L2 to L3 larvae (default 201.6 year-1)
#' @param L2.in vector of L2 from main matrix for each individual
#' @param a.H proportion of L3 shed per bite (on any blood host); default is 0.8
#' @param g length of the gonotrophic cycle (0.0096)
#' @param mu.v per-capita mortality rate of blackfly vectors (default 26 year-1)
#' @param sigma.L0 per-capita mortality rate of L3 in blackfly vectors
#'
#' @returns vector of mean L3 per individual
calc.L3 <- function(nutwo, L2.in, a.H, g, mu.v, sigma.L0)

{
  out <- (nutwo * L2.in) / ((a.H / g) + mu.v + sigma.L0)

  return(out)
}



#### Current file: R/mf_dynamics_functions.R 

#' @title
#' rotate matrix
#' @description
#' function to rotate matrix, used in mf function.
#'
#' @param x matrix to rotate.
#'
#' @returns rotated matrix (?)
rotate <- function(x) {
  x_rotated <- t(apply(x, 2, rev))

  return(x_rotated)
}


#' @title
#' change microfilariae number per human
#' @description
#' function calculates change in the number of microfilariae (mf) (offspring of adult worms) for each mf age compartment in each human using RK4 method
#' this is called in ep.equi.sim for each mf age class.
#' @param dat main matrix with mf inputs for each individual
#' @param num.comps number of age classes in adult worms (fecundity rate is a function of adult worm age), default = 21 age classes (specified by c_max in Hamley et al. 2019).
#' @param mf.cpt age class under consideration.
#' @param num.mf.comps number of age classes in mf (same as adult worm), default = 21 age classes (specified by c_max in Hamley et al. 2019).
#' @param ws column in main matrix where worm age compartments start for each individual (row in the main matrix).
#' @param DT timestep
#' @param time.each.comp duration fo each age class in mf (default value ins 0.125 year; q_M in Hamley et al. 2019).
#' @param mu.rates.mf mf mortality rate (modelled as a weibull distribution, with mu.mf1 (y_M) and mu.mf2 (d_M) parameters) as a function of age.
#' @param fec.rates fecundity rate in adult female worms (m(a) in Hamley et al. 2019); given by epsilon (1.158), fec.w.1 (F), fec.w.2 (fec.w.2).
#' @param mf.move.rate determines mf aging (moving rate to the next age class); given by 1 / duration spent in each mf age class.
#' @param up constant to allow for very large yet finite mf effect upon treatment (mu in Hamley et al. 2019; default value of 9.6 x 10-3) & determines extent of IVM induced mortality (with kap).
#' @param kap shape parameter for excess mf mortality following IVM treatment (kappa in Hamley et al. 2019; default value of 1.25).
#' @param iteration iteration moving through each timepoint (based on timestep).
#' @param treat.vec vector contains how long since each person treated, mortality rate due to ivermectin decays with time since treatment.
#' @param give.treat input (if 1 then treatment given).
#' @param treat.start treatment start time.
#'
#' @returns vector for each individual in population (to replace specific mf age compartment column for all individuals in main matrix)
change.micro <- function(dat, num.comps, mf.cpt, num.mf.comps, ws, DT, time.each.comp, mu.rates.mf, fec.rates, mf.move.rate,
                         up, kap, iteration, treat.vec, give.treat, treat.start)

{
  N <- length(dat[,1])

  # indexes for fertile worms (to use in production of mf)
  fert.worms.start <-  ws + num.comps * 2
  fert.worms.end <-  (ws - 1) + num.comps * 3

  # indexes to check if there are males (males start is just 'ws')
  # there must be >= 1 male worm for females to produce microfilariae
  mal.worms.end <- (ws - 1) + num.comps
  mf.mu <- rep(mu.rates.mf[mf.cpt], N)
  fert.worms <- dat[, fert.worms.start:fert.worms.end] #number of fertile females worms

  # increases microfilarial mortality if treatment has started
  if(give.treat == 1 & iteration >= treat.start)
  {
    tao <- ((iteration - 1) * DT) - treat.vec # tao is zero if treatment has been given at this timestep

    mu.mf.prime <- ((tao + up) ^ (- kap)) # additional mortality due to ivermectin treatment

    mu.mf.prime[which(is.na(mu.mf.prime) == TRUE)] <- 0

    mf.mu <- mf.mu + mu.mf.prime

  }

  # if the first age class of microfilariae
  if(mf.cpt == 1)
  {
    mp <- rep(0, N)

    inds.fec <- which(rowSums(dat[, ws : mal.worms.end]) > 0); mp[inds.fec] <- 1 # need to check there is more than one male

    k1 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = 0)  #fert worms and epin are vectors
    k2 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k1/2)
    k3 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k2/2)
    k4 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k3)

    out <- dat[, 6 + mf.cpt] + DT/6 * (k1 + 2 * k2 + 2* k3 + k4)

  }

  # if age class of microfilariae is >1
  if(mf.cpt > 1)
  {
    k1 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = 0)
    k2 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k1/2)
    k3 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k2/2)
    k4 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k3)

    out <- dat[, 6 + mf.cpt] + DT/6 * (k1 + 2 * k2 + 2 * k3 + k4)

  }

  return(out)
}

#' @title
#' microfilariae age class 1 derivatives
#' @description
#' function called during RK4 for first age class of microfilariae (mf), which is a function of fertile adult female worms.
#' @param fert.worms matrix (or vector?) of fertile females worms in all age-classes per individual.
#' @param mf.in matrix of mf in first age-class per individual.
#' @param ep.in fecundity of female worms (based on fec.rates in adult female worms (m(a) in Hamley et al. 2019))
#' @param mf.mort mf.mu rate specified for first age class in each human
#' @param mf.move determines mf aging (moving rate to the next age class).
#' @param mp vector (?) length of human population
#' @param k.in previous k (for k1 = 0) for RDK4 method.
#'
#' @returns vector for each individual
derivmf.one <- function(fert.worms, mf.in, ep.in, mf.mort, mf.move, mp, k.in)  #fert worms and epin are vectors
{
  new.in <- (rotate(fert.worms) * ep.in) #need to rotate matrix to each column is multiplied by respective fecundity rate, not each row
  new.in <- rotate(rotate(rotate(new.in)))
  new.in <- rowSums(new.in)

  mort.temp <- mf.mort * (mf.in + k.in)
  move.temp <- mf.move * (mf.in + k.in)

  mort.temp[which(mort.temp < 0)] <- 0
  move.temp [which(move.temp < 0)] <- 0

  if(length(which(mf.mort * (mf.in + k.in) < 0)) > 0) {print('MF NEGATIVE1')}
  if(length(which(mf.move * (mf.in + k.in) < 0)) > 0) {print('MF NEGATIVE2')}


  out <- mp * new.in - mort.temp - move.temp

  return(out)
}


#' @title
#' microfilariae age classes > 1 derivatives
#' @description
#' function called during RK4 for > 1 age classes of microfilariae (mf).
#' @param mf.in matrix of mf in > 1 age classes per individual.
#' @param mf.mort mf.mu rate specified for > 1 age classes in each human.
#' @param mf.move determines mf aging (moving rate to the next age class).
#' @param mf.comp.minus.one previous mf age class (compartment).
#' @param k.in previous k (for k1 = 0) for RDK4 method.
#'
#' @returns vector for each individual
derivmf.rest <- function(mf.in, mf.mort, mf.move, mf.comp.minus.one, k.in)
{
  move.last <- mf.comp.minus.one * mf.move
  mort.temp <- mf.mort * (mf.in + k.in)
  move.temp <- mf.move * (mf.in + k.in)

  move.last[which(move.last < 0)] <- 0
  mort.temp[which(mort.temp < 0)] <- 0
  move.temp [which(move.temp < 0)] <- 0

  if(length(which(mf.mort * (mf.in + k.in) < 0)) > 0) {print('WARNING MF NEGATIVE3')}
  if(length(which(mf.move * (mf.in + k.in) < 0)) > 0) {print('WARNING MF NEGATIVE4')}


  out <- move.last - mort.temp - move.temp

  return(out)
}


#' @title
#' calculate mf per skin snip
#' @description
#' function calculates number of mf in skin snip for all people.
#' people are tested for the presence of mf using a skin snip, we assume mf are overdispersed in the skin.
#' @param ss.wt weight of the skin snip
#' @param num.ss number of skin snips taken (default set to 2)
#' @param slope.kmf slope value governing linear relationship between degree of mf overdispersion and adult female worms
#' @param int.kMf initial value governing linear relationship between degree of mf overdispersion and adult female worms
#' @param data data is the matrix tracking age compartments of mf and W per individual
#' @param nfw.start column (first age compartment) in matrix where non-fertile female worms begin
#' @param fw.end column (last age compartment) in matrix where fertile female worms end
#' @param mf.start column (first age compartment) in matrix where mf begin
#' @param mf.end column (last age compartment) in matrix where mf ends
#' @param pop.size human population size
#' @param kM.const.toggle if set to YES then kM is a constant (default = 15)
#'
#' @returns element (1) in list is mean of mf per skin snip; element (2) contains all mf per skin snip for each individual
mf.per.skin.snip <- function(ss.wt, num.ss, slope.kmf, int.kMf, data, nfw.start, fw.end,  ###check vectorization
                             mf.start, mf.end, pop.size, kM.const.toggle)

{

  all.mfobs <- c()

  if(isTRUE(kM.const.toggle)){
    kmf <- 0 * (rowSums(data[,nfw.start:fw.end])) + 15}
  else {
     kmf <- slope.kmf * (rowSums(data[,nfw.start:fw.end])) + int.kMf #rowSums(da... sums up adult worms for all individuals giving a vector of kmfs
  }

  mfobs <- rnbinom(pop.size, size = kmf, mu = ss.wt * (rowSums(data[,mf.start:mf.end])))

  nans <- which(mfobs == 'NaN'); mfobs[nans] <- 0

  if(num.ss > 1)

  {

    tot.ss.mf <- matrix(, nrow = length(data[,1]), ncol = num.ss) # error?
    tot.ss.mf[,1] <- mfobs

    for(j in 2 : (num.ss)) #could be vectorized

    {

      temp <- rnbinom(pop.size, size = kmf, mu = ss.wt * (rowSums(data[,mf.start:mf.end])))

      nans <- which(temp == 'NaN'); temp[nans] <- 0

      tot.ss.mf[,j] <- temp

    }

    mfobs <- rowSums(tot.ss.mf)

  }

  mfobs <- mfobs / (ss.wt * num.ss)

  list(mean(mfobs), mfobs)

}

#' @title
#' mf prevalence
#' @description
#' calculates mf prevalence in people based on a skin snip
#' @param age minimum age for giving a skin snip (therefore age from which prevalence calculated)
#' @param ss.in takes mf per skin snip count object for each individual to convert to binary variable for prevalence
#' @param main.dat main matrix contain age of each individual (to ensure only calculate prevalence based o age from which skin snips are taken)
#'
#' @returns value for prevalence
prevalence.for.age <- function(age, ss.in, main.dat)

{
  inds <- which(main.dat[,2] >= age)

  out <- length(which(ss.in[[2]][inds] > 0)) / length(inds)

  return(out)
}



#' @title
#' mf prevalence by strata
#' @description
#' calculates mf prevalence in people based on a skin snip for different population strata (age and sex)
#' @param ss.in takes mf per skin snip count object for each individual to convert to binary variable for prevalence
#' @param main.dat main matrix contain age of each individual (to ensure only calculate prevalence based o age from which skin snips are taken)
#' @param lwr_age lower age to measure prevalence from
#' @param upr_age lower age to measure prevalence to
#' @param sex  sex of strata to measure prevalence in
#'
#' @returns value for prevalence
prevalence.for.age_sex.strata <- function(ss.in, main.dat, lwr_age, upr_age, sex)

{
  if(sex == "male"){
    sex_ind <- 1
  } else {
    sex_ind <- 0
  }

  inds <- which(main.dat[,2] >= lwr_age & main.dat[,2] < upr_age & main.dat[,3] == sex_ind)

  out <- length(which(ss.in[[2]][inds] > 0)) / length(inds)

  return(out)
}


#' @title
#' mf prevalence by strata (including compliance)
#' @description
#' calculates mf prevalence in people based on a skin snip for different population strata (age, sex and compliance)
#' @param ss.in takes mf per skin snip count object for each individual to convert to binary variable for prevalence
#' @param main.dat main matrix contain age of each individual (to ensure only calculate prevalence based o age from which skin snips are taken)
#' @param lwr_age lower age to measure prevalence from
#' @param upr_age lower age to measure prevalence to
#' @param sex  sex of strata to measure prevalence in
#'
#' @returns value for prevalence
prevalence.for.age_sex_compl.strata <- function(ss.in, main.dat, lwr_age, upr_age, sex, compliance)

{
  if(sex == "male"){
    sex_ind <- 1
  } else {
    sex_ind <- 0
  }

  inds <- which(main.dat[,2] >= lwr_age & main.dat[,2] < upr_age)

  out <- length(which(ss.in[[2]][inds] > 0)) / length(inds)

  return(out)
}





#### Current file: R/Model_wrappers.R 

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
#' @param Q Ratio of Male Exposure vs Female Exposure
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
                        N.in = 400,
                        calc_ov16=FALSE,
                        ov16_equilibrium=NA,
                        ov16_store_times = c(),
                        no_prev_run=FALSE,
                        custom_treat_params = list(),
                        seroreversion = "none",
                        Q = 1.2)


{
  # ====================== #
  # Set-up time parameters #

  DT <- 1/366
  time.its <- round(time.its / (DT))
  year_its <- seq(0, time.its, 366)
  # if(give.treat == 1) #calculate timesteps at which treatment is given
  # {times.of.treat.in <- seq(treat.start, treat.stop - (treat.int / DT), treat.int / DT)}
  # else {times.of.treat.in <- 0}

  if(length(ov16_store_times) > 0) {
    ov16_store_times <- round(ov16_store_times / (DT))
  }
  
  if(give.treat == 1)
  {
    treat.stop <- round(treat.stop / (DT))
    if(treat.start >= 1) {treat.start <-  round( (treat.start) / (DT)) + 1}
    if(treat.start == 0) {treat.start <-  1}
    if(length(ov16_store_times) == 0) {
      ov16_store_times <- c(treat.start-1, treat.stop, treat.stop+(1/DT)+1)
    }
  } else {
    if(length(ov16_store_times) == 0) {
      ov16_store_times <- c(time.its-1)
    }
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
  N = N.in #human population size
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

  # delete E0 and q
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

      Ov16_Seropositive_serorevert <- rep(0, N)
      Ov16_Seropositive_L3_serorevert <- 
      Ov16_Seropositive_L4_serorevert <- rep(0, N)
      Ov16_Seropositive_mating_no_mf_serorevert <- rep(0, N)
      Ov16_Seropositive_mating_detectable_mf_serorevert <- rep(0, N)
      Ov16_Seropositive_mating_any_mf_serorevert <- rep(0, N)


      mf_indv_prev <- rep(0, N)


      Ov16_Seropositive_matrix <- matrix(0, nrow=N, ncol=length(ov16_store_times)*9)
      Ov16_Seropositive_Serorevert_matrix <- matrix(0, nrow=N, ncol=length(ov16_store_times)*9)
      matrix_index <- 1

      prev_Ov16 <- 0

      # 80% of pop is able to mount Antibody response to Ov16
      all.mats.temp[,num.cols+ov16.col] <- 1#sample(rep(c(1,1,1,1,1,1,1,1,0,0), N*0.1), N)
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

      Ov16_Seropositive_serorevert <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive_serorevert
      Ov16_Seropositive_L3_serorevert <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive_l3_serorevert
      Ov16_Seropositive_L4_serorevert <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive_l4_serorevert
      Ov16_Seropositive_mating_no_mf_serorevert <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive_mating_no_mf_serorevert
      Ov16_Seropositive_mating_detectable_mf_serorevert <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive_mating_detectable_mf_serorevert
      Ov16_Seropositive_mating_any_mf_serorevert <- ov16_equilibrium$ov16_seropos_outputs$ov16_seropositive_mating_any_mf_serorevert

      mf_indv_prev <- ov16_equilibrium$mf_indv_prev

      Ov16_Seropositive_matrix <- matrix(0, nrow=N, ncol=length(ov16_store_times)*9)
      Ov16_Seropositive_Serorevert_matrix <- matrix(0, nrow=N, ncol=length(ov16_store_times)*9)
      matrix_index <- 1

      prev_Ov16 <- sum(Ov16_Seropositive_mating_any_mf_serorevert)/N
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

    e.f <- 1/(sex.rat*(Q-1) + 1) # 0.9
    e.m <- Q * e.f # 1.08

    s.a.exp <- rep(0, N)

    alpha.m <- exp(-age.exp.m * (all.mats.cur[mls, 2]))
    gamma.m <- 1/mean(alpha.m)
    s.a.exp[mls] <- alpha.m * e.m * gamma.m


    alpha.f <- exp(-age.exp.f * (all.mats.cur[fmls, 2]))
    gamma.f <- 1/mean(alpha.f)
    s.a.exp[fmls] <- alpha.f * e.f * gamma.f

    norm.ex.vec <- ex.vec * (1 / mean(ex.vec)) #normalize so mean = 1 (matt: normalising the indvidual-specific exposure from line 565)

    tot.ex.ai <- s.a.exp * norm.ex.vec # matt: combine sex/age specific exposure + individual specific exposure (total exposure to blackfly bites)
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

    if(calc_ov16) {
      # exposure checks
      any_juvy_worms <- (rowSums(all.mats.temp[, worms.start:num.cols]) > 0)
      any_l3_exposure <- (l.extras[,1] > 0)
      any_larvae <- (rowSums(l.extras) > 0)
      l4_development <- (l.extras[,floor(length(l.extras[1,])/2)] > 0)
      any_worms <- rowSums(all.mats.temp[,worms.start:fw.end])
      mating_worm <- ((rowSums(all.mats.temp[,worms.start:nfw.start])) > 0 & (rowSums(all.mats.temp[, fw.start:fw.end]) > 0))
      mating_worm_detectable_mf <- (mating_worm & temp.mf[[2]] > 0)
      mating_worm_any_mf <- (mating_worm & (rowSums(all.mats.temp[,mf.start:mf.end]) > 0))

      indv_antibody_response <- all.mats.temp[,91]

      findPositives <- function(exposure_array, curr_array, antibody_resp, doSerorevert=FALSE) {
        curr_array[which(exposure_array == TRUE & curr_array == 0 & antibody_resp == 1)] <- 1
        # hard seroreversion
        if(doSerorevert & seroreversion == "no_infection") {
          curr_array[which(curr_array == 1 & any_larvae == FALSE & exposure_array == FALSE & any_worms == FALSE & rowSums(all.mats.temp[,mf.start:mf.end]) == 0)] <- 0
        }
        if(doSerorevert & seroreversion == "absence_of_trigger") {
          curr_array[which(curr_array == 1 & exposure_array == FALSE)] <- 0
        }
        return(curr_array)
      }

      Ov16_Seropositive <- findPositives(any_juvy_worms, Ov16_Seropositive, indv_antibody_response)
      Ov16_Seropositive_L3 <- findPositives(any_l3_exposure, Ov16_Seropositive_L3, indv_antibody_response)
      Ov16_Seropositive_L4 <- findPositives(l4_development, Ov16_Seropositive_L4, indv_antibody_response)
      Ov16_Seropositive_mating_no_mf <- findPositives(mating_worm, Ov16_Seropositive_mating_no_mf, indv_antibody_response)
      Ov16_Seropositive_mating_detectable_mf <- findPositives(mating_worm_detectable_mf, Ov16_Seropositive_mating_detectable_mf, indv_antibody_response)
      Ov16_Seropositive_mating_any_mf <- findPositives(mating_worm_any_mf, Ov16_Seropositive_mating_any_mf, indv_antibody_response)

      Ov16_Seropositive_serorevert <- findPositives(any_juvy_worms, Ov16_Seropositive_serorevert, indv_antibody_response, doSerorevert=TRUE)
      Ov16_Seropositive_L3_serorevert <- findPositives(any_l3_exposure, Ov16_Seropositive_L3_serorevert, indv_antibody_response, doSerorevert=TRUE)
      Ov16_Seropositive_L4_serorevert <- findPositives(l4_development, Ov16_Seropositive_L4_serorevert, indv_antibody_response, doSerorevert=TRUE)
      Ov16_Seropositive_mating_no_mf_serorevert <- findPositives(mating_worm, Ov16_Seropositive_mating_no_mf_serorevert, indv_antibody_response, doSerorevert=TRUE)
      Ov16_Seropositive_mating_detectable_mf_serorevert <- findPositives(mating_worm_detectable_mf, Ov16_Seropositive_mating_detectable_mf_serorevert, indv_antibody_response, doSerorevert=TRUE)
      Ov16_Seropositive_mating_any_mf_serorevert <- findPositives(mating_worm_any_mf, Ov16_Seropositive_mating_any_mf_serorevert, indv_antibody_response, doSerorevert=TRUE)

      mf_indv_prev <- as.integer(temp.mf[[2]] > 0)
      prev_Ov16 <- c(prev_Ov16, sum(Ov16_Seropositive_mating_any_mf_serorevert)/N)
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

        Ov16_Seropositive_serorevert[to.die] <- 0
        Ov16_Seropositive_L3_serorevert[to.die] <- 0
        Ov16_Seropositive_L4_serorevert[to.die] <- 0
        Ov16_Seropositive_mating_no_mf_serorevert[to.die] <- 0
        Ov16_Seropositive_mating_detectable_mf_serorevert[to.die] <- 0
        Ov16_Seropositive_mating_any_mf_serorevert[to.die] <- 0
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
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-5] <- Ov16_Seropositive_serorevert
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-4] <- Ov16_Seropositive_L3_serorevert
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-3] <- Ov16_Seropositive_L4_serorevert
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-2] <- Ov16_Seropositive_mating_no_mf_serorevert
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index-1] <- Ov16_Seropositive_mating_detectable_mf_serorevert
      Ov16_Seropositive_Serorevert_matrix[,9*matrix_index] <- Ov16_Seropositive_mating_any_mf_serorevert

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
        ov16_seropos_outputs <- list(prev_Ov16, mf_indv_prev, Ov16_Seropositive_matrix, Ov16_Seropositive_Serorevert_matrix)
        names(ov16_seropos_outputs) <- c('ov16_seroprevalence', 'mf_indv_prevalence', 'ov16_seropositive_matrix', 'ov16_seropositive_matrix_serorevert')
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
        ov16_output <- list(prev_Ov16, mf_indv_prev, Ov16_Seropositive_matrix, Ov16_Seropositive_Serorevert_matrix)
        names(ov16_output) <- c('ov16_seroprevalence', 'mf_indv_prevalence', 'ov16_seropositive_matrix', 'ov16_seropositive_matrix_serorevert')
        outp <- append(outp, ov16_output)
      }
      return(outp)
    }
  }


}


#### Current file: runModelRCSTogo.R 

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
