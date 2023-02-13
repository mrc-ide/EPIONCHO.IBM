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
#' @param lam.m embryostatic effects of ivermectin; lam.m is the max rate of treatment-induced sterility (lambda_max in Hamley et al. 2019 supp): default is 32.4 year-1
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


