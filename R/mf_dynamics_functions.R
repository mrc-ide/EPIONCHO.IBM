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
#' @param dat
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
#'
#' @returns element (1) in list is mean of mf per skin snip; element (2) contains all mf per skin snip for each individual
mf.per.skin.snip <- function(ss.wt, num.ss, slope.kmf, int.kMf, data, nfw.start, fw.end,  ###check vectorization
                             mf.start, mf.end, pop.size)

{

  all.mfobs <- c()

  kmf <- slope.kmf * (rowSums(data[,nfw.start:fw.end])) + int.kMf #rowSums(da... sums up adult worms for all individuals giving a vector of kmfs

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
