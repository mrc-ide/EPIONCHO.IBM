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

