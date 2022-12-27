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
