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


# ========================================================================================================== #
#                               New compliance structure (based on beta-binomial)                            #
# ========================================================================================================== #

#' @title
#' generateNeverTreat function
#' @description
#' function to generate a vector of individuals that will never be treated not matter how many rounds of treatment performed
#' this is an additional proportion on top of the possibility that some proportion of individuals who will never have been
#' treated after n rounds with the new compliance structure
#'
#' @param N human population size (default = 400).
#' @param probNeverTreat proportion of individuals that will never be treated
#'
#' @return vector of length of pop size with individuals never treated (1)
generateNeverTreat <- function(N, probNeverTreat){
  return(rbinom(N, size = 1, prob = probNeverTreat))
}


# generate never Treat individuals
# probneverTreat = 0.1
# neverTreat = generateNeverTreat(N, probneverTreat)

#' @title
#' initializePTreat function
#' @description
#' function to generate a vector of probability of treatment for each individual using the beta-binomial distribution
#' before the first treatment round
#'
#' @param N human population size (default = 400).Matt G: Number of elements.
#' @param cov Coverage parameter (0-1)
#' @param rho rho parameter (0-1) which is the correlation between treatment rounds,
#' where rho = 0 means random participation and rho = 1 means each person will attend round k if they attended first round
#' which is equivalent to the systematic compliance scheme (same individuals at each round)
#'
#' @return vector of length of pop size with individuals probability of treatment (0-1)
initializePTreat <- function(N, cov, rho) {
  # Calculate alpha and beta parameters for the beta distribution
  alpha = cov * (1 - rho) / rho
  beta = (1 - cov) * (1 - rho) / rho

  # Draw N probabilities from the beta distribution
  pTreat = rbeta(N, alpha, beta)

  # Return the generated pTreat values
  return(pTreat)
}

# cov = 0.6
# rho = 0.2
# pTreat = initializePTreat(N, cov, rho)


#' @title
#' checkForZeroPTreat function
#' @description
#' function to check for any introduced zeros in the pTreat column (due to introduce of new individuals) and then
#' redrawing pTreat value for these individuals
#'
#' @param pTreat probability of treatment column from compliance.mat (to be checked for zeros)
#' @param cov Coverage parameter (0-1)
#' @param rho rho parameter (0-1) which is the correlation between treatment rounds,
#'
#'
#' @return vector of length of pop size with individuals probability of treatment (0-1)
checkForZeroPTreat <- function(pTreat, cov, rho) {
  # Find indices where pTreat values are zero
  zeros <- (pTreat == 0)

  # If there are any zero values, regenerate them using the beta distribution
  if (any(zeros)) {
    # Calculate alpha and beta parameters for the beta distribution
    alpha = cov * (1 - rho) / rho
    beta = (1 - cov) * (1 - rho) / rho

    # Generate new pTreat values only for the zero positions
    pTreat[zeros] = rbeta(sum(zeros), alpha, beta)
  }

  # Return the updated pTreat vector
  return(pTreat)
}


#' @title
#' editPTreat function
#' @description
#' update the pTreat for each individual when the cov and/or rho value changes between MDA round.
#' draw new probabilities (newPTreats) from beta-binomal and sort (ascending) them,
#' rank the original pTreat values (indices) for sorting i.e., 440 index is highest prob,
#' and finally assign values from newPTreats to correct location in pTreat to update this
#'
#' @param pTreat probability of treatment column from compliance.mat (to be updated)
#' @param cov Coverage parameter (0-1)
#' @param rho rho parameter (0-1) which is the correlation between treatment rounds,
#'
#' @return vector of length of pop size with individuals probability of treatment (0-1)
editPTreat <- function(pTreat, cov, rho) {
  # Define the parameters for the probability of treatment
  alpha <- cov * (1 - rho) / rho
  beta <- (1 - cov) * (1 - rho) / rho

  # Get the number of elements in pTreat
  N = length(pTreat)

  # Draw new probabilities from the beta distribution and sort them
  newPTreats <- sort(rbeta(N, alpha, beta))

  # Rank the original pTreat values to get indices for sorting
  indices <- rank(pTreat)

  # Assign the values from newPTreats to the appropriate places in pTreat vector
  pTreat <- newPTreats[indices]
  return(pTreat)
}


#' @title
#' check eligibility of individuals for treatment function
#' @description
#' Eligible individuals which satisfy conditions of correct age (=> 5 yrs),
#' are not considered a never_treated individual, and probability of treatment
#' (pTreat) > than random number drawn from uniform distribution
#'
#' @param comp.mat compliance matrix
#' @param all.dt current main matrix with information on individuals (age)
#' @param minAgeMDA minimum age for MDA (set to 5 yrs)
#' @param maxAgeMDA maximum age for MDA (set to 80 yrs)
#'
#' @return updated complicance matrix with column for treatment (1) or not (0),
#' and a vector with positions for individuals to be treated
check_eligibility <- function(comp.mat, all.dt, minAgeMDA, maxAgeMDA) {

  # set up current age for eligibility criteria check for MDA (current age col)
  comp.mat[,1] <- all.dt[,2]

  # Identify hosts in eligible age range (hostsEligibleAge col)
  comp.mat[,4] <- (comp.mat[,1] >= minAgeMDA) & (comp.mat[,1] < maxAgeMDA)

  # Generate all random numbers at once (random number col)
  comp.mat[,5] <- runif(length(comp.mat[,1]))

  # find eligible hosts in matrix (assign as 1) in 6th column (eligible host col)
  comp.mat[,6] <- ifelse(comp.mat[,4] == 1 & comp.mat[,5] < comp.mat[,3] & comp.mat[,2] == 0, 1, 0)

  # Find indices of hosts meeting conditions (if needed in treatment funcs)
  eligible_hosts <- which(comp.mat[,6] == 1)

  return(list(comp.mat, eligible_hosts))

}
