# ============================ #
# at set-up / before first MDA #

# treat.start = 1
# treat.stop = 16
# DT = 1/366
# treat.int = 1
# treat.stop <- round(treat.stop / (DT))
# if(treat.start >= 1) {treat.start <-  round( (treat.start) / (DT)) + 1}
# if(treat.start == 0) {treat.start <-  1}
#
# times.of.treat.in <- seq(treat.start, treat.stop, treat.int / DT)
#
# i <- 5125
#
# if(i == times.of.treat.in[1]){
#   cat(i, "is equal to the first value in the vector")
# }
#
# if (i %in% times.of.treat.in[-1]) {
#   cat(i, "is equal to a value in the vector, except the first one.\n")
# }

N = 440
probneverTreat = 0.1
coverage = 0.6
correlation = 0.2

compliance.mat <- matrix(nrow=N, ncol=4) # col 1 = age, col 2 = never_treat,
                                         # col 3 = probability of treatment, col 4 = to be treated in this round
compliance.mat[,2] = generateNeverTreat(N, probneverTreat) # never treat col (mat[,2])

# individual probability of treatment values
cov = coverage # whatever the coverage of this MDA is
rho = correlation # whatever the correlation of this MDA is
compliance.mat[,3] = initializePTreat(N, cov, rho) # initialize pTreat - correlation for each individual (mat[,3])

# previous coverage #
prevCov = cov # set prevCov to coverage value used
prevRho = rho # set prevRho to correlation value used


# ========================================== #
# at each subsequent MDA (associated iter i) #

# logic if cov or rho updated for next MDA round

#if((prevCov != coverage) | (prevRho != correlation)){

  # 1) check and update/redraw any zero values introduced in pTreat for individuals since last MDA round
  #compliance.mat[20,2] = 0 # test = introduce 1 zero at p20 in pTreat (compliance[,3])
  compliance.mat[,3] = checkForZeroPTreat(pTreat = compliance.mat[,3], prevCov, prevRho)

  # 2) assign everyone a new/updated pTreat value for the next MDA round if cov and/or rho have changed
  cov = coverage
  rho = correlation
  compliance.mat[,3] = editPTreat(pTreat = compliance.mat[,3], cov = 0.75, rho) # update pTreat = test with coverage going up to 0.75
  prevCov = cov
  prevRho = rho

#}







# =================================================== #
#       Functions to call                             #

generateNeverTreat <- function(N, probNeverTreat){
  return(rbinom(N, size = 1, prob = probNeverTreat))
}


initializePTreat <- function(N, cov, rho) {
  # Calculate alpha and beta parameters for the beta distribution
  alpha = cov * (1 - rho) / rho
  beta = (1 - cov) * (1 - rho) / rho

  # Draw N probabilities from the beta distribution
  pTreat = rbeta(N, alpha, beta)

  # Return the generated pTreat values
  return(pTreat)
}



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
