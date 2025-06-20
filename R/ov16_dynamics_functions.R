#' @title
#' determine_serostatus
#'
#' @description
#' Determine new Ov16 positive cases, and new individuals who serorevert
#'
#' @param exposure_array the yearly probability
#' @param curr_array the number of days in a year (355 or 366)
#' @param seroreversion the number of days in a year (355 or 366)
#' @param doSerorevert Calculate seroreversion
#' @param seroreversion_arrays list of arrays to use to see if seroreversion occurs
#'
#' @returns a list of all individuals and their ov16 serostatus (1 if positive, 0 if negative)
determine_serostatus <- function(
  exposure_array, curr_array, do_serorevert = FALSE, seroreversion_arrays = list()
) {
  curr_array[which(exposure_array == TRUE & curr_array == 0)] <- 1

  if (do_serorevert) {
    curr_array[
      which(
        curr_array == 1 &
          seroreversion_arrays[["any_larvae_arr"]] == FALSE &
          seroreversion_arrays[["any_worms_arr"]] == FALSE
      )
    ] <- 0
  }

  return(curr_array)
}

seroprevalence_for_age <- function(
  serostatus, ages, lower_age = 5, upper_age = 81, sensitivity = 1, specificity = 1
) {
  subset_inds <- which(ages >= lower_age & ages < upper_age)
  total_inds <- length(subset_inds)

  prob <- runif(total_inds)
  new_serostatus <- rep(0, total_inds)
  pos <- which(serostatus == 1)
  neg <- which(serostatus == 0)
  if (length(pos) > 0) {
    new_serostatus[pos] <- as.numeric(prob[pos] <= sensitivity)
  }
  if (length(neg) > 0) {
    new_serostatus[neg] <- as.numeric(prob[neg] > specificity)
  }

  seropositives <- which(new_serostatus[subset_inds] == 1)
  return(length(seropositives) / total_inds)
}

#' @title
#' Determine seroprevalence
#'
#' @description
#' Determine seroprevalence for provided age groups, with or without diagnostic adjustment
#'
#' @param serostatus_no_seroreversion the serostatus without seroreversion
#' @param serostatus_finite_seroreversion the serostatus with finite seroreversion
#' @param ages the ages of the population
#' @param age_groups a list of age groups list(c(5, 80), c(..), ..)
#' @param diagnostic_adjustments a vector of length 2 containing sens, spec (0.8, 0.99)
#'
#' @returns vector of seroprevalence for all provided age groups in order
calculate_seroprevalence_across_age_groups <- function(
  serostatus_no_seroreversion, serostatus_finite_seroreversion,
  ages, age_groups, diagnostic_adjustments = c(1, 1)
) {
  output_ov16_data <- rep(NA, length(age_groups) * 2)
  for (age_group_index in 1:length(age_groups)) {
    age_group <- age_groups[[age_group_index]]
    output_ov16_data[age_group_index] <- seroprevalence_for_age(
      serostatus_no_seroreversion, ages, age_group[1], age_group[2],
      diagnostic_adjustments[1], diagnostic_adjustments[2]
    )
    output_ov16_data[age_group_index + length(age_groups)] <- seroprevalence_for_age(
      serostatus_finite_seroreversion, ages, age_group[1], age_group[2],
      diagnostic_adjustments[1], diagnostic_adjustments[2]
    )
  }

  return(output_ov16_data)
}