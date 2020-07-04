#' {so17_eqn27}: Approximate leaf temperature using equation 27 of Schymanski and Or 2017
#'
#' @inheritParams tleaves
#'
#' @export

so17_eqn27 <- function(leaf_par, enviro_par, constants, progress = TRUE,
                       quiet = FALSE, set_units = TRUE) {

  #this makes sure that the inputs are in correct units

  if (set_units) {
    leaf_par %<>% tealeaves::leaf_par()
    enviro_par %<>% tealeaves::enviro_par()
    constants %<>% tealeaves::constants()
  } else {
    if (!quiet) warning("tleaves: units have not been checked prior to solving")
  } #this checks everything is in correct units

  pars <- c(leaf_par, enviro_par)
  par_units <- purrr::map(pars, units) %>% #map units to parameters
    magrittr::set_names(names(pars))

  #this drops units to increase processing speed

  pars %<>%
    purrr::map_if(~ inherits(.x, "units"), drop_units) %>% #drop units to improve processing speed
    make_parameter_sets()

  constants %<>%
    purrr::map_if(~ inherits(.x, "units"), drop_units) #drop units to improve processing speed

  c_H <- .get_ch(pars, constants)
  c_E <- .get_ce(pars, constants)

  ###Below is SO_27 equation

  # delta_Ta: Slope of saturation vapour pressure at air temperature (Pa/K). Eqn B27 (approx)
  # See notes on M_w and lambda_E in .get_ce function
  M_w <- 0.018
  lambda_E <- 2.45e6
  delta_Ta <- (611 * lambda_E * M_w * exp(lambda_E * M_w / constants$R * (1 / 273 - 1 / pars$T_air))) / (constants$R * pars$T_air ^ 2)

  # p_air = P_wa: vapour pressure in the atmosphere (Pa)
  # p_sat = P_was: saturation vapour pressure in the atmosphere (Pa)
  # .get_ps return pressure in kPa, hence conversion by 1000x to Pa
  p_sat <- 1000 * tealeaves:::.get_ps(pars$T_air, pars$P, TRUE)
  p_air <- p_sat * pars$RH

  # Original equation 27 (pg 690)
  # (R_s + c_H * T_air + c_E * (delta_Ta * T_air + p_air - p_sat) + a_sh * e_l * s * (3 * T_air ^ 4 + T_w ^ 4) / (c_H + c_E * delta_Ta + 4 * a_sh * abs_l * s * T_air ^ 3)
  #
  # I assume:
  # a_sh = 2
  # T_air = T_w
  R_s <- 0 # calculated from input, this ultimately should be calculated from the input


  T_leaf <- (R_s + c_H * pars$T_air + c_E * (delta_Ta * pars$T_air + p_air - p_sat) + 2 * pars$abs_l * constants$s * (4 * pars$T_air ^ 4)) / (c_H + c_E * delta_Ta + 4 * 2 * pars$abs_l * constants$s * pars$T_air ^ 3)

  T_leaf # in K

  #' @return Value in degrees Kelvin (K) of class \code{units}

  #' @inheritParams .get_H #inherit the parameters as mucn as possible
  #'
  #' @details
  #'
  #' Schymanski and Or function can be used to obtain analytical expressions for  H,  L,  and  S_r that satisfy the energy balance (S_r).
  #' Alternatively, the value of  T_leaf  obtained from the Schymanski and Or function for specific conditions could be used to calculate any of the energy balance components using the fundamental equations.
  #' In this case, bias in  T_leaf  due to simplifying assumptions included in the derivation of Schymanski and Or function could lead to a mismatch in the leaf energy balance calculation.
  #'
  #' Below is the same equation as Schymanski and Or 2017 but with variable nomenclature consistent with tealeaves:
  #'
  #' \deqn{T_\mathrm{leaf} = (S_\mathrm{r} + H(T_\mathrm{air}) + ?(? + ? - p_\mathrm{sat}) + ?(alpha_\mathrm{l})sigma(3T_\mathrm{air})^4 + (?_\mathrm)^4)) / (H + L? + 4?(alpha_\mathrm{l})(sigmaT_\mathrm{air})^3)}{T_leaf = (S_r + H (T_air) + ? (? + ? - p_sat) + ?(alpha_l)sigma(3T_air^4 = ?^4)) / (H + L? + 4?alpha_l sigma T_air^3)}
  #'
  #' Below is the original Schymanski and Or 2017 equation:
  #'
  #' \deqn{T_\mathrm{l} = (R_\mathrm{s} + C_\mathrm{H}(T_\mathrm{a}) + C_\mathrm{E}(delta_\mathrm{eTa} + P_\mathrm{wa} - P_\mathrm{was}) + a_\mathrm{sh}(epsilon_\mathrm{l})sigma(3T_\mathrm{a}^4 + T_\mathrm{w}^4)) / (C_\mathrm{H} + C_\mathrm{E}(delta_\mathrm{eTa}) + 4a_\mathrm{sh}(epsilon_\mathrm{l})sigma(T_\mathrm{a}^3))}{T_l = (R_s + C_H(T_a) + C_E(delta_eTa + P_wa - P_was) + a_sh(epsilon_l)sigma(3T_a^4 = T_w^4)) / (C_H + C_E(delta_eTa) + 4a_sh(epsilon_l)sigma(T_a^3))}
  #'
  #' \tabular{lllll}{
  #' \emph{Symbol} \tab \emph{S & O 2017} \tab \emph{R} \tab \emph{Description} \tab \emph{Units} \cr
  #' \eqn{T_leaf} \tab \eqn{T_l} \tab \code{T_leaf} \tab leaf temperature \tab K \cr
  #' \eqn{S_r} \tab \eqn{R_s} \tab \code{S_r} \tab absorbed short-wave radiation \tab W m^-2 \cr
  #' \eqn{?} \tab \eqn{C_H} \tab \code{?} \tab transfer coefficients for latent sensible heat \tab unitless \cr
  #' \eqn{T_air} \tab \eqn{T_a} \tab \code{T_air} \tab leaf temperature \tab K \cr
  #' \eqn{?} \tab \eqn{C_E} \tab \code{?} \tab transfer coefficients for latent heat \tab unitless \cr
  #' \eqn{?} \tab \eqn{delta_eTa} \tab \code{?} \tab slope of saturation vapour pressure at air temperature \tab unitless \cr
  #' \eqn{?} \tab \eqn{P_wa} \tab \code{?} \tab vapour pressure in the atmosphere \tab kPa \cr
  #' \eqn{p_sat} \tab \eqn{P_was} \tab \code{p_sat} \tab saturating water vapour pressure \tab kPa \cr
  #' \eqn{?} \tab \eqn{a_sh} \tab \code{?} \tab fraction of projected area exchanging sensible heat with the air \tab 1 \cr
  #' \eqn{alpha_l} \tab \eqn{epsilon_l} \tab \code{abs_l} \tab absorbtivity of long-wave radiation \tab unitless \cr
  #' \eqn{signma} \tab \eqn{s} \tab \code{signma} \tab Stefan–Boltzmann constant \tab W K^-1 m^-2 \cr
  #' \eqn{T_air} \tab \eqn{T_air} \tab \code{T_a} \tab Stefan–Boltzmann constant \tab K \cr
  #' \eqn{?} \tab \eqn{T_w} \tab \code{?} \tab surroundings temperature \tab K \cr
  #'
  #' }
  #'
  #'
  #'
  #'
  #'
  #' @examples
  #' library(tealeaves)
  #'
  #' leaf_par <- make_leafpar()
  #' enviro_par <- make_enviropar()
  #' constants <- make_constants()
  #' tleaf(leaf_par, enviro_par, constants)
  #'
  #' tealeaves:::so17_eqn27(leaf_par, enviro_par, constants)

}

