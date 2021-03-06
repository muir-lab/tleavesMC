% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Eq_27.R
\name{so17_eqn27}
\alias{so17_eqn27}
\title{{so17_eqn27}: Approximate leaf temperature using equation 27 of Schymanski and Or 2017}
\usage{
so17_eqn27(
  leaf_par,
  enviro_par,
  constants,
  progress = TRUE,
  quiet = FALSE,
  set_units = TRUE
)
}
\arguments{
\item{leaf_par}{A list of leaf parameters. This can be generated using the \code{make_leafpar} function.}

\item{enviro_par}{A list of environmental parameters. This can be generated using the \code{make_enviropar} function.}

\item{constants}{A list of physical constants. This can be generated using the \code{make_constants} function.}

\item{progress}{Logical. Should a progress bar be displayed?}

\item{quiet}{Logical. Should messages be displayed?}

\item{set_units}{Logical. Should \code{units} be set? The function is faster when FALSE, but input must be in correct units or else results will be incorrect without any warning.}
}
\value{
Value in degrees Kelvin (K) of class \code{units}
}
\description{
{so17_eqn27}: Approximate leaf temperature using equation 27 of Schymanski and Or 2017
}
\details{
Schymanski and Or function can be used to obtain analytical expressions for  H,  L,  and  S_r that satisfy the energy balance (S_r).
Alternatively, the value of  T_leaf  obtained from the Schymanski and Or function for specific conditions could be used to calculate any of the energy balance components using the fundamental equations.
In this case, bias in  T_leaf  due to simplifying assumptions included in the derivation of Schymanski and Or function could lead to a mismatch in the leaf energy balance calculation.

Below is the same equation as Schymanski and Or 2017 but with variable nomenclature consistent with tealeaves:

\deqn{T_\mathrm{leaf} = (S_\mathrm{r} + H(T_\mathrm{air}) + ?(? + ? - p_\mathrm{sat}) + ?(alpha_\mathrm{l})sigma(3T_\mathrm{air})^4 + (?_\mathrm)^4)) / (H + L? + 4?(alpha_\mathrm{l})(sigmaT_\mathrm{air})^3)}{T_leaf = (S_r + H (T_air) + ? (? + ? - p_sat) + ?(alpha_l)sigma(3T_air^4 = ?^4)) / (H + L? + 4?alpha_l sigma T_air^3)}

Below is the original Schymanski and Or 2017 equation:

\deqn{T_\mathrm{l} = (R_\mathrm{s} + C_\mathrm{H}(T_\mathrm{a}) + C_\mathrm{E}(delta_\mathrm{eTa} + P_\mathrm{wa} - P_\mathrm{was}) + a_\mathrm{sh}(epsilon_\mathrm{l})sigma(3T_\mathrm{a}^4 + T_\mathrm{w}^4)) / (C_\mathrm{H} + C_\mathrm{E}(delta_\mathrm{eTa}) + 4a_\mathrm{sh}(epsilon_\mathrm{l})sigma(T_\mathrm{a}^3))}{T_l = (R_s + C_H(T_a) + C_E(delta_eTa + P_wa - P_was) + a_sh(epsilon_l)sigma(3T_a^4 = T_w^4)) / (C_H + C_E(delta_eTa) + 4a_sh(epsilon_l)sigma(T_a^3))}

\tabular{lllll}{
\emph{Symbol} \tab \emph{S & O 2017} \tab \emph{R} \tab \emph{Description} \tab \emph{Units} \cr
\eqn{T_leaf} \tab \eqn{T_l} \tab \code{T_leaf} \tab leaf temperature \tab K \cr
\eqn{S_r} \tab \eqn{R_s} \tab \code{S_r} \tab absorbed short-wave radiation \tab W m^-2 \cr
\eqn{?} \tab \eqn{C_H} \tab \code{?} \tab transfer coefficients for latent sensible heat \tab unitless \cr
\eqn{T_air} \tab \eqn{T_a} \tab \code{T_air} \tab leaf temperature \tab K \cr
\eqn{?} \tab \eqn{C_E} \tab \code{?} \tab transfer coefficients for latent heat \tab unitless \cr
\eqn{?} \tab \eqn{delta_eTa} \tab \code{?} \tab slope of saturation vapour pressure at air temperature \tab unitless \cr
\eqn{?} \tab \eqn{P_wa} \tab \code{?} \tab vapour pressure in the atmosphere \tab kPa \cr
\eqn{p_sat} \tab \eqn{P_was} \tab \code{p_sat} \tab saturating water vapour pressure \tab kPa \cr
\eqn{?} \tab \eqn{a_sh} \tab \code{?} \tab fraction of projected area exchanging sensible heat with the air \tab 1 \cr
\eqn{alpha_l} \tab \eqn{epsilon_l} \tab \code{abs_l} \tab absorbtivity of long-wave radiation \tab unitless \cr
\eqn{signma} \tab \eqn{s} \tab \code{signma} \tab Stefan–Boltzmann constant \tab W K^-1 m^-2 \cr
\eqn{T_air} \tab \eqn{T_air} \tab \code{T_a} \tab Stefan–Boltzmann constant \tab K \cr
\eqn{?} \tab \eqn{T_w} \tab \code{?} \tab surroundings temperature \tab K \cr

}
}
\examples{
library(tealeaves)

leaf_par <- make_leafpar()
enviro_par <- make_enviropar()
constants <- make_constants()
tleaf(leaf_par, enviro_par, constants)

tealeaves:::so17_eqn27(leaf_par, enviro_par, constants)
}
