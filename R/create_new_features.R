is_missing <- function(x) x <= -999

#' Create New Features
#'
#' To improve classification performance, we can calculate these extra features using the provided primitive features.
#' Each of these new features has a physical meaning that will hopefully be useful to the classification algorithms.
#'
#' @param df A dataframe containing data from the HiggsML challenge. Must contain all primitive (PRI_) features as columns.
#'
#' @return The same dataframe given as an argument but with columns added for each of the extra features.
#' @export
#'
#' @examples
#' \dontrun{
#'   create_new_features(training_dataset)
#' }
create_new_features <- function(df) {
  df[is_missing(df)] <- NA # useful to use NA in calculations below

  # abs
  df$SEB_abs_eta_tau <- abs(df$PRI_tau_eta)
  df$SEB_abs_eta_lep <- abs(df$PRI_lep_eta)
  df$SEB_abs_eta_jet1 <- abs(df$PRI_jet_leading_eta)
  df$SEB_abs_eta_jet2 <- abs(df$PRI_jet_subleading_eta)
  df$SEB_deltaeta_tau_lep <- abs(df$PRI_tau_eta-df$PRI_lep_eta)
  df$SEB_deltaeta_tau_jet1 <- abs(df$PRI_tau_eta-df$PRI_jet_leading_eta)
  df$SEB_deltaeta_tau_jet2 <- abs(df$PRI_tau_eta-df$PRI_jet_subleading_eta)
  df$SEB_deltaeta_lep_jet1 <- abs(df$PRI_lep_eta-df$PRI_jet_leading_eta)
  df$SEB_deltaeta_lep_jet2 <- abs(df$PRI_lep_eta-df$PRI_jet_subleading_eta)
  df$SEB_deltaeta_jet_jet <- abs(df$PRI_jet_leading_eta-df$PRI_jet_subleading_eta)

  # prodeta
  df$SEB_prodeta_tau_lep <- df$PRI_tau_eta*df$PRI_lep_eta
  df$SEB_prodeta_tau_jet1 <- df$PRI_tau_eta*df$PRI_jet_leading_eta
  df$SEB_prodeta_tau_jet2 <- df$PRI_tau_eta*df$PRI_jet_subleading_eta
  df$SEB_prodeta_lep_jet1 <- df$PRI_lep_eta*df$PRI_jet_leading_eta
  df$SEB_prodeta_lep_jet2 <- df$PRI_lep_eta*df$PRI_jet_subleading_eta
  df$SEB_prodeta_jet_jet <- df$PRI_jet_leading_eta*df$PRI_jet_subleading_eta

  #aux
  SEB_deltaphi_tau_lep <- abs(df$PRI_tau_phi-df$PRI_lep_phi)
  SEB_deltaphi_tau_lep <- ifelse(SEB_deltaphi_tau_lep>pi,2*pi-SEB_deltaphi_tau_lep, SEB_deltaphi_tau_lep)
  SEB_deltaphi_tau_jet1 <- abs(df$PRI_tau_phi-df$PRI_jet_leading_phi)
  SEB_deltaphi_tau_jet1 <- ifelse(SEB_deltaphi_tau_jet1>pi,2*pi-SEB_deltaphi_tau_jet1,SEB_deltaphi_tau_jet1)
  SEB_deltaphi_tau_jet2 <- abs(df$PRI_tau_phi-df$PRI_jet_subleading_phi)
  SEB_deltaphi_tau_jet2 <- ifelse(SEB_deltaphi_tau_jet2>pi,2*pi-SEB_deltaphi_tau_jet2,SEB_deltaphi_tau_jet2)
  SEB_deltaphi_lep_jet1 <- abs(df$PRI_lep_phi-df$PRI_jet_leading_phi)
  SEB_deltaphi_lep_jet1 <- ifelse(SEB_deltaphi_lep_jet1>pi,2*pi-SEB_deltaphi_lep_jet1,SEB_deltaphi_lep_jet1)
  SEB_deltaphi_lep_jet2 <- abs(df$PRI_lep_phi-df$PRI_jet_subleading_phi)
  SEB_deltaphi_lep_jet2 <- ifelse(SEB_deltaphi_lep_jet2>pi,2*pi-SEB_deltaphi_lep_jet2,SEB_deltaphi_lep_jet2)
  SEB_deltaphi_jet_jet <- abs(df$PRI_jet_leading_phi-df$PRI_jet_subleading_phi)
  SEB_deltaphi_jet_jet <- ifelse(SEB_deltaphi_jet_jet>pi, 2*pi-SEB_deltaphi_jet_jet, SEB_deltaphi_jet_jet)

  # deltar
  df$SEB_deltar_tau_lep <- sqrt(df$SEB_deltaeta_tau_lep^2+SEB_deltaphi_tau_lep^2)
  df$SEB_deltar_tau_jet1 <- sqrt(df$SEB_deltaeta_tau_jet1^2+SEB_deltaphi_tau_jet1^2)
  df$SEB_deltar_tau_jet2 <- sqrt(df$SEB_deltaeta_tau_jet2^2+SEB_deltaphi_tau_jet2^2)
  df$SEB_deltar_lep_jet1 <- sqrt(df$SEB_deltaeta_lep_jet1^2+SEB_deltaphi_lep_jet1^2)
  df$SEB_deltar_lep_jet2 <- sqrt(df$SEB_deltaeta_lep_jet2^2+SEB_deltaphi_lep_jet2^2)
  df$SEB_deltar_jet_jet <- sqrt(df$SEB_deltaeta_jet_jet^2+SEB_deltaphi_jet_jet^2)

  #aux
  d <- df$PRI_tau_phi - df$PRI_lep_phi
  d <- 1.0 - 2.0*((d>pi)|((d<0) & (d>-pi)))
  a <- sin(df$PRI_met_phi-df$PRI_lep_phi)
  b <- sin(df$PRI_tau_phi-df$PRI_met_phi)
  df$TIM_met_phi_centrality <- d*(a+b)/sqrt(a^2+b^2)

  # centrality
  df$SEB_lep_eta_centrality <- exp(-4.0*(df$PRI_lep_eta-(df$PRI_jet_leading_eta+df$PRI_jet_subleading_eta)/2)^2/(df$PRI_jet_leading_eta-df$PRI_jet_subleading_eta)^2)
  df$SEB_tau_eta_centrality <- exp(-4.0*(df$PRI_tau_eta-(df$PRI_jet_leading_eta+df$PRI_jet_subleading_eta)/2)^2/(df$PRI_jet_leading_eta-df$PRI_jet_subleading_eta)^2)

  # pt2
  df$SEB_pt2_met_tau <- ((df$PRI_met*cos(df$PRI_met_phi) + df$PRI_tau_pt*cos(df$PRI_tau_phi))^2
                         +(df$PRI_met*sin(df$PRI_met_phi) + df$PRI_tau_pt*sin(df$PRI_tau_phi))^2)
  df$SEB_pt2_met_lep <- ((df$PRI_met*cos(df$PRI_met_phi) + df$PRI_lep_pt*cos(df$PRI_lep_phi))^2
                         +(df$PRI_met*sin(df$PRI_met_phi) + df$PRI_lep_pt*sin(df$PRI_lep_phi))^2)
  df$SEB_pt2_met_jet1 <- ((df$PRI_met*cos(df$PRI_met_phi) + df$PRI_jet_leading_pt*cos(df$PRI_jet_leading_phi))^2
                          +(df$PRI_met*sin(df$PRI_met_phi) + df$PRI_jet_leading_pt*sin(df$PRI_jet_leading_phi))^2)
  df$SEB_pt2_met_jet2 <- ((df$PRI_met*cos(df$PRI_met_phi) + df$PRI_jet_subleading_pt*cos(df$PRI_jet_subleading_phi))^2
                          +(df$PRI_met*sin(df$PRI_met_phi) + df$PRI_jet_subleading_pt*sin(df$PRI_jet_subleading_phi))^2)
  df$SEB_pt2_tau_lep <- ((df$PRI_tau_pt*cos(df$PRI_tau_phi) + df$PRI_lep_pt*cos(df$PRI_lep_phi))^2
                         +(df$PRI_tau_pt*sin(df$PRI_tau_phi) + df$PRI_lep_pt*sin(df$PRI_lep_phi))^2)
  df$SEB_pt2_tau_jet1 <- ((df$PRI_tau_pt*cos(df$PRI_tau_phi) + df$PRI_jet_leading_pt*cos(df$PRI_jet_leading_phi))^2
                          +(df$PRI_tau_pt*sin(df$PRI_tau_phi) + df$PRI_jet_leading_pt*sin(df$PRI_jet_leading_phi))^2)
  df$SEB_pt2_tau_jet2 <- ((df$PRI_tau_pt*cos(df$PRI_tau_phi) + df$PRI_jet_subleading_pt*cos(df$PRI_jet_subleading_phi))^2
                          +(df$PRI_tau_pt*sin(df$PRI_tau_phi) + df$PRI_jet_subleading_pt*sin(df$PRI_jet_subleading_phi))^2)
  df$SEB_pt2_lep_jet1 <- ((df$PRI_lep_pt*cos(df$PRI_lep_phi) + df$PRI_jet_leading_pt*cos(df$PRI_jet_leading_phi))^2
                          +(df$PRI_lep_pt*sin(df$PRI_lep_phi) + df$PRI_jet_leading_pt*sin(df$PRI_jet_leading_phi))^2)
  df$SEB_pt2_lep_jet2 <- ((df$PRI_lep_pt*cos(df$PRI_lep_phi) + df$PRI_jet_subleading_pt*cos(df$PRI_jet_subleading_phi))^2
                          +(df$PRI_lep_pt*sin(df$PRI_lep_phi) + df$PRI_jet_subleading_pt*sin(df$PRI_jet_subleading_phi))^2)
  df$SEB_pt2_jet_jet <- ((df$PRI_jet_leading_pt*cos(df$PRI_jet_leading_phi) + df$PRI_jet_subleading_pt*cos(df$PRI_jet_subleading_phi))^2
                         +(df$PRI_jet_leading_pt*sin(df$PRI_jet_leading_phi) + df$PRI_jet_subleading_pt*sin(df$PRI_jet_subleading_phi))^2)

  # trans_mass
  # WARNING some square roots below cause NAN errors
  # TODO maybe want to handle the sqrt errors more carefully / add indicator
  df$SEB_trans_mass_met_tau <- sqrt((df$PRI_met+df$PRI_tau_pt)^2-df$SEB_pt2_met_tau)
  df$SEB_trans_mass_met_lep <- sqrt((df$PRI_met+df$PRI_lep_pt)^2-df$SEB_pt2_met_lep)
  df$SEB_trans_mass_met_jet1 <- sqrt((df$PRI_met+df$PRI_jet_leading_pt)^2-df$SEB_pt2_met_jet1)
  df$SEB_trans_mass_met_jet2 <- sqrt((df$PRI_met+df$PRI_jet_subleading_pt)^2-df$SEB_pt2_met_jet2)
  df$SEB_trans_mass_tau_lep <- sqrt((df$PRI_tau_pt+df$PRI_lep_pt)^2-df$SEB_pt2_tau_lep)
  df$SEB_trans_mass_tau_jet1 <- sqrt((df$PRI_tau_pt+df$PRI_jet_leading_pt)^2-df$SEB_pt2_tau_jet1)
  df$SEB_trans_mass_tau_jet2 <- sqrt((df$PRI_tau_pt+df$PRI_jet_subleading_pt)^2-df$SEB_pt2_tau_jet2)
  df$SEB_trans_mass_lep_jet1 <- sqrt((df$PRI_lep_pt+df$PRI_jet_leading_pt)^2-df$SEB_pt2_lep_jet1)
  df$SEB_trans_mass_lep_jet2 <- sqrt((df$PRI_lep_pt+df$PRI_jet_subleading_pt)^2-df$SEB_pt2_lep_jet2)
  df$SEB_trans_mass_jet_jet <- sqrt((df$PRI_jet_leading_pt+df$PRI_jet_subleading_pt)^2-df$SEB_pt2_jet_jet)

  # p2
  df$SEB_p2_tau_lep <- df$SEB_pt2_tau_lep +(df$PRI_tau_pt*sinh(df$PRI_tau_eta) + df$PRI_lep_pt*sinh(df$PRI_lep_eta))^2
  df$SEB_p2_tau_jet1 <- df$SEB_pt2_tau_jet1 +(df$PRI_tau_pt*sinh(df$PRI_tau_eta) + df$PRI_jet_leading_pt*sinh(df$PRI_jet_leading_eta))^2
  df$SEB_p2_tau_jet2 <- df$SEB_pt2_tau_jet2 +(df$PRI_tau_pt*sinh(df$PRI_tau_eta) + df$PRI_jet_subleading_pt*sinh(df$PRI_jet_subleading_eta))^2
  df$SEB_p2_lep_jet1 <- df$SEB_pt2_lep_jet1 +(df$PRI_lep_pt*sinh(df$PRI_lep_eta) + df$PRI_jet_leading_pt*sinh(df$PRI_jet_leading_eta))^2
  df$SEB_p2_lep_jet2 <- df$SEB_pt2_lep_jet2 +(df$PRI_lep_pt*sinh(df$PRI_lep_eta) + df$PRI_jet_subleading_pt*sinh(df$PRI_jet_subleading_eta))^2
  df$SEB_p2_jet_jet <- df$SEB_pt2_jet_jet +(df$PRI_jet_leading_pt*sinh(df$PRI_jet_leading_eta) + df$PRI_jet_subleading_pt*sinh(df$PRI_jet_subleading_eta))^2

  # E
  df$E_tau <- df$PRI_tau_pt*cosh(df$PRI_tau_eta)
  df$E_lep <- df$PRI_lep_pt*cosh(df$PRI_lep_eta)
  df$E_jet1 <- df$PRI_jet_leading_pt*cosh(df$PRI_jet_leading_eta)
  df$E_jet2 <- df$PRI_jet_subleading_pt*cosh(df$PRI_jet_subleading_eta)

  # mass
  df$SEB_mass_tau_lep <- sqrt((df$E_tau+df$E_lep)^2-df$SEB_p2_tau_lep)
  df$SEB_mass_tau_jet1 <- sqrt((df$E_tau+df$E_jet1)^2-df$SEB_p2_tau_jet1)
  df$SEB_mass_tau_jet2 <- sqrt((df$E_tau+df$E_jet2)^2-df$SEB_p2_tau_jet2)
  df$SEB_mass_lep_jet1 <- sqrt((df$E_lep+df$E_jet1)^2-df$SEB_p2_lep_jet1)
  df$SEB_mass_lep_jet2 <- sqrt((df$E_lep+df$E_jet2)^2-df$SEB_p2_lep_jet2)
  df$SEB_mass_jet_jet <- sqrt((df$E_jet1+df$E_jet2)^2-df$SEB_p2_jet_jet)

  # aux
  sum_px <- df$PRI_met*cos(df$PRI_met_phi) + df$PRI_tau_pt*cos(df$PRI_tau_phi) + df$PRI_lep_pt*cos(df$PRI_lep_phi)
  sum_py <- df$PRI_met*sin(df$PRI_met_phi) + df$PRI_tau_pt*sin(df$PRI_tau_phi) + df$PRI_lep_pt*sin(df$PRI_lep_phi)
  df$SEB_pt_met_tau_lep <- sqrt((sum_px)^2 +(sum_py)^2)

  sum_px_2 <- sum_px + fillna(df$PRI_jet_leading_pt*cos(df$PRI_jet_leading_phi),0.)
  sum_py_2 <- sum_py + fillna(df$PRI_jet_leading_pt*sin(df$PRI_jet_leading_phi),0.)
  df$SEB_pt_met_tau_lep_jet1 <- sqrt((sum_px_2)^2 +(sum_py_2)^2)

  sum_px_3 <- sum_px_2 + fillna(df$PRI_jet_subleading_pt*cos(df$PRI_jet_subleading_phi),0.)
  sum_py_3 <- sum_py_2 + fillna(df$PRI_jet_subleading_pt*sin(df$PRI_jet_subleading_phi),0.)
  df$SEB_pt_met_tau_lep_jet1_jet2 <- sqrt((sum_px_3)^2 +(sum_py_3)^2)

  df$SEB_sum_pt_met_tau_lep <- df$PRI_met + df$PRI_tau_pt + df$PRI_lep_pt
  df$SEB_sum_pt_met_tau_lep_jet1 <- df$SEB_sum_pt_met_tau_lep + fillna(df$PRI_jet_leading_pt,0.)
  df$SEB_sum_pt_met_tau_lep_jet1_jet2 <- df$SEB_sum_pt_met_tau_lep_jet1 + fillna(df$PRI_jet_subleading_pt,0.)
  df$SEB_sum_pt_met_tau_lep_jet_all <- df$SEB_sum_pt_met_tau_lep_jet1 + df$PRI_jet_all_pt

  df$SEB_sum_pt <- df$PRI_tau_pt + df$PRI_lep_pt + df$PRI_jet_all_pt

  df$SEB_pt_ratio_lep_tau <- df$PRI_lep_pt/df$PRI_tau_pt


  df[is.na(df)] <- -999.0 # move back to using NA
  # we shouldn't have to add any more missingness
  # indicators because hopefullly missingness in
  # any of the new features is indicated by
  # existing missingness indicators
  df
}
