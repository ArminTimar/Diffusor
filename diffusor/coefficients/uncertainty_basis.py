"""Where each sampled log D scatter (``sigma_logD``) comes from.

Every law with a ``sigma_logD`` must have an entry here, and each entry starts
with one of three words:

``published``  the source states this scatter (location given);
``derived``    computed by Diffusor from a statement in the source (how is given);
``assumed``    Diffusor's choice, with no basis in the source beyond what is said.

Audited 2 October 2026 against the rendered source pages.
The registry attaches these strings to the laws and the test suite fails if a
law with a ``sigma_logD`` has none. The basis is shown on the coefficient page
and in every exported methods paragraph that samples the coefficient.
"""

_SIEVWRIGHT = ("assumed: 0.2 log units for every element. Sievwright et al. (2020, p. 10) state that "
               "uncertainties on individual log D values are typically below 0.2 log units (1 sigma) and "
               "give no uncertainty for the Table 5 fit parameters. Their Supplementary Table S1 lists the 1 sigma "
               "of the 164 individual log D values that have one (median 0.12, 84 % at or below 0.2, largest "
               "0.51) and the fits were weighted by it; refitting eq. 5 to those values (Diffusor's reconstruction, not a "
               "published number) gives a 1 sigma of log D along the fit of 0.02-0.32 (median 0.08) depending on "
               "the element, so 0.2 is a representative, conservative width. For the T-scaled entries the "
               "borrowed activation energies add an uncertainty that is not included")
_AGGARWAL = ("assumed: 0.3 log units. Van Orman & Crispin (2010) Tables 10 and 11 give no uncertainties "
             "for these buffer-path Arrhenius parameters")
_VOC12 = ("assumed: 0.3 log units. Van Orman & Crispin (2010) Table 12 gives no uncertainties. For orientation, "
          "at 1200 and 1300 C the Fe, Co, Mn and Ti rows reproduce the 65 raw tracer values that Aggarwal & Dieckmann (2002b) "
          "measured at x_Ti = 0 and 0.2 with an rms of 0.13 log units, and the interpolation to their "
          "x_Ti = 0.1 values has an rms of 0.30; below 1100 C no data exist")

SIGMA_LOGD_BASIS = {
    "ol_FeMg_dohmen_chakraborty2007_tamed": (
        "assumed: 0.21 log units. Dohmen & Chakraborty (2007, abstract and p. 425) state only that the "
        "equations reproduce all 113 experimental data points with a maximum deviation of 0.5 log units, "
        "which is not a standard deviation. 0.21 treats that maximum as about 2.4 sigma (for 113 normally "
        "distributed points the expected maximum would be about 2.8 sigma, i.e. sigma = 0.18)"),
    "ol_FeMg_dohmen_chakraborty2007_ped": (
        "assumed: 0.21 log units, as for eq. 27. The paper states a maximum deviation of 0.5 log units "
        "(p. 425) and no standard deviation"),
    "ol_FeMg_chakraborty1997": (
        "assumed: 0.3 log units. The paper gives errors on D0 and Q only (the level of those errors is "
        "not stated)"),
    "opx_FeMg_dohmen2016": (
        "assumed: 0.2 log units, the overall experimental uncertainty of individual D values that "
        "Dohmen et al. (2016, p. 2214) state from reproducibility and time-series experiments (level not "
        "stated). They give errors on log D0, Q and n (level not stated) but no scatter of log D about "
        "the fit"),
    "opx_FeMg_dohmen2016_fs1": (
        "assumed: 0.2 log units, the same overall experimental uncertainty of individual D values that "
        "Dohmen et al. (2016, p. 2214) state for the study; no scatter about the Opx8 fit is given"),
    "opx_FeMg_dias2025": (
        "assumed: 0.34 log units, the figure Dias, Dohmen & Behrens (2025, p. 208) give for how well the "
        "equations reproduce the experimental data 'on average' (the statistic is not defined in the main "
        "text; it refers to the supplement). Taking it as the Gaussian 1 sigma of log D is Diffusor's choice; "
        "if it were a mean absolute deviation the sigma would be about 0.43"),
    "opx_FeMg_dias_dohmen2024": "assumed: 0.2 log units. The paper gives the error of Q only",
    "opx_Lu_dias2025": "assumed: 0.3 log units. Dias et al. (2025) give the error of Q only",
    "opx_Ce_dias2025": "assumed: 0.3 log units. Dias et al. (2025) give the error of Q only",
    "opx_Eu_dias2025": "assumed: 0.3 log units. Dias et al. (2025) give the error of Q only",
    "opx_FeMg_ganguly_tazzoli1994_nofo2": (
        "published: Ganguly & Tazzoli (1994, p. 934), the standard error in the estimated value of D is "
        "approximately one order of magnitude"),
    "opx_FeMg_ganguly_tazzoli1994": (
        "published: Ganguly & Tazzoli (1994, p. 934), the standard error in the estimated value of D is "
        "approximately one order of magnitude"),
    "cpx_FeMg_muller2013": (
        "derived: 0.5 log units, half of the 'within 1 log unit' agreement of the Arrhenius fit with "
        "the measurements stated by Mueller et al. (2013, p. 1569, envelope of Fig. 5b). Reading that "
        "envelope as about two standard deviations is Diffusor's conversion, not a statement of the paper "
        "(which also gives 0.4 log units as the reproducibility of individual D values)"),
    "cpx_FeMg_dimanov_sautter2000": (
        "assumed: 0.5 log units. The paper gives 1 sigma for log D0 (0.32) and Q (64 kJ/mol) without "
        "their covariance"),
    "cpx_CaMg_brady1983": (
        "derived: 0.3 log units = log10(2), from the 'uncertainty of a factor of 2' in the abstract of "
        "Brady & McCallister (1983)"),
    "plag_Mg_vanorman2014": (
        "assumed: 0.25 log units. Van Orman et al. (2014) give 2 sigma errors of the three fit parameters "
        "without their covariance"),
    "plag_Mg_costa2003": ("assumed: 0.3 log units. Costa et al. (2003) give no uncertainty for eq. 8, an "
                          "estimate rather than a fit; Van Orman et al. (2014) find it over-predicts D by up "
                          "to a factor of 10 at low temperature and low An"),
    "plag_Mg_audetat2026": (
        "assumed: 0.25 log units. Audetat et al. give errors of the four fit parameters without their "
        "covariance or confidence level"),
    "plag_Sr_grocolas2025": (
        "assumed: 0.25 log units. The Monte Carlo uses the covariance mode for this law by default; this "
        "scatter applies only if log D at T is chosen"),
    "plag_Ba_grocolas2025": (
        "assumed: 0.25 log units. The Monte Carlo uses the covariance mode for this law by default; this "
        "scatter applies only if log D at T is chosen"),
    "plag_Sr_giletti_casserly1994": (
        "assumed: 0.3 log units. Giletti & Casserly (1994, p. 3791) state that eq. 1 and the individual "
        "fits differ by less than a factor of two over the measured range"),
    "plag_Sr_cherniak_watson1994": (
        "assumed: 0.3 log units. Grocolas et al. (2025) eq. 13 gives parameter errors without covariance"),
    "plag_Ba_cherniak2002": (
        "assumed: 0.3 log units. Grocolas et al. (2025) eq. 14 gives parameter errors without covariance"),
    "plag_Li_pohl2024_interstitial": (
        "assumed: 0.3 log units. Pohl et al. (2024) give errors of log D0 and Q without their covariance"),
    "plag_Li_pohl2024_vacancy": (
        "assumed: 0.3 log units. Pohl et al. (2024) give errors of log D0 and Q without their covariance"),
    "plag_NaSiCaAl_grove1984": (
        "derived: 0.7 log units = log10(5), the upper end of the statement of Grove et al. (1984, p. 2116) "
        "that each D is correct to within a factor of 2 to 5, applied by Diffusor as 1 sigma (their own "
        "estimate for each D is about +/-0.5 in ln D, or 0.22 log10 units)"),
    "kfs_Sr_cherniak1996": (
        "published (secondary): 0.03 log units, the uncertainty Chamberlain et al. (2014, p. 14) state for D0 "
        "and E of this law; the derivation is in their Electronic Appendix 7, which Diffusor does not have, "
        "and no confidence level is given"),
    "kfs_Ba_cherniak2002": (
        "published (secondary): 0.12 log units, the uncertainty Chamberlain et al. (2014, p. 14) state for D0 "
        "and E of this law; the derivation is in their Electronic Appendix 7, which Diffusor does not have, "
        "and no confidence level is given"),
    "mt_FeTi_freer_hauptman1978": (
        "assumed: 0.5 log units. Freer & Hauptman (1978) give D0 as 3.85 (+1.68/-1.11) x 1e-3 cm2/s and "
        "the activation energy as 2.23 +/- 0.04 eV at 3 mol% Ti, and beta = 13.3 +/- 4.8"),
    "mt_FeTi_aragon1984": (
        "assumed: 0.5 log units. Van Orman & Crispin (2010) Table 11 gives no uncertainty"),
}
for _sp in ("Ti", "Fe", "Mn", "Co", "Cr", "Al"):
    SIGMA_LOGD_BASIS[f"mt_{_sp}_vanorman_crispin2010"] = _VOC12
for _sp in ("Mn", "Co", "Ni", "Mg", "Zn", "Sc", "Al", "Ga", "In", "Y", "Cr", "Lu", "V3+", "Ti", "V4+",
            "Zr", "Hf", "U", "Nb", "Ta", "Mo"):
    SIGMA_LOGD_BASIS[f"mt_{_sp}_sievwright2020"] = _SIEVWRIGHT
for _sp in ("Mn", "Co", "Al", "Cr", "Ti"):
    SIGMA_LOGD_BASIS[f"mt_{_sp}_sievwright2020_1150"] = _SIEVWRIGHT
for _k in ("Fe_aggarwal2002_WM", "Fe_aggarwal2002_MH", "Ti_aggarwal2002_WM", "Ti_aggarwal2002_MH",
           "Ti_aggarwal2002_WM_xti02", "Ti_aggarwal2002_MH_xti02", "Fe_aggarwal2002_WM_xti02",
           "Fe_aggarwal2002_MH_xti02"):
    SIGMA_LOGD_BASIS[f"mt_{_k}"] = _AGGARWAL


# Errors the sources quote for laws whose coefficient is held fixed (no covariance, no scatter).
QUOTED_ERRORS = {
    # The papers print these +/- values without saying what they are (the one-sigma statements in the
    # Cherniak papers refer to the individual D values of the data tables).
    'rt_Hf_cherniak2007_c': 'Cherniak et al. (2007) abstract: Q = 169 +/- 36 kJ/mol (confidence level not stated)',
    'rt_Hf_cherniak2007_a': ('Cherniak et al. (2007) abstract: Q = 227 +/- 62 kJ/mol (normal to c; '
                             'confidence level not stated)'),
    'ap_La_cherniak2000_in': ('Cherniak (2000) abstract: Q = 324 +/- 9 kJ/mol; Table 5: log D0 = -6.58 +/- 0.38 '
                              '(confidence level not stated)'),
    'ap_Nd_cherniak2000_in': ('Cherniak (2000) abstract: Q = 348 +/- 13 kJ/mol; Table 5: log D0 = -5.61 +/- 0.41 '
                              '(confidence level not stated)'),
    'ap_Dy_cherniak2000_in': ('Cherniak (2000) abstract: Q = 340 +/- 11 kJ/mol; Table 5: log D0 = -6.01 +/- 0.46 '
                              '(confidence level not stated)'),
    'ap_Yb_cherniak2000_in': ('Cherniak (2000) abstract: Q = 292 +/- 23 kJ/mol; Table 5: log D0 = -7.88 +/- 0.95 '
                              '(confidence level not stated)'),
    'ap_Sm_cherniak2000_implant': ('Cherniak (2000) abstract: Q = 298 +/- 17 kJ/mol; Table 5: log D0 = -6.20 +/- 0.74 '
                                   '(confidence level not stated)'),
    'ap_Nd_cherniak2000_out': ('Cherniak (2000) abstract: Q = 392 +/- 31 kJ/mol; Table 5: log D0 = -5.03 +/- 1.11 '
                               '(confidence level not stated)'),
    'zrn_Sm_cherniak1997': ('Cherniak et al. (1997) abstract: log D0 = 8.46 +/- 1.61, Q = 841 +/- 57 kJ/mol '
                            '(confidence level not stated)'),
    'zrn_Dy_cherniak1997': ('Cherniak et al. (1997) abstract: log D0 = 5.36 +/- 0.21, Q = 734 +/- 35 kJ/mol '
                            '(confidence level not stated)'),
    'zrn_Yb_cherniak1997': ('Cherniak et al. (1997) abstract: log D0 = 7.40 +/- 1.15, Q = 769 +/- 34 kJ/mol '
                            '(confidence level not stated)'),
    'zrn_Ti_bloch2022_c': ('Bloch et al. (2022) abstract: log D0 = 1.34 +/- 1.44, Q = 555,425 +/- 44,820 J/mol '
                           '(confidence level not stated)'),
    'zrn_Ti_cherniak2007_perp_c': ('Cherniak & Watson (2007) abstract: Q = 754 +/- 56 kJ/mol; results text: '
                                   'log D0 = 2.52 +/- 1.72 (confidence level not stated)'),
    'mnz_Pb_cherniak2004': ('Cherniak et al. (2004) abstract: Q = 592 +/- 39 kJ/mol (the results text and Fig. 3 '
                            'caption print 594 +/- 39; the log D0 error is printed as +/- 0.014 in the text and '
                            '+/- 1.350 in the Fig. 3 caption; confidence level not stated)'),
    'xtm_Sm_cherniak2006': ('Cherniak (2006) abstract: Q = 441 +/- 12 kJ/mol; results text: log D0 = -3.824 +/- 0.420 '
                            '(confidence level not stated)'),
    'xtm_Dy_cherniak2006': ('Cherniak (2006) abstract: Q = 349 +/- 16 kJ/mol; results text: log D0 = -7.046 +/- 0.585 '
                            '(confidence level not stated)'),
    'xtm_Yb_cherniak2006': ('Cherniak (2006) abstract: Q = 362 +/- 13 kJ/mol; results text: log D0 = -6.407 +/- 0.460 '
                            '(confidence level not stated)'),
    'xtm_Pb_cherniak2006': ('Cherniak (2006) abstract: Q = 382 +/- 64 kJ/mol; results text: log D0 = -8.529 +/- 2.064 '
                            '(confidence level not stated)'),
    'ol_Ni_petry2004': ('Petry et al. (2004, Fig. 6, p. 4184) print E_D = 216 +/- 6 kJ/mol and '
                        'D0 = 3.84 (+2.5/-1.5) x 1e-9 m2/s for the Fo90 fixed-fO2 fit used here (confidence level '
                        'not stated); the text puts the overall uncertainty of D at about a factor of two'),
}
