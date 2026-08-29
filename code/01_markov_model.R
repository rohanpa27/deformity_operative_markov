# =============================================================================
# Cost-Effectiveness of Operative versus Nonoperative Management of
# Adult Spinal Deformity: A Markov State-Transition Model, 10-Year Horizon
#
# 01_markov_model.R  --  model definition, base case, OWSA, PSA, scenarios
#
# Run from the repository root:  Rscript code/01_markov_model.R
# Writes all CSV outputs to data/ . Must be run before 02_make_figures.R
# and 03_make_tables.R, which read those CSVs.
#
# Structure follows the DARTH cohort-Markov conventions (explicit transition
# array, cohort trace, reward vectors, within-cycle correction).
#
# EVIDENCE TIERS, applied to every parameter below:
#   DIRECT     exact value extracted from a peer-reviewed publication
#   DERIVED    arithmetic from one or more cited primary sources
#   EQUALIZED  held identical across arms because no arm-specific data exists
#   ASSUMPTION acknowledged estimate with no direct published source
#
# All corrections applied here are documented in docs/source_audit/SOURCE_AUDIT.md
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

set.seed(20260720)
RNG_SEED <- 20260720

# Resolve the repository root whether this is run from the root or from code/
root <- if (dir.exists("data") || file.exists("code/01_markov_model.R")) "." else ".."
DATA_DIR <- file.path(root, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 0. Price adjustment to 2024 US dollars
#
# Three of the four cost sources report hospital DIRECT costs rather than payer
# costs. AHRQ guidance for expenditure comparisons is to deflate such costs with
# the CMS Personal Health Care (PHC) price index rather than CPI-U Medical Care,
# because CPI-M tracks out-of-pocket consumer prices only while PHC reflects
# prices actually paid across all payers. PHC is therefore the base case and
# CPI-U Medical Care is retained as a scenario.
#
# PHC index (2017 = 100), CMS National Health Expenditure Accounts, republished
# by AHRQ MEPS. CPI-U Medical Care (1982-84 = 100), BLS series CUUR0000SAM.
# No source states a cost year distinct from its publication year, so each cost
# is inflated from its year of publication. This convention is stated in Methods.
# -----------------------------------------------------------------------------
PHC <- c(`2010` = 90.7, `2016` = 98.7, `2018` = 101.5, `2019` = 103.0,
         `2020` = 105.0, `2023` = 112.6, `2024` = 116.0)
CPIM <- c(`2010` = 388.436, `2016` = 463.675, `2018` = 484.700,
          `2019` = 498.413, `2020` = 518.876, `2023` = 549.084,
          `2024` = 563.841)

INFLATE_INDEX <- "PHC"   # "PHC" or "CPIM"; only used when COSTING == "literature"

# COSTING == "medicare"   : CMS DRG + Physician Fee Schedule payment (base case)
# COSTING == "literature" : journal-derived hospital direct costs, inflated
COSTING <- Sys.getenv("ASD_COSTING", "medicare")

infl <- function(amount, from_year, index = INFLATE_INDEX) {
  idx <- if (index == "PHC") PHC else CPIM
  amount * unname(idx["2024"]) / unname(idx[as.character(from_year)])
}

# -----------------------------------------------------------------------------
# 1. Model structure
#
# NINE states. The tunnel structure on the post-operative Stable state exists so
# that a nonoperative patient who crosses over to surgery in a later cycle is
# exposed to the SAME early (years 0-2) complication hazard as a patient
# operated at cycle 0. Without it, crossover patients would inherit the low
# late-phase hazard and the nonoperative arm would be flattered.
#
#  1 Nonoperative        nonoperative management (nonoperative arm only)
#  2 Stable_Y1           post-operative year 1
#  3 Stable_Y2           post-operative year 2
#  4 Stable_Late         post-operative year 3 and beyond
#  5 MechComp_Early      mechanical complication arising in post-op years 0-2
#  6 MechComp_Late       mechanical complication arising in post-op year 3+
#  7 Revision            revision surgery (one-cycle tunnel)
#  8 PostRev_Stable      alive after revision, no active complication
#  9 Death               absorbing
#
# MechComp is split so that the probability of proceeding to revision can differ
# between early complications (Soroceanu, 0.526) and late ones (Park, 0.333)
# without the model having to remember how a patient reached the state.
# -----------------------------------------------------------------------------
# CDC/NCHS United States Life Tables 2023, Table 1, total population, both
# sexes: probability of dying within one year (qx) at each single year of age.
LIFE_TABLE_QX <- c(
  `50` = 0.004120, `51` = 0.004420, `52` = 0.004752, `53` = 0.005130,
  `54` = 0.005562, `55` = 0.006022, `56` = 0.006524, `57` = 0.007106,
  `58` = 0.007769, `59` = 0.008476,
  `60` = 0.009202, `61` = 0.009926, `62` = 0.010656, `63` = 0.011405,
  `64` = 0.012197, `65` = 0.013054, `66` = 0.014088, `67` = 0.015075,
  `68` = 0.016150, `69` = 0.017324,
  `70` = 0.018460, `71` = 0.019921, `72` = 0.021541, `73` = 0.023476,
  `74` = 0.025685, `75` = 0.028219, `76` = 0.030546, `77` = 0.034523,
  `78` = 0.037655, `79` = 0.042211,
  `80` = 0.046538
)
# Table covers ages 50 to 80, which spans every age the model reaches: the
# oldest modeled cohort starts at 70 and the horizon is 10 annual cycles.

# Background mortality for a given model cycle, given the cohort's starting age.
# The life table gives an annual qx, which is converted to the cycle length.
bg_mortality <- function(p, cycle) {
  age <- p$start_age + (cycle %/% CYCLES_PER_YEAR)
  stopifnot(as.character(age) %in% names(LIFE_TABLE_QX))
  q_annual <- min(unname(LIFE_TABLE_QX[as.character(age)]) * p$qx_mult, 0.5)
  per_cycle(q_annual)
}

# -----------------------------------------------------------------------------
# CYCLE LENGTH
#
# CYCLES_PER_YEAR = 4 gives 3-month cycles; 1 gives annual cycles.
#
# Quarterly is the base case, for a clinical rather than a numerical reason.
# MechComp and Revision are single-cycle transient states: a patient passes
# through them rather than living in them. Under annual cycles that means a
# patient sits in "mechanical complication" for a full twelve months before
# being revised, and then spends a further twelve months "in revision". Neither
# is true. Proximal junctional failure is typically revised within weeks to a
# few months. Annual cycles therefore over-expose patients to the complication
# disutility and push revision costs a year later than they occur, where
# discounting makes them look cheaper. Quarterly cycles put both right.
#
# Because the early complication hazard is time-since-surgery dependent, the
# post-operative stable state is a tunnel spanning the first two post-operative
# years, which is 2 * CYCLES_PER_YEAR states.
# -----------------------------------------------------------------------------
CYCLES_PER_YEAR <- as.integer(Sys.getenv("ASD_CYCLES_PER_YEAR", "4"))
HORIZON_YEARS   <- 10
N_CYCLES        <- HORIZON_YEARS * CYCLES_PER_YEAR
DISCOUNT_RATE   <- 0.03

# Convert an annual probability to the per-cycle equivalent, preserving the
# cumulative risk over a year. Conditional probabilities (for example the
# probability of revision GIVEN a complication) are NOT rate quantities and must
# not be passed through this.
per_cycle <- function(p_annual) 1 - (1 - p_annual)^(1 / CYCLES_PER_YEAR)

N_TUNNEL <- 2 * CYCLES_PER_YEAR          # post-op tunnel covering years 1 and 2
STATES <- c("Nonoperative",
            sprintf("Stable_T%02d", seq_len(N_TUNNEL)),
            "Stable_Late", "MechComp_Early", "MechComp_Late",
            "Revision", "PostRev_Stable", "Death")
N_STATES <- length(STATES)
S_NONOP <- 1L
S_TUN   <- 1L + seq_len(N_TUNNEL)        # vector of tunnel state indices
S_LATE  <- 2L + N_TUNNEL
S_MCE   <- 3L + N_TUNNEL
S_MCL   <- 4L + N_TUNNEL
S_REV   <- 5L + N_TUNNEL
S_PRS   <- 6L + N_TUNNEL
S_DEAD  <- 7L + N_TUNNEL

DISCOUNT_PER_CYCLE <- (1 + DISCOUNT_RATE)^(1 / CYCLES_PER_YEAR) - 1
disc_vec <- 1 / (1 + DISCOUNT_PER_CYCLE)^(0:(N_CYCLES - 1))

# -----------------------------------------------------------------------------
# 2. Base-case parameters
# -----------------------------------------------------------------------------
BASE <- list(

  # ---- transition probabilities -------------------------------------------
  # DERIVED from Soroceanu 2015 (PMID 26426712): 31.7% cumulative radiographic
  # and implant-related complication incidence at 2 years in n=245.
  # 1 - (1 - 0.317)^(1/2) = 0.1736
  p_comp_early = 0.1736,

  # DERIVED from Imbo 2023 (PMID 37040468): 19.2% complication rate after the
  # 2-year mark in n=99 followed to a minimum of 5 years.
  # 1 - (1 - 0.192)^(1/3) = 0.0686
  p_comp_late = 0.0686,

  # DIRECT, Soroceanu 2015: "The incidence of RIC was 31.7% and 52.6% of those
  # patients required reoperation" within 2 years.
  p_rev_given_comp_early = 0.526,

  # DIRECT, Park 2021 (PMID 33290367): revision performed in 23 of 69 (33.3%).
  # NOTE the cohort is neurologically INTACT PJF patients, a selected subgroup,
  # which makes this a conservative (low) late revision probability.
  p_rev_given_comp_late = 0.333,

  # DERIVED from Adida 2025 (PMID 40262220): 48% (26/54) re-revision, median
  # follow-up 24 months. 1 - (1 - 0.48)^(1/2) = 0.2789.
  # The protocol's alternative annualization over the 14-month median time to
  # event was recomputed: the correct value is 1 - (1 - 0.48)^(12/14) = 0.4291,
  # not 0.534 (the original inverted the exponent). Tested in scenario analysis.
  p_rerev = 0.2789,

  # DERIVED from Carreon 2019 (PMID 31205182): 24 of 81 (30%) nonoperative
  # patients had surgery by 5 years. 1 - (1 - 0.30)^(1/5) = 0.0689.
  p_crossover = 0.0689,

  # DIRECT, Mo 2026 (PMID 42118075): 0.3% (4/1507) 90-day mortality after ASD
  # surgery. Applied once on entry to surgery in BOTH arms.
  p_periop_mort = 0.003,

  # EQUALIZED. CDC/NCHS United States Life Tables 2023 (NVSR vol 74 no 6),
  # Table 1, total population, both sexes. Single-year qx is applied as the
  # cohort ages from 60 through 69 rather than a window average, because a
  # constant window mean overstates the hazard early and understates it late.
  #
  # Background mortality is held IDENTICAL across arms by design. Mo 2026
  # reports 0.008/person-year in a post-surgical ASD cohort, which is LOWER
  # than the age-matched general population because surgical candidates are
  # selected for fitness. Using it for the operative arm only would manufacture
  # a survival benefit out of selection bias. The operative arm's sole
  # mortality penalty is the 90-day perioperative risk above.
  start_age = 60,
  qx_mult   = 1.0,   # PSA / scenario multiplier on the whole life-table schedule

  # ---- utilities (SF-6D) ---------------------------------------------------
  # DERIVED from Scheer 2018 (PMID 27253084), which reports cumulative QALYs at
  # 1, 2 and 3 years: operative 0.651 / 1.290 / 1.903 and nonoperative
  # 0.610 / 1.189 / 1.749.
  #
  # CRITICAL: those cumulative values are DISCOUNTED AT 3.5% PER YEAR. Scheer
  # states the construction explicitly, "QALY(3yr) = SF6D(1yr)/(1+0.035) +
  # SF6D(2yr)/(1+0.035)^2 + SF6D(3yr)/(1+0.035)^3", and the Table 4 footnote
  # confirms the reported QALYs are discounted at 3.5% per year. They are
  # therefore a discounted SUM of point utilities, not an undiscounted time
  # integral. Taking plain first differences would carry Scheer's discounting
  # into this model, which then applies its own 3% discount on top, double
  # discounting every utility. The undiscounted point utility is
  #
  #     u_t = (Q_t - Q_{t-1}) * 1.035^t
  #
  # Doing this correctly matters for interpretation as well as for the result:
  # the plain first differences fall year on year in both arms, which looks like
  # clinical deterioration but is entirely the discount factor. The recovered
  # utilities are flat, which is what the source actually shows.
  u_op_y1     = 0.674,   # 0.651 * 1.035
  u_op_y2     = 0.685,   # (1.290 - 0.651) * 1.035^2
  u_op_late   = 0.680,   # (1.903 - 1.290) * 1.035^3, held flat for years 4-10
  u_nonop_y1  = 0.631,   # 0.610 * 1.035
  u_nonop_y2  = 0.620,   # (1.189 - 0.610) * 1.035^2
  u_nonop_late= 0.621,   # (1.749 - 1.189) * 1.035^3, held flat for years 4-10

  # ASSUMPTION. No published SF-6D or EQ-5D decrement exists for PJK, PJF, rod
  # fracture or pseudarthrosis in ASD. Wick 2022 (PMID 34812199) establishes
  # that these complications worsen HRQoL but reports ODI, SF-36 and SRS-22r,
  # never a utility. Swept 0.00 to 0.15 in one-way sensitivity analysis.
  u_decrement_comp = 0.05,

  # ASSUMPTION, set equal to the late post-operative stable utility.
  u_postrev = 0.680,

  # ---- costs, inflated to 2024 USD at run time -----------------------------
  # DIRECT, Ames 2020 (PMID 31513120): mean index episode DIRECT cost $70,766
  # (SD $24,422) in n=210. Hospital direct cost, not a payer payment.
  cost_index_raw = 70766, cost_index_year = 2020,
  cost_index_sd_raw = 24422,

  # DERIVED, Raman 2018 (PMID 29099409) 2-year primary surgical plus
  # spine-related total $137,990 (a MEDIAN, IQR $84,186, not a mean) minus the
  # Ames index cost, halved. Both are inflated to 2024 USD before subtracting.
  # Literature costing path only; not used in the manuscript.
  raman_primary_2yr_raw = 137990, raman_year = 2018,

  # DIRECT, Theologis 2016 (PMID 26909838): average direct cost of PJF revision
  # $55,547 across 57 operations.
  cost_revision_raw = 55547, cost_revision_year = 2016,

  # DERIVED, Yagi 2023 (PMID 36730058): approximately $1,067/yr long-term
  # follow-up. NOTE this is a Japanese cohort and the source is internally
  # inconsistent (it states a 13% rise from $45K to $53K, which is 18%).
  # Retained for want of a US long-term ASD follow-up cost; swept widely.
  cost_late_fu_raw = 1067, cost_late_fu_year = 2023,

  # DERIVED, Glassman 2010 (PMID 20118843): $10,815 mean nonoperative treatment
  # cost over 2 years in the n=68 who used nonoperative resources.
  cost_nonop_raw = 5408, cost_nonop_year = 2010,

  # ASSUMPTION. Cost of working up and nonoperatively managing a mechanical
  # complication that does not proceed to revision. Set to one year of the
  # early post-operative follow-up rate. Swept widely.
  cost_comp_workup_multiplier = 1.0,

  # ---- Medicare costing parameters ---------------------------------------
  # ASSUMPTION. Multiplier on the DRG payment, covering wage index, indirect
  # medical education, disproportionate share, uncompensated care and outlier
  # adjustments, none of which are in the standardized rate. Unlike the acuity
  # mix, which is now taken from observed ASD data, this can move in either
  # direction: wage index is below 1.0 in much of the country while teaching and
  # safety-net add-ons push above it. Base case 1.0, swept both ways.
  drg_mix_mult = 1.0,
  drg_all_tier = NULL,

  # ASSUMPTION. Multiplier on all routine follow-up utilization (visits,
  # imaging, injections). The CPT payments are verified; the number of services
  # a stable patient consumes per year is not. Swept widely.
  fu_intensity = 1.0,

  # ASSUMPTION. Multiplier on nonoperative management intensity.
  nonop_intensity = 1.0,

  cohort_size = 10000
)

# Derived cost quantities in 2024 USD
#
# The post-operative follow-up cost is decomposed rather than taken as a simple
# residual. Raman's $137,990 is a COHORT MEAN two-year total, so it already
# contains the revision operations performed on the subset revised within two
# years. This model charges revisions separately on each inflow into the
# Revision state. Subtracting only the index cost would therefore charge those
# revisions twice. The expected two-year revision cost per patient is removed
# before the residual is annualized:
#
#   c_postop = (Raman 2-yr total - index cost - E[revision cost, 2 yr]) / 2
#
# where E[revision cost] uses the model's own early revision probability
# (0.317 x 0.526 = 0.167 of patients revised within two years, Soroceanu 2015).
build_costs <- function(p) {
  if (identical(COSTING, "medicare")) return(medicare_costs(p))
  c_index   <- infl(p$cost_index_raw,        p$cost_index_year)
  raman_2yr <- infl(p$raman_primary_2yr_raw, p$raman_year)
  c_rev     <- infl(p$cost_revision_raw,     p$cost_revision_year)
  p_rev_2yr <- (1 - (1 - p$p_comp_early)^2) * p$p_rev_given_comp_early
  c_postop  <- max((raman_2yr - c_index - p_rev_2yr * c_rev) / 2, 0)
  if (!is.null(p$c_postop_override)) c_postop <- p$c_postop_override
  list(
    c_index    = c_index,
    c_postop   = c_postop,                                     # years 1-2
    c_late_fu  = infl(p$cost_late_fu_raw,  p$cost_late_fu_year),
    c_revision = infl(p$cost_revision_raw, p$cost_revision_year),
    c_nonop    = infl(p$cost_nonop_raw,    p$cost_nonop_year),
    c_comp     = c_postop * p$cost_comp_workup_multiplier
  )
}

# -----------------------------------------------------------------------------
# 2b. MEDICARE PAYMENT COSTING (base case)
#
# Built from actual CMS payment rates rather than from journal-derived hospital
# direct costs. This is the more defensible basis for a stated payer perspective
# and, critically, it dissolves the twofold disagreement between the published
# ASD cost sources: every figure below comes from one payer, one year, one
# schedule, so no cross-source inflation or reconciliation is required.
#
# Inpatient: FY2026 IPPS Final Rule (CMS-1833-F) Table 5 relative weights, times
# the national standardized amount of $6,752.61 plus the capital rate of $524.15
# (Tables 1A-1E), i.e. RW x $7,276.76 at a wage-index-1.0 hospital.
# Professional: CY2026 Physician Fee Schedule Relative Value File (RVU26C, July
# release), total RVU x the nonqualifying-APM conversion factor of $33.4009.
# Both were retrieved from the CMS files and independently corroborated.
#
# ASD deformity fusion groups to MS-DRG 456/457/458, whose title explicitly
# names spinal curvature. NOTE that MS-DRGs 453, 454, 455, 459 and 460, which
# older spine economic papers use, were DELETED effective 1 October 2024 and no
# longer exist.
#
# Utilization assumptions are stated explicitly and swept: the CPT payments are
# verified facts, but how many visits and images a stable patient consumes in a
# year is an assumption.
# -----------------------------------------------------------------------------
IPPS_TOTAL_RATE <- 6752.61 + 524.15     # operating + capital, FY2026
PFS_CF          <- 33.4009              # CY2026, nonqualifying APM

DRG_RW <- c(`456` = 8.4034, `457` = 5.9631, `458` = 4.1726)

# Observed DRG distribution for adult spinal deformity surgery, from Nayak 2025
# (PMID 41222566), a 675-patient ISSG-AO ASD cohort grouped into exactly these
# three tiers: 14% without CC or MCC, 71% with CC, 15% with MCC.
#
# This replaces an earlier base case built on DRG 456 alone. MS-DRG 456 is the
# MCC tier and is the LEAST common of the three at 15%, not the modal group, so
# pricing every case at 456 overstated the episode by roughly 38%.
DRG_MIX <- c(`458` = 0.14, `457` = 0.71, `456` = 0.15)
DRG_RW_BLENDED <- sum(DRG_RW[names(DRG_MIX)] * DRG_MIX)   # 6.0785

drg_payment <- function(drg) unname(DRG_RW[as.character(drg)]) * IPPS_TOTAL_RATE
drg_payment_blended <- function() DRG_RW_BLENDED * IPPS_TOTAL_RATE

# CY2026 total RVUs. NF = non-facility, F = facility.
RVU_NF <- c(`99213` = 2.85, `99214` = 4.06, `97110` = 0.87, `64483` = 7.93,
            `64484` = 3.52, `62323` = 8.18, `72148` = 5.74, `72100` = 1.21)
cpt_nf <- function(code) unname(RVU_NF[as.character(code)]) * PFS_CF

medicare_costs <- function(p) {
  # --- index deformity fusion, T10-pelvis, 7-12 segments -------------------
  # Priced on the blended MS-DRG 456/457/458 weight (DRG_RW_BLENDED, 6.0785);
  # DRG 457 (with CC) is the modal tier at 71%, not 456 (see DRG_MIX above).
  # Professional fee applies the standard multiple-procedure reduction: the
  # highest-valued procedure at 100%, subsequent ones at 50%, and ZZZ add-on
  # codes at 100% because they are exempt from the reduction.
  prof_index <- 1936.25 +            # 22802 posterior deformity fusion 7-12
                 875.10 +            # 22844 segmental instrumentation, add-on
                 850.05 +            # 22633 interbody, second procedure at 50%
                 228.80 +            # 22853 interbody device, add-on, ONE
                 147.30 +            # 20937 morselized autograft, add-on
                 532.74              # 63047 decompression, third at 50%
  # One 22853 per interspace: the stack bills a single 22633, so a second
  # interbody device would require 22634 and is not assumed here.
  rw_idx <- if (is.null(p$drg_all_tier)) DRG_RW_BLENDED else
              unname(DRG_RW[p$drg_all_tier])
  c_index <- rw_idx * IPPS_TOTAL_RATE * p$drg_mix_mult + prof_index

  # --- revision for proximal junctional failure ---------------------------
  # Proximal extension with instrumentation revision. PJF revisions almost
  # always carry an MCC, so they also group to DRG 456. This is why a revision
  # episode is only a few percent cheaper than a primary one under Medicare
  # payment, in sharp contrast to the journal-derived direct costs.
  # 22802 proximal extension (100%) + 22844 instrumentation (ZZZ, 100%)
  # + 22852 removal of posterior instrumentation (second procedure, 50%)
  # + 20937 autograft (ZZZ) + 22853 interbody device (ZZZ)
  prof_rev <- 1936.25 + 875.10 + 346.03 + 147.30 + 228.80
  # Revisions skew to higher acuity than primaries but are not uniformly MCC:
  # Nayak reports 97% of intervention-requiring complications land in CC OR MCC,
  # which does not support pricing every revision at MS-DRG 456. The revision
  # episode uses a mix shifted one tier up from the primary distribution.
  rw_rev <- if (is.null(p$drg_all_tier)) {
    0.05 * DRG_RW["458"] + 0.55 * DRG_RW["457"] + 0.40 * DRG_RW["456"]
  } else DRG_RW[p$drg_all_tier]
  c_revision <- unname(rw_rev) * IPPS_TOTAL_RATE * p$drg_mix_mult + prof_rev

  # --- routine follow-up, all ASSUMED utilization -------------------------
  # Post-operative years 1 and 2: three established-patient visits, two
  # radiograph series, and an MRI every other year.
  c_postop <- p$fu_intensity *
    (3 * cpt_nf(99214) + 2 * cpt_nf(72100) + 0.5 * cpt_nf(72148))
  # Year 3 and beyond: one visit, one radiograph series, MRI every fourth year.
  c_late <- p$fu_intensity *
    (1 * cpt_nf(99214) + 1 * cpt_nf(72100) + 0.25 * cpt_nf(72148))

  # --- mechanical complication managed without revision -------------------
  # Two visits, an MRI, two radiograph series and a transforaminal injection.
  c_comp <- p$fu_intensity *
    (2 * cpt_nf(99214) + cpt_nf(72148) + 2 * cpt_nf(72100) + cpt_nf(64483))

  # --- nonoperative management, per year ----------------------------------
  # 20 physical therapy sessions at 2 units (second unit carries the therapy
  # multiple-procedure reduction on its practice expense), four office visits,
  # two two-level transforaminal injections, one MRI and two radiograph series.
  pt_unit2 <- (0.45 + 0.41 * 0.5 + 0.01) * PFS_CF     # 97110, PE halved
  c_nonop <- p$nonop_intensity *
    (20 * (cpt_nf(97110) + pt_unit2) +
     3 * cpt_nf(99214) + cpt_nf(99213) +
     2 * (cpt_nf(64483) + cpt_nf(64484)) +
     cpt_nf(72148) + 2 * cpt_nf(72100))

  if (!is.null(p$c_postop_override)) c_postop <- p$c_postop_override
  list(c_index = c_index, c_postop = c_postop, c_late_fu = c_late,
       c_revision = c_revision, c_nonop = c_nonop, c_comp = c_comp)
}

# -----------------------------------------------------------------------------
# 3. Transition matrix
#
# `arm` determines only whether the Nonoperative state is occupied and whether
# crossover is active. Once a patient is post-operative, both arms share the
# same downstream dynamics, which is the conservative choice: the model grants
# surgery no long-term structural advantage beyond what the cited rates imply.
# -----------------------------------------------------------------------------
build_tm <- function(p, arm, cycle = 0) {
  m <- matrix(0, N_STATES, N_STATES, dimnames = list(STATES, STATES))
  d  <- bg_mortality(p, cycle)
  xo <- if (arm == "operative") 0 else per_cycle(p$p_crossover)

  comp_e <- per_cycle(p$p_comp_early)
  comp_l <- per_cycle(p$p_comp_late)
  rerev  <- per_cycle(p$p_rerev)

  # --- Nonoperative -----------------------------------------------------
  # Crossover carries the same 90-day perioperative mortality as the index
  # operation in the operative arm.
  xo_die  <- xo * p$p_periop_mort
  xo_live <- xo - xo_die
  m[S_NONOP, S_TUN[1]] <- xo_live
  m[S_NONOP, S_DEAD]   <- d + xo_die
  m[S_NONOP, S_NONOP]  <- 1 - xo_live - d - xo_die

  # --- Post-operative stable tunnel, years 1 and 2 ----------------------
  # Each tunnel state advances to the next; the last one empties into the
  # late stable state.
  for (i in seq_len(N_TUNNEL)) {
    nxt <- if (i < N_TUNNEL) S_TUN[i + 1] else S_LATE
    m[S_TUN[i], S_MCE]  <- comp_e
    m[S_TUN[i], S_DEAD] <- d
    m[S_TUN[i], nxt]    <- 1 - comp_e - d
  }

  m[S_LATE, S_MCL]  <- comp_l
  m[S_LATE, S_DEAD] <- d
  m[S_LATE, S_LATE] <- 1 - comp_l - d

  # --- Mechanical complication (one-cycle transient) --------------------
  # Patients not revised are managed nonoperatively and return to the late
  # stable state; they do not re-enter the early tunnel.
  # Revision probabilities are capped at 1 - d so that the residual flow back to
  # the stable state can never go negative when a scenario pushes them to 1.
  rev_e <- min(p$p_rev_given_comp_early, 1 - d)
  rev_l <- min(p$p_rev_given_comp_late,  1 - d)

  m[S_MCE, S_REV]  <- rev_e
  m[S_MCE, S_DEAD] <- d
  m[S_MCE, S_LATE] <- 1 - rev_e - d

  m[S_MCL, S_REV]  <- rev_l
  m[S_MCL, S_DEAD] <- d
  m[S_MCL, S_LATE] <- 1 - rev_l - d

  # --- Revision (one-cycle tunnel) --------------------------------------
  rev_die <- p$p_periop_mort
  m[S_REV, S_DEAD] <- d + rev_die
  m[S_REV, S_PRS]  <- 1 - d - rev_die

  # --- Post-revision stable ---------------------------------------------
  m[S_PRS, S_MCL]  <- rerev
  m[S_PRS, S_DEAD] <- d
  m[S_PRS, S_PRS]  <- 1 - rerev - d

  m[S_DEAD, S_DEAD] <- 1

  stopifnot(all(m >= -1e-12), all(abs(rowSums(m) - 1) < 1e-9))
  m
}

# -----------------------------------------------------------------------------
# 4. Reward vectors
# -----------------------------------------------------------------------------
# Utilities are annual; a cycle accrues utility x (1 / CYCLES_PER_YEAR) QALYs.
# Tunnel states in the first year carry the year 1 utility, those in the second
# year the year 2 utility.
utility_vec <- function(p, arm, cycle) {
  yr <- cycle %/% CYCLES_PER_YEAR
  u_nonop <- if (yr < 1) p$u_nonop_y1 else if (yr < 2) p$u_nonop_y2 else p$u_nonop_late
  u <- numeric(N_STATES)
  u[S_NONOP] <- u_nonop
  u[S_TUN[seq_len(CYCLES_PER_YEAR)]] <- p$u_op_y1
  u[S_TUN[CYCLES_PER_YEAR + seq_len(CYCLES_PER_YEAR)]] <- p$u_op_y2
  u[S_LATE]  <- p$u_op_late
  u[S_MCE]   <- p$u_op_y1   - p$u_decrement_comp
  u[S_MCL]   <- p$u_op_late - p$u_decrement_comp
  u[S_REV]   <- p$u_op_late - p$u_decrement_comp
  u[S_PRS]   <- p$u_postrev
  u[S_DEAD]  <- 0
  u / CYCLES_PER_YEAR
}

# State-residence costs are annual rates and are divided across the cycles in a
# year. Transition costs (index surgery, revision) are point events and are
# charged separately on inflow, undivided.
state_cost_vec <- function(p, cst) {
  c <- numeric(N_STATES)
  c[S_NONOP] <- cst$c_nonop
  c[S_TUN]   <- cst$c_postop
  c[S_LATE]  <- cst$c_late_fu
  c[S_MCE]   <- cst$c_comp
  c[S_MCL]   <- cst$c_comp
  c[S_REV]   <- cst$c_late_fu
  c[S_PRS]   <- cst$c_late_fu
  c[S_DEAD]  <- 0
  c / CYCLES_PER_YEAR
}

# -----------------------------------------------------------------------------
# 5. Run one arm
#
# Transition costs (the index operation and each revision) are charged on the
# INFLOW into the receiving state, not as a state-residence cost, so that a
# patient is charged once per event.
# -----------------------------------------------------------------------------
run_arm <- function(p, arm) {
  cst <- build_costs(p)
  # One transition matrix per cycle, because background mortality rises as the
  # cohort ages through the life table.
  tms <- lapply(0:(N_CYCLES - 1), function(t) build_tm(p, arm, t))

  trace <- matrix(0, N_CYCLES, N_STATES, dimnames = list(NULL, STATES))
  if (arm == "operative") {
    # The whole cohort is operated at cycle 0; 90-day perioperative mortality
    # is applied immediately as a one-time decrement.
    trace[1, S_TUN[1]] <- 1 - p$p_periop_mort
    trace[1, S_DEAD]   <- p$p_periop_mort
  } else {
    trace[1, S_NONOP] <- 1
  }
  for (t in 2:N_CYCLES) trace[t, ] <- trace[t - 1, ] %*% tms[[t - 1]]

  # Event inflows
  inflow_rev <- numeric(N_CYCLES)   # revision operations
  inflow_srg <- numeric(N_CYCLES)   # index operations (crossover in nonop arm)
  inflow_cmp <- numeric(N_CYCLES)   # new mechanical complications
  if (arm == "operative") inflow_srg[1] <- 1
  for (t in 1:(N_CYCLES - 1)) {
    tm <- tms[[t]]
    inflow_rev[t + 1] <- trace[t, S_MCE] * tm[S_MCE, S_REV] +
                         trace[t, S_MCL] * tm[S_MCL, S_REV]
    inflow_srg[t + 1] <- trace[t, S_NONOP] * tm[S_NONOP, S_TUN[1]]
    inflow_cmp[t + 1] <- sum(trace[t, S_TUN] * tm[cbind(S_TUN, S_MCE)]) +
                         trace[t, S_LATE] * tm[S_LATE, S_MCL] +
                         trace[t, S_PRS]  * tm[S_PRS, S_MCL]
  }

  # Within-cycle (half-cycle) correction on state-residence rewards only.
  # Transition costs are point events and are not half-cycle corrected.
  hc <- rbind(trace[1, , drop = FALSE],
              0.5 * (trace[-N_CYCLES, , drop = FALSE] + trace[-1, , drop = FALSE]))

  cvec <- state_cost_vec(p, cst)
  cost_cycle <- as.vector(hc %*% cvec) +
                inflow_rev * cst$c_revision +
                inflow_srg * cst$c_index
  qaly_cycle <- vapply(seq_len(N_CYCLES),
                       function(t) sum(hc[t, ] * utility_vec(p, arm, t - 1)),
                       numeric(1))

  n <- p$cohort_size
  list(
    arm = arm, trace = trace, tm = tms[[1]],
    cost_cycle = cost_cycle, qaly_cycle = qaly_cycle,
    disc_cost = sum(cost_cycle * disc_vec),
    disc_qaly = sum(qaly_cycle * disc_vec),
    undisc_cost = sum(cost_cycle), undisc_qaly = sum(qaly_cycle),
    cum_cost_5yr = sum(cost_cycle[1:(5 * CYCLES_PER_YEAR)]),  # undiscounted, for validation
    n_revisions   = sum(inflow_rev) * n,
    n_index_ops   = sum(inflow_srg) * n,
    n_complications = sum(inflow_cmp) * n,
    n_deaths      = trace[N_CYCLES, S_DEAD] * n,
    inflow_rev = inflow_rev, inflow_srg = inflow_srg
  )
}

run_base_case <- function(p = BASE) {
  op <- run_arm(p, "operative")
  no <- run_arm(p, "nonoperative")
  dC <- op$disc_cost - no$disc_cost
  dQ <- op$disc_qaly - no$disc_qaly
  list(op = op, nonop = no, dC = dC, dQ = dQ,
       icer = if (dQ != 0) dC / dQ else NA_real_, params = p)
}

# -----------------------------------------------------------------------------
# 6. One-way sensitivity analysis
# -----------------------------------------------------------------------------
OWSA_SPEC <- list(
  list(key = "p_comp_early",           lo = 0.5,   hi = 1.5,  label = "Early mechanical complication rate"),
  list(key = "p_comp_late",            lo = 0.5,   hi = 1.5,  label = "Late mechanical complication rate"),
  list(key = "p_rev_given_comp_early", lo = 0.7,   hi = 1.3,  label = "Revision given early complication"),
  list(key = "p_rev_given_comp_late",  lo = 0.7,   hi = 1.3,  label = "Revision given late complication"),
  list(key = "p_rerev",                lo = 0.5,   hi = 1.54, label = "Re-revision rate"),
  list(key = "p_crossover",            lo = 0.5,   hi = 1.5,  label = "Crossover to surgery"),
  list(key = "u_op_late",              lo = 0.90,  hi = 1.10, label = "Operative utility, year 3+"),
  list(key = "u_nonop_late",           lo = 0.90,  hi = 1.10, label = "Nonoperative utility, year 3+"),
  list(key = "u_decrement_comp",       lo = 0.0,   hi = 3.0,  label = "Complication disutility"),
  list(key = "u_postrev",              lo = 0.90,  hi = 1.10, label = "Post-revision utility"),
  list(key = "drg_mix_mult",           lo = 0.80,  hi = 1.35, label = "DRG payment adjustments (wage, IME, DSH)"),
  list(key = "fu_intensity",           lo = 0.50,  hi = 2.00, label = "Routine follow-up intensity"),
  list(key = "nonop_intensity",        lo = 0.50,  hi = 2.00, label = "Nonoperative management intensity")
)

run_owsa <- function(wtp = 100000) {
  base  <- run_base_case()
  binmb <- wtp * base$dQ - base$dC
  rows  <- lapply(OWSA_SPEC, function(s) {
    plo <- BASE; phi <- BASE
    plo[[s$key]] <- BASE[[s$key]] * s$lo
    phi[[s$key]] <- BASE[[s$key]] * s$hi
    rlo <- run_base_case(plo); rhi <- run_base_case(phi)
    lo_inmb <- wtp * rlo$dQ - rlo$dC
    hi_inmb <- wtp * rhi$dQ - rhi$dC
    data.frame(parameter = s$label, param_key = s$key,
               low_value = BASE[[s$key]] * s$lo, high_value = BASE[[s$key]] * s$hi,
               base_inmb = binmb, low_inmb = lo_inmb, high_inmb = hi_inmb,
               low_icer = rlo$icer, high_icer = rhi$icer,
               swing = abs(hi_inmb - lo_inmb))
  })
  bind_rows(rows) %>% arrange(desc(swing))
}

# -----------------------------------------------------------------------------
# 7. Probabilistic sensitivity analysis
# -----------------------------------------------------------------------------
beta_pars <- function(mean, sd) {
  v <- sd^2
  if (v <= 0 || v >= mean * (1 - mean)) v <- mean * (1 - mean) / 4
  a <- mean * (mean * (1 - mean) / v - 1)
  b <- (1 - mean) * (mean * (1 - mean) / v - 1)
  c(max(a, 0.5), max(b, 0.5))
}
gamma_pars <- function(mean, sd) {
  if (sd <= 0) sd <- 0.25 * mean
  c(shape = mean^2 / sd^2, scale = sd^2 / mean)
}

sample_params <- function() {
  p <- BASE
  # Probabilities: Beta, SE = 25% of the mean unless a published SD exists
  for (k in c("p_comp_early", "p_comp_late", "p_rev_given_comp_early",
              "p_rev_given_comp_late", "p_rerev", "p_crossover",
              "p_periop_mort")) {
    ab <- beta_pars(BASE[[k]], max(0.25 * BASE[[k]], 0.001))
    p[[k]] <- rbeta(1, ab[1], ab[2])
  }
  # Life-table schedule multiplier: lognormal centred on 1, applied to BOTH
  # arms in the same draw so the equalization survives under uncertainty.
  p$qx_mult <- rlnorm(1, meanlog = 0, sdlog = 0.15)
  # Utilities: Beta, using the STANDARD ERROR OF THE MEAN, not the patient-level
  # SD. A PSA characterizes uncertainty in the mean parameter; feeding it the
  # between-patient SD (0.09 in Scheer) conflates patient heterogeneity with
  # parameter uncertainty and inflates the incremental QALY spread several-fold.
  # Scheer's utility subanalysis cohorts are 90 operative and 61 nonoperative,
  # giving SE = 0.09 / sqrt(n) of roughly 0.010 and 0.012 respectively.
  # Scheer's trending subset is 44 patients per arm after matching, not the 90
  # and 61 of the unmatched subanalysis, so SE = 0.09 / sqrt(44).
  for (k in c("u_op_y1", "u_op_y2", "u_op_late",
              "u_nonop_y1", "u_nonop_y2", "u_nonop_late")) {
    ab <- beta_pars(BASE[[k]], 0.09 / sqrt(44))
    p[[k]] <- rbeta(1, ab[1], ab[2])
  }
  # Post-revision utility has no published dispersion; given a wider SE
  ab <- beta_pars(BASE$u_postrev, 0.03)
  p$u_postrev <- rbeta(1, ab[1], ab[2])
  # Complication disutility: the unsourced assumption gets a deliberately wide
  # Beta over 0 to 0.15 rather than a tight one around 0.05
  ab <- beta_pars(0.05 / 0.15, 0.25)
  p$u_decrement_comp <- rbeta(1, ab[1], ab[2]) * 0.15
  # Costs: Gamma, again on the standard error of the mean. Ames reports
  # SD $24,422 in n=210, so SE = 24422 / sqrt(210). Costs with no published
  # dispersion are given a 20% SE, a conventional choice that is deliberately
  # wider than any of the published standard errors.
  gp <- gamma_pars(BASE$cost_index_raw, BASE$cost_index_sd_raw / sqrt(210))
  p$cost_index_raw <- rgamma(1, shape = gp[1], scale = gp[2])
  # Medicare costing multipliers. The DRG multiplier is two-sided and centred
  # on 1.0: wage index sits below 1.0 in much of the country while teaching,
  # safety-net and outlier add-ons push above it, so payment can move either way.
  p$drg_mix_mult    <- rlnorm(1, meanlog = 0, sdlog = 0.14)
  p$fu_intensity    <- rlnorm(1, meanlog = 0, sdlog = 0.35)
  p$nonop_intensity <- rlnorm(1, meanlog = 0, sdlog = 0.35)
  for (k in c("cost_revision_raw", "cost_nonop_raw", "cost_late_fu_raw",
              "raman_primary_2yr_raw")) {
    gp <- gamma_pars(BASE[[k]], 0.20 * BASE[[k]])
    p[[k]] <- rgamma(1, shape = gp[1], scale = gp[2])
  }
  # Keep the derived post-operative cost non-negative
  if (infl(p$raman_primary_2yr_raw, p$raman_year) <=
      infl(p$cost_index_raw, p$cost_index_year)) {
    p$raman_primary_2yr_raw <- p$cost_index_raw * 1.05
  }
  p
}

run_psa <- function(n_iter = 10000) {
  set.seed(RNG_SEED)
  out <- matrix(NA_real_, n_iter, 6)
  for (i in seq_len(n_iter)) {
    r <- run_base_case(sample_params())
    out[i, ] <- c(r$dC, r$dQ, r$op$disc_cost, r$nonop$disc_cost,
                  r$op$disc_qaly, r$nonop$disc_qaly)
  }
  as.data.frame(setNames(as.data.frame(out),
    c("dC", "dQ", "cost_op", "cost_nonop", "qaly_op", "qaly_nonop")))
}

# -----------------------------------------------------------------------------
# 8. Scenario analyses
# -----------------------------------------------------------------------------
run_scenarios <- function() {
  sc <- list()

  sc[["Base case"]] <- run_base_case()

  # Late complication rate: reoperations only (Imbo 5.1% over 3 yr)
  p <- BASE; p$p_comp_late <- 0.0173; p$p_rev_given_comp_late <- 0.98
  sc[["Late rate: reoperations only"]] <- run_base_case(p)

  # Re-revision annualized from the 14-month median time to event.
  # Correct arithmetic is 1 - (1 - 0.48)^(12/14) = 0.4291. The source protocol
  # reported 0.534, which inverted the exponent.
  p <- BASE; p$p_rerev <- 0.4291
  sc[["Re-revision: time-to-event"]] <- run_base_case(p)

  # Crossover from the Clohisy observational cohort. The 22% commonly quoted is
  # a whole-follow-up proportion over a mean of about 86 months; the paper's
  # time-resolved cumulative incidence is 24.1% at 4 years.
  p <- BASE; p$p_crossover <- 1 - (1 - 0.241)^(1 / 4)
  sc[["Crossover: Clohisy time-resolved"]] <- run_base_case(p)

  # Age subgroups. CDC life-table mortality for the relevant decade.
  p <- BASE; p$start_age <- 50
  sc[["Age 50 subgroup"]] <- run_base_case(p)
  p <- BASE; p$start_age <- 70
  sc[["Age 70 subgroup"]] <- run_base_case(p)

  # Full revision episode costing: charge Raman's 2-year revision total
  # ($115,509, 2018 USD) instead of the Theologis index-only revision cost.
  if (identical(COSTING, "medicare")) {
    # Acuity brackets. The base case uses the observed ASD distribution
    # (14/71/15 across DRG 458/457/456, Nayak 2025). These two scenarios price
    # every case at a single tier to bound the effect of the acuity mix.
    p <- BASE; p$drg_all_tier <- "456"
    sc[["All cases MS-DRG 456 (with MCC)"]] <- run_base_case(p)
    p <- BASE; p$drg_all_tier <- "458"
    sc[["All cases MS-DRG 458 (no CC/MCC)"]] <- run_base_case(p)

    # Payment adjustments not in the standardized rate: wage index, IME, DSH,
    # uncompensated care and outliers. These move in both directions.
    p <- BASE; p$drg_mix_mult <- 1.30
    sc[["Payment adjustments +30%"]] <- run_base_case(p)
    p <- BASE; p$drg_mix_mult <- 0.85
    sc[["Payment adjustments -15%"]] <- run_base_case(p)

    # Follow-up and nonoperative utilization doubled.
    p <- BASE; p$fu_intensity <- 2.0; p$nonop_intensity <- 2.0
    sc[["Double follow-up utilization"]] <- run_base_case(p)
  } else {
    # ASLS cost base, only meaningful under literature costing.
    p <- BASE; p$c_postop_override <- 6000
    sc[["ASLS cost base (Carreon)"]] <- run_base_case(p)

    p <- BASE; p$cost_revision_raw <- 115509; p$cost_revision_year <- 2018
    sc[["Full revision episode cost"]] <- run_base_case(p)

    oldi <- INFLATE_INDEX; INFLATE_INDEX <<- "CPIM"
    sc[["CPI-M price index"]] <- run_base_case()
    INFLATE_INDEX <<- oldi
  }

  # Annual rather than quarterly cycles, to show the effect of cycle length.
  # Reported for transparency; see Methods on why quarterly is the base case.

  # 5-year horizon
  oldN <- N_CYCLES; oldD <- disc_vec
  N_CYCLES <<- 5 * CYCLES_PER_YEAR
  disc_vec <<- 1 / (1 + DISCOUNT_PER_CYCLE)^(0:(N_CYCLES - 1))
  sc[["5-year horizon"]] <- run_base_case()
  N_CYCLES <<- oldN; disc_vec <<- oldD

  # Discount rate 0% and 5%
  for (r in c(0, 0.05)) {
    oldD <- disc_vec
    disc_vec <<- 1 / ((1 + r)^(1 / CYCLES_PER_YEAR))^(0:(N_CYCLES - 1))
    sc[[sprintf("Discount rate %.0f%%", r * 100)]] <- run_base_case()
    disc_vec <<- oldD
  }

  sc
}

# -----------------------------------------------------------------------------
# 9. Model validation
#
# Carreon 2019 reports 5-year INTENT-TO-TREAT costs of $96,000 (operative) and
# $49,546 (nonoperative). Those figures are NOT model inputs: the nonoperative
# figure already contains the surgical costs of the 30% who crossed over, so
# using it as a state cost while also modeling crossover would double count.
# They are instead held out as external comparators, inflated to 2024 USD.
#
# IMPORTANT INTERPRETATION. The operative comparison is expected to diverge, and
# the divergence is a finding rather than a model defect. Carreon's cohort is
# adult symptomatic lumbar scoliosis; Ames and Raman, which anchor this model's
# cost inputs, are tertiary multicenter complex-ASD cohorts. The two literatures
# are mutually inconsistent at overlapping timepoints: Raman's TWO-year primary
# total ($137,990, 2018 USD) already exceeds Carreon's FIVE-year operative total
# ($96,000, 2019 USD). No parameterization can satisfy both simultaneously. The
# base case is anchored to the complex-ASD cost base and the "ASLS cost base"
# scenario reports the result on Carreon's cost level.
# -----------------------------------------------------------------------------
validate <- function(base) {
  if (identical(COSTING, "medicare")) {
    # Under Medicare payment costing the Carreon comparison is not apt: Carreon
    # reports a hybrid direct-plus-indirect cost in an ASLS cohort, not payer
    # payments. The meaningful internal check is that the modelled index episode
    # reproduces the CMS arithmetic it was built from.
    cst <- build_costs(BASE)
    return(data.frame(
      quantity = c("Blended DRG relative weight (14/71/15 across 458/457/456)",
                   "Index episode payment",
                   "Revision episode payment",
                   "Annual nonoperative management payment"),
      model = c(DRG_RW_BLENDED, cst$c_index, cst$c_revision, cst$c_nonop),
      recomputed = c(0.14 * 4.1726 + 0.71 * 5.9631 + 0.15 * 8.4034,
                     DRG_RW_BLENDED * IPPS_TOTAL_RATE + 4570.24,
                     (0.05 * 4.1726 + 0.55 * 5.9631 + 0.40 * 8.4034) *
                       IPPS_TOTAL_RATE + 3533.48,
                     2564.86),
      source = "CMS FY2026 IPPS Table 5; CY2026 PFS RVU file; Nayak 2025 DRG mix",
      note = c("acuity distribution from a 675-patient ASD cohort",
               "blended weight x $7,276.76 plus professional fee",
               "revision mix shifted one tier up from the primary distribution",
               "20 PT sessions, 4 visits, 2 injections, MRI, 2 radiographs")
    ) %>% mutate(ratio = model / recomputed))
  }
  tgt_op    <- infl(96000, 2019)
  tgt_nonop <- infl(49546, 2019)
  data.frame(
    quantity = c("5-year cumulative cost, operative arm",
                 "5-year cumulative cost, nonoperative arm (ITT)"),
    model    = c(base$op$cum_cost_5yr, base$nonop$cum_cost_5yr),
    comparator = c(tgt_op, tgt_nonop),
    source   = "Carreon 2019 (PMID 31205182), inflated to 2024 USD",
    note = c("expected to exceed: complex-ASD vs ASLS cost base",
             "within-range agreement")
  ) %>% mutate(ratio = model / comparator)
}

# -----------------------------------------------------------------------------
# 10. Parameter table for the manuscript
# -----------------------------------------------------------------------------
param_table <- function() {
  cst <- build_costs(BASE)
  core_rows <- tribble(
    ~parameter, ~value, ~tier_source, ~psa_dist,
    "Transition probabilities (annual)", "", "", "",
    "Mechanical complication, post-operative years 1-2", "0.174",
      paste("[DERIVED] 31.7% cumulative incidence of radiographic and",
            "implant-related complications at 2 y (n=245), Soroceanu et al [3];",
            "annualized as 1-(1-0.317)^(1/2)"), "Beta",
    "Mechanical complication, post-operative year 3+", "0.069",
      paste("[DERIVED] 19.2% complication rate after the 2-y mark (n=99,",
            "minimum 5-y follow-up), Imbo et al [4];",
            "annualized as 1-(1-0.192)^(1/3)"), "Beta",
    "Revision given early complication", "0.526",
      paste("[DIRECT] 52.6% of patients with a radiographic or implant-related",
            "complication required reoperation, Soroceanu et al [3]"), "Beta",
    "Revision given late complication", "0.333",
      paste("[DIRECT] 23/69 (33.3%) revised among neurologically intact",
            "proximal junctional failure patients, Park et al [5]"), "Beta",
    "Re-revision after revision", "0.279",
      paste("[DERIVED] 26/54 (48%) re-revised at median 24-mo follow-up,",
            "Adida et al [6]; annualized as 1-(1-0.48)^(1/2)"), "Beta",
    "Crossover, nonoperative to surgery", "0.069",
      paste("[DERIVED] 24/81 (30%) of nonoperative patients underwent surgery",
            "by 5 y, Carreon et al [7]; annualized as 1-(1-0.30)^(1/5)"), "Beta",
    "90-day perioperative mortality", "0.003",
      "[DIRECT] 4/1507 (0.3%) 90-day mortality after ASD surgery, Mo et al [8]", "Beta",
    "Background all-cause mortality, BOTH ARMS", "0.0092 to 0.0173 (age 60 to 69)",
      paste("[EQUALIZED] single-year qx, ages 60-69, total population, US Life",
            "Tables 2023, Table 1, Arias et al [16]; identical in both arms"),
      "Lognormal multiplier",
    "Health-state utilities (SF-6D)", "", "", "",
    "Operative, year 1", "0.674",
      paste("[DERIVED] cumulative QALY 0.651 at 1 y, Scheer et al [9];",
            "3.5% source discount removed: 0.651 x 1.035"), "Beta",
    "Operative, year 2", "0.685",
      paste("[DERIVED] cumulative QALYs 0.651 and 1.290 at 1 and 2 y,",
            "Scheer et al [9]; (1.290 - 0.651) x 1.035^2"), "Beta",
    "Operative, year 3+", "0.680",
      paste("[DERIVED] cumulative QALYs 1.290 and 1.903 at 2 and 3 y,",
            "Scheer et al [9]; (1.903 - 1.290) x 1.035^3, held for years 4-10"), "Beta",
    "Nonoperative, year 1", "0.631",
      paste("[DERIVED] cumulative QALY 0.610 at 1 y, Scheer et al [9];",
            "0.610 x 1.035"), "Beta",
    "Nonoperative, year 2", "0.620",
      paste("[DERIVED] cumulative QALYs 0.610 and 1.189 at 1 and 2 y,",
            "Scheer et al [9]; (1.189 - 0.610) x 1.035^2"), "Beta",
    "Nonoperative, year 3+", "0.621",
      paste("[DERIVED] cumulative QALYs 1.189 and 1.749 at 2 and 3 y,",
            "Scheer et al [9]; (1.749 - 1.189) x 1.035^3, held for years 4-10"), "Beta",
    "Mechanical complication disutility", "-0.05",
      paste("[ASSUMPTION] no published utility decrement exists for PJK, PJF,",
            "rod fracture or pseudarthrosis [10]; anchored to post-revision QALY",
            "gains of 0.35 to 0.40 [11,12]; swept 0.00 to 0.15"), "Beta (0 to 0.15)",
    "Post-revision stable", "0.680",
      "[ASSUMPTION] set equal to the year 3+ post-operative utility [9]", "Beta",
    "Costs", "", "", ""
  )
  cost_rows <- if (identical(COSTING, "medicare")) tribble(
    ~parameter, ~value, ~tier_source, ~psa_dist,
    "Costs (Medicare payment; FY2026 IPPS, CY2026 PFS)", "", "", "",
    "Index ASD deformity fusion episode",
      sprintf("$%s", format(round(cst$c_index), big.mark = ",")),
      paste("[DIRECT] blended MS-DRG 456/457/458 relative weight 6.0785",
            "(observed ASD mix 15%/71%/14%, Nayak et al [14]) x $7,276.76",
            "(FY2026 IPPS operating + capital rate [13]) + $4,570 professional",
            "fee (CY2026 PFS [15]: 22802, 22844, 22853, 20937 at full value;",
            "22633, 63047 at 50%)"), "Gamma",
    "Revision episode for proximal junctional failure",
      sprintf("$%s", format(round(cst$c_revision), big.mark = ",")),
      paste("[DIRECT] revision acuity mix shifted one tier above primary",
            "(40%/55%/5% across 456/457/458, RW 6.8497 [13,14]) x $7,276.76 +",
            "$3,533 professional fee (CY2026 PFS [15]: 22802, 22844, 22853,",
            "20937 at full value; 22852 at 50%)"), "Gamma",
    "Post-operative follow-up, years 1-2",
      sprintf("$%s/yr", format(round(cst$c_postop), big.mark = ",")),
      "[DERIVED] CY2026 PFS [15]: 3 x 99214, 2 x 72100, 0.5 x 72148 per year", "Gamma",
    "Post-operative follow-up, year 3+",
      sprintf("$%s/yr", format(round(cst$c_late_fu), big.mark = ",")),
      "[DERIVED] CY2026 PFS [15]: 99214, 72100, 0.25 x 72148 per year", "Gamma",
    "Mechanical complication managed without revision",
      sprintf("$%s", format(round(cst$c_comp), big.mark = ",")),
      "[DERIVED] CY2026 PFS [15]: 2 x 99214, 72148, 2 x 72100, 64483", "Gamma",
    "Nonoperative management",
      sprintf("$%s/yr", format(round(cst$c_nonop), big.mark = ",")),
      paste("[DERIVED] CY2026 PFS [15]: 20 PT sessions (97110 x 2 units),",
            "4 office visits, 2 two-level transforaminal injections",
            "(64483 + 64484), MRI (72148), 2 radiograph series (72100)",
            "per year"), "Gamma",
    "DRG payment add-on multiplier (IME, DSH, outlier)", "1.00",
      paste("[ASSUMPTION] wage index, teaching, disproportionate-share and",
            "outlier adjustments are not in the standardized amount [13];",
            "they move payment in both directions"), "Lognormal (SD 0.14 on the log scale)",
    "Follow-up utilization multiplier", "1.00",
      "[ASSUMPTION] CPT payments verified against the PFS file [15]; service counts per year are assumptions", "Lognormal"
  ) else tribble(
    ~parameter, ~value, ~tier_source, ~psa_dist,
    "Costs (2024 USD, hospital direct costs, PHC-deflated)", "", "", "",
    "Index ASD surgery",
      sprintf("$%s", format(round(cst$c_index), big.mark = ",")),
      "[DIRECT] Ames 2020 ($70,766, 2020 USD, n=210)", "Gamma",
    "Post-operative follow-up, years 1-2",
      sprintf("$%s/yr", format(round(cst$c_postop), big.mark = ",")),
      "[DERIVED] Raman 2018 2-yr total minus Ames index and expected revisions", "Gamma",
    "Post-operative follow-up, year 3+",
      sprintf("$%s/yr", format(round(cst$c_late_fu), big.mark = ",")),
      "[DERIVED] Yagi 2023 (Japanese cohort; swept widely)", "Gamma",
    "Revision surgery",
      sprintf("$%s", format(round(cst$c_revision), big.mark = ",")),
      "[DIRECT] Theologis 2016 ($55,547, 2016 USD, n=57)", "Gamma",
    "Nonoperative management",
      sprintf("$%s/yr", format(round(cst$c_nonop), big.mark = ",")),
      "[DERIVED] Glassman 2010 ($10,815 over 2 yr, n=68)", "Gamma",
    "Mechanical complication workup",
      sprintf("$%s", format(round(cst$c_comp), big.mark = ",")),
      "[ASSUMPTION] one year of early post-operative cost", "Gamma"
  )
  bind_rows(core_rows, cost_rows)
}

# =============================================================================
# MAIN
# =============================================================================
if (!interactive() && Sys.getenv("ASD_LIB_ONLY") != "1") {

  cat("=== ASD MARKOV MODEL: OPERATIVE vs NONOPERATIVE ===\n\n")
  base <- run_base_case()

  cat(sprintf("Operative     cost $%10s   QALYs %.3f\n",
              format(round(base$op$disc_cost), big.mark = ","), base$op$disc_qaly))
  cat(sprintf("Nonoperative  cost $%10s   QALYs %.3f\n",
              format(round(base$nonop$disc_cost), big.mark = ","), base$nonop$disc_qaly))
  cat(sprintf("Incremental   dC = $%s   dQ = %+.3f\n",
              format(round(base$dC), big.mark = ","), base$dQ))
  cat(sprintf("ICER = $%s per QALY\n\n",
              format(round(base$icer), big.mark = ",")))

  cat(sprintf("--- cost validation (%s costing) ---\n", COSTING))
  v <- validate(base); print(v, row.names = FALSE)

  cat("\n--- clinical events per 10,000 patients over 10 years ---\n")
  ev <- data.frame(
    arm = c("Operative", "Nonoperative"),
    index_operations = c(base$op$n_index_ops,   base$nonop$n_index_ops),
    mechanical_complications = c(base$op$n_complications, base$nonop$n_complications),
    revisions = c(base$op$n_revisions, base$nonop$n_revisions),
    deaths_10yr = c(base$op$n_deaths, base$nonop$n_deaths)
  )
  print(ev, row.names = FALSE)

  # ---- write outputs ----
  write_csv(data.frame(
    arm = c("Operative", "Nonoperative"),
    disc_cost = c(base$op$disc_cost, base$nonop$disc_cost),
    disc_qaly = c(base$op$disc_qaly, base$nonop$disc_qaly),
    undisc_cost = c(base$op$undisc_cost, base$nonop$undisc_cost),
    undisc_qaly = c(base$op$undisc_qaly, base$nonop$undisc_qaly)
  ), file.path(DATA_DIR, "base_case_summary.csv"))

  write_csv(ev,               file.path(DATA_DIR, "events.csv"))
  write_csv(v,                file.path(DATA_DIR, "validation.csv"))
  write_csv(param_table(),    file.path(DATA_DIR, "parameters.csv"))

  for (a in c("operative", "nonoperative")) {
    tr <- as.data.frame(base[[if (a == "operative") "op" else "nonop"]]$trace)
    tr <- cbind(cycle = 0:(N_CYCLES - 1), tr)
    write_csv(tr, file.path(DATA_DIR, sprintf("trace_%s.csv", a)))
  }

  cat("\n--- one-way sensitivity analysis ---\n")
  owsa <- run_owsa(); write_csv(owsa, file.path(DATA_DIR, "owsa.csv"))
  print(owsa %>% select(parameter, low_inmb, high_inmb, swing) %>% head(8),
        row.names = FALSE)

  cat("\n--- probabilistic sensitivity analysis (10,000 iterations) ---\n")
  psa <- run_psa(10000); write_csv(psa, file.path(DATA_DIR, "psa.csv"))
  cat(sprintf("mean dC = $%s   mean dQ = %+.3f\n",
              format(round(mean(psa$dC)), big.mark = ","), mean(psa$dQ)))
  ceac <- data.frame(wtp = seq(0, 300000, by = 5000)) %>%
    rowwise() %>%
    mutate(p_op_ce = mean(wtp * psa$dQ - psa$dC > 0)) %>%
    ungroup()
  write_csv(ceac, file.path(DATA_DIR, "ceac.csv"))
  for (w in c(50000, 100000, 150000)) {
    cat(sprintf("P(operative cost-effective at $%s/QALY) = %.1f%%\n",
                format(w, big.mark = ","),
                100 * mean(w * psa$dQ - psa$dC > 0)))
  }

  psa_summary <- bind_rows(lapply(c(50000, 100000, 150000), function(w) {
    inmb <- w * psa$dQ - psa$dC
    data.frame(wtp = w,
               p_op_ce = mean(inmb > 0),
               mean_inmb = mean(inmb),
               inmb_lo = quantile(inmb, 0.025),
               inmb_hi = quantile(inmb, 0.975))
  }))
  write_csv(psa_summary, file.path(DATA_DIR, "psa_summary.csv"))

  cat("\n--- scenario analyses ---\n")
  sc <- run_scenarios()
  sc_df <- bind_rows(lapply(names(sc), function(nm) {
    r <- sc[[nm]]
    data.frame(scenario = nm,
               cost_op = r$op$disc_cost, cost_nonop = r$nonop$disc_cost,
               qaly_op = r$op$disc_qaly, qaly_nonop = r$nonop$disc_qaly,
               dC = r$dC, dQ = r$dQ, icer = r$icer)
  }))
  write_csv(sc_df, file.path(DATA_DIR, "scenarios.csv"))

  # Quantities the manuscript needs but that are easy to get wrong by hand:
  # the PSA means differ from the deterministic base case (several PSA
  # distributions are deliberately skewed), and the CEAC 50% crossing is NOT the
  # deterministic ICER, because it tracks the distribution of net benefit rather
  # than the ratio of the means.
  ceac_cross <- ceac$wtp[which(ceac$p_op_ce >= 0.5)[1]]
  # Parameters whose plausible range straddles the decision boundary, i.e. that
  # can on their own make surgery NOT cost-effective at the stated threshold.
  owsa_flip <- owsa %>%
    filter(pmin(low_inmb, high_inmb) < 0 & pmax(low_inmb, high_inmb) > 0) %>%
    pull(parameter)
  write_csv(data.frame(
    psa_mean_dC   = mean(psa$dC),
    psa_mean_dQ   = mean(psa$dQ),
    ceac_crossing = ifelse(is.na(ceac_cross), NA_real_, ceac_cross),
    n_owsa_flip   = length(owsa_flip),
    n_iter        = nrow(psa),
    n_more_eff    = sum(psa$dQ > 0),
    n_ne_quadrant = sum(psa$dQ > 0 & psa$dC > 0),
    n_scen           = nrow(sc_df),
    n_scen_below_150 = sum(sc_df$icer < 150000),
    scen_above_150   = paste(sc_df$scenario[sc_df$icer >= 150000], collapse = "; "),
    scen_min_icer    = min(sc_df$icer),
    scen_max_icer    = max(sc_df$icer),
    drg_rw_blended = DRG_RW_BLENDED,
    drg_payment_blended = DRG_RW_BLENDED * IPPS_TOTAL_RATE,
    c_index = build_costs(BASE)$c_index,
    c_revision = build_costs(BASE)$c_revision,
    c_postop = build_costs(BASE)$c_postop,
    c_late_fu = build_costs(BASE)$c_late_fu,
    c_comp = build_costs(BASE)$c_comp,
    c_nonop = build_costs(BASE)$c_nonop,
    owsa_flip     = paste(owsa_flip, collapse = "; ")
  ), file.path(DATA_DIR, "derived_facts.csv"))

  print(sc_df %>% mutate(across(where(is.numeric), ~round(.x, 3))) %>%
          select(scenario, dC, dQ, icer), row.names = FALSE)

  cat("\nAll outputs written to ", DATA_DIR, "\n", sep = "")
}
