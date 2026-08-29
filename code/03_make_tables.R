# =============================================================================
# 03_make_tables.R -- manuscript tables, written as CSV to tables/
#
# Run AFTER 01_markov_model.R.  Rscript code/03_make_tables.R
# 04_make_manuscript.py reads these CSVs and renders them into the .docx.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(scales)
})

root <- if (dir.exists("data")) "." else ".."
DATA <- file.path(root, "data")
TABS <- file.path(root, "tables")
dir.create(TABS, showWarnings = FALSE, recursive = TRUE)

usd  <- function(x) dollar(round(x), accuracy = 1)
q3   <- function(x) sprintf("%.3f", x)

# ---- Table 1: model input parameters ---------------------------------------
read_csv(file.path(DATA, "parameters.csv"), show_col_types = FALSE) %>%
  rename(Parameter = parameter, `Base value` = value,
         `Evidence tier and source` = tier_source,
         `PSA distribution` = psa_dist) %>%
  write_csv(file.path(TABS, "Table_1_Parameters.csv"))

# ---- Table 2: base-case cost-effectiveness ---------------------------------
bcs <- read_csv(file.path(DATA, "base_case_summary.csv"), show_col_types = FALSE)
op  <- bcs %>% filter(arm == "Operative")
no  <- bcs %>% filter(arm == "Nonoperative")
dC  <- op$disc_cost - no$disc_cost
dQ  <- op$disc_qaly - no$disc_qaly

tibble(
  Strategy = c("Nonoperative management", "Operative management",
               "Undiscounted total (nonoperative)", "Undiscounted total (operative)"),
  `Discounted cost (USD)` = c(usd(no$disc_cost), usd(op$disc_cost),
                                   usd(no$undisc_cost), usd(op$undisc_cost)),
  `Discounted QALYs` = c(q3(no$disc_qaly), q3(op$disc_qaly),
                         q3(no$undisc_qaly), q3(op$undisc_qaly)),
  `Incremental cost` = c("Reference", usd(dC), "", ""),
  `Incremental QALYs` = c("Reference", sprintf("+%s", q3(dQ)), "", ""),
  `ICER (USD per QALY)` = c("Reference", usd(dC / dQ), "", "")
) %>% write_csv(file.path(TABS, "Table_2_Base_Case.csv"))

# ---- Table 3: clinical events ----------------------------------------------
read_csv(file.path(DATA, "events.csv"), show_col_types = FALSE) %>%
  transmute(
    `Outcome (per 10,000 patients over 10 years)` = arm,
    `Index operations`         = comma(round(index_operations)),
    `Mechanical complications` = comma(round(mechanical_complications)),
    `Revision operations`      = comma(round(revisions)),
    `Deaths at 10 years`       = comma(round(deaths_10yr))
  ) %>%
  pivot_longer(-1, names_to = "Outcome", values_to = "v") %>%
  pivot_wider(names_from = 1, values_from = v) %>%
  write_csv(file.path(TABS, "Table_3_Clinical_Events.csv"))

# ---- Table 4: probabilistic sensitivity analysis ---------------------------
read_csv(file.path(DATA, "psa_summary.csv"), show_col_types = FALSE) %>%
  transmute(
    `Willingness-to-pay threshold` = paste0(usd(wtp), " per QALY"),
    `P(operative cost-effective)`    = percent(p_op_ce, accuracy = 0.1),
    `P(nonoperative cost-effective)` = percent(1 - p_op_ce, accuracy = 0.1),
    `Mean incremental NMB (USD)` = usd(mean_inmb),
    `95% credible interval` = paste0(usd(inmb_lo), " to ", usd(inmb_hi))
  ) %>% write_csv(file.path(TABS, "Table_4_PSA.csv"))

# ---- Table 5: scenario and subgroup analyses -------------------------------
read_csv(file.path(DATA, "scenarios.csv"), show_col_types = FALSE) %>%
  transmute(
    `Scenario or subgroup` = scenario,
    `Incremental cost (USD)` = usd(dC),
    `Incremental QALYs` = sprintf("+%s", q3(dQ)),
    `ICER (USD per QALY)` = usd(icer),
    `Cost-effective at $150,000 per QALY` =
      ifelse(icer < 150000, "Yes", "No")
  ) %>% write_csv(file.path(TABS, "Table_5_Scenarios.csv"))

# ---- Table 6: model validation ---------------------------------------------
read_csv(file.path(DATA, "validation.csv"), show_col_types = FALSE) %>%
  transmute(Quantity = quantity,
            # the first row is a relative weight, the rest are dollars
            `Model value` = ifelse(model < 100, sprintf("%.4f", model), usd(model)),
            `Recomputed from CMS files` = ifelse(recomputed < 100,
                                                 sprintf("%.4f", recomputed),
                                                 usd(recomputed)),
            `Ratio` = sprintf("%.3f", ratio),
            Basis = note) %>%
  write_csv(file.path(TABS, "Table_6_Validation.csv"))

cat("Tables written to ", TABS, "\n", sep = "")
for (f in list.files(TABS)) cat("  ", f, "\n", sep = "")
