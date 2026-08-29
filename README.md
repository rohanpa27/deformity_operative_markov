# deformity_operative_markov

Markov state-transition model comparing operative with nonoperative management of
adult spinal deformity (ASD) over a 10-year horizon, from the US Medicare payer
perspective.

This repository holds the analytic code only. Every model input is drawn from
published literature or from public CMS payment files; the model contains no
patient-level data.

## Model

A 15-state cohort Markov model with 3-month cycles, a 10-year horizon, half-cycle
correction on state-residence rewards, and 3% annual discounting of costs and QALYs.
A hypothetical cohort of 10,000 patients aged 60 enters in either arm.

States are nonoperative management, a post-operative stable tunnel spanning the first
two post-operative years (one state per cycle, eight states), late post-operative
stable, early mechanical complication, late mechanical complication, revision surgery,
post-revision stable, and death. The tunnel exists so that a patient who crosses over
from nonoperative management in a later year faces the early complication hazard
appropriate to their own first post-operative year rather than the model's.

Costs are Medicare payments. The index admission is priced on the blended MS-DRG
456/457/458 relative weight (6.0785), weighted by the acuity mix observed in a
675-patient ASD cohort, times the FY2026 IPPS combined operating and capital rate.
Professional fees, routine follow-up, complication workup and nonoperative management
are built code by code from the CY2026 Physician Fee Schedule. Background mortality
comes from CDC/NCHS United States Life Tables 2023 and is held identical across arms.

## Requirements

R 4.5.3 with `dplyr`, `tidyr`, `readr`, `ggplot2`, `scales`, `forcats` and `patchwork`.

```r
install.packages(c("dplyr", "tidyr", "readr", "ggplot2", "scales", "forcats", "patchwork"))
```

## Running

From the repository root, in order:

```bash
Rscript code/01_markov_model.R   # model, PSA, one-way sensitivity, scenarios -> data/*.csv
Rscript code/02_make_figures.R   # five figures at 300 dpi -> figures/*.png
Rscript code/03_make_tables.R    # six tables -> tables/*.csv
```

Scripts 02 and 03 read the CSVs written by 01, so 01 must run first. Each script
resolves paths relative to the repository root and can be run from the root or from
`code/`. The `data/`, `figures/` and `tables/` directories are created on first run and
are not tracked here.

The probabilistic analysis is seeded (`set.seed(20260720)`), so a clean run reproduces
every reported figure exactly.

## Scripts

| Script | Purpose |
| --- | --- |
| `code/01_markov_model.R` | Single source of truth for every parameter. Builds the transition matrices, runs both arms, and performs the one-way sensitivity analysis, the 10,000-iteration probabilistic sensitivity analysis and 14 prespecified scenario and subgroup analyses. |
| `code/02_make_figures.R` | Model diagram, tornado plot, cost-effectiveness plane, acceptability curves, cohort trace and event counts. |
| `code/03_make_tables.R` | Parameter table, base case, clinical events, probabilistic results, scenarios and the cost-arithmetic reproduction. |

Every parameter in `01_markov_model.R` carries an evidence tier in its comment:
DIRECT (exact published value), DERIVED (arithmetic from cited primary sources),
EQUALIZED (held identical across arms because no arm-specific evidence exists), or
ASSUMPTION (acknowledged estimate with no direct published source).

## Configuration

Two environment variables change the analysis:

- `ASD_COSTING` selects `medicare` (default, CMS payment schedules) or `literature`
  (journal-derived hospital direct costs, retained only as a cross-check).
- `ASD_CYCLES_PER_YEAR` sets cycle length; `4` (default) gives 3-month cycles.

## License

MIT. See `LICENSE`.
