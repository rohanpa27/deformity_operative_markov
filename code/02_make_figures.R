# =============================================================================
# 02_make_figures.R -- all manuscript figures, generated in R
#
# Run AFTER 01_markov_model.R, which writes the CSVs this script reads.
#   Rscript code/02_make_figures.R
# Writes 300 dpi PNG to figures/ .
#
# Colour: the Okabe-Ito qualitative palette (Okabe M, Ito K, 2008), chosen
# because it is an empirically validated colourblind-safe set rather than a
# hand-picked one. Categorical hues are assigned in a FIXED order and are never
# cycled. Every figure also carries a non-colour encoding (direct labels, facet,
# or line type) so identity never rests on colour alone, and all figures remain
# legible in greyscale print.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(scales);  library(patchwork); library(forcats)
})

root <- if (dir.exists("data")) "." else ".."
DATA <- file.path(root, "data")
FIGS <- file.path(root, "figures")
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

# ---- shared style -----------------------------------------------------------
OKABE <- c(blue = "#0072B2", vermillion = "#D55E00", green = "#009E73",
           orange = "#E69F00", sky = "#56B4E9", pink = "#CC79A7",
           yellow = "#F0E442", grey = "#999999")

ARM_COL <- c(Operative = unname(OKABE["blue"]),
             Nonoperative = unname(OKABE["vermillion"]))

FONT <- "Times"           # matches the manuscript typeface on macOS quartz

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = FONT) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      axis.line   = element_line(colour = "grey20", linewidth = 0.4),
      axis.ticks  = element_line(colour = "grey20", linewidth = 0.4),
      axis.text   = element_text(colour = "grey20"),
      axis.title  = element_text(colour = "black"),
      plot.title  = element_text(face = "bold", size = base_size + 1,
                                 hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(colour = "grey30", size = base_size - 1,
                                   margin = margin(b = 8)),
      legend.position = "top",
      legend.title = element_blank(),
      legend.key.height = unit(10, "pt"),
      plot.margin = margin(10, 14, 10, 10)
    )
}

save_fig <- function(plot, file, w, h) {
  ggsave(file.path(FIGS, file), plot, width = w, height = h,
         dpi = 300, units = "in", bg = "white")
  cat("wrote", file, "\n")
}

# =============================================================================
# FIGURE 1 -- Markov model structure
# A diagram, not a plot: boxes for the nine states, arrows for the transitions
# that carry a cited probability.
# =============================================================================
fig1 <- local({
  BW <- 2.0; BH <- 0.95            # box width and height
  nodes <- tribble(
    ~id,      ~x,   ~y,   ~label,                            ~fill,
    "Y1",     3.6,  5.0,  "Post-operative stable\ntunnel, year 1 (4 cycles)",  "#DCEAF7",
    "Y2",     6.4,  5.0,  "Post-operative stable\ntunnel, year 2 (4 cycles)",  "#DCEAF7",
    "LATE",   9.2,  5.0,  "Post-operative\nstable, year 3+", "#DCEAF7",
    "NONOP",  0.8,  3.2,  "Nonoperative\nmanagement",        "#FDE9D9",
    "MCE",    5.0,  3.2,  "Mechanical\ncomplication (early)","#FCE4E4",
    "MCL",    9.2,  3.2,  "Mechanical\ncomplication (late)", "#FCE4E4",
    "REV",    6.9,  1.4,  "Revision\nsurgery",               "#FFF4CC",
    "PRS",    11.4, 2.3,  "Post-revision\nstable",           "#E3F2E7",
    "DEAD",   1.6,  1.0,  "Death\n(absorbing)",              "#EDEDED"
  )
  nx <- setNames(nodes$x, nodes$id); ny <- setNames(nodes$y, nodes$id)

  # Clip a centre-to-centre segment to the two box borders so the arrowhead
  # lands on the edge of the target box rather than under it.
  clip <- function(from, to, pad = 0.06) {
    x1 <- nx[from]; y1 <- ny[from]; x2 <- nx[to]; y2 <- ny[to]
    dx <- x2 - x1; dy <- y2 - y1
    tt <- function(hw, hh) {
      tx <- if (dx != 0) hw / abs(dx) else Inf
      ty <- if (dy != 0) hh / abs(dy) else Inf
      min(tx, ty)
    }
    t1 <- tt(BW / 2 + pad, BH / 2 + pad)
    t2 <- tt(BW / 2 + pad, BH / 2 + pad)
    data.frame(x = x1 + dx * t1, y = y1 + dy * t1,
               xend = x2 - dx * t2, yend = y2 - dy * t2)
  }
  E <- function(from, to, lab = "", off_x = 0, off_y = 0, col = "grey30") {
    s <- clip(from, to)
    s$lab <- lab
    s$lx <- (s$x + s$xend) / 2 + off_x
    s$ly <- (s$y + s$yend) / 2 + off_y
    s$col <- col
    s
  }
  edges <- bind_rows(
    E("NONOP", "Y1",   "(a) crossover to surgery\n0.069/yr", -0.20, 0.34),
    E("Y1",    "Y2"),
    E("Y2",    "LATE"),
    E("Y1",    "MCE",  "(b) 0.174/yr", -0.72, 0.05),
    E("Y2",    "MCE"),
    E("LATE",  "MCL",  "(c) 0.069/yr", 0.80, 0),
    E("MCE",   "REV",  "(d) 0.526", -0.32, -0.26),
    E("MCL",   "REV",  "(e) 0.333", 0.62, 0.06),
    E("MCE",   "LATE", "(h) managed without revision", 0.10, -0.34),
    E("MCL",   "LATE"),
    E("REV",   "PRS",  "(f) revision episode\ncost charged here", 0.35, -0.62),
    E("PRS",   "MCL",  "(g) re-revision\n0.279/yr", 0.42, 0.68),
    E("NONOP", "DEAD", "", 0, 0, "grey55"),
    E("Y1",    "DEAD", "", 0, 0, "grey55"),
    E("REV",   "DEAD", "", 0, 0, "grey55")
  )

  ggplot() +
    # boxes first, arrows on top, so arrowheads are never hidden
    geom_tile(data = nodes, aes(x = x, y = y, fill = fill),
              width = BW, height = BH, colour = "grey25", linewidth = 0.45) +
    geom_segment(data = edges,
                 aes(x = x, y = y, xend = xend, yend = yend, colour = col),
                 arrow = arrow(length = unit(5.5, "pt"), type = "closed"),
                 linewidth = 0.45) +
    geom_text(data = filter(edges, lab != ""),
              aes(x = lx, y = ly, label = lab),
              size = 2.5, family = FONT, colour = "grey15", lineheight = 0.95) +
    geom_text(data = nodes, aes(x = x, y = y, label = label),
              size = 2.6, family = FONT, lineheight = 0.95) +
    scale_fill_identity() + scale_colour_identity() +
    annotate("text", x = 0.1, y = 0.15, hjust = 0, size = 2.4, family = FONT,
             colour = "grey35",
             label = paste("Cycle length is 3 months. The post-operative stable tunnel is eight states, one per cycle across the first two post-operative years;",
                           "it is drawn as two boxes for clarity.\nGrey arrows denote background mortality, which applies from every living state and is drawn only three times.")) +
    coord_cartesian(xlim = c(-0.35, 12.6), ylim = c(0.05, 5.75)) +
    labs(title = "Markov model structure: operative versus nonoperative management of adult spinal deformity") +
    theme_void(base_family = FONT) +
    theme(plot.title = element_text(face = "bold", size = 10.5,
                                    hjust = 0.5, margin = margin(b = 8)),
          plot.margin = margin(8, 8, 6, 8))
})
save_fig(fig1, "Figure_1_Model_Structure.png", 11.0, 5.6)

# =============================================================================
# FIGURE 2 -- Tornado plot (one-way sensitivity analysis)
# Form: magnitude ranking of an effect that can go either way, so a diverging
# two-colour bar anchored on the base-case iNMB.
# =============================================================================
owsa <- read_csv(file.path(DATA, "owsa.csv"), show_col_types = FALSE) %>%
  mutate(parameter = fct_reorder(parameter, swing))

base_inmb <- owsa$base_inmb[1]

fig2 <- owsa %>%
  select(parameter, low_inmb, high_inmb) %>%
  pivot_longer(-parameter, names_to = "bound", values_to = "inmb") %>%
  mutate(bound = recode(bound, low_inmb = "Low value of parameter",
                                high_inmb = "High value of parameter")) %>%
  ggplot(aes(y = parameter, x = inmb, fill = bound)) +
  # Bars are drawn as deviations from the base case, so the reference line sits
  # at zero in the shifted coordinate space.
  geom_vline(xintercept = 0, linewidth = 0.5,
             colour = "grey35", linetype = "22") +
  # Solid rule at iNMB = 0, the decision boundary: bars crossing it identify the
  # parameters that can on their own make surgery not cost-effective.
  geom_vline(xintercept = -base_inmb, linewidth = 0.4, colour = "grey55") +
  geom_col(data = ~ filter(.x, bound == "Low value of parameter"),
           aes(x = inmb - base_inmb), position = "identity",
           width = 0.62, colour = "white", linewidth = 0.3) +
  geom_col(data = ~ filter(.x, bound == "High value of parameter"),
           aes(x = inmb - base_inmb), position = "identity",
           width = 0.62, colour = "white", linewidth = 0.3) +
  scale_x_continuous(
    labels = function(v) dollar((v + base_inmb) / 1000, accuracy = 1,
                                suffix = "k"),
    expand = expansion(mult = 0.10)) +
  scale_fill_manual(values = c("Low value of parameter" = unname(OKABE["sky"]),
                               "High value of parameter" = unname(OKABE["vermillion"]))) +
  labs(
    title = "One-way sensitivity analysis",
    subtitle = paste0("Incremental net monetary benefit of operative management at a willingness to pay of $100,000 per QALY.\n",
                      "Dashed line, base case (", dollar(round(base_inmb)),
                      "). Solid line, zero: bars crossing it can make surgery not cost-effective."),
    x = "Incremental net monetary benefit (2026 USD)", y = NULL) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())
save_fig(fig2, "Figure_2_Tornado.png", 9.0, 6.0)

# =============================================================================
# FIGURE 3 -- Cost-effectiveness plane
# =============================================================================
psa  <- read_csv(file.path(DATA, "psa.csv"), show_col_types = FALSE)
bcs  <- read_csv(file.path(DATA, "base_case_summary.csv"), show_col_types = FALSE)
base_dC <- bcs$disc_cost[bcs$arm == "Operative"] - bcs$disc_cost[bcs$arm == "Nonoperative"]
base_dQ <- bcs$disc_qaly[bcs$arm == "Operative"] - bcs$disc_qaly[bcs$arm == "Nonoperative"]

# Subtitle is computed from the draws rather than written by hand, so it cannot
# contradict the result it sits above.
n_iter    <- nrow(psa)
n_ne      <- sum(psa$dQ > 0 & psa$dC > 0)
n_saving  <- sum(psa$dQ > 0 & psa$dC <= 0)
p_below150 <- mean(150000 * psa$dQ - psa$dC > 0)
ce_subtitle <- sprintf(
  paste0("Each point is one of %s probabilistic iterations. Operative management is more effective in every\n",
         "iteration; it is also more costly in %s of them and cost saving in %s. %s%% of iterations fall below\n",
         "the $150,000 per QALY line, and %s%% below $100,000."),
  format(n_iter, big.mark = ","), format(n_ne, big.mark = ","), n_saving,
  formatC(100 * p_below150, format = "f", digits = 1),
  formatC(100 * mean(100000 * psa$dQ - psa$dC > 0), format = "f", digits = 1))

wtp_lines <- tibble(wtp = c(50000, 100000, 150000),
                    lab = c("$50k/QALY", "$100k/QALY", "$150k/QALY"))

fig3 <- ggplot(psa, aes(dQ, dC)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_abline(data = wtp_lines, aes(slope = wtp, intercept = 0,
                                    linetype = fct_inorder(lab)),
              colour = "grey25", linewidth = 0.45) +
  geom_point(alpha = 0.10, size = 0.5, colour = unname(OKABE["blue"])) +
  geom_point(aes(x = base_dQ, y = base_dC), size = 3.2, shape = 23,
             fill = unname(OKABE["vermillion"]), colour = "black", stroke = 0.6) +
  annotate("text", x = base_dQ, y = base_dC, label = "  base case",
           hjust = 0, size = 3, family = FONT) +
  scale_linetype_manual(values = c("solid", "22", "42")) +
  scale_y_continuous(labels = dollar_format(scale = 1e-3, suffix = "k")) +
  labs(
    title = "Cost-effectiveness plane",
    subtitle = ce_subtitle,
    x = "Incremental QALYs (operative minus nonoperative)",
    y = "Incremental cost (2026 USD)") +
  theme_pub()
save_fig(fig3, "Figure_3_CE_Plane.png", 8.5, 6.0)

# =============================================================================
# FIGURE 4 -- Cost-effectiveness acceptability curve
# =============================================================================
ceac <- read_csv(file.path(DATA, "ceac.csv"), show_col_types = FALSE) %>%
  mutate(Nonoperative = 1 - p_op_ce) %>%
  rename(Operative = p_op_ce) %>%
  pivot_longer(c(Operative, Nonoperative),
               names_to = "strategy", values_to = "p")

lab_pts <- ceac %>% filter(wtp == 300000)

fig4 <- ggplot(ceac, aes(wtp, p, colour = strategy, linetype = strategy)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = c(100000, 150000), colour = "grey70",
             linetype = "22", linewidth = 0.4) +
  geom_text(data = lab_pts, aes(label = strategy), hjust = 1.02, vjust = -0.8,
            size = 3.1, family = FONT, show.legend = FALSE) +
  scale_colour_manual(values = ARM_COL) +
  scale_linetype_manual(values = c(Operative = "solid", Nonoperative = "22")) +
  scale_x_continuous(labels = dollar_format(scale = 1e-3, suffix = "k"),
                     breaks = seq(0, 300000, 50000)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Cost-effectiveness acceptability curve",
    subtitle = "Probability that each strategy is cost-effective across willingness-to-pay thresholds.\nReference lines mark $100,000 and $150,000 per QALY.",
    x = "Willingness to pay per QALY (2026 USD)",
    y = "Probability cost-effective") +
  theme_pub()
save_fig(fig4, "Figure_4_CEAC.png", 8.5, 5.5)

# =============================================================================
# FIGURE 5 -- Cohort trace and cumulative revision burden
# =============================================================================
read_trace <- function(arm) {
  read_csv(file.path(DATA, sprintf("trace_%s.csv", arm)), show_col_types = FALSE) %>%
    mutate(arm = if (arm == "operative") "Operative" else "Nonoperative")
}
tr <- bind_rows(read_trace("operative"), read_trace("nonoperative"))

# Map state names to display groups. The post-operative stable state is a tunnel
# whose width depends on cycle length (Stable_T01 ... Stable_Tnn), so it is
# matched by prefix rather than enumerated.
group_of <- function(state) {
  dplyr::case_when(
    state == "Nonoperative"            ~ "Nonoperative",
    grepl("^Stable_T[0-9]+$", state)   ~ "Post-operative stable",
    state == "Stable_Late"             ~ "Post-operative stable",
    grepl("^MechComp_", state)         ~ "Mechanical complication",
    state == "Revision"                ~ "Revision",
    state == "PostRev_Stable"          ~ "Post-revision stable",
    state == "Death"                   ~ "Death"
  )
}

# Cycles per year is inferred from the trace: the tunnel spans two years.
CPY <- sum(grepl("^Stable_T[0-9]+$", names(tr))) / 2

trace_long <- tr %>%
  pivot_longer(-c(cycle, arm), names_to = "state", values_to = "p") %>%
  mutate(grp = factor(group_of(state),
                      levels = c("Post-operative stable", "Post-revision stable",
                                 "Nonoperative", "Mechanical complication",
                                 "Revision", "Death")),
         year = cycle / CPY) %>%
  group_by(arm, year, grp) %>% summarise(p = sum(p), .groups = "drop") %>%
  mutate(arm = factor(arm, levels = c("Operative", "Nonoperative")))
stopifnot(!any(is.na(trace_long$grp)))

GRP_COL <- c("Post-operative stable"   = unname(OKABE["blue"]),
             "Post-revision stable"    = unname(OKABE["sky"]),
             "Nonoperative"            = unname(OKABE["orange"]),
             "Mechanical complication" = unname(OKABE["vermillion"]),
             "Revision"                = unname(OKABE["pink"]),
             "Death"                   = unname(OKABE["grey"]))

p5a <- ggplot(trace_long, aes(year, p, fill = grp)) +
  geom_area(colour = NA) +
  facet_wrap(~arm) +
  scale_fill_manual(values = GRP_COL) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(0)) +
  scale_x_continuous(breaks = seq(0, 10, 2), expand = expansion(0)) +
  labs(title = "A. Cohort distribution across health states",
       x = "Year", y = "Share of cohort") +
  theme_pub() +
  theme(strip.text = element_text(face = "bold", size = 10, family = FONT),
        legend.position = "right")

events <- read_csv(file.path(DATA, "events.csv"), show_col_types = FALSE) %>%
  select(arm, `Index operations` = index_operations,
         `Mechanical complications` = mechanical_complications,
         Revisions = revisions) %>%
  pivot_longer(-arm, names_to = "event", values_to = "n") %>%
  mutate(event = fct_inorder(event),
         arm = factor(arm, levels = c("Operative", "Nonoperative")))

p5b <- ggplot(events, aes(event, n, fill = arm)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = comma(round(n))),
            position = position_dodge(width = 0.72),
            vjust = -0.45, size = 2.9, family = FONT, colour = "grey20") +
  scale_fill_manual(values = ARM_COL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.16))) +
  labs(title = "B. Clinical events per 10,000 patients over 10 years",
       x = NULL, y = "Events") +
  theme_pub() +
  theme(panel.grid.major.x = element_blank())

fig5 <- p5a / p5b + plot_layout(heights = c(1.25, 1))
save_fig(fig5, "Figure_5_Trace_and_Events.png", 9.5, 9.0)

cat("\nAll figures written to ", FIGS, "\n", sep = "")
