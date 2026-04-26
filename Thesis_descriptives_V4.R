# ===============================================================
# THESIS DESCRIPTIVES SCRIPT
# Retirement x Pre-retirement Volunteering x Self-Rated Health
# Scientifically clean version aligned with final thesis design
# ===============================================================

# -------------------------
# Packages
# -------------------------
library(tidyverse)
library(gtsummary)
library(gt)
library(ggplot2)
library(scales)
library(survey)
library(haven)

options(dplyr.summarise.inform = FALSE)

# -------------------------
# Paths
# -------------------------
root <- "C:/Users/Samuel/Desktop/Hertie/Thesis"
out  <- file.path(root, "output")
fig_dir <- file.path(out, "figures")
tab_dir <- file.path(out, "tables")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Helper functions
# -------------------------
save_plot_dual <- function(plot_obj, filename_base, width = 8, height = 5) {
  ggsave(
    filename = file.path(fig_dir, paste0(filename_base, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height
  )
  ggsave(
    filename = file.path(fig_dir, paste0(filename_base, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
}

save_gt_html <- function(gt_obj, filename) {
  gtsave(data = gt_obj, filename = file.path(tab_dir, filename))
}

strip_labelled <- function(df) {
  df %>%
    mutate(across(where(haven::is.labelled), haven::zap_labels))
}

# -------------------------
# Load data
# -------------------------
analytic_pw <- readRDS(file.path(out, "soep_analytic_personwaves_with_weights.rds"))
analytic_person <- readRDS(file.path(out, "soep_analytic_personlevel_with_weights.rds"))
flow <- read_csv(file.path(out, "sample_flow_builder.csv"), show_col_types = FALSE)
retention_by_sex <- read_csv(file.path(out, "retention_by_sex.csv"), show_col_types = FALSE)
timing_gap <- read_csv(file.path(out, "volunteering_gap_distribution.csv"), show_col_types = FALSE)

# ===============================================================
# 0) Basic checks
# ===============================================================
stopifnot(all(!is.na(analytic_person$vol_pre)))
stopifnot(all(!is.na(analytic_person$vol_pre_binary)))
stopifnot(all(!is.na(analytic_person$vol_pre_cat)))
stopifnot(all(!is.na(analytic_person$vol_gap_years)))

# ===============================================================
# 1) Labels and factor handling
# Keep raw ISCED categories unless verified mapping is available
# ===============================================================
analytic_person <- analytic_person %>%
  strip_labelled() %>%
  mutate(
    pid = as.character(pid),
    sex = factor(as.character(sex), levels = c("Male", "Female")),
    migback = factor(as.character(migback), levels = c("No", "Yes")),
    isced_clean = factor(isced_clean),
    vol_pre_binary = factor(
      as.character(vol_pre_binary),
      levels = c("No volunteering", "Any volunteering")
    ),
    vol_pre_cat = factor(
      as.character(vol_pre_cat),
      levels = c("Never", "Rare", "Occasional", "Frequent", "Very frequent")
    ),
    vol_gap_years = factor(as.numeric(vol_gap_years), levels = c(1, 2, 3)),
    age = as.numeric(age),
    hinc = as.numeric(hinc),
    srh_rev = as.numeric(srh_rev),
    phrf = suppressWarnings(as.numeric(phrf))
  )

analytic_pw <- analytic_pw %>%
  strip_labelled() %>%
  mutate(
    pid = as.character(pid),
    sex = factor(as.character(sex), levels = c("Male", "Female")),
    migback = factor(as.character(migback), levels = c("No", "Yes")),
    vol_pre_binary = factor(
      as.character(vol_pre_binary),
      levels = c("No volunteering", "Any volunteering")
    ),
    vol_pre_cat = factor(
      as.character(vol_pre_cat),
      levels = c("Never", "Rare", "Occasional", "Frequent", "Very frequent")
    ),
    srh_rev = as.numeric(srh_rev),
    event_time = as.integer(event_time)
  )

# ===============================================================
# 2) Sample restriction flow (appendix)
# ===============================================================
flow_gt <- flow %>%
  gt() %>%
  tab_header(
    title = md("**Table A1. Sample restriction flow**"),
    subtitle = "SOEP 2003-2023; final analytic window spans 2 years before to 5 years after retirement"
  ) %>%
  fmt_number(columns = N_individuals, decimals = 0)

save_gt_html(flow_gt, "A1_sample_flow.html")

# ===============================================================
# 3) Main descriptive table
# Main text: keep lean and aligned with manuscript
# ===============================================================
tbl_main <- analytic_person %>%
  select(
    age, sex, migback, isced_clean, hinc, srh_rev, vol_pre_binary
  ) %>%
  tbl_summary(
    missing = "no",
    type = list(
      age ~ "continuous",
      hinc ~ "continuous",
      srh_rev ~ "continuous"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      age ~ "Age",
      sex ~ "Sex",
      migback ~ "Migration background",
      isced_clean ~ "Highest observed ISCED category",
      hinc ~ "Monthly household income",
      srh_rev ~ "Self-rated health (reversed)",
      vol_pre_binary ~ "Any pre-retirement volunteering"
    )
  ) %>%
  modify_caption("**Table 1. Descriptive characteristics of the final analytic sample**")

save_gt_html(as_gt(tbl_main), "01_table1_descriptives_final_sample.html")

# ===============================================================
# 4) Weighted descriptive table (appendix)
# phrf = German population proportion weights
# This is descriptive weighting only, not selection correction
# ===============================================================
if ("phrf" %in% names(analytic_person) && any(!is.na(analytic_person$phrf))) {
  analytic_person_w <- analytic_person %>%
    filter(!is.na(phrf), phrf > 0)
  
  if (nrow(analytic_person_w) > 0) {
    svy_person <- svydesign(ids = ~1, weights = ~phrf, data = analytic_person_w)
    
    tbl_main_weighted <- tbl_svysummary(
      svy_person,
      include = c(
        age, sex, migback, isced_clean, hinc, srh_rev,
        vol_pre_binary, vol_pre_cat, vol_gap_years
      ),
      type = list(
        age ~ "continuous",
        hinc ~ "continuous",
        srh_rev ~ "continuous"
      ),
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      missing = "no"
    ) %>%
      modify_caption("**Table A2. phrf-weighted descriptive characteristics of the final analytic sample**")
    
    save_gt_html(as_gt(tbl_main_weighted), "A2_weighted_descriptives_phrf.html")
  }
}

# ===============================================================
# 5) Volunteering frequency and timing (appendix)
# ===============================================================
tbl_vol_details <- analytic_person %>%
  select(vol_pre_cat, vol_gap_years) %>%
  tbl_summary(
    missing = "no",
    statistic = all_categorical() ~ "{n} ({p}%)",
    label = list(
      vol_pre_cat ~ "Pre-retirement volunteering frequency",
      vol_gap_years ~ "Years between volunteering measure and retirement"
    )
  ) %>%
  modify_caption("**Table A3. Distribution of volunteering frequency and timing**")

save_gt_html(as_gt(tbl_vol_details), "A3_volunteering_frequency_timing.html")

# ===============================================================
# 6) Descriptives by sex (appendix)
# No p-values: descriptive, not inferential
# ===============================================================
tbl_by_sex <- analytic_person %>%
  select(
    sex, age, migback, isced_clean, hinc, srh_rev,
    vol_pre_binary, vol_pre_cat, vol_gap_years
  ) %>%
  tbl_summary(
    by = sex,
    missing = "no",
    type = list(
      age ~ "continuous",
      hinc ~ "continuous",
      srh_rev ~ "continuous"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) %>%
  modify_caption("**Table A4. Final analytic sample by sex**")

save_gt_html(as_gt(tbl_by_sex), "A4_descriptives_by_sex.html")

# ===============================================================
# 7) Retention by sex (appendix)
# ===============================================================
retention_gt <- retention_by_sex %>%
  gt() %>%
  tab_header(
    title = md("**Table A5. Retention from eligible retirement-event sample to final analytic sample, by sex**")
  ) %>%
  fmt_number(columns = c(eligible, included), decimals = 0) %>%
  fmt_percent(columns = retained_share, decimals = 1)

save_gt_html(retention_gt, "A5_retention_by_sex.html")

p_retention <- retention_by_sex %>%
  ggplot(aes(x = sex, y = retained_share)) +
  geom_col(width = 0.6) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Retention into final analytic sample by sex",
    x = NULL,
    y = "Retention rate"
  ) +
  theme_minimal(base_size = 12)

save_plot_dual(p_retention, "A1_retention_by_sex", width = 7, height = 5)

# ===============================================================
# 8) Timing of volunteering measure (appendix)
# ===============================================================
timing_gap <- timing_gap %>%
  mutate(vol_gap_years = factor(vol_gap_years, levels = c(1, 2, 3)))

timing_gap_gt <- timing_gap %>%
  gt() %>%
  tab_header(
    title = md("**Table A6. Timing of pre-retirement volunteering measure**"),
    subtitle = "Gap between observed volunteering measure and retirement year"
  ) %>%
  fmt_percent(columns = share, decimals = 1)

save_gt_html(timing_gap_gt, "A6_timing_gap_table.html")

p_gap <- timing_gap %>%
  ggplot(aes(x = vol_gap_years, y = share)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Timing of volunteering measure relative to retirement",
    x = "Years before retirement",
    y = "Share of final sample"
  ) +
  theme_minimal(base_size = 12)

save_plot_dual(p_gap, "A2_timing_gap_distribution", width = 7, height = 5)

# ===============================================================
# 9) Event-time composition (appendix)
# ===============================================================
event_dist <- analytic_pw %>%
  count(event_time, name = "n_py") %>%
  complete(event_time = -2:5, fill = list(n_py = 0)) %>%
  mutate(share_py = n_py / sum(n_py))

event_dist_gt <- event_dist %>%
  gt() %>%
  tab_header(
    title = md("**Table A7. Distribution of person-years by event time**"),
    subtitle = "Restricted event-time window: -2 to +5 years"
  ) %>%
  fmt_percent(columns = share_py, decimals = 1)

save_gt_html(event_dist_gt, "A7_event_time_distribution.html")

# ===============================================================
# 10) Distribution of pre-retirement volunteering (appendix)
# ===============================================================
vol_dist <- analytic_person %>%
  count(vol_pre_cat, name = "n") %>%
  mutate(share = n / sum(n))

vol_dist_gt <- vol_dist %>%
  gt() %>%
  tab_header(
    title = md("**Table A8. Distribution of pre-retirement volunteering**")
  ) %>%
  fmt_percent(columns = share, decimals = 1)

save_gt_html(vol_dist_gt, "A8_volunteering_distribution.html")

p_vol <- analytic_person %>%
  ggplot(aes(x = vol_pre_cat)) +
  geom_bar() +
  labs(
    title = "Pre-retirement volunteering frequency",
    x = NULL,
    y = "Number of individuals"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

save_plot_dual(p_vol, "A3_volunteering_distribution", width = 8, height = 5)

# ===============================================================
# 11) Main descriptive figure for thesis
# Descriptive only: no naive confidence intervals for panel data
# ===============================================================
plot_df <- analytic_pw %>%
  group_by(event_time, vol_pre_binary) %>%
  summarise(
    mean_srh = mean(srh_rev, na.rm = TRUE),
    n = sum(!is.na(srh_rev)),
    .groups = "drop"
  )

p_srh_vol <- ggplot(
  plot_df,
  aes(x = event_time, y = mean_srh, group = vol_pre_binary, linetype = vol_pre_binary)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = -2:5) +
  labs(
    title = "Self-rated health around retirement",
    subtitle = "Unadjusted mean profile by pre-retirement volunteering",
    x = "Years relative to retirement",
    y = "Mean self-rated health (reversed; higher = better)",
    linetype = "Pre-retirement volunteering"
  ) +
  theme_minimal(base_size = 12)

save_plot_dual(p_srh_vol, "01_mean_srh_eventtime_by_volunteering", width = 8, height = 5)

# ===============================================================
# 12) Mean SRH around retirement by sex (appendix)
# ===============================================================
plot_sex <- analytic_pw %>%
  group_by(event_time, sex) %>%
  summarise(
    mean_srh = mean(srh_rev, na.rm = TRUE),
    n = sum(!is.na(srh_rev)),
    .groups = "drop"
  )

p_srh_sex <- ggplot(
  plot_sex,
  aes(x = event_time, y = mean_srh, group = sex, linetype = sex)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = -2:5) +
  labs(
    title = "Mean self-rated health around retirement by sex",
    x = "Years relative to retirement",
    y = "Mean self-rated health (reversed; higher = better)",
    linetype = "Sex"
  ) +
  theme_minimal(base_size = 12)

save_plot_dual(p_srh_sex, "A4_mean_srh_eventtime_by_sex", width = 8, height = 5)

# ===============================================================
# 13) Retirement age distribution by sex (appendix)
# ===============================================================
retirement_age_df <- analytic_pw %>%
  filter(event_time == 0) %>%
  distinct(pid, .keep_all = TRUE) %>%
  mutate(retirement_age = age)

p_ret_age <- ggplot(retirement_age_df, aes(x = retirement_age, linetype = sex)) +
  geom_density(linewidth = 0.9, na.rm = TRUE) +
  labs(
    title = "Distribution of retirement age by sex",
    x = "Retirement age",
    y = "Density",
    linetype = "Sex"
  ) +
  theme_minimal(base_size = 12)

save_plot_dual(p_ret_age, "A5_retirement_age_by_sex", width = 8, height = 5)

# ===============================================================
# 14) Small thesis-writeup numbers file
# ===============================================================
writeup_numbers <- tibble(
  n_persons = n_distinct(analytic_person$pid),
  n_person_years = nrow(analytic_pw),
  pct_female = mean(analytic_person$sex == "Female", na.rm = TRUE),
  mean_age = mean(analytic_person$age, na.rm = TRUE),
  mean_srh = mean(analytic_person$srh_rev, na.rm = TRUE),
  pct_any_volunteering = mean(analytic_person$vol_pre_binary == "Any volunteering", na.rm = TRUE),
  median_vol_gap = median(as.numeric(as.character(analytic_person$vol_gap_years)), na.rm = TRUE)
)

write_csv(writeup_numbers, file.path(out, "10_writeup_numbers_descriptives.csv"))


# ===============================================================
# CONCEPTUAL DIAGRAM: Retirement, volunteering, and health
# Thesis-ready version
# Packages: DiagrammeR, DiagrammeRsvg, rsvg
# ===============================================================

# Packages
if (!requireNamespace("DiagrammeR", quietly = TRUE)) install.packages("DiagrammeR")
if (!requireNamespace("DiagrammeRsvg", quietly = TRUE)) install.packages("DiagrammeRsvg")
if (!requireNamespace("rsvg", quietly = TRUE)) install.packages("rsvg")
if (!requireNamespace("htmlwidgets", quietly = TRUE)) install.packages("htmlwidgets")

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(htmlwidgets)



concept_diagram <- grViz("
digraph thesis_concept {

  graph [
    layout = dot,
    rankdir = LR,
    splines = spline,
    overlap = false,
    outputorder = edgesfirst,
    bgcolor = 'white',
    pad = 0.35,
    nodesep = 0.45,
    ranksep = 0.9
  ]

  node [
    shape = box,
    style = 'rounded,filled',
    fontname = Helvetica,
    fontsize = 18,
    color = '#5B6574',
    penwidth = 1.4,
    margin = '0.18,0.12'
  ]

  edge [
    color = '#6F7B8A',
    penwidth = 1.5,
    arrowsize = 0.8
  ]

  # -------------------------------------------------------------
  # Left block: background context
  # -------------------------------------------------------------
  subgraph cluster_context {
    label = 'Context and sources of heterogeneity'
    fontsize = 20
    fontname = Helvetica
    color = '#C9D0D8'
    style = 'rounded,dashed'
    pencolor = '#C9D0D8'
    penwidth = 1.4
    bgcolor = 'transparent'

    preadv [
      label = 'Pre-existing advantage\\n(health, personality,\\nprior social integration)',
      fillcolor = '#F4F6F8'
    ]

    ses [
      label = 'Socioeconomic resources\\nand access to engagement',
      fillcolor = '#F4F6F8'
    ]

    gender [
      label = 'Gendered life-course\\npatterns',
      fillcolor = '#F4F6F8'
    ]
  }

  # -------------------------------------------------------------
  # Core transition + moderator
  # -------------------------------------------------------------
  volunteering [
    label = 'Pre-retirement\\nvolunteering',
    fillcolor = '#DCE8F4',
    color = '#56779A',
    penwidth = 1.6
  ]

  retirement [
    label = 'Retirement\\ntransition',
    fillcolor = '#EEF1F4'
  ]

  # -------------------------------------------------------------
  # Mechanisms
  # -------------------------------------------------------------
  social [
    label = 'Sustained social\\nintegration',
    fillcolor = '#F8F8F8'
  ]

  purpose [
    label = 'Preserved routine\\nand purpose',
    fillcolor = '#F8F8F8'
  ]

  loneliness [
    label = 'Reduced loneliness\\nand isolation',
    fillcolor = '#F8F8F8'
  ]

  stress [
    label = 'Less work strain\\nand more rest',
    fillcolor = '#F8F8F8'
  ]

  roleloss [
    label = 'Loss of work role\\nand identity change',
    fillcolor = '#F8F8F8'
  ]

  # -------------------------------------------------------------
  # Outcome
  # -------------------------------------------------------------
  health [
    label = 'Self-rated health\\ntrajectory around retirement',
    fillcolor = '#E7EFE8',
    color = '#607565',
    penwidth = 1.7
  ]

  # -------------------------------------------------------------
  # Structural ranks
  # -------------------------------------------------------------
  { rank = same; preadv; ses; gender }
  { rank = same; volunteering; retirement }
  { rank = same; social; purpose; loneliness; stress; roleloss }

  # -------------------------------------------------------------
  # Main causal pathways
  # -------------------------------------------------------------
  volunteering -> social
  volunteering -> purpose
  volunteering -> loneliness

  retirement -> stress
  retirement -> roleloss
  retirement -> purpose [color = '#9AA6B5', penwidth = 1.2]

  social -> health
  purpose -> health
  loneliness -> health
  stress -> health
  roleloss -> health

  # -------------------------------------------------------------
  # Background factors shaping both engagement and health
  # Dashed arrows indicate selection / confounding structure
  # -------------------------------------------------------------
  edge [style = dashed, color = '#8A95A3', penwidth = 1.3, arrowsize = 0.7]

  preadv -> volunteering
  preadv -> health

  ses -> volunteering
  ses -> health

  gender -> volunteering
  gender -> health

  # Optional curved background paths
  preadv -> health [constraint = false]
  ses -> health [constraint = false]
  gender -> health [constraint = false]
}
")

concept_diagram

svg_txt <- export_svg(concept_diagram)
writeLines(svg_txt, "conceptual_diagram.svg")

# High-resolution PNG for Word
rsvg_png(charToRaw(svg_txt), file = "conceptual_diagram.png", width = 2400, height = 1200)

# PDF version
rsvg_pdf(charToRaw(svg_txt), file = "conceptual_diagram.pdf")
# Save interactive HTML
htmlwidgets::saveWidget(
  concept_diagram,
  file = file.path(fig_dir, "00_conceptual_diagram.html"),
  selfcontained = TRUE
)


# Export to PNG
rsvg::rsvg_png(
  charToRaw(svg_txt),
  file = file.path(fig_dir, "00_conceptual_diagram.png"),
  width = 2400,
  height = 1300
)



message("Conceptual diagram saved to: ", fig_dir)

message("Descriptive script completed successfully.")
