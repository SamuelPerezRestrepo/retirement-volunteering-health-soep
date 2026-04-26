# ===============================================================
# SOEP THESIS DATASET BUILDER
# Calendar-year event time, explicit pre-retirement volunteering
# Author: Samuel Pérez Restrepo
# Purpose: Build panel for Retirement x Pre-retirement Volunteering x SRH
# ===============================================================

# -------------------------
# Packages
# -------------------------
library(tidyverse)
library(haven)
library(gtsummary)
library(gt)
library(broom)
library(sandwich)
library(lmtest)
library(ggplot2)
library(survey)

options(dplyr.summarise.inform = FALSE)

# -------------------------
# Paths
# -------------------------
root <- "C:/Users/Samuel/Desktop/Hertie/Thesis"
out  <- file.path(root, "output")
fig_dir <- file.path(out, "figures")
tab_dir <- file.path(out, "tables")

dir.create(out, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Helper: convert SOEP negative codes to NA
# -------------------------
clean_vars <- function(df) {
  df %>%
    mutate(across(where(is.numeric), ~ ifelse(. < 0, NA, .)))
}

# -------------------------
# Load raw files
# -------------------------
ppathl <- readRDS(file.path(root, "ppathl.rds"))
pl     <- readRDS(file.path(root, "pl.rds"))
hl     <- readRDS(file.path(root, "hl.rds"))
pgen   <- readRDS(file.path(root, "pgen.rds"))

if (file.exists(file.path(root, "phrf.rds"))) {
  phrf_extra <- readRDS(file.path(root, "phrf.rds"))
} else {
  phrf_extra <- NULL
}

# ===============================================================
# 1) Backbone: person-year panel skeleton
# ===============================================================
backbone <- ppathl %>%
  select(pid, syear, gebjahr, sex, migback, netto, any_of("phrf")) %>%
  clean_vars() %>%
  mutate(
    age = syear - gebjahr,
    sex = case_when(
      sex == 1 ~ "Male",
      sex == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    migback = case_when(
      migback == 1 ~ "No",
      migback %in% c(2, 3) ~ "Yes",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(
    pl %>% select(pid, syear, hid),
    by = c("pid", "syear")
  )

if (!"phrf" %in% names(backbone) && !is.null(phrf_extra)) {
  backbone <- backbone %>%
    left_join(
      phrf_extra %>% select(pid, syear, phrf),
      by = c("pid", "syear")
    )
}

if ("phrf" %in% names(backbone)) {
  backbone <- backbone %>%
    mutate(phrf = as.numeric(phrf))
}

# ===============================================================
# 2) Civic engagement: volunteering only
# Harmonize pli0096_* into a single volunteering measure
# Recoding:
#   4 = very frequent
#   3 = frequent
#   2 = occasional
#   1 = rare
#   0 = never
# ===============================================================
ce_cols <- c("pli0096_h", "pli0096_v1", "pli0096_v2", "pli0096_v3")

ce_df <- pl %>%
  select(pid, syear, any_of(ce_cols)) %>%
  clean_vars() %>%
  rowwise() %>%
  mutate(
    ce_voluntary_raw = coalesce(!!!syms(ce_cols)),
    ce_voluntary = case_when(
      ce_voluntary_raw == 1 ~ 4,
      ce_voluntary_raw == 2 ~ 3,
      ce_voluntary_raw == 3 ~ 2,
      ce_voluntary_raw == 4 ~ 1,
      ce_voluntary_raw == 5 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup() %>%
  select(pid, syear, ce_voluntary)

# ===============================================================
# 3) Health: self-rated health
# SOEP original coding: 1 = very good ... 5 = bad
# Reverse-coded later as srh_rev = 6 - srh
# ===============================================================
health_df <- pl %>%
  select(pid, syear, ple0008) %>%
  clean_vars() %>%
  transmute(
    pid,
    syear,
    srh = ple0008
  )

# ===============================================================
# 4) Income: household net monthly income
# ===============================================================
income_df <- hl %>%
  select(hid, syear, hlc0005_h, hlc0005_v1, hlc0005_v2) %>%
  clean_vars() %>%
  mutate(
    hinc = coalesce(hlc0005_h, hlc0005_v1, hlc0005_v2)
  ) %>%
  select(hid, syear, hinc)

# ===============================================================
# 5) Retirement: absorbing status and first observed transition
# ===============================================================
ret_df <- pl %>%
  select(
    pid, syear,
    retired_now   = plc0311,
    retired_prev1 = plb0282_h,
    retired_prev2 = plb0282_v1,
    retired_prev3 = pab0005
  ) %>%
  clean_vars() %>%
  mutate(
    retired_status_raw = case_when(
      retired_now   == 1 ~ 1L,
      retired_prev1 == 1 ~ 1L,
      retired_prev2 == 1 ~ 1L,
      retired_prev3 == 1 ~ 1L,
      TRUE ~ 0L
    )
  ) %>%
  arrange(pid, syear) %>%
  group_by(pid) %>%
  mutate(
    retired_status = cummax(retired_status_raw)
  ) %>%
  ungroup() %>%
  select(pid, syear, retired_status_raw, retired_status)

# ===============================================================
# 6) Education: yearly ISCED from pgen
# Official SOEP education variable: pgisced97 in pgen
# We will later define person-level education as the highest observed
# ISCED up to (and including) the retirement year.
# ===============================================================
edu_yearly <- pgen %>%
  select(pid, syear, pgisced97) %>%
  clean_vars() %>%
  transmute(
    pid,
    syear,
    isced_year = as.numeric(pgisced97)
  ) %>%
  mutate(
    isced_year = if_else(isced_year %in% 1:6, isced_year, NA_real_)
  )

# ===============================================================
# 7) Merge core panel
# ===============================================================
panel0 <- backbone %>%
  left_join(ce_df,      by = c("pid", "syear")) %>%
  left_join(health_df,  by = c("pid", "syear")) %>%
  left_join(income_df,  by = c("hid", "syear")) %>%
  left_join(ret_df,     by = c("pid", "syear")) %>%
  left_join(edu_yearly, by = c("pid", "syear")) %>%
  filter(syear >= 2003, syear <= 2023) %>%
  arrange(pid, syear)

# ===============================================================
# 8) Stepwise sample restriction
# IMPORTANT:
# Define observed retirement transition AFTER age restriction
# so the event is observed within the analytic age frame
# ===============================================================
step0 <- panel0
n0 <- n_distinct(step0$pid)

step1 <- step0 %>%
  filter(age >= 50, age <= 70) %>%
  arrange(pid, syear) %>%
  group_by(pid) %>%
  mutate(
    prev_retired = lag(retired_status_raw),
    retirement_event = if_else(
      !is.na(prev_retired) & prev_retired == 0L & retired_status_raw == 1L,
      1L, 0L
    ),
    retire_year = if_else(
      any(retirement_event == 1L),
      syear[which(retirement_event == 1L)[1]],
      NA_integer_
    ),
    event_time = if_else(!is.na(retire_year), syear - retire_year, NA_integer_),
    retired_status = cummax(retired_status_raw)
  ) %>%
  ungroup()
n1 <- n_distinct(step1$pid)

step2 <- step1 %>%
  filter(!is.na(retire_year)) %>%
  group_by(pid) %>%
  filter(any(event_time == 0)) %>%
  ungroup()
n2 <- n_distinct(step2$pid)

# ===============================================================
# 8b) Education: highest observed ISCED up to retirement
# This avoids using post-retirement education information.
# ===============================================================
edu_person <- step2 %>%
  filter(!is.na(isced_year), syear <= retire_year) %>%
  group_by(pid) %>%
  summarise(
    isced_clean = max(isced_year, na.rm = TRUE),
    .groups = "drop"
  )

step2 <- step2 %>%
  left_join(edu_person, by = "pid")

# ===============================================================
# 8c) Pre-retirement volunteering moderator
# Most recent non-missing volunteering observation in the
# 1- to 3-year window before retirement
# ===============================================================
vol_pre_df <- step2 %>%
  filter(!is.na(ce_voluntary)) %>%
  filter(syear < retire_year, syear >= retire_year - 3) %>%
  arrange(pid, desc(syear)) %>%
  group_by(pid) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    pid,
    vol_pre = ce_voluntary,
    vol_pre_year = syear,
    vol_gap_years = retire_year - vol_pre_year
  )

step3 <- step2 %>%
  filter(event_time >= -2, event_time <= 5) %>%
  group_by(pid) %>%
  filter(any(event_time == 0)) %>%
  ungroup()
n3 <- n_distinct(step3$pid)

# Save eligible retirement-event sample for diagnostics
saveRDS(step3, file.path(out, "soep_step3_pre_final_restriction.rds"))

# Final inclusion IDs: must have valid pre-retirement volunteering
step4_ids <- step3 %>%
  distinct(pid) %>%
  inner_join(vol_pre_df, by = "pid") %>%
  filter(!is.na(vol_pre))

step4 <- step3 %>%
  inner_join(step4_ids, by = "pid")

n4 <- n_distinct(step4$pid)

# Safety checks
stopifnot(all(!is.na(step4$vol_pre)))
stopifnot(all(!is.na(step4$vol_gap_years)))
stopifnot(all(
  step4 %>%
    distinct(pid, retire_year) %>%
    count(pid) %>%
    pull(n) == 1
))

# ===============================================================
# 9) Flow table
# ===============================================================
flow <- tibble(
  Step = c(
    "Initial SOEP sample (2003-2023)",
    "Age restriction: 50-70",
    "Observed retirement transition (first event)",
    "Event-time window: -2 to +5 years",
    "Valid pre-retirement volunteering measure",
    "Final analytical sample"
  ),
  N_individuals = c(n0, n1, n2, n3, n4, n4)
)

flow %>%
  gt() %>%
  tab_header(title = "Sample Restriction Flow - SOEP 2003-2023") %>%
  fmt_number(columns = N_individuals, decimals = 0) %>%
  gtsave(file.path(tab_dir, "sample_flow_builder.png"))

write_csv(flow, file.path(out, "sample_flow_builder.csv"))

# ===============================================================
# 10) Build final person-year panel
# ===============================================================
analytic_pw <- step4 %>%
  arrange(pid, syear) %>%
  mutate(
    srh_rev = if_else(!is.na(srh), 6 - srh, NA_real_),
    vol_pre_binary = case_when(
      vol_pre == 0 ~ "No volunteering",
      vol_pre >= 1 ~ "Any volunteering",
      TRUE ~ NA_character_
    ),
    vol_pre_cat = case_when(
      vol_pre == 0 ~ "Never",
      vol_pre == 1 ~ "Rare",
      vol_pre == 2 ~ "Occasional",
      vol_pre == 3 ~ "Frequent",
      vol_pre == 4 ~ "Very frequent",
      TRUE ~ NA_character_
    )
  )

# Safety checks after coding
stopifnot(all(!is.na(analytic_pw$vol_pre_binary)))
stopifnot(all(!is.na(analytic_pw$vol_pre_cat)))

# ===============================================================
# 11) Person-level descriptive snapshot
# Prefer event_time = -1; if unavailable, fallback to -2
# This gives a consistent pre-retirement descriptive baseline
# ===============================================================
analytic_person <- analytic_pw %>%
  group_by(pid) %>%
  mutate(
    snapshot_priority = case_when(
      event_time == -1 ~ 1L,
      event_time == -2 ~ 2L,
      TRUE ~ 99L
    )
  ) %>%
  arrange(pid, snapshot_priority, syear) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    pid, syear, event_time, age, sex, migback, isced_clean,
    hinc, srh, srh_rev, retire_year,
    vol_pre, vol_pre_year, vol_gap_years,
    vol_pre_binary, vol_pre_cat, phrf
  )

# Safety checks
stopifnot(n_distinct(analytic_person$pid) == n4)
stopifnot(all(!is.na(analytic_person$vol_pre)))
stopifnot(all(!is.na(analytic_person$vol_gap_years)))

# ===============================================================
# 12) Diagnostics: included vs excluded from eligible step3 sample
# Snapshot rule:
# prefer event_time = -1, then -2, then 0 if needed
# ===============================================================
step3_snapshot <- step3 %>%
  group_by(pid) %>%
  mutate(
    snapshot_priority = case_when(
      event_time == -1 ~ 1L,
      event_time == -2 ~ 2L,
      event_time == 0  ~ 3L,
      TRUE ~ 99L
    )
  ) %>%
  arrange(pid, snapshot_priority, syear) %>%
  slice(1) %>%
  ungroup() %>%
  select(pid, age, sex, migback, isced_clean, hinc, srh, event_time, phrf)

saveRDS(step3_snapshot, file.path(out, "soep_step3_snapshot.rds"))

included_pids <- analytic_person$pid

excluded_persons <- step3_snapshot %>%
  filter(!pid %in% included_pids) %>%
  mutate(group = "Excluded")

included_persons <- step3_snapshot %>%
  filter(pid %in% included_pids) %>%
  mutate(group = "Included")

compare_df <- bind_rows(included_persons, excluded_persons) %>%
  transmute(
    group = group,
    age = as.numeric(age),
    sex = as.character(sex),
    migback = as.character(migback),
    isced_clean = as.integer(isced_clean),
    hinc = as.numeric(hinc),
    srh = as.numeric(srh)
  )

compare_tbl <- compare_df %>%
  tbl_summary(
    by = group,
    missing = "ifany",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      age ~ "Age",
      sex ~ "Sex",
      migback ~ "Migration background",
      isced_clean ~ "ISCED-1997",
      hinc ~ "Monthly household income",
      srh ~ "Self-rated health"
    )
  ) %>%
  add_p() %>%
  modify_header(label = "**Variable**") %>%
  modify_caption("Included vs excluded individuals among retirement-event cases in the event-time window")

gtsave(as_gt(compare_tbl), file.path(tab_dir, "included_vs_excluded.png"))

# ===============================================================
# 13) Retention by sex
# ===============================================================
retention_by_sex <- step3 %>%
  distinct(pid, sex, .keep_all = FALSE) %>%
  mutate(included = if_else(pid %in% included_pids, 1L, 0L)) %>%
  group_by(sex) %>%
  summarise(
    eligible = n(),
    included = sum(included),
    retained_share = included / eligible,
    .groups = "drop"
  )

write_csv(retention_by_sex, file.path(out, "retention_by_sex.csv"))

# ===============================================================
# 14) Distribution of timing gap
# ===============================================================
timing_gap <- analytic_person %>%
  count(vol_gap_years, name = "n") %>%
  mutate(share = n / sum(n))

write_csv(timing_gap, file.path(out, "volunteering_gap_distribution.csv"))

# ===============================================================
# 15) Weighted descriptives using phrf if available
# ===============================================================
if ("phrf" %in% names(analytic_person) && any(!is.na(analytic_person$phrf))) {
  analytic_person_w <- analytic_person %>%
    filter(!is.na(phrf), phrf > 0)
  
  svy_person <- svydesign(ids = ~1, weights = ~phrf, data = analytic_person_w)
  
  w_age <- svymean(~age, design = svy_person, na.rm = TRUE)
  w_srh <- svymean(~srh_rev, design = svy_person, na.rm = TRUE)
  
  weighted_means <- tibble(
    variable = c("age", "srh_rev"),
    weighted_mean = c(as.numeric(coef(w_age)), as.numeric(coef(w_srh))),
    se = c(as.numeric(SE(w_age)), as.numeric(SE(w_srh)))
  )
  
  write_csv(weighted_means, file.path(out, "weighted_means.csv"))
  
  w_vol <- svytable(~vol_pre_cat, design = svy_person)
  w_vol_df <- as.data.frame(w_vol, stringsAsFactors = FALSE)
  names(w_vol_df) <- c("vol_pre_cat", "weighted_count")
  
  w_vol_df <- w_vol_df %>%
    mutate(
      vol_pre_cat = as.character(vol_pre_cat),
      weighted_count = as.numeric(weighted_count),
      weighted_prop = weighted_count / sum(weighted_count)
    )
  
  write_csv(w_vol_df, file.path(out, "weighted_volunteering_distribution.csv"))
}

# ===============================================================
# 16) Save final outputs
# ===============================================================
saveRDS(analytic_pw, file.path(out, "soep_analytic_personwaves_with_weights.rds"))
saveRDS(analytic_person, file.path(out, "soep_analytic_personlevel_with_weights.rds"))

# Optional legacy filenames
saveRDS(analytic_pw, file.path(out, "soep_retirement_voluntary_srh_personwaves_2003_2023.rds"))
saveRDS(analytic_person, file.path(out, "soep_retirement_voluntary_srh_personlevel_2003_2023.rds"))

message("Builder finished. Final persons: ", n4, " ; person-years: ", nrow(analytic_pw))