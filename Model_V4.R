# ===============================================================
# THESIS EMPIRICAL STRATEGY
# Retirement x Pre-retirement Volunteering x Self-Rated Health
# Final polished version for thesis production
# ===============================================================

# Main strategy:
# 1) Pooled OLS with event-time interactions, clustered SE by pid
# 2) Event-time model with individual fixed effects
# 3) Robustness: pooled ordered logit
# 4) Supplementary heterogeneity: gender and education-based SES
# 5) Diagnostics: multicollinearity, rank, heteroskedasticity,
#    cell counts, weight distribution, sample structure, optional
#    panel serial-correlation check
# 6) Appendix robustness: survey-year fixed effects
#
# Data requirement:
# - output/soep_analytic_personwaves_with_weights.rds
# - output/soep_step3_snapshot.rds
# ===============================================================

# -------------------------
# Packages
# -------------------------
library(tidyverse)
library(broom)
library(ggplot2)
library(fixest)
library(gt)
library(modelsummary)
library(ordinal)
library(haven)
library(lmtest)

options(dplyr.summarise.inform = FALSE)

# -------------------------
# Paths
# -------------------------
root <- "C:/Users/Samuel/Desktop/Hertie/Thesis"
out  <- file.path(root, "output")
fig_dir <- file.path(out, "figures")
tab_dir <- file.path(out, "tables")
mod_dir <- file.path(out, "models")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Load data
# -------------------------
analytic_pw <- readRDS(file.path(out, "soep_analytic_personwaves_with_weights.rds"))
step3_snapshot <- readRDS(file.path(out, "soep_step3_snapshot.rds"))

# -------------------------
# Sanity checks: builder alignment
# -------------------------
required_cols <- c(
  "pid", "syear", "event_time", "srh", "srh_rev",
  "vol_pre", "vol_pre_year", "vol_gap_years",
  "vol_pre_binary", "vol_pre_cat",
  "age", "sex", "migback", "isced_clean", "hinc", "phrf"
)

missing_cols <- setdiff(required_cols, names(analytic_pw))
if (length(missing_cols) > 0) {
  stop(
    "The analysis file is missing required builder variables: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ===============================================================
# HELPERS
# ===============================================================

strip_labelled <- function(df) {
  df %>%
    mutate(across(where(haven::is.labelled), haven::zap_labels))
}

mode_value <- function(x) {
  ux <- na.omit(x)
  if (length(ux) == 0) return(NA_character_)
  names(sort(table(ux), decreasing = TRUE))[1]
}

save_plot_dual <- function(plot_obj, filename_base, width = 8, height = 5) {
  pdf_path <- file.path(fig_dir, paste0(filename_base, ".pdf"))
  png_path <- file.path(fig_dir, paste0(filename_base, ".png"))
  
  ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
  
  tryCatch(
    {
      ggsave(
        filename = pdf_path,
        plot = plot_obj,
        width = width,
        height = height,
        device = grDevices::pdf
      )
    },
    error = function(e) {
      message("PDF export failed for ", filename_base, ": ", e$message)
      message("PNG version was still saved successfully.")
    }
  )
}

save_gt_html <- function(gt_obj, filename) {
  gt::gtsave(data = gt_obj, filename = file.path(tab_dir, filename))
}

save_df_html <- function(df, title, file, notes = NULL) {
  tab <- modelsummary::datasummary_df(
    df,
    output = "gt"
  ) |>
    gt::tab_header(title = title)
  
  if (!is.null(notes)) {
    tab <- tab |> gt::tab_source_note(source_note = notes)
  }
  
  gt::gtsave(tab, file)
}

save_model_html <- function(models,
                            title,
                            file,
                            notes = NULL,
                            estimate = "{estimate}{stars}",
                            statistic = "({std.error})",
                            stars = TRUE,
                            coef_omit = NULL,
                            gof_omit = "AIC|BIC|Log.Lik|Adj|RMSE|Std.Errors|F|FE") {
  tab <- modelsummary::modelsummary(
    models,
    output = "gt",
    estimate = estimate,
    statistic = statistic,
    stars = stars,
    coef_omit = coef_omit,
    gof_omit = gof_omit
  ) |>
    gt::tab_header(title = title)
  
  if (!is.null(notes)) {
    tab <- tab |> gt::tab_source_note(source_note = notes)
  }
  
  gt::gtsave(tab, file)
}

thesis_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey20"),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

vol_colors <- c(
  "No volunteering" = "#1A1A1A",
  "Any volunteering" = "#6F6F6F"
)

vol_linetypes <- c(
  "No volunteering" = "solid",
  "Any volunteering" = "dashed"
)

compute_vif_from_lm <- function(model) {
  mm <- model.matrix(model)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  
  if (ncol(mm) < 2) {
    return(tibble(term = colnames(mm), vif = NA_real_))
  }
  
  vif_vals <- sapply(seq_len(ncol(mm)), function(j) {
    y_j <- mm[, j]
    x_j <- mm[, -j, drop = FALSE]
    fit_j <- lm(y_j ~ x_j)
    r2_j <- summary(fit_j)$r.squared
    1 / (1 - r2_j)
  })
  
  tibble(
    term = colnames(mm),
    vif = as.numeric(vif_vals)
  ) %>%
    mutate(
      flag = case_when(
        vif >= 10 ~ "High (>=10)",
        vif >= 5 ~ "Moderate (>=5)",
        TRUE ~ "Low"
      )
    ) %>%
    arrange(desc(vif))
}

extract_wald_row <- function(wobj, label) {
  vals <- unlist(wobj)
  nms <- names(vals)
  
  stat <- if ("stat" %in% nms) vals["stat"] else vals[1]
  pval <- if ("p" %in% nms) vals["p"] else vals[length(vals)]
  df1  <- if ("df1" %in% nms) vals["df1"] else NA_real_
  df2  <- if ("df2" %in% nms) vals["df2"] else NA_real_
  
  tibble(
    test = label,
    statistic = as.numeric(stat),
    df1 = as.numeric(df1),
    df2 = as.numeric(df2),
    p_value = as.numeric(pval)
  )
}

# ===============================================================
# 0) PREPARE ANALYSIS DATA
# ===============================================================

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
    event_time = as.integer(event_time),
    srh = as.numeric(srh),
    srh_rev = as.numeric(srh_rev),
    age = as.numeric(age),
    hinc = as.numeric(hinc),
    hinc_1k = as.numeric(hinc) / 1000,
    phrf = suppressWarnings(as.numeric(phrf)),
    isced_num = suppressWarnings(as.numeric(as.character(isced_clean)))
  )

step3_snapshot <- step3_snapshot %>%
  strip_labelled() %>%
  mutate(
    pid = as.character(pid),
    sex = as.character(sex),
    migback = as.character(migback),
    age = as.numeric(age),
    hinc = as.numeric(hinc),
    srh = as.numeric(srh),
    isced_clean = suppressWarnings(as.numeric(as.character(isced_clean))),
    phrf = suppressWarnings(as.numeric(phrf))
  )

analysis_pw <- analytic_pw %>%
  filter(event_time >= -2, event_time <= 5) %>%
  mutate(
    event_time_f = factor(event_time, levels = -2:5),
    event_time_f = relevel(event_time_f, ref = "-1"),
    
    vol_any = if_else(vol_pre_binary == "Any volunteering", 1L, 0L),
    
    isced_clean_f = factor(
      if_else(is.na(isced_num), "Missing", as.character(isced_num))
    ),
    
    ses_edu = case_when(
      !is.na(isced_num) & isced_num <= 4 ~ "Lower / medium education",
      !is.na(isced_num) & isced_num >= 5 ~ "Higher education",
      TRUE ~ NA_character_
    ),
    ses_edu = factor(
      ses_edu,
      levels = c("Lower / medium education", "Higher education")
    ),
    
    age_c = age - mean(age, na.rm = TRUE),
    hinc_1k_c = hinc_1k - mean(hinc_1k, na.rm = TRUE),
    
    srh_rev_ord = case_when(
      srh == 5 ~ "Bad",
      srh == 4 ~ "Poor",
      srh == 3 ~ "Satisfactory",
      srh == 2 ~ "Good",
      srh == 1 ~ "Very good",
      TRUE ~ NA_character_
    ),
    srh_rev_ord = ordered(
      srh_rev_ord,
      levels = c("Bad", "Poor", "Satisfactory", "Good", "Very good")
    )
  )

main_df <- analysis_pw %>%
  select(
    pid, syear,
    srh, srh_rev, srh_rev_ord,
    event_time, event_time_f,
    vol_pre, vol_pre_year, vol_gap_years,
    vol_pre_binary, vol_pre_cat, vol_any,
    age, age_c, sex, migback,
    hinc_1k, hinc_1k_c, phrf,
    isced_num, isced_clean_f, ses_edu
  ) %>%
  filter(
    !is.na(pid),
    !is.na(srh_rev),
    !is.na(event_time),
    !is.na(event_time_f),
    !is.na(vol_pre_binary),
    !is.na(vol_pre_cat),
    !is.na(age_c),
    !is.na(sex),
    !is.na(migback),
    !is.na(hinc_1k_c)
  )

# ===============================================================
# 0b) Selection model for IPW
# ===============================================================
included_model_pids <- main_df %>%
  distinct(pid) %>%
  pull(pid)

sel_df_person <- step3_snapshot %>%
  mutate(
    included_model = if_else(pid %in% included_model_pids, 1L, 0L),
    sex = as.character(sex),
    migback = as.character(migback),
    isced_miss = if_else(is.na(isced_clean), 1L, 0L),
    hinc_miss  = if_else(is.na(hinc), 1L, 0L),
    srh_miss   = if_else(is.na(srh), 1L, 0L),
    isced_clean = if_else(is.na(isced_clean), -99, isced_clean),
    hinc = if_else(is.na(hinc), median(hinc, na.rm = TRUE), hinc),
    srh  = if_else(is.na(srh), median(srh, na.rm = TRUE), srh)
  )

sel_mod <- glm(
  included_model ~ age + sex + migback + isced_clean + hinc + srh +
    isced_miss + hinc_miss + srh_miss,
  data = sel_df_person,
  family = binomial(),
  na.action = na.exclude
)

sel_df_person <- sel_df_person %>%
  mutate(
    p_incl = predict(sel_mod, newdata = sel_df_person, type = "response"),
    p_incl_trim = pmin(pmax(p_incl, 0.01), 0.99),
    ipw_stab = mean(included_model, na.rm = TRUE) / p_incl_trim
  )

write_csv(
  sel_df_person %>% select(pid, p_incl, p_incl_trim, ipw_stab),
  file.path(out, "selection_probabilities_person_modelsample.csv")
)

main_df <- main_df %>%
  left_join(sel_df_person %>% select(pid, ipw_stab), by = "pid") %>%
  mutate(
    analysis_weight = if_else(!is.na(ipw_stab), ipw_stab, 1),
    phrf = as.numeric(phrf)
  )

ord_df <- main_df %>%
  filter(!is.na(srh_rev_ord))

# Small data summary for write-up
analysis_numbers <- tibble(
  n_persons = dplyr::n_distinct(main_df$pid),
  n_person_years = nrow(main_df),
  mean_srh_rev = mean(main_df$srh_rev, na.rm = TRUE),
  share_any_volunteering = mean(main_df$vol_pre_binary == "Any volunteering", na.rm = TRUE),
  share_female = mean(main_df$sex == "Female", na.rm = TRUE)
)

save_df_html(
  analysis_numbers,
  title = "Analysis Sample Summary",
  file = file.path(tab_dir, "01_analysis_sample_numbers.html"),
  notes = "Sample restricted to observed retirement transitions with valid pre-retirement volunteering and complete covariates for the pooled models."
)

saveRDS(main_df, file.path(mod_dir, "main_df_analysis.rds"))
saveRDS(ord_df, file.path(mod_dir, "ord_df_analysis.rds"))

# ===============================================================
# 1) MAIN MODEL: POOLED OLS
# ===============================================================
pooled_formula_bin <- srh_rev ~ event_time_f * vol_pre_binary +
  age_c + sex + migback + isced_clean_f + hinc_1k_c

pooled_formula_cat <- srh_rev ~ event_time_f * vol_pre_cat +
  age_c + sex + migback + isced_clean_f + hinc_1k_c

m_ols_bin <- fixest::feols(
  pooled_formula_bin,
  data = main_df,
  cluster = ~pid
)

m_ols_cat <- fixest::feols(
  pooled_formula_cat,
  data = main_df,
  cluster = ~pid
)

saveRDS(m_ols_bin, file.path(mod_dir, "m_ols_bin.rds"))
saveRDS(m_ols_cat, file.path(mod_dir, "m_ols_cat.rds"))

save_model_html(
  models = list(
    "Pooled OLS: Binary volunteering" = m_ols_bin,
    "Pooled OLS: Volunteering frequency" = m_ols_cat
  ),
  title = "Pooled OLS Models with Event-Time Interactions",
  file = file.path(tab_dir, "02_pooled_ols_models.html"),
  notes = "Standard errors clustered at the individual level. Outcome: reverse-coded self-rated health. Continuous covariates centered; household income entered in thousands."
)

# ===============================================================
# 2) FE EVENT-TIME MODEL
# ===============================================================
m_fe_bin <- fixest::feols(
  srh_rev ~ i(event_time, ref = -1) + i(event_time, vol_any, ref = -1) | pid,
  data = main_df,
  cluster = ~pid
)

m_fe_cat <- fixest::feols(
  srh_rev ~ i(event_time, ref = -1) + i(event_time, vol_pre_cat, ref = -1) | pid,
  data = main_df,
  cluster = ~pid
)

saveRDS(m_fe_bin, file.path(mod_dir, "m_fe_bin.rds"))
saveRDS(m_fe_cat, file.path(mod_dir, "m_fe_cat.rds"))

save_model_html(
  models = list(
    "FE event-time: Binary volunteering" = m_fe_bin,
    "FE event-time: Volunteering frequency" = m_fe_cat
  ),
  title = "Individual Fixed-Effects Event-Time Models",
  file = file.path(tab_dir, "03_fe_eventtime_models.html"),
  notes = "Standard errors clustered at the individual level. Fixed effects absorb time-invariant individual characteristics."
)

# ===============================================================
# 3) ORDERED LOGIT ROBUSTNESS
# ===============================================================
m_clm_bin <- ordinal::clm(
  srh_rev_ord ~ event_time_f * vol_pre_binary +
    age_c + sex + migback + isced_clean_f + hinc_1k_c,
  data = ord_df,
  link = "logit",
  Hess = TRUE
)

saveRDS(m_clm_bin, file.path(mod_dir, "m_clm_bin.rds"))

m_clm_cat <- NULL
try({
  m_clm_cat <- ordinal::clm(
    srh_rev_ord ~ event_time_f * vol_pre_cat +
      age_c + sex + migback + isced_clean_f + hinc_1k_c,
    data = ord_df,
    link = "logit",
    Hess = TRUE
  )
  saveRDS(m_clm_cat, file.path(mod_dir, "m_clm_cat.rds"))
}, silent = FALSE)

if (!is.null(m_clm_cat)) {
  save_model_html(
    models = list(
      "Ordered logit: Binary volunteering" = m_clm_bin,
      "Ordered logit: Volunteering frequency" = m_clm_cat
    ),
    title = "Ordinal Robustness Models",
    file = file.path(tab_dir, "04_ordered_logit_models.html"),
    notes = "Ordered logit cumulative-link models. Model-based standard errors."
  )
} else {
  save_model_html(
    models = list(
      "Ordered logit: Binary volunteering" = m_clm_bin
    ),
    title = "Ordinal Robustness Model",
    file = file.path(tab_dir, "04_ordered_logit_models.html"),
    notes = "Binary volunteering model only. Frequency model omitted due to instability in sparse categories."
  )
}

# ===============================================================
# 4) GENDER HETEROGENEITY
# ===============================================================
main_male <- main_df %>% filter(sex == "Male")
main_female <- main_df %>% filter(sex == "Female")

m_ols_male <- fixest::feols(
  srh_rev ~ event_time_f * vol_pre_binary +
    age_c + migback + isced_clean_f + hinc_1k_c,
  data = main_male,
  cluster = ~pid
)

m_ols_female <- fixest::feols(
  srh_rev ~ event_time_f * vol_pre_binary +
    age_c + migback + isced_clean_f + hinc_1k_c,
  data = main_female,
  cluster = ~pid
)

m_ols_sex_interaction <- fixest::feols(
  srh_rev ~ event_time_f * vol_pre_binary * sex +
    age_c + migback + isced_clean_f + hinc_1k_c,
  data = main_df,
  cluster = ~pid
)

saveRDS(m_ols_male, file.path(mod_dir, "m_ols_male.rds"))
saveRDS(m_ols_female, file.path(mod_dir, "m_ols_female.rds"))
saveRDS(m_ols_sex_interaction, file.path(mod_dir, "m_ols_sex_interaction.rds"))

save_model_html(
  models = list(
    "Men" = m_ols_male,
    "Women" = m_ols_female,
    "Sex interaction" = m_ols_sex_interaction
  ),
  title = "Gender Heterogeneity Models",
  file = file.path(tab_dir, "05_gender_heterogeneity_models.html"),
  notes = "Pooled OLS models with standard errors clustered at the individual level."
)

# ===============================================================
# 5) SES HETEROGENEITY
# ===============================================================
ses_df <- main_df %>%
  filter(!is.na(ses_edu))

m_ols_ses_interaction <- fixest::feols(
  srh_rev ~ event_time_f * vol_pre_binary * ses_edu +
    age_c + sex + migback + hinc_1k_c,
  data = ses_df,
  cluster = ~pid
)

saveRDS(m_ols_ses_interaction, file.path(mod_dir, "m_ols_ses_interaction.rds"))

save_model_html(
  models = list(
    "SES interaction" = m_ols_ses_interaction
  ),
  title = "Education-Based SES Heterogeneity Model",
  file = file.path(tab_dir, "06_ses_heterogeneity_model.html"),
  notes = "Pooled OLS model with standard errors clustered at the individual level."
)

# ===============================================================
# 6) IPW ROBUSTNESS
# ===============================================================
m_ols_bin_ipw <- fixest::feols(
  pooled_formula_bin,
  data = main_df,
  weights = ~analysis_weight,
  cluster = ~pid
)

saveRDS(m_ols_bin_ipw, file.path(mod_dir, "m_ols_bin_ipw.rds"))

save_model_html(
  models = list(
    "IPW pooled OLS" = m_ols_bin_ipw
  ),
  title = "IPW Robustness Model",
  file = file.path(tab_dir, "07_weighted_robustness_model.html"),
  notes = "Pooled OLS model weighted by stabilized inverse-probability weights; standard errors clustered at the individual level."
)

# ===============================================================
# 7) phrf-WEIGHTED ROBUSTNESS
# ===============================================================
main_df_phrf <- main_df %>%
  filter(!is.na(phrf), phrf > 0)

if (nrow(main_df_phrf) > 0) {
  m_ols_bin_phrf <- fixest::feols(
    pooled_formula_bin,
    data = main_df_phrf,
    weights = ~phrf,
    cluster = ~pid
  )
  
  saveRDS(m_ols_bin_phrf, file.path(mod_dir, "m_ols_bin_phrf.rds"))
  
  save_model_html(
    models = list(
      "Pooled OLS (phrf-weighted)" = m_ols_bin_phrf
    ),
    title = "Population-Weighted Robustness Model",
    file = file.path(tab_dir, "07b_phrf_weighted_robustness_model.html"),
    notes = "Pooled OLS model weighted by phrf survey/population weights; standard errors clustered at the individual level. Appendix-only robustness."
  )
}

# ===============================================================
# 8) JOINT WALD TESTS FOR DIFFERENTIAL TRAJECTORIES
# ===============================================================
pooled_inter_terms <- names(coef(m_ols_bin))[grepl(
  "event_time_f.*:vol_pre_binaryAny volunteering",
  names(coef(m_ols_bin))
)]

fe_inter_terms <- names(coef(m_fe_bin))[grepl(
  "event_time::.*:vol_any",
  names(coef(m_fe_bin))
)]

pooled_post_terms <- pooled_inter_terms[!grepl("event_time_f-2", pooled_inter_terms, fixed = TRUE)]
fe_post_terms <- fe_inter_terms[!grepl("event_time::\\-2", fe_inter_terms)]

wt_pooled_all  <- fixest::wald(m_ols_bin, pooled_inter_terms, print = FALSE)
wt_pooled_post <- fixest::wald(m_ols_bin, pooled_post_terms, print = FALSE)
wt_fe_all      <- fixest::wald(m_fe_bin, fe_inter_terms, print = FALSE)
wt_fe_post     <- fixest::wald(m_fe_bin, fe_post_terms, print = FALSE)

joint_tests <- bind_rows(
  extract_wald_row(wt_pooled_all,  "Pooled OLS: all volunteering x event-time interactions"),
  extract_wald_row(wt_pooled_post, "Pooled OLS: post-retirement volunteering x event-time interactions"),
  extract_wald_row(wt_fe_all,      "FE model: all volunteering x event-time interactions"),
  extract_wald_row(wt_fe_post,     "FE model: post-retirement volunteering x event-time interactions")
)

save_df_html(
  joint_tests,
  title = "Joint Tests of Differential Health Trajectories",
  file = file.path(tab_dir, "08_joint_tests_eventtime_interactions.html"),
  notes = "Wald tests for the volunteering-by-event-time interaction terms."
)

# ===============================================================
# 9) SURVEY-YEAR FE ROBUSTNESS
# ===============================================================
m_ols_bin_syear <- fixest::feols(
  srh_rev ~ event_time_f * vol_pre_binary +
    age_c + sex + migback + isced_clean_f + hinc_1k_c | syear,
  data = main_df,
  cluster = ~pid
)

m_fe_bin_syear <- fixest::feols(
  srh_rev ~ i(event_time, ref = -1) + i(event_time, vol_any, ref = -1) | pid + syear,
  data = main_df,
  cluster = ~pid
)

saveRDS(m_ols_bin_syear, file.path(mod_dir, "m_ols_bin_syear.rds"))
saveRDS(m_fe_bin_syear, file.path(mod_dir, "m_fe_bin_syear.rds"))

save_model_html(
  models = list(
    "Pooled OLS + survey-year FE" = m_ols_bin_syear,
    "Individual FE + survey-year FE" = m_fe_bin_syear
  ),
  title = "Survey-Year Fixed Effects Robustness Models",
  file = file.path(tab_dir, "09_syear_fe_robustness.html"),
  notes = "Appendix robustness including survey-year fixed effects."
)

# ===============================================================
# 10) DIAGNOSTICS
# ===============================================================

# 10a) VIFs / rank / exact collinearity
lm_diag <- lm(pooled_formula_bin, data = main_df)

vif_df <- compute_vif_from_lm(lm_diag)
save_df_html(
  vif_df,
  title = "Multicollinearity Diagnostics (Pooled OLS Design)",
  file = file.path(tab_dir, "10_diagnostics_vif.html"),
  notes = "Column-level VIFs from the pooled OLS design matrix. Factor variables are expanded into indicator columns."
)

design_mm <- model.matrix(lm_diag)
if ("(Intercept)" %in% colnames(design_mm)) {
  design_mm_noint <- design_mm[, colnames(design_mm) != "(Intercept)", drop = FALSE]
} else {
  design_mm_noint <- design_mm
}

condition_number <- kappa(design_mm_noint)
alias_complete <- alias(lm_diag)$Complete

diag_rank <- tibble(
  diagnostic = c(
    "Condition number",
    "Design matrix columns",
    "Design matrix rank",
    "Exact collinearity detected"
  ),
  value = c(
    round(condition_number, 2),
    ncol(design_mm_noint),
    qr(design_mm_noint)$rank,
    ifelse(is.null(alias_complete), "No", "Yes")
  )
)

save_df_html(
  diag_rank,
  title = "Design Matrix Diagnostics",
  file = file.path(tab_dir, "11_diagnostics_design_matrix.html"),
  notes = "Condition number and rank diagnostics for the pooled OLS specification."
)

# 10b) Heteroskedasticity diagnostic
bp_test <- lmtest::bptest(lm_diag)

bp_df <- tibble(
  statistic = as.numeric(bp_test$statistic),
  df = as.numeric(bp_test$parameter),
  p_value = as.numeric(bp_test$p.value)
)

save_df_html(
  bp_df,
  title = "Breusch-Pagan Heteroskedasticity Diagnostic",
  file = file.path(tab_dir, "12_diagnostics_bptest.html"),
  notes = "Diagnostic only. Main inference uses individual-level clustered standard errors."
)

# 10c) Event-time by volunteering cell counts
cell_counts <- main_df %>%
  count(event_time, vol_pre_binary, name = "n") %>%
  complete(
    event_time = -2:5,
    vol_pre_binary = factor(
      c("No volunteering", "Any volunteering"),
      levels = c("No volunteering", "Any volunteering")
    ),
    fill = list(n = 0)
  ) %>%
  group_by(event_time) %>%
  mutate(share_within_event_time = n / sum(n)) %>%
  ungroup()

save_df_html(
  cell_counts,
  title = "Cell Counts by Event Time and Volunteering Group",
  file = file.path(tab_dir, "13_diagnostics_eventtime_cells.html"),
  notes = "Used to inspect support and sparsity across event-time × volunteering cells."
)

# 10d) IPW distribution
ipw_summary <- sel_df_person %>%
  summarise(
    n_persons = n(),
    min_ipw = min(ipw_stab, na.rm = TRUE),
    p1_ipw = quantile(ipw_stab, 0.01, na.rm = TRUE),
    p5_ipw = quantile(ipw_stab, 0.05, na.rm = TRUE),
    median_ipw = median(ipw_stab, na.rm = TRUE),
    mean_ipw = mean(ipw_stab, na.rm = TRUE),
    p95_ipw = quantile(ipw_stab, 0.95, na.rm = TRUE),
    p99_ipw = quantile(ipw_stab, 0.99, na.rm = TRUE),
    max_ipw = max(ipw_stab, na.rm = TRUE)
  )

save_df_html(
  ipw_summary,
  title = "IPW Distribution Summary",
  file = file.path(tab_dir, "14_diagnostics_ipw_summary.html"),
  notes = "Stabilized inverse-probability weights for selection into the estimation sample."
)

# 10e) Person-wave distribution
wave_dist <- main_df %>%
  count(pid, name = "n_waves") %>%
  count(n_waves, name = "n_persons") %>%
  arrange(n_waves)

save_df_html(
  wave_dist,
  title = "Person-Wave Distribution in Estimation Sample",
  file = file.path(tab_dir, "15_diagnostics_personwave_distribution.html"),
  notes = "Number of observed event-time person-years per individual in the estimation sample."
)

# 10f) Optional panel serial-correlation diagnostic
if (requireNamespace("plm", quietly = TRUE)) {
  serial_note <- tryCatch(
    {
      pdata <- plm::pdata.frame(
        main_df %>% select(pid, syear, srh_rev, event_time),
        index = c("pid", "syear")
      )
      fe_serial_model <- plm::plm(
        srh_rev ~ factor(event_time),
        data = pdata,
        model = "within"
      )
      serial_test <- plm::pbgtest(fe_serial_model, order = 1)
      
      tibble(
        statistic = as.numeric(serial_test$statistic),
        df = as.numeric(serial_test$parameter),
        p_value = as.numeric(serial_test$p.value)
      )
    },
    error = function(e) {
      tibble(
        statistic = NA_real_,
        df = NA_real_,
        p_value = NA_real_
      )
    }
  )
  
  save_df_html(
    serial_note,
    title = "Optional Panel Serial-Correlation Diagnostic",
    file = file.path(tab_dir, "16_diagnostics_panel_serial_correlation.html"),
    notes = "Breusch-Godfrey/Wooldridge-style diagnostic using an analogous within model. Main inference uses individual-level clustered standard errors."
  )
}

# ===============================================================
# 11) MAIN-TEXT PLOT: POOLED OLS PREDICTIONS
# ===============================================================
pred_grid_bin <- tidyr::expand_grid(
  event_time_f = factor(-2:5, levels = -2:5),
  vol_pre_binary = factor(
    c("No volunteering", "Any volunteering"),
    levels = c("No volunteering", "Any volunteering")
  ),
  age_c = 0,
  sex = factor(mode_value(main_df$sex), levels = levels(main_df$sex)),
  migback = factor(mode_value(main_df$migback), levels = levels(main_df$migback)),
  isced_clean_f = factor(
    mode_value(main_df$isced_clean_f),
    levels = levels(main_df$isced_clean_f)
  ),
  hinc_1k_c = 0
) %>%
  mutate(
    event_time = as.numeric(as.character(event_time_f))
  )

pred_vals <- predict(m_ols_bin, newdata = pred_grid_bin, se.fit = TRUE)

pred_plot_df <- pred_grid_bin %>%
  mutate(
    fit = pred_vals$fit,
    se = pred_vals$se.fit,
    ci_low = fit - 1.96 * se,
    ci_high = fit + 1.96 * se
  )

saveRDS(pred_plot_df, file.path(mod_dir, "predicted_values_pooled_ols_binary.rds"))

p_pred_bin <- ggplot(
  pred_plot_df,
  aes(
    x = event_time, y = fit,
    group = vol_pre_binary,
    color = vol_pre_binary,
    linetype = vol_pre_binary
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = -2:5) +
  scale_color_manual(values = vol_colors, drop = FALSE) +
  scale_linetype_manual(values = vol_linetypes, drop = FALSE) +
  labs(
    title = "Predicted self-rated health around retirement",
    subtitle = "Adjusted predictions from pooled OLS with event-time interactions",
    x = "Years relative to retirement",
    y = "Predicted self-rated health (reversed; higher = better)",
    color = "Pre-retirement volunteering",
    linetype = "Pre-retirement volunteering"
  ) +
  thesis_theme()

save_plot_dual(p_pred_bin, "02_predicted_srh_pooled_ols_binary", width = 8, height = 5)

# ===============================================================
# 12) MAIN-TEXT PLOT: FE EVENT-TIME FIGURE
# ===============================================================
build_fe_plot_binary <- function(model_obj) {
  b <- coef(model_obj)
  V <- vcov(model_obj)
  coef_names <- names(b)
  
  get_base_name <- function(k) {
    nm <- grep(paste0("event_time::", k), coef_names, value = TRUE)
    nm[!grepl(":vol_any", nm)][1]
  }
  
  get_int_name <- function(k) {
    nm <- grep(paste0("event_time::", k), coef_names, value = TRUE)
    nm[grepl(":vol_any", nm)][1]
  }
  
  times <- -2:5
  
  out_no <- purrr::map_dfr(times, function(k) {
    if (k == -1) {
      tibble(
        vol_group = "No volunteering",
        event_time = k,
        estimate = 0,
        conf.low = 0,
        conf.high = 0
      )
    } else {
      base_name <- get_base_name(k)
      est <- unname(b[base_name])
      se <- sqrt(V[base_name, base_name])
      
      tibble(
        vol_group = "No volunteering",
        event_time = k,
        estimate = est,
        conf.low = est - 1.96 * se,
        conf.high = est + 1.96 * se
      )
    }
  })
  
  out_yes <- purrr::map_dfr(times, function(k) {
    if (k == -1) {
      tibble(
        vol_group = "Any volunteering",
        event_time = k,
        estimate = 0,
        conf.low = 0,
        conf.high = 0
      )
    } else {
      base_name <- get_base_name(k)
      int_name  <- get_int_name(k)
      
      if (is.na(int_name) || is.null(int_name)) {
        est <- unname(b[base_name])
        se <- sqrt(V[base_name, base_name])
      } else {
        est <- unname(b[base_name] + b[int_name])
        var_est <- V[base_name, base_name] +
          V[int_name, int_name] +
          2 * V[base_name, int_name]
        se <- sqrt(var_est)
      }
      
      tibble(
        vol_group = "Any volunteering",
        event_time = k,
        estimate = est,
        conf.low = est - 1.96 * se,
        conf.high = est + 1.96 * se
      )
    }
  })
  
  bind_rows(out_no, out_yes) %>%
    mutate(
      vol_group = factor(
        vol_group,
        levels = c("No volunteering", "Any volunteering")
      )
    ) %>%
    arrange(vol_group, event_time)
}

fe_plot_df <- build_fe_plot_binary(m_fe_bin)
saveRDS(fe_plot_df, file.path(mod_dir, "eventstudy_fe_binary_plotdata.rds"))

p_fe <- ggplot(
  fe_plot_df,
  aes(
    x = event_time, y = estimate,
    group = vol_group,
    color = vol_group,
    linetype = vol_group
  )
) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.12, linewidth = 0.5) +
  scale_x_continuous(breaks = -2:5) +
  scale_color_manual(values = vol_colors, drop = FALSE) +
  scale_linetype_manual(values = vol_linetypes, drop = FALSE) +
  labs(
    title = "Within-person self-rated health around retirement",
    subtitle = "Event-time model with individual fixed effects; reference year = -1",
    x = "Years relative to retirement",
    y = "Change in self-rated health",
    color = "Pre-retirement volunteering",
    linetype = "Pre-retirement volunteering"
  ) +
  thesis_theme()

save_plot_dual(p_fe, "03_eventstudy_fe_binary", width = 8, height = 5)

# ===============================================================
# 13) APPENDIX HETEROGENEITY PLOTS
# ===============================================================

# 13a) Sex heterogeneity
pred_grid_sex <- tidyr::expand_grid(
  event_time_f = factor(-2:5, levels = -2:5),
  vol_pre_binary = factor(
    c("No volunteering", "Any volunteering"),
    levels = c("No volunteering", "Any volunteering")
  ),
  sex = factor(c("Male", "Female"), levels = levels(main_df$sex)),
  age_c = 0,
  migback = factor(mode_value(main_df$migback), levels = levels(main_df$migback)),
  isced_clean_f = factor(
    mode_value(main_df$isced_clean_f),
    levels = levels(main_df$isced_clean_f)
  ),
  hinc_1k_c = 0
) %>%
  mutate(
    event_time = as.numeric(as.character(event_time_f))
  )

pred_vals_sex <- predict(m_ols_sex_interaction, newdata = pred_grid_sex, se.fit = TRUE)

pred_plot_sex <- pred_grid_sex %>%
  mutate(
    fit = pred_vals_sex$fit,
    se = pred_vals_sex$se.fit,
    ci_low = fit - 1.96 * se,
    ci_high = fit + 1.96 * se
  )

saveRDS(pred_plot_sex, file.path(mod_dir, "predicted_values_gender_heterogeneity.rds"))

p_het_sex <- ggplot(
  pred_plot_sex,
  aes(
    x = event_time, y = fit,
    group = vol_pre_binary,
    color = vol_pre_binary,
    linetype = vol_pre_binary
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.10, linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~sex) +
  scale_x_continuous(breaks = -2:5) +
  scale_color_manual(values = vol_colors, drop = FALSE) +
  scale_linetype_manual(values = vol_linetypes, drop = FALSE) +
  labs(
    title = "Predicted self-rated health around retirement by sex",
    subtitle = "Adjusted pooled OLS predictions from the sex interaction model",
    x = "Years relative to retirement",
    y = "Predicted self-rated health (reversed; higher = better)",
    color = "Pre-retirement volunteering",
    linetype = "Pre-retirement volunteering"
  ) +
  thesis_theme()

save_plot_dual(p_het_sex, "A9_heterogeneity_by_sex", width = 9, height = 5)

# 13b) SES heterogeneity
pred_grid_ses <- tidyr::expand_grid(
  event_time_f = factor(-2:5, levels = -2:5),
  vol_pre_binary = factor(
    c("No volunteering", "Any volunteering"),
    levels = c("No volunteering", "Any volunteering")
  ),
  ses_edu = factor(
    c("Lower / medium education", "Higher education"),
    levels = levels(ses_df$ses_edu)
  ),
  age_c = 0,
  sex = factor(mode_value(ses_df$sex), levels = levels(ses_df$sex)),
  migback = factor(mode_value(ses_df$migback), levels = levels(ses_df$migback)),
  hinc_1k_c = 0
) %>%
  mutate(
    event_time = as.numeric(as.character(event_time_f))
  )

pred_vals_ses <- predict(m_ols_ses_interaction, newdata = pred_grid_ses, se.fit = TRUE)

pred_plot_ses <- pred_grid_ses %>%
  mutate(
    fit = pred_vals_ses$fit,
    se = pred_vals_ses$se.fit,
    ci_low = fit - 1.96 * se,
    ci_high = fit + 1.96 * se
  )

saveRDS(pred_plot_ses, file.path(mod_dir, "predicted_values_ses_heterogeneity.rds"))

p_het_ses <- ggplot(
  pred_plot_ses,
  aes(
    x = event_time, y = fit,
    group = vol_pre_binary,
    color = vol_pre_binary,
    linetype = vol_pre_binary
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.10, linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ses_edu) +
  scale_x_continuous(breaks = -2:5) +
  scale_color_manual(values = vol_colors, drop = FALSE) +
  scale_linetype_manual(values = vol_linetypes, drop = FALSE) +
  labs(
    title = "Predicted self-rated health around retirement by education-based SES",
    subtitle = "Adjusted pooled OLS predictions from the SES interaction model",
    x = "Years relative to retirement",
    y = "Predicted self-rated health (reversed; higher = better)",
    color = "Pre-retirement volunteering",
    linetype = "Pre-retirement volunteering"
  ) +
  thesis_theme()

save_plot_dual(p_het_ses, "A10_heterogeneity_by_ses", width = 9, height = 5)

# ===============================================================
# 14) THESIS-READY SUMMARY TABLES
# ===============================================================
save_model_html(
  models = list(
    "Pooled OLS: Binary volunteering" = m_ols_bin,
    "FE event-time: Binary volunteering" = m_fe_bin,
    "Ordered logit: Binary volunteering" = m_clm_bin
  ),
  title = "Core Binary-Specification Models",
  file = file.path(tab_dir, "17_core_binary_models.html"),
  notes = "Core specification table. Full interaction details are easier to interpret in the figures."
)

model_overview <- tibble(
  model_id = c("M1", "M2", "M3", "M4", "M5", "M6", "M7"),
  description = c(
    "Pooled OLS, binary volunteering, clustered SE by pid",
    "Pooled OLS, volunteering frequency, clustered SE by pid",
    "Individual FE event-time model, binary volunteering, clustered SE by pid",
    "Individual FE event-time model, volunteering frequency, clustered SE by pid",
    "Ordered logit robustness, binary volunteering",
    "Sex heterogeneity: pooled OLS with three-way interaction",
    "SES heterogeneity: pooled OLS with three-way interaction"
  )
)

save_df_html(
  model_overview,
  title = "Empirical Strategy Overview",
  file = file.path(tab_dir, "18_model_overview.html"),
  notes = "Main and supplementary models aligned with the thesis design."
)

writeup_models <- tibble(
  main_sample_persons = dplyr::n_distinct(main_df$pid),
  main_sample_personyears = nrow(main_df),
  ordered_logit_personyears = nrow(ord_df),
  male_stratified_persons = dplyr::n_distinct(main_male$pid),
  female_stratified_persons = dplyr::n_distinct(main_female$pid),
  ses_analysis_persons = dplyr::n_distinct(ses_df$pid)
)

save_df_html(
  writeup_models,
  title = "Model Sample Sizes",
  file = file.path(tab_dir, "19_model_sample_sizes.html"),
  notes = "Counts correspond to estimation samples used in the analysis script."
)

# ===============================================================
# 15) BASIC REPRODUCIBILITY FOOTPRINT
# ===============================================================
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_analysis.txt"))

message("Empirical strategy script completed successfully.")