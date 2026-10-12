# =============================================================================
# Helper functions for "Predicting the age of Alzheimer's symptom onset from a
# blood test: evaluating a prognostic claim" (Insel & Donohue, medRxiv
# 10.64898/2026.03.26.26349230).
#
# All functions take the analysis data set `dat` (one row per participant per
# clock method; see the data-preparation chunk in index.qmd) or its long
# version `dat_long` (one row per participant x method x predictor).
#
# Key variables (all in years):
#   baseline_plasma_age    age at blood draw (first plasma collection)
#   baseline_plasma        raw %p-tau217 at that blood draw
#   time_from_positivity   clock-estimated time since %p-tau217 positivity at
#                          the blood draw (negative = not yet positive)
#   ptau217_positivity_age clock-estimated age at %p-tau217 positivity
#                          = baseline_plasma_age - time_from_positivity
#   age_of_symptom_onset   observed age at symptom onset (the outcome)
# =============================================================================

# Labels for the three predictors in Figure 1 / Supplementary Figures S1-S2.
# Defined once so panel factors, segment data and constraint ribbons all use
# identical strings (facet matching in ggplot2 is by exact factor level).
PANEL_AGE  <- "Age at blood draw (yrs)"
PANEL_TIME <- "Estimated time from %p-tau217 positivity (yrs)"
PANEL_POS  <- "Estimated age at %p-tau217 positivity (yrs)"
PANEL_LEVELS <- c(PANEL_AGE, PANEL_TIME, PANEL_POS)

# ----------------------------------------------------------------------------
# Regression summary used in every panel label: slope, intercept, adjusted R2
# from ordinary least squares of y on x.
# ----------------------------------------------------------------------------
lm_summary <- function(df) {
  model <- lm(y ~ x, data = df)
  tibble(
    beta = coef(model)[2],
    alpha = coef(model)[1],
    adj_r2 = summary(model)$adj.r.squared
  )
}

# HTML label (rendered by ggtext::geom_richtext) reporting the fitted line,
# adjusted R2, Spearman's rho and N. `adj_r2_extra` optionally appends a
# second adjusted R2 in parentheses (used for the restricted-sample fits in
# Figure 2), colored to match the restricted-sample regression line.
stats_label <- function(beta, alpha, adj_r2, rho, n, adj_r2_extra = NULL,
  extra_color = "darkred") {
  r2_txt <- sprintf("%.2f", adj_r2)
  if (!is.null(adj_r2_extra)) {
    r2_txt <- sprintf("%s (<span style='color:%s'>%.2f</span>)",
      r2_txt, extra_color, adj_r2_extra)
  }
  sprintf(
    "y = %.2fx + %.1f <br>Adjusted R<sup>2</sup>=%s<br>Spearman ρ = %.2f<br>N = %d",
    beta, alpha, r2_txt, rho, n)
}

# ----------------------------------------------------------------------------
# Figure 1 / S1 / S2: per-panel regression statistics and label positions.
# One row per predictor (panel) x clock method within a cohort.
# ----------------------------------------------------------------------------
compute_panel_stats <- function(dat_long, source_name) {
  dat_long %>%
    filter(Source == source_name, !is.na(panel)) %>%
    group_by(panel, Method) %>%
    summarise(
      n = sum(complete.cases(x, y)),
      rho = cor(x, y, method = "spearman", use = "complete.obs"),
      stats = list(lm_summary(pick(x, y))),
      x_min = min(x, na.rm = TRUE),
      x_max = max(x, na.rm = TRUE),
      y_min = min(y, na.rm = TRUE),
      y_max = max(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    unnest(stats) %>%
    mutate(
      label = stats_label(beta, alpha, adj_r2, rho, n),
      # default position: upper left of each panel (can be overridden in
      # plot_panels() so all rows share a common height)
      x_pos = x_min + 0.02 * (x_max - x_min),
      y_pos = y_max - 0.10 * (y_max - y_min)
    )
}

# ----------------------------------------------------------------------------
# Figure 1 / S1 / S2: horizontal segments in the right-hand column joining
# each participant's age at blood draw to their estimated age at positivity.
# The segment length is the estimated time from positivity, i.e. the only
# part of the clock predictor that comes from the blood test.
# ----------------------------------------------------------------------------
make_seg_data <- function(dat, source_name) {
  dat %>%
    filter(Source == source_name) %>%
    mutate(
      y = age_of_symptom_onset,
      x_start = baseline_plasma_age,
      x_end = ptau217_positivity_age,
      panel = factor(PANEL_POS, levels = PANEL_LEVELS)
    )
}

# ----------------------------------------------------------------------------
# Figure 1 / S1 / S2: the bound imposed by the follow-up window (pink band in
# the left-hand column). Onset was only observed between min(onset - draw age)
# and max(onset - draw age) years from the blood draw, so for any age at blood
# draw x, age at onset must lie in [x + min, x + max]. This band is fixed by
# the design and contains no biomarker information.
# ----------------------------------------------------------------------------
make_constraint <- function(dat_long, source_name) {
  base <- dat_long %>%
    filter(Source == source_name, panel == PANEL_AGE)

  age_range <- range(base$x) %>% round()
  window <- range(base$y - base$x)   # observed onset window around the draw

  expand_grid(
    x = seq(age_range[1], age_range[2]),
    panel = factor(PANEL_AGE, levels = PANEL_LEVELS),
    Method = unique(dat_long$Method)
  ) %>%
    mutate(
      y_min = x + window[1],
      y_max = x + window[2],
      y = NA_real_
    )
}

# ----------------------------------------------------------------------------
# Figure 1 / S1 / S2: rows = clock method (TIRA / SILA / Random), columns =
# the three predictors of age at symptom onset.
#   y_pos_override: common vertical position for the statistics labels
# ----------------------------------------------------------------------------
plot_panels <- function(dat_long, source_name, panel_stats, seg_data,
  constraint_data, y_pos_override = NULL) {

  if (!is.null(y_pos_override)) {
    panel_stats$y_pos <- y_pos_override
  }

  ggplot(dat_long %>% filter(Source == source_name, !is.na(panel)),
    aes(x = x, y = y)) +
    geom_point(size = 2.5, alpha = 0.4) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_grid(rows = vars(Method), cols = vars(panel), scales = "free_x",
      switch = "x") +
    geom_segment(
      data = seg_data,
      aes(x = x_start, xend = x_end, y = y, yend = y),
      color = "grey30",
      alpha = 0.4,
      inherit.aes = FALSE
    ) +
    geom_richtext(
      data = panel_stats,
      aes(x = x_pos, y = y_pos, label = label),
      hjust = 0, vjust = 1,
      size = 3.5,
      fill = NA, label.color = NA,
      show.legend = FALSE
    ) +
    geom_ribbon(
      data = constraint_data,
      aes(x = x, ymin = y_min, ymax = y_max),
      alpha = 0.3, fill = "pink"
    ) +
    labs(
      x = "",
      y = "Age at symptom onset (yrs)"
    ) +
    theme_jama() +
    theme(
      strip.placement = "outside",
      strip.text.y = element_text(size = 15),
      strip.clip = "off",
      panel.spacing = unit(2, "lines")
    )
}

# ----------------------------------------------------------------------------
# Figure 3 and the sensitivity analysis: the four OLS models of age at
# symptom onset, fit separately in each cohort x clock-method sample.
#   1. Age at blood draw alone
#   2. Clock-estimated age at %p-tau217 positivity (the original report's
#      predictor)
#   3. "Random" age at positivity: age at blood draw minus a random time from
#      positivity (same participants as model 1-2)
#   4. Age at blood draw + raw %p-tau217
# Returns adjusted R2 and AIC (with delta AIC relative to model 1).
# ----------------------------------------------------------------------------
MODEL_LEVELS <- c(
  "Age at blood draw",
  "Age of %p-tau217 positivity",
  "Random age of %p-tau217 positivity",
  "Age at blood draw + %p-tau217")

fit_model_comparison <- function(dat) {
  map_dfr(c("Knight ADRC", "ADNI"), function(src) {
    map_dfr(c("SILA", "TIRA"), function(meth) {

      df_sub <- dat %>% filter(Source == src, Method == meth)
      # random comparator restricted to the same participants
      df_sub_rand <- dat %>%
        filter(Source == src, Method == "Random", id %in% df_sub$id)
      if (nrow(df_sub) != nrow(df_sub_rand)) {
        stop("random comparator sample differs from ", src, "/", meth)
      }

      mods <- list(
        lm(age_of_symptom_onset ~ baseline_plasma_age, data = df_sub),
        lm(age_of_symptom_onset ~ ptau217_positivity_age, data = df_sub),
        lm(age_of_symptom_onset ~ ptau217_positivity_age, data = df_sub_rand),
        lm(age_of_symptom_onset ~ baseline_plasma_age + baseline_plasma,
          data = df_sub))

      aic <- map_dbl(mods, AIC)
      tibble(
        Source = src,
        Method = meth,
        N = nrow(df_sub),
        Model = MODEL_LEVELS,
        R2 = map_dbl(mods, ~ summary(.x)$r.squared),
        Adj_R2 = map_dbl(mods, ~ summary(.x)$adj.r.squared),
        AIC = aic,
        delta_AIC = aic - aic[1])
    })
  })
}
