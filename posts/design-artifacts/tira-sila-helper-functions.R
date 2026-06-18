# ----------------------------
# helper: regression summary
# ----------------------------
lm_summary <- function(df) {
  model <- lm(y ~ x, data = df)
  tibble(
    beta = coef(model)[2],
    alpha = coef(model)[1],
    adj_r2 = summary(model)$adj.r.squared
  )
}

# ----------------------------
# helper: panel stats
# ----------------------------
compute_panel_stats <- function(dat_long, source_name) {
  dat_long %>%
    filter(Source == source_name, !is.na(panel)) %>%
    group_by(panel, Method) %>%
    summarise(
      n = sum(complete.cases(x, y)),
      rho = cor(x, y, method = "spearman", use = "complete.obs"),
      stats = list(lm_summary(cur_data())),
      x_min = min(x, na.rm = TRUE),
      x_max = max(x, na.rm = TRUE),
      y_min = min(y, na.rm = TRUE),
      y_max = max(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    unnest(stats) %>%
    mutate(
      label  = sprintf(
        "y = %.2fx + %.1f <br>Adjusted R<sup>2</sup>=%.2f<br>Spearman \u03c1 = %.2f<br>N = %d",
        beta, alpha, adj_r2, rho, n),
      x_range = x_max - x_min,
      y_range = y_max - y_min,
      x_pos = x_min + 0.02 * x_range,
      y_pos = y_max - 0.10 * y_range
    )
}

# ----------------------------
# helper: segment data
# ----------------------------
make_seg_data <- function(dat, source_name, panel_levels) {
  dat %>%
    filter(Source == source_name) %>%
    mutate(
      y = age_of_symptom_onset,
      x_start = baseline_plasma_age,
      x_end = ptau217_positivity_age,
      panel = factor(
        "Estimated age at %p-tau217 positivity (yrs)",
        levels = panel_levels
      )
    )
}

# ----------------------------
# helper: constraint ribbon
# ----------------------------
make_constraint <- function(dat_long, source_name) {
  base <- dat_long %>%
    filter(Source == source_name, panel == "Baseline age (yrs)")
  
  age_range <- range(base$x) %>% round()
  obs_range <- range(base$y - base$x)
  
  expand_grid(
    x = seq(age_range[1], age_range[2]),
    panel = factor("Baseline age (yrs)", levels = levels(dat_long$panel)),
    Method = unique(dat_long$Method)
  ) %>%
    mutate(
      y_min = x + obs_range[1],
      y_max = x + obs_range[2],
      y = NA_real_
    )
}

# ----------------------------
# generic plotting function
# ----------------------------
plot_panels <- function(dat_long, source_name, panel_stats, seg_data, constraint_data, y_pos_override = NULL) {
  
  if (!is.null(y_pos_override)) {
    panel_stats$y_pos <- y_pos_override
  }
  
  ggplot(dat_long %>% filter(Source == source_name, !is.na(panel)),
    aes(x = x, y = y)) +
    geom_point(size = 2.5, alpha = 0.4) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_grid(rows = vars(Method), cols = vars(panel), scales = "free_x", switch = "x") +
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

# ----------------------------
# helper: commonality table
# ----------------------------
compute_commonality <- function(df, dv, ivlist) {
  res <- commonalityCoefficients(
    df,
    dv = dv,
    ivlist = ivlist
  )
  
  adj_r2 <- summary(lm(
    age_of_symptom_onset ~ baseline_plasma_age + time_from_positivity,
    data = df
  ))$adj.r.squared
  
  res$CC %>%
    as.data.frame() %>%
    rownames_to_column("Effect") %>%
    bind_rows(tibble(Effect = "Adj.", Coefficient = adj_r2))
}
