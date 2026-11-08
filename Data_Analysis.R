# Jan Blackburn
# Math 257: Multivariate Analysis Project
# Data Analysis File that imports the data analysis file created 
# and peforms the multivariate analysis for the project
# Initial code was generated using Claude and was modified by me as needed

#
# Outcomes (per House district, 2024): turnout and Democratic two-party share
#
# Input:   analysis_2024.csv   (from build_analysis_data.R)
# Output:  output/*.csv        tables for the report
#          figures/*.png       figures for the report
#          output/results.rds  everything, for the Quarto report to read
#
# Steps
#   1. Prepare the modeling sample
#   2. Descriptives and collinearity check
#   3. Univariate models: two separate multiple regressions
#   4. Multivariate model: one regression with both outcomes
#   5. Significance tests: Bonferroni t-tests vs. Pillai's trace
#   6. Cross-outcome tests: same effect on turnout and share?
#   7. Residual correlation between the outcomes
#   8. Joint confidence ellipses vs. Bonferroni rectangles
#   9. Missing-data cost of the joint model
#  10. Assumption diagnostics
#  11. Sensitivity check: poverty rate in place of the Gini index
#  12. Report exhibits: the draft report's tables, with plain-English names,
#      written to output/exhibits/ (figures are labeled in steps 8 and 10)
#
# Data update (35 columns): the analysis file now also carries pct_female,
# pct_low_inc, pct_high_inc, gini and pct_poverty. The main model adds gini and
# pct_female; the other three are too collinear with median income to enter
# together (see the candidate VIF table in step 2 and the check in step 11).
# =============================================================================

library(tidyverse)   # data handling and ggplot2
library(car)         # vif(), Anova(), linearHypothesis()
library(sandwich)    # vcovCL(): state-clustered standard errors
library(lmtest)      # coeftest()
library(broom)       # tidy() model output into data frames

rm(list=ls())

dir.create("output",  showWarnings = FALSE)
dir.create("output/exhibits", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
# -----------------------------------------------------------------------------
# 1. Prepare the modeling sample
# -----------------------------------------------------------------------------
# Shares are rescaled from proportions (0-1) to percentage points (0-100) so a
# coefficient reads as "points of outcome per 1 point of predictor".

# All pct_* columns (including the new pct_female, pct_low_inc, pct_high_inc
# and pct_poverty) are rescaled by starts_with("pct_"). The Gini index is a 0-1
# index, not a share, so it is rescaled separately: one unit = 0.01 of the index.
# key is read as text so codes like "0101" keep their leading zero.
# totalvotes uses -1 as a missing code (FL-20, OK-3); it is reference only and
# not used in the models, but is set to NA so it cannot leak into summaries.

analysis <- read_csv("analysis_2024.csv", show_col_types = FALSE,
                     col_types = cols(key = col_character(), state_po = col_character(),
                                      contested = col_logical(), .default = col_guess()),
                     na = c("", "NA")) %>%
  mutate(totalvotes = na_if(totalvotes, -1),
         across(c(turnout, dshare, starts_with("pct_")), ~ .x * 100),
         gini = gini * 100)

stopifnot(nrow(analysis) == 435,
          all(c("pct_female", "pct_low_inc", "pct_high_inc", "gini", "pct_poverty") %in% names(analysis)))

model_data <- analysis %>%
  filter(contested, !is.na(turnout)) %>%          # 396 districts
  mutate(z_turnout = as.numeric(scale(turnout)),  # standardized outcomes for
         z_dshare  = as.numeric(scale(dshare)))   # the cross-outcome tests

predictors <- c("pct_black", "pct_hisp", "pct_asian", "pct_ba", "inc_10k",
                "pct_18_29", "pct_65plus", "log_density", "gini", "pct_female",
                "inc_dem", "inc_rep")
# Open seats are the incumbency baseline. pct_white_nh is left out because
# it is nearly a linear combination of the other race shares.
# New fields: gini (inequality) and pct_female (share of CVAP that is female)
# are added. pct_female has little spread across districts (SD about 1.2
# points), so expect an imprecise coefficient. pct_low_inc, pct_high_inc and
# pct_poverty are left out: pct_high_inc correlates 0.98 with inc_10k, and
# pct_low_inc correlates 0.94 with pct_poverty, so together they push VIFs
# past 25-100. pct_poverty is tried in place of gini in step 11.
excluded_income <- c("pct_low_inc", "pct_high_inc", "pct_poverty")

# Plain-English names for every variable that appears in an exhibit.
#   label:      used in tables (includes units)
#   short:      used in figure panels, where space is tight
#   definition: used in the variable-definitions table
# CVAP = citizen voting-age population.
var_labels <- tribble(
  ~variable,      ~label,                                   ~short,                        ~definition,
  "turnout",      "Turnout (%)",                            "Turnout",                     "House votes cast as a share of CVAP",
  "dshare",       "Democratic two-party share (%)",         "Democratic share",            "Democratic votes / (Democratic + Republican votes)",
  "pct_black",    "Black share (%)",                        "Black share",                 "Share of CVAP that is Black",
  "pct_hisp",     "Hispanic share (%)",                     "Hispanic share",              "Share of CVAP that is Hispanic",
  "pct_asian",    "Asian share (%)",                        "Asian share",                 "Share of CVAP that is Asian",
  "pct_ba",       "Bachelor's degree or higher (%)",        "Bachelor's degree or higher", "Share of adults with a bachelor's degree or higher",
  "inc_10k",      "Median household income ($10,000s)",     "Median household income",     "Median household income, in $10,000s",
  "gini",         "Income inequality (Gini index x 100)",   "Income inequality",           "Gini index of household income, times 100",
  "pct_18_29",    "Age 18-29 (%)",                          "Age 18-29",                   "Share of CVAP aged 18-29",
  "pct_65plus",   "Age 65 and over (%)",                    "Age 65 and over",             "Share of CVAP aged 65 and over",
  "pct_female",   "Female share (%)",                       "Female share",                "Share of CVAP that is female",
  "log_density",  "Population density (log10)",            "Population density",          "log10 of persons per square mile",
  "inc_dem",      "Democratic incumbent",                   "Democratic incumbent",        "1 if a Democratic incumbent ran (open seat = baseline)",
  "inc_rep",      "Republican incumbent",                   "Republican incumbent",        "1 if a Republican incumbent ran (open seat = baseline)",
  "pct_poverty",  "Poverty rate (%)",                       "Poverty rate",                "Share of people below the poverty line",
  "pct_low_inc",  "Low-income households (%)",              "Low-income households",       "Share of households with low income",
  "pct_high_inc", "High-income households (%)",             "High-income households",      "Share of households with high income"
)
# Look up a label for a vector of variable codes, keeping their order
lab <- function(x, type = "label") var_labels[[type]][match(x, var_labels$variable)]

rhs <- paste(predictors, collapse = " + ")
f_turnout <- as.formula(paste("turnout ~", rhs))
f_dshare  <- as.formula(paste("dshare ~",  rhs))
f_mv      <- as.formula(paste("cbind(turnout, dshare) ~", rhs))
f_mv_z    <- as.formula(paste("cbind(z_turnout, z_dshare) ~", rhs))

cat("Modeling sample:", nrow(model_data), "districts\n")

# -----------------------------------------------------------------------------
# 2. Descriptives and collinearity check
# -----------------------------------------------------------------------------

descriptives <- model_data %>%
  select(turnout, dshare, all_of(predictors)) %>%
  pivot_longer(everything(), names_to = "variable") %>%
  group_by(variable) %>%
  summarise(mean = mean(value), sd = sd(value), min = min(value),
            median = median(value), max = max(value), .groups = "drop")

outcome_cor <- cor(model_data$turnout, model_data$dshare)
cat("Raw correlation, turnout vs. Democratic share:", round(outcome_cor, 3), "\n")

# VIF depends only on the predictors, so one model covers both outcomes.
vif_table <- tibble(predictor = predictors,
                    vif = unname(vif(lm(f_turnout, data = model_data))))
print(vif_table)   # above ~5 deserves a sentence in the report; above 10, rethink
# Expect inc_10k and pct_ba near 6-7 with gini in the model; all below 10.

# Candidate model with every new field, to document why three were left out
f_candidate <- as.formula(paste("turnout ~",
                                paste(c(predictors, excluded_income), collapse = " + ")))
vif_candidate <- tibble(predictor = c(predictors, excluded_income),
                        vif = unname(vif(lm(f_candidate, data = model_data))),
                        in_model = predictor %in% predictors)
print(vif_candidate)

# Correlations among the income measures
income_cor <- cor(model_data[, c("inc_10k", "pct_ba", "gini", excluded_income)])
print(round(income_cor, 2))

# -----------------------------------------------------------------------------
# 3. Univariate models: two separate multiple regressions
# -----------------------------------------------------------------------------

fit_turnout <- lm(f_turnout, data = model_data)
fit_dshare  <- lm(f_dshare,  data = model_data)

# Classical and state-clustered standard errors. Districts in the same state
# share election rules and statewide races, so their errors may be correlated.
tidy_uni <- function(fit, outcome) {
  cl <- coeftest(fit, vcov = vcovCL(fit, cluster = ~ state_po))
  tidy(fit) %>%
    mutate(outcome = outcome,
           se_cluster = cl[, "Std. Error"],
           p_cluster  = cl[, "Pr(>|t|)"])
}

uni_coefs <- bind_rows(tidy_uni(fit_turnout, "turnout"),
                       tidy_uni(fit_dshare,  "dshare"))

uni_fit <- bind_rows(glance(fit_turnout) %>% mutate(outcome = "turnout"),
                     glance(fit_dshare)  %>% mutate(outcome = "dshare")) %>%
  select(outcome, r.squared, adj.r.squared, F = statistic, p.value, df, df.residual, nobs)

# -----------------------------------------------------------------------------
# 4. Multivariate model: one regression with both outcomes
# -----------------------------------------------------------------------------

fit_mv <- lm(f_mv, data = model_data)

# Check: with the same predictors in both equations, the coefficients are
# identical to the separate regressions. Any difference between the approaches
# comes from inference, not estimation.
same_coefs <- isTRUE(all.equal(unname(coef(fit_mv)[, "turnout"]), unname(coef(fit_turnout)))) &&
  isTRUE(all.equal(unname(coef(fit_mv)[, "dshare"]),  unname(coef(fit_dshare))))
cat("Multivariate coefficients equal separate OLS coefficients:", same_coefs, "\n")

# -----------------------------------------------------------------------------
# 5. Significance tests: Bonferroni t-tests vs. Pillai's trace
# -----------------------------------------------------------------------------

# Pillai's trace and its F approximation from a car::linearHypothesis() result
pillai_from_lh <- function(lh) {
  H <- lh$SSPH; E <- lh$SSPE
  p <- ncol(E); q <- lh$df; dfe <- lh$df.residual
  V <- sum(diag(H %*% solve(H + E)))
  s <- min(p, q); m <- (abs(p - q) - 1) / 2; n <- (dfe - p - 1) / 2
  F <- ((2 * n + s + 1) / (2 * m + s + 1)) * V / (s - V)
  df1 <- s * (2 * m + s + 1); df2 <- s * (2 * n + s + 1)
  tibble(pillai = V, F = F, df1 = df1, df2 = df2,
         p_value = pf(F, df1, df2, lower.tail = FALSE))
}

mv_tests <- map_dfr(predictors, ~ pillai_from_lh(linearHypothesis(fit_mv, .x)) %>%
                      mutate(term = .x, .before = 1))
# Cross-check against car's own table: these match Anova(fit_mv, test = "Pillai")
mv_anova <- Anova(fit_mv, type = "II", test.statistic = "Pillai")

alpha <- 0.05
test_comparison <- uni_coefs %>%
  filter(term %in% predictors) %>%
  select(term, outcome, p.value) %>%
  pivot_wider(names_from = outcome, values_from = p.value, names_prefix = "p_") %>%
  left_join(mv_tests %>% select(term, pillai, p_pillai = p_value), by = "term") %>%
  mutate(
    sig_any_uni_bonf = p_turnout < alpha / 2 | p_dshare < alpha / 2,  # 2 tests per predictor
    sig_mv           = p_pillai < alpha,
    agreement = case_when(
      sig_any_uni_bonf & sig_mv    ~ "Both significant",
      !sig_any_uni_bonf & !sig_mv  ~ "Neither significant",
      sig_mv                       ~ "Joint only",
      TRUE                         ~ "Separate only"))
print(test_comparison)

# Overall test of all predictors: two univariate F-tests vs. one Pillai test
overall_mv <- pillai_from_lh(linearHypothesis(fit_mv, predictors))

# -----------------------------------------------------------------------------
# 6. Cross-outcome tests: same effect on turnout and share?
# -----------------------------------------------------------------------------
# On standardized outcomes, test H0: b_turnout = b_dshare for each predictor
# (response contrast P = (1, -1)). Rejecting means the predictor moves turnout
# and Democratic share by different amounts, in SD units. To compare with
# Republican share instead, use P = (1, 1).

fit_mv_z <- lm(f_mv_z, data = model_data)

cross_test <- function(term, contrast) {
  lh <- linearHypothesis(fit_mv_z, term, P = matrix(contrast, ncol = 1))
  F  <- (lh$SSPH[1, 1] / lh$df) / (lh$SSPE[1, 1] / lh$df.residual)
  tibble(F = F, p_value = pf(F, lh$df, lh$df.residual, lower.tail = FALSE))
}

cross_outcome <- tibble(term = predictors) %>%
  mutate(b_z_turnout = coef(fit_mv_z)[term, "z_turnout"],
         b_z_dshare  = coef(fit_mv_z)[term, "z_dshare"],
         equal = map(term, cross_test, contrast = c(1, -1))) %>%
  unnest(equal)
print(cross_outcome)

# -----------------------------------------------------------------------------
# 7. Residual correlation between the outcomes
# -----------------------------------------------------------------------------
# How turnout and share move together after the predictors are accounted for.
# Its size and sign govern how much the joint tests differ from separate ones.

res <- resid(fit_mv)
resid_cor_test <- cor.test(res[, "turnout"], res[, "dshare"])
resid_cor <- unname(resid_cor_test$estimate)
cat("Residual correlation:", round(resid_cor, 3),
    "| p =", signif(resid_cor_test$p.value, 3), "\n")

# -----------------------------------------------------------------------------
# 8. Joint confidence ellipses vs. Bonferroni rectangles
# -----------------------------------------------------------------------------
# For selected predictors, the 95% joint region for (b_turnout, b_dshare) from
# the multivariate model vs. the Bonferroni rectangle from two separate 97.5%
# intervals. Standardized outcomes put both axes on the same scale.

ellipse_terms <- c("pct_ba", "pct_hisp", "pct_65plus")
V_all <- vcov(fit_mv_z)
dfe   <- df.residual(fit_mv_z)

joint_region <- function(term, level = 0.95, n_pts = 200) {
  idx <- paste0(c("z_turnout:", "z_dshare:"), term)
  b <- coef(fit_mv_z)[term, ]
  S <- V_all[idx, idx]
  radius <- sqrt(2 * dfe / (dfe - 1) * qf(level, 2, dfe - 1))   # Hotelling T-squared, p = 2
  theta <- seq(0, 2 * pi, length.out = n_pts)
  pts <- t(b + radius * t(chol(S)) %*% rbind(cos(theta), sin(theta)))
  tibble(term = term, x = pts[, 1], y = pts[, 2])
}

bonf_rect <- map_dfr(ellipse_terms, function(term) {
  idx <- paste0(c("z_turnout:", "z_dshare:"), term)
  b <- coef(fit_mv_z)[term, ]; se <- sqrt(diag(V_all[idx, idx]))
  tq <- qt(1 - 0.05 / 4, dfe)        # 97.5% two-sided interval each
  tibble(term = term, b_t = b[1], b_s = b[2],
         xmin = b[1] - tq * se[1], xmax = b[1] + tq * se[1],
         ymin = b[2] - tq * se[2], ymax = b[2] + tq * se[2])
})

ellipses <- map_dfr(ellipse_terms, joint_region)

# Plain-English panel names, in the same panel order as the draft report
panel_levels <- lab(sort(ellipse_terms), "short")
ellipses  <- ellipses  %>% mutate(term = factor(lab(term, "short"), levels = panel_levels))
bonf_rect <- bonf_rect %>% mutate(term = factor(lab(term, "short"), levels = panel_levels))

p_ellipse <- ggplot() +
  geom_rect(data = bonf_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "grey45", linetype = "dashed") +
  geom_path(data = ellipses, aes(x, y), linewidth = 0.8) +
  geom_point(data = bonf_rect, aes(b_t, b_s)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  facet_wrap(~ term, scales = "free") +
  labs(x = "Effect on turnout (SD units)", y = "Effect on Democratic share (SD units)",
       title = "95% joint confidence region (solid) vs. Bonferroni rectangle (dashed)") +
  theme_minimal()
ggsave("figures/confidence_ellipses.png", p_ellipse, width = 9, height = 3.5, dpi = 300)

# -----------------------------------------------------------------------------
# 9. Missing-data cost of the joint model
# -----------------------------------------------------------------------------
# The joint model needs both outcomes, so it drops uncontested races. A
# turnout-only model can keep them (433 districts with valid vote counts).

turnout_full_data <- analysis %>% filter(!is.na(turnout))
fit_turnout_full  <- lm(f_turnout, data = turnout_full_data)

missing_cost <- full_join(
  tidy(fit_turnout)      %>% select(term, est_396 = estimate, se_396 = std.error),
  tidy(fit_turnout_full) %>% select(term, est_full = estimate, se_full = std.error),
  by = "term") %>%
  mutate(change = est_full - est_396)
cat("Turnout model N: contested only =", nobs(fit_turnout),
    "| all valid =", nobs(fit_turnout_full), "\n")

# -----------------------------------------------------------------------------
# 10. Assumption diagnostics
# -----------------------------------------------------------------------------

diag_data <- model_data %>%
  select(key, state_po, district) %>%
  mutate(fitted_turnout = fitted(fit_mv)[, "turnout"], res_turnout = res[, "turnout"],
         fitted_dshare  = fitted(fit_mv)[, "dshare"],  res_dshare  = res[, "dshare"],
         rstd_turnout = rstandard(fit_turnout), rstd_dshare = rstandard(fit_dshare),
         mahal = mahalanobis(res, center = colMeans(res), cov = cov(res)))

# Univariate: residuals vs. fitted and normal Q-Q, per outcome
p_resid <- diag_data %>%
  pivot_longer(c(fitted_turnout, fitted_dshare, res_turnout, res_dshare),
               names_to = c(".value", "outcome"), names_sep = "_") %>%
  mutate(outcome = factor(lab(outcome, "short"), levels = lab(c("turnout", "dshare"), "short"))) %>%
  ggplot(aes(fitted, res)) + geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~ outcome, scales = "free") +
  labs(x = "Fitted value (percentage points)", y = "Residual (percentage points)") +
  theme_minimal()
ggsave("figures/residuals_vs_fitted.png", p_resid, width = 8, height = 3.5, dpi = 300)

p_qq <- diag_data %>%
  pivot_longer(c(rstd_turnout, rstd_dshare), names_to = "outcome", values_to = "r") %>%
  mutate(outcome = factor(lab(sub("rstd_", "", outcome), "short"),
                          levels = lab(c("turnout", "dshare"), "short"))) %>%
  ggplot(aes(sample = r)) + stat_qq(alpha = 0.5) + stat_qq_line(colour = "red") +
  facet_wrap(~ outcome) + labs(x = "Theoretical quantile", y = "Standardized residual") +
  theme_minimal()
ggsave("figures/qq_univariate.png", p_qq, width = 8, height = 3.5, dpi = 300)

# Multivariate: squared Mahalanobis distances of the residual pairs should
# follow a chi-square distribution with 2 df if the errors are bivariate normal.
p_mahal <- diag_data %>%
  arrange(mahal) %>%
  mutate(theoretical = qchisq(ppoints(n()), df = 2)) %>%
  ggplot(aes(theoretical, mahal)) + geom_point(alpha = 0.5) +
  geom_abline(colour = "red") +
  labs(x = expression(chi^2 * "(2) quantile"), y = "Squared Mahalanobis distance",
       title = "Bivariate normality check of residuals") + theme_minimal()
ggsave("figures/qq_mahalanobis.png", p_mahal, width = 5, height = 4, dpi = 300)

# Districts unusual jointly but not on either outcome alone
outliers <- diag_data %>%
  mutate(mv_outlier  = mahal > qchisq(0.975, 2),
         uni_outlier = abs(rstd_turnout) > 2.5 | abs(rstd_dshare) > 2.5) %>%
  filter(mv_outlier | uni_outlier) %>%
  arrange(desc(mahal)) %>%
  select(key, state_po, district, mahal, rstd_turnout, rstd_dshare, mv_outlier, uni_outlier)

# -----------------------------------------------------------------------------
# 11. Sensitivity check: poverty rate in place of the Gini index
# -----------------------------------------------------------------------------
# Both describe the economic side of a district that median income misses.
# If the Pillai results for the other predictors hold up when pct_poverty
# replaces gini, the conclusions don't hinge on which one was chosen.

predictors_alt <- c(setdiff(predictors, "gini"), "pct_poverty")
f_mv_alt <- as.formula(paste("cbind(turnout, dshare) ~",
                             paste(predictors_alt, collapse = " + ")))
fit_mv_alt <- lm(f_mv_alt, data = model_data)

vif_alt <- tibble(predictor = predictors_alt,
                  vif = unname(vif(lm(as.formula(paste("turnout ~", paste(predictors_alt, collapse = " + "))),
                                      data = model_data))))
print(vif_alt)   # expect pct_poverty, inc_10k and pct_ba around 7-8

sensitivity <- map_dfr(predictors_alt, ~ pillai_from_lh(linearHypothesis(fit_mv_alt, .x)) %>%
                         mutate(term = .x, .before = 1)) %>%
  select(term, pillai_alt = pillai, p_alt = p_value) %>%
  full_join(mv_tests %>% select(term, pillai_main = pillai, p_main = p_value), by = "term") %>%
  mutate(sig_main = p_main < alpha, sig_alt = p_alt < alpha,
         same_conclusion = sig_main == sig_alt)
print(sensitivity)

# -----------------------------------------------------------------------------
# 12. Report exhibits
# -----------------------------------------------------------------------------
# Each table below is the finished version of an exhibit in the draft report:
# rows in report order, plain-English variable names, numbers already rounded
# and formatted as text. The column names are the report's column headers, so
# in Quarto read them with check.names = FALSE and pass straight to kable().

fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
fmt_b_se <- function(b, se) sprintf("%.3f (%.3f)", b, se)
fmt_F <- function(F) ifelse(F < 1, sprintf("%.2f", F), sprintf("%.1f", F))
fmt_num <- function(x, digits) sprintf(paste0("%.", digits, "f"), x)  # digits may vary by row

# Predictor order used in the report's tables (differs from the model formula)
report_order <- c("pct_black", "pct_hisp", "pct_asian", "pct_ba", "inc_10k", "gini",
                  "pct_18_29", "pct_65plus", "pct_female", "log_density",
                  "inc_dem", "inc_rep")
stopifnot(setequal(report_order, predictors))

# Variables shown in the descriptives table (incumbency is summarized as counts)
desc_vars <- c("turnout", "dshare", "pct_black", "pct_hisp", "pct_asian", "pct_ba",
               "inc_10k", "gini", "pct_18_29", "pct_65plus", "pct_female", "log_density")

# Section 2: variable definitions
ex_definitions <- var_labels %>%
  filter(variable %in% c("turnout", "dshare", predictors)) %>%
  select(Variable = label, Definition = definition)

# Section 2 / Table 1: descriptive statistics (n = 396)
ex_descriptives <- descriptives %>%
  filter(variable %in% desc_vars) %>%
  mutate(digits = if_else(variable %in% c("inc_10k", "log_density"), 2, 1),
         across(c(mean, sd, min, median, max), ~ fmt_num(.x, digits))) %>%
  arrange(match(variable, desc_vars)) %>%
  transmute(Variable = lab(variable), Mean = mean, SD = sd, Min = min,
            Median = median, Max = max)

# Incumbency counts quoted in Section 2
incumbency_counts <- model_data %>%
  summarise(`Democratic incumbent` = sum(inc_dem), `Republican incumbent` = sum(inc_rep),
            `Open seat` = sum(inc_dem == 0 & inc_rep == 0))
print(incumbency_counts)

# Section 4.1: overall fit of each regression and the overall Pillai test
ex_model_fit <- uni_fit %>%
  mutate(Outcome = lab(outcome)) %>%
  transmute(Outcome, `R-squared` = sprintf("%.3f", r.squared),
            `Adjusted R-squared` = sprintf("%.3f", adj.r.squared),
            F = sprintf("F(%d, %d) = %.1f", as.integer(df), as.integer(df.residual), F),
            p = fmt_p(p.value)) %>%
  bind_rows(tibble(Outcome = "Both outcomes jointly (Pillai's trace)",
                   `R-squared` = "", `Adjusted R-squared` = "",
                   F = sprintf("F(%d, %d) = %.1f", as.integer(overall_mv$df1),
                               as.integer(overall_mv$df2), overall_mv$F),
                   p = fmt_p(overall_mv$p_value)))
ex_model_fit$`Pillai's trace` <- c("", "", sprintf("%.2f", overall_mv$pillai))

# Section 4.2: coefficients (classical SEs) and joint Pillai tests
ex_coefficients <- uni_coefs %>%
  filter(term %in% predictors) %>%
  transmute(term, outcome, b_se = fmt_b_se(estimate, std.error), p = fmt_p(p.value)) %>%
  pivot_wider(names_from = outcome, values_from = c(b_se, p)) %>%
  left_join(mv_tests %>% transmute(term, pillai = sprintf("%.3f", pillai),
                                   joint_p = fmt_p(p_value)), by = "term") %>%
  arrange(match(term, report_order)) %>%
  transmute(Predictor = lab(term),
            `Turnout b (SE)` = b_se_turnout, `Turnout p` = p_turnout,
            `Dem. share b (SE)` = b_se_dshare, `Dem. share p` = p_dshare,
            Pillai = pillai, `Joint p` = joint_p)

# Section 4.3: joint vs. separate (Bonferroni) conclusions for each predictor
ex_test_comparison <- test_comparison %>%
  arrange(match(term, report_order)) %>%
  transmute(Predictor = lab(term),
            `Turnout p` = fmt_p(p_turnout), `Dem. share p` = fmt_p(p_dshare),
            `Significant separately (Bonferroni 0.025)` = if_else(sig_any_uni_bonf, "Yes", "No"),
            `Joint p` = fmt_p(p_pillai),
            `Significant jointly (0.05)` = if_else(sig_mv, "Yes", "No"),
            Agreement = agreement)

# Section 4.4: equal-effect tests on standardized outcomes
cross_df <- as.integer(df.residual(fit_mv_z))
ex_cross_outcome <- cross_outcome %>%
  arrange(match(term, report_order)) %>%
  transmute(Predictor = lab(term),
            `Turnout (SD)` = sprintf("%.3f", b_z_turnout),
            `Dem. share (SD)` = sprintf("%.3f", b_z_dshare),
            F = fmt_F(F), p = fmt_p(p_value)) %>%
  rename(!!sprintf("F(1, %d)", cross_df) := F)

# Section 4.6: state-clustered standard errors
ex_clustered <- uni_coefs %>%
  filter(term %in% predictors) %>%
  transmute(term, outcome, b_se = fmt_b_se(estimate, se_cluster), p = fmt_p(p_cluster)) %>%
  pivot_wider(names_from = outcome, values_from = c(b_se, p)) %>%
  arrange(match(term, report_order)) %>%
  transmute(Predictor = lab(term),
            `Turnout b (clustered SE)` = b_se_turnout, `Turnout p` = p_turnout,
            `Dem. share b (clustered SE)` = b_se_dshare, `Dem. share p` = p_dshare)

# Section 4.6: poverty rate in place of the Gini index
ex_sensitivity <- sensitivity %>%
  arrange(match(term, c(report_order, "pct_poverty"))) %>%
  transmute(Predictor = lab(term),
            `Pillai (main model)` = if_else(is.na(pillai_main), "", sprintf("%.3f", pillai_main)),
            `Joint p (main model)` = if_else(is.na(p_main), "", fmt_p(p_main)),
            `Pillai (poverty model)` = if_else(is.na(pillai_alt), "", sprintf("%.3f", pillai_alt)),
            `Joint p (poverty model)` = if_else(is.na(p_alt), "", fmt_p(p_alt)))

# Section 4.6: missing-data cost (turnout model, 396 vs. 433 districts)
ex_missing_cost <- missing_cost %>%
  filter(term %in% predictors) %>%
  arrange(match(term, report_order)) %>%
  transmute(Predictor = lab(term),
            `b, contested (n = 396)` = sprintf("%.3f", est_396),
            `b, all valid (n = 433)` = sprintf("%.3f", est_full),
            Change = sprintf("%+.3f", change))

# Section 5: districts beyond the 97.5% chi-square(2) cutoff
outlier_note <- function(rt, rs) {
  part <- function(r, name) case_when(abs(r) > 3.5 & r > 0 ~ paste(name, "far above prediction"),
                                      abs(r) > 3.5         ~ paste(name, "far below prediction"),
                                      r >  2.5             ~ paste(name, "above prediction"),
                                      r < -2.5             ~ paste(name, "below prediction"),
                                      TRUE                 ~ NA_character_)
  notes <- map2_chr(part(rt, "Turnout"), part(rs, "Share"),
                    ~ paste(na.omit(c(.x, .y)), collapse = "; "))
  if_else(notes == "", "Joint outlier only", notes)
}

ex_outliers <- outliers %>%
  filter(mv_outlier) %>%
  arrange(desc(mahal)) %>%
  transmute(District = paste0(state_po, "-", if_else(district == 0, "AL", as.character(district))),
            `Mahalanobis D²` = sprintf("%.1f", mahal),
            `Turnout resid.` = sprintf("%.2f", rstd_turnout),
            `Dem. share resid.` = sprintf("%.2f", rstd_dshare),
            Note = outlier_note(rstd_turnout, rstd_dshare))

# Section 5: collinearity (final model and the candidate model with all income measures)
ex_vif <- vif_candidate %>%
  arrange(desc(in_model), match(predictor, c(report_order, excluded_income))) %>%
  left_join(vif_table %>% rename(vif_final = vif), by = "predictor") %>%
  transmute(Predictor = lab(predictor),
            `In final model` = if_else(in_model, "Yes", "No"),
            `VIF, final model` = if_else(is.na(vif_final), "", sprintf("%.1f", vif_final)),
            `VIF, all income measures` = sprintf("%.1f", vif))

exhibits <- list(
  definitions      = ex_definitions,
  descriptives     = ex_descriptives,
  model_fit        = ex_model_fit,
  coefficients     = ex_coefficients,
  test_comparison  = ex_test_comparison,
  cross_outcome    = ex_cross_outcome,
  clustered_se     = ex_clustered,
  sensitivity      = ex_sensitivity,
  missing_cost     = ex_missing_cost,
  outliers         = ex_outliers,
  vif              = ex_vif
)
iwalk(exhibits, ~ write_csv(.x, file.path("output/exhibits", paste0(.y, ".csv"))))

# -----------------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------------

write_csv(descriptives,    "output/descriptives.csv")
write_csv(vif_table,       "output/vif.csv")
write_csv(uni_coefs,       "output/univariate_coefficients.csv")
write_csv(uni_fit,         "output/univariate_fit.csv")
write_csv(mv_tests,        "output/multivariate_pillai_tests.csv")
write_csv(test_comparison, "output/test_comparison.csv")
write_csv(cross_outcome,   "output/cross_outcome_tests.csv")
write_csv(missing_cost,    "output/missing_data_cost.csv")
write_csv(outliers,        "output/outliers.csv")
write_csv(vif_candidate,   "output/vif_candidate.csv")
write_csv(sensitivity,     "output/sensitivity_poverty.csv")

saveRDS(list(
  n = nrow(model_data), predictors = predictors, excluded_income = excluded_income,
  descriptives = descriptives, outcome_cor = outcome_cor, vif = vif_table,
  vif_candidate = vif_candidate, income_cor = income_cor,
  fit_mv_alt = fit_mv_alt, vif_alt = vif_alt, sensitivity = sensitivity,
  fit_turnout = fit_turnout, fit_dshare = fit_dshare, fit_mv = fit_mv, fit_mv_z = fit_mv_z,
  uni_coefs = uni_coefs, uni_fit = uni_fit, same_coefs = same_coefs,
  mv_tests = mv_tests, mv_anova = mv_anova, overall_mv = overall_mv,
  test_comparison = test_comparison, cross_outcome = cross_outcome,
  resid_cor = resid_cor, resid_cor_test = resid_cor_test,
  missing_cost = missing_cost, outliers = outliers,
  var_labels = var_labels, incumbency_counts = incumbency_counts,
  exhibits = exhibits
), "output/results.rds")

cat("Done. Tables in output/, report exhibits in output/exhibits/, figures in figures/.\n")