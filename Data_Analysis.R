# Jan Blackburn
# Math 257: Multivariate Analysis Project
# Data analysis file: reads the analysis data set and produces every table,
# figure and number used in Report_Draft_v3.qmd
# Initial code was generated using Claude and was modified by me as needed
#
# Outcomes (per House district, 2024): turnout and Democratic two-party share
#
# Input:   analysis_2024.csv        (from build_analysis_data.R)
# Output:  output/exhibits/*.csv    the report's three tables, formatted
#          output/*.csv             unformatted results behind the tables and text
#          figures/*.png            the report's three Appendix B figures
#          console                  the other numbers quoted in the report text
#
# Steps
#   1. Prepare the modeling sample
#   2. Descriptives and the collinearity check behind the excluded measures
#   3. Separate regressions and the multivariate regression
#   4. Joint tests: Pillai's trace vs. Bonferroni t-tests
#   5. Cross-outcome tests: same effect on turnout and share?
#   6. Residual correlation between the outcomes
#   7. Assumption diagnostics (Appendix B figures)
#   8. Report tables
#   9. Numbers quoted in the report text
# =============================================================================

library(tidyverse)   # data handling and ggplot2
library(car)         # vif(), linearHypothesis()
library(broom)       # tidy() / glance() model output into data frames

rm(list = ls())

dir.create("output/exhibits", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Prepare the modeling sample
# -----------------------------------------------------------------------------
# Shares are rescaled from proportions (0-1) to percentage points (0-100) so a
# coefficient reads as "points of outcome per 1 point of predictor". The Gini
# index is rescaled the same way: one unit = 0.01 of the index.
# key is read as text so codes like "0101" keep their leading zero.

analysis <- read_csv("analysis_2024.csv", show_col_types = FALSE,
                     col_types = cols(key = col_character(), state_po = col_character(),
                                      contested = col_logical(), .default = col_guess()),
                     na = c("", "NA")) %>%
  mutate(across(c(turnout, dshare, starts_with("pct_")), ~ .x * 100),
         gini = gini * 100)

stopifnot(nrow(analysis) == 435)

model_data <- analysis %>%
  filter(contested, !is.na(turnout)) %>%          # 396 districts
  mutate(z_turnout = as.numeric(scale(turnout)),  # standardized outcomes for
         z_dshare  = as.numeric(scale(dshare)))   # the cross-outcome tests

# Open seats are the incumbency baseline. pct_white_nh is left out because it
# is nearly a linear combination of the other race shares; the three income
# measures below are left out for collinearity with median income (step 2).
predictors <- c("pct_black", "pct_hisp", "pct_asian", "pct_ba", "inc_10k",
                "pct_18_29", "pct_65plus", "log_density", "gini", "pct_female",
                "inc_dem", "inc_rep")
excluded_income <- c("pct_low_inc", "pct_high_inc", "pct_poverty")

# Predictor order used in the report's tables (differs from the model formula)
report_order <- c("pct_black", "pct_hisp", "pct_asian", "pct_ba", "inc_10k", "gini",
                  "pct_18_29", "pct_65plus", "pct_female", "log_density",
                  "inc_dem", "inc_rep")
stopifnot(setequal(report_order, predictors))

# Plain-English names used in the report's tables
var_labels <- c(
  turnout     = "Turnout (%)",
  dshare      = "Democratic two-party share (%)",
  pct_black   = "Black share (%)",
  pct_hisp    = "Hispanic share (%)",
  pct_asian   = "Asian share (%)",
  pct_ba      = "Bachelor's degree or higher (%)",
  inc_10k     = "Median household income ($10,000s)",
  gini        = "Income inequality (Gini index x 100)",
  pct_18_29   = "Age 18-29 (%)",
  pct_65plus  = "Age 65 and over (%)",
  pct_female  = "Female share (%)",
  log_density = "Population density (log10)",
  inc_dem     = "Democratic incumbent",
  inc_rep     = "Republican incumbent"
)
lab <- function(x) unname(var_labels[x])

# Short outcome names for figure panels, in panel order
outcome_names <- c(turnout = "Turnout", dshare = "Democratic share")

rhs <- paste(predictors, collapse = " + ")
f_turnout <- as.formula(paste("turnout ~", rhs))
f_dshare  <- as.formula(paste("dshare ~",  rhs))
f_mv      <- as.formula(paste("cbind(turnout, dshare) ~", rhs))
f_mv_z    <- as.formula(paste("cbind(z_turnout, z_dshare) ~", rhs))

# -----------------------------------------------------------------------------
# 2. Descriptives and collinearity check
# -----------------------------------------------------------------------------

descriptives <- model_data %>%
  select(turnout, dshare, all_of(predictors)) %>%
  pivot_longer(everything(), names_to = "variable") %>%
  group_by(variable) %>%
  summarise(mean = mean(value), sd = sd(value), min = min(value),
            median = median(value), max = max(value), .groups = "drop")

incumbency_counts <- model_data %>%
  summarise(`Democratic incumbent` = sum(inc_dem), `Republican incumbent` = sum(inc_rep),
            `Open seat` = sum(inc_dem == 0 & inc_rep == 0))

outcome_cor <- cor(model_data$turnout, model_data$dshare)

# Supports the Section 2 statement that the three income measures were left
# out: with them added, VIFs climb well past 10.
# VIF depends only on the predictors, so one model covers both outcomes.
vif_table <- tibble(predictor = predictors,
                    vif = unname(vif(lm(f_turnout, data = model_data))))

f_candidate <- as.formula(paste("turnout ~",
                                paste(c(predictors, excluded_income), collapse = " + ")))
vif_candidate <- tibble(predictor = c(predictors, excluded_income),
                        vif = unname(vif(lm(f_candidate, data = model_data))),
                        in_model = predictor %in% predictors)

# -----------------------------------------------------------------------------
# 3. Separate regressions and the multivariate regression
# -----------------------------------------------------------------------------

fit_turnout <- lm(f_turnout, data = model_data)
fit_dshare  <- lm(f_dshare,  data = model_data)
fit_mv      <- lm(f_mv,      data = model_data)

# Section 3.2: with the same predictors in both equations, the multivariate
# coefficients are identical to the separate OLS coefficients.
stopifnot(isTRUE(all.equal(unname(coef(fit_mv)[, "turnout"]), unname(coef(fit_turnout)))),
          isTRUE(all.equal(unname(coef(fit_mv)[, "dshare"]),  unname(coef(fit_dshare)))))

uni_coefs <- bind_rows(tidy(fit_turnout) %>% mutate(outcome = "turnout"),
                       tidy(fit_dshare)  %>% mutate(outcome = "dshare"))

uni_fit <- bind_rows(glance(fit_turnout) %>% mutate(outcome = "turnout"),
                     glance(fit_dshare)  %>% mutate(outcome = "dshare")) %>%
  select(outcome, r.squared, adj.r.squared, F = statistic, p.value, df, df.residual, nobs)

# -----------------------------------------------------------------------------
# 4. Joint tests: Pillai's trace vs. Bonferroni t-tests
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

# Overall test: all twelve predictors have no effect on either outcome
overall_mv <- pillai_from_lh(linearHypothesis(fit_mv, predictors))

# One joint test per predictor: H0 (b_turnout, b_dshare) = (0, 0)
mv_tests <- map_dfr(predictors, ~ pillai_from_lh(linearHypothesis(fit_mv, .x)) %>%
                      mutate(term = .x, .before = 1))

# Joint test at 0.05 vs. two separate t-tests at the Bonferroni 0.025
alpha <- 0.05
test_comparison <- uni_coefs %>%
  filter(term %in% predictors) %>%
  select(term, outcome, p.value) %>%
  pivot_wider(names_from = outcome, values_from = p.value, names_prefix = "p_") %>%
  left_join(mv_tests %>% select(term, pillai, p_pillai = p_value), by = "term") %>%
  mutate(sig_any_uni_bonf = p_turnout < alpha / 2 | p_dshare < alpha / 2,  # 2 tests per predictor
         sig_mv           = p_pillai < alpha,
         agreement = case_when(
           sig_any_uni_bonf & sig_mv   ~ "Both significant",
           !sig_any_uni_bonf & !sig_mv ~ "Neither significant",
           sig_mv                      ~ "Joint only",
           TRUE                        ~ "Separate only"))

# -----------------------------------------------------------------------------
# 5. Cross-outcome tests: same effect on turnout and share?
# -----------------------------------------------------------------------------
# On standardized outcomes, test H0: b*_turnout = b*_dshare for each predictor
# with the response contrast c = (1, -1).

fit_mv_z <- lm(f_mv_z, data = model_data)

cross_test <- function(term) {
  lh <- linearHypothesis(fit_mv_z, term, P = matrix(c(1, -1), ncol = 1))
  F  <- (lh$SSPH[1, 1] / lh$df) / (lh$SSPE[1, 1] / lh$df.residual)
  tibble(F = F, p_value = pf(F, lh$df, lh$df.residual, lower.tail = FALSE))
}

cross_outcome <- tibble(term = predictors) %>%
  mutate(b_z_turnout = coef(fit_mv_z)[term, "z_turnout"],
         b_z_dshare  = coef(fit_mv_z)[term, "z_dshare"],
         test = map(term, cross_test)) %>%
  unnest(test)

# -----------------------------------------------------------------------------
# 6. Residual correlation between the outcomes
# -----------------------------------------------------------------------------

res <- resid(fit_mv)
resid_cor_test <- cor.test(res[, "turnout"], res[, "dshare"])

# -----------------------------------------------------------------------------
# 7. Assumption diagnostics (Appendix B figures)
# -----------------------------------------------------------------------------

diag_data <- tibble(
  fitted_turnout = fitted(fit_mv)[, "turnout"], res_turnout = res[, "turnout"],
  fitted_dshare  = fitted(fit_mv)[, "dshare"],  res_dshare  = res[, "dshare"],
  rstd_turnout = rstandard(fit_turnout), rstd_dshare = rstandard(fit_dshare),
  mahal = mahalanobis(res, center = colMeans(res), cov = cov(res)))

# Figure B1: residuals vs. fitted values
p_resid <- diag_data %>%
  pivot_longer(c(fitted_turnout, fitted_dshare, res_turnout, res_dshare),
               names_to = c(".value", "outcome"), names_sep = "_") %>%
  mutate(outcome = factor(outcome_names[outcome], levels = outcome_names)) %>%
  ggplot(aes(fitted, res)) + geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~ outcome, scales = "free") +
  labs(x = "Fitted value (percentage points)", y = "Residual (percentage points)") +
  theme_minimal()
ggsave("figures/residuals_vs_fitted.png", p_resid, width = 8, height = 3.5, dpi = 300)

# Figure B2: normal Q-Q plots of standardized residuals
p_qq <- diag_data %>%
  pivot_longer(c(rstd_turnout, rstd_dshare), names_to = "outcome", values_to = "r",
               names_prefix = "rstd_") %>%
  mutate(outcome = factor(outcome_names[outcome], levels = outcome_names)) %>%
  ggplot(aes(sample = r)) + stat_qq(alpha = 0.5) + stat_qq_line(colour = "red") +
  facet_wrap(~ outcome) + labs(x = "Theoretical quantile", y = "Standardized residual") +
  theme_minimal()
ggsave("figures/qq_univariate.png", p_qq, width = 8, height = 3.5, dpi = 300)

# Figure B3: squared Mahalanobis distances of the residual pairs should follow
# a chi-square distribution with 2 df if the errors are bivariate normal.
p_mahal <- diag_data %>%
  arrange(mahal) %>%
  mutate(theoretical = qchisq(ppoints(n()), df = 2)) %>%
  ggplot(aes(theoretical, mahal)) + geom_point(alpha = 0.5) +
  geom_abline(colour = "red") +
  labs(x = expression(chi^2 * "(2) quantile"), y = "Squared Mahalanobis distance",
       title = "Bivariate normality check of residuals") + theme_minimal()
ggsave("figures/qq_mahalanobis.png", p_mahal, width = 5, height = 4, dpi = 300)

# -----------------------------------------------------------------------------
# 8. Report tables
# -----------------------------------------------------------------------------
# Rows in report order, plain-English names, numbers formatted as text. The
# column names are the report's column headers; Quarto reads them with
# check.names = FALSE and passes them straight to kable().

fmt_p    <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
fmt_b_se <- function(b, se) sprintf("%.3f (%.3f)", b, se)
fmt_F    <- function(F) ifelse(F < 1, sprintf("%.2f", F), sprintf("%.1f", F))
fmt_num  <- function(x, digits) sprintf(paste0("%.", digits, "f"), x)  # digits may vary by row

# Table 1 (@tbl-descriptives): descriptive statistics (incumbency is reported
# as counts in the text instead)
desc_vars <- c("turnout", "dshare", setdiff(report_order, c("inc_dem", "inc_rep")))
ex_descriptives <- descriptives %>%
  filter(variable %in% desc_vars) %>%
  arrange(match(variable, desc_vars)) %>%
  mutate(digits = if_else(variable %in% c("inc_10k", "log_density"), 2, 1),
         across(c(mean, sd, min, median, max), ~ fmt_num(.x, digits))) %>%
  transmute(Variable = lab(variable), Mean = mean, SD = sd, Min = min,
            Median = median, Max = max)

# Table 2 (@tbl-coefficients): coefficients and joint Pillai tests
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

# Table 3 (@tbl-outcome-test): equal-effect tests on standardized outcomes
ex_cross_outcome <- cross_outcome %>%
  arrange(match(term, report_order)) %>%
  transmute(Predictor = lab(term),
            `Turnout (SD)` = sprintf("%.3f", b_z_turnout),
            `Dem. share (SD)` = sprintf("%.3f", b_z_dshare),
            F = fmt_F(F), p = fmt_p(p_value)) %>%
  rename(!!sprintf("F(1, %d)", as.integer(df.residual(fit_mv_z))) := F)

write_csv(ex_descriptives,  "output/exhibits/descriptives.csv")
write_csv(ex_coefficients,  "output/exhibits/coefficients.csv")
write_csv(ex_cross_outcome, "output/exhibits/cross_outcome.csv")

# Unformatted results behind the tables and text, same file names as before
write_csv(descriptives,    "output/descriptives.csv")
write_csv(vif_table,       "output/vif.csv")
write_csv(vif_candidate,   "output/vif_candidate.csv")
write_csv(uni_coefs,       "output/univariate_coefficients.csv")
write_csv(uni_fit,         "output/univariate_fit.csv")
write_csv(mv_tests,        "output/multivariate_pillai_tests.csv")
write_csv(test_comparison, "output/test_comparison.csv")
write_csv(cross_outcome,   "output/cross_outcome_tests.csv")

# -----------------------------------------------------------------------------
# 9. Numbers quoted in the report text
# -----------------------------------------------------------------------------
# The report types these numbers into the prose; check them here after any
# change to the data or model.

pct <- function(x) sprintf("%.1f%%", 100 * x)
fit_line <- function(o) {
  g <- uni_fit %>% filter(outcome == o)
  sprintf("%-16s R2 = %s, adj. R2 = %.3f, F(%d, %d) = %.1f, p = %s",
          outcome_names[o], pct(g$r.squared), g$adj.r.squared,
          as.integer(g$df), as.integer(g$df.residual), g$F, fmt_p(g$p.value))
}

cat("\n=== Section 2: Data ===\n")
cat("Modeling sample:", nrow(model_data), "districts;",
    nrow(analysis) - nrow(model_data), "excluded\n")
cat("Excluded districts by state (most first):\n")
print(analysis %>% anti_join(model_data, by = "key") %>% count(state_po, sort = TRUE))
print(incumbency_counts)
cat("Raw correlation, turnout vs. Democratic share:", sprintf("%.2f", outcome_cor), "\n")
cat("VIFs (final model, and with the three income measures added):\n")
print(vif_candidate %>% left_join(vif_table %>% rename(vif_final = vif), by = "predictor") %>%
        select(predictor, in_model, vif_final, vif_all = vif))

cat("\n=== Section 4.1: Model fit and joint tests ===\n")
cat(fit_line("turnout"), "\n")
cat(fit_line("dshare"), "\n")
cat(sprintf("Overall Pillai's trace = %.2f, F(%d, %d) = %.1f, p = %s\n",
            overall_mv$pillai, as.integer(overall_mv$df1), as.integer(overall_mv$df2),
            overall_mv$F, fmt_p(overall_mv$p_value)))
cat("Joint vs. separate (Bonferroni) conclusions:\n")
print(count(test_comparison, agreement))
print(test_comparison %>% filter(agreement %in% c("Joint only", "Separate only",
                                                  "Neither significant")))
cat(sprintf("Residual correlation: r = %.3f, p = %.2f\n",
            resid_cor_test$estimate, resid_cor_test$p.value))

cat("\n=== Section 4.2: Cross-outcome tests ===\n")
cat("Predictors with significantly different effects: p < 0.05:",
    sum(cross_outcome$p_value < 0.05), "| p < 0.001:", sum(cross_outcome$p_value < 0.001), "\n")
cat("Not significantly different:",
    paste(lab(cross_outcome$term[cross_outcome$p_value >= 0.05]), collapse = "; "), "\n")
opposite <- cross_outcome %>%
  left_join(test_comparison, by = "term") %>%
  filter(p_turnout < alpha, p_dshare < alpha, sign(b_z_turnout) != sign(b_z_dshare))
cat("Significant for both outcomes, opposite signs:", paste(lab(opposite$term), collapse = "; "), "\n")

cat("\nDone. Results in output/, report tables in output/exhibits/, figures in figures/.\n")
