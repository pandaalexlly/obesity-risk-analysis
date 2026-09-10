# ==============================================================================
# 全国肥胖风险因素分析：非参数检验与稳健推断
# ==============================================================================
# 数据：obesity_level.csv（20,758 条记录，18 个原始字段）
# 目标：分析人口属性与生活方式因素和 BMI 的统计关联，并在 OLS 假设
#       不充分时使用 HC3 与 case bootstrap 进行稳健性检验。
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 环境准备
# ------------------------------------------------------------------------------
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "rstatix",
  "nortest", "lmtest", "sandwich", "car"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "缺少 R 包：", paste(missing_packages, collapse = ", "),
      "\n请先运行：\ninstall.packages(c(",
      paste0('"', missing_packages, '"', collapse = ", "),
      "))"
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(rstatix)
  library(nortest)
  library(lmtest)
  library(sandwich)
  library(car)
})

SEED <- 42
BOOT_B <- 1000
set.seed(SEED)

# ------------------------------------------------------------------------------
# 1. 路径与数据读取
# ------------------------------------------------------------------------------
data_candidates <- c(
  file.path(getwd(), "data", "obesity_level.csv"),
  file.path(getwd(), "obesity_level.csv")
)

existing_data <- data_candidates[file.exists(data_candidates)]

if (length(existing_data) == 0) {
  stop("未找到 obesity_level.csv。请放在项目根目录或 data/ 文件夹。")
}

data_path <- existing_data[1]
output_dir <- file.path(getwd(), "outputs")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

obesity <- read.csv(
  data_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expected_cols <- c(
  "Gender", "Age", "Height", "Weight",
  "family_history_with_overweight", "FAVC",
  "FCVC", "NCP", "CAEC", "SMOKE",
  "CH2O", "SCC", "FAF", "TUE",
  "CALC", "MTRANS"
)

missing_cols <- setdiff(expected_cols, names(obesity))
if (length(missing_cols) > 0) {
  stop("数据缺少字段：", paste(missing_cols, collapse = ", "))
}

obesity <- obesity %>%
  filter(complete.cases(across(all_of(expected_cols)))) %>%
  mutate(
    BMI = Weight / Height^2,
    Age_group = cut(
      Age,
      breaks = c(-Inf, 20, 30, 45, Inf),
      labels = c("20岁以下", "20-30岁", "30-45岁", "45岁以上"),
      right = FALSE
    ),
    Gender = factor(Gender, levels = c("Female", "Male")),
    family_history_with_overweight = factor(
      family_history_with_overweight,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),
    FAVC = factor(FAVC, levels = c(0, 1), labels = c("No", "Yes")),
    SMOKE = factor(SMOKE, levels = c(0, 1), labels = c("No", "Yes")),
    SCC = factor(SCC, levels = c(0, 1), labels = c("No", "Yes")),
    CAEC = factor(
      CAEC,
      levels = c("0", "Sometimes", "Frequently", "Always")
    ),
    CALC = factor(
      CALC,
      levels = c("0", "Sometimes", "Frequently")
    ),
    MTRANS = factor(
      MTRANS,
      levels = c(
        "Automobile", "Public_Transportation",
        "Walking", "Bike", "Motorbike"
      )
    )
  )

cat("数据维度：", nrow(obesity), "×", ncol(obesity), "\n")
cat(sprintf("BMI 均值：%.2f；中位数：%.2f\n",
            mean(obesity$BMI), median(obesity$BMI)))

# ------------------------------------------------------------------------------
# 2. 描述性统计与分布
# ------------------------------------------------------------------------------
group_summary <- obesity %>%
  group_by(Gender, Age_group) %>%
  summarise(
    N = n(),
    Mean_BMI = mean(BMI),
    Median_BMI = median(BMI),
    SD_BMI = sd(BMI),
    Q1 = quantile(BMI, 0.25),
    Q3 = quantile(BMI, 0.75),
    .groups = "drop"
  )

write.csv(
  group_summary,
  file.path(output_dir, "01_BMI_group_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

p_box <- ggplot(obesity, aes(x = Age_group, y = BMI, fill = Gender)) +
  geom_boxplot(
    alpha = 0.75,
    outlier.alpha = 0.15,
    position = position_dodge(width = 0.8)
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "不同性别与年龄段的 BMI 分布",
    x = "年龄段",
    y = "BMI",
    fill = "性别"
  )

ggsave(
  file.path(output_dir, "01_BMI_性别年龄段分布.png"),
  p_box, width = 9, height = 6, dpi = 180
)

# ------------------------------------------------------------------------------
# 3. 非参数组间比较
# ------------------------------------------------------------------------------
gender_test <- wilcox.test(
  BMI ~ Gender,
  data = obesity,
  conf.int = TRUE,
  exact = FALSE
)

age_test <- kruskal.test(
  BMI ~ Age_group,
  data = obesity
)

age_dunn <- obesity %>%
  dunn_test(
    BMI ~ Age_group,
    p.adjust.method = "bonferroni"
  )

nonparam_summary <- data.frame(
  Test = c("Wilcoxon: Gender", "Kruskal-Wallis: Age group"),
  Statistic = c(
    unname(gender_test$statistic),
    unname(age_test$statistic)
  ),
  P_value = c(gender_test$p.value, age_test$p.value)
)

write.csv(
  nonparam_summary,
  file.path(output_dir, "02_nonparametric_tests.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  age_dunn,
  file.path(output_dir, "03_age_Dunn_Bonferroni.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------------------------
# 4. 控制年龄和性别后的秩变换分析
# ------------------------------------------------------------------------------
# 对 BMI 进行秩变换后，以嵌套模型比较目标因素是否提供额外解释信息。
# 该分析用于稳健的组间关联比较，不作因果解释。
rank_adjusted_test <- function(data, factor_var) {
  tmp <- data
  tmp$BMI_rank <- rank(tmp$BMI, ties.method = "average")

  reduced_formula <- BMI_rank ~ Age + Gender
  full_formula <- as.formula(
    paste("BMI_rank ~ Age + Gender +", factor_var)
  )

  reduced_model <- lm(reduced_formula, data = tmp)
  full_model <- lm(full_formula, data = tmp)

  cmp <- anova(reduced_model, full_model)

  data.frame(
    Variable = factor_var,
    Df = cmp$Df[2],
    F_value = cmp$F[2],
    P_value = cmp$`Pr(>F)`[2]
  )
}

rank_variables <- c(
  "FAVC",
  "family_history_with_overweight",
  "CAEC",
  "FAF",
  "MTRANS"
)

rank_results <- bind_rows(
  lapply(
    rank_variables,
    function(v) rank_adjusted_test(obesity, v)
  )
)

write.csv(
  rank_results,
  file.path(output_dir, "04_rank_adjusted_tests.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------------------------
# 5. 多变量 OLS 模型
# ------------------------------------------------------------------------------
# 数值型变量标准化，使不同量纲的连续变量系数更便于比较。
model_data <- obesity %>%
  select(
    BMI, Gender, Age,
    family_history_with_overweight, FAVC,
    FCVC, NCP, CAEC, SMOKE,
    CH2O, SCC, FAF, TUE,
    CALC, MTRANS
  ) %>%
  mutate(
    across(
      c(Age, FCVC, NCP, CH2O, FAF, TUE),
      ~ as.numeric(scale(.x))
    )
  )

model_formula <- BMI ~
  Gender + Age +
  family_history_with_overweight + FAVC +
  FCVC + NCP + CAEC + SMOKE +
  CH2O + SCC + FAF + TUE +
  CALC + MTRANS

ols_model <- lm(model_formula, data = model_data)

model_fit <- data.frame(
  N = nobs(ols_model),
  R_squared = summary(ols_model)$r.squared,
  Adjusted_R_squared = summary(ols_model)$adj.r.squared,
  Residual_SE = summary(ols_model)$sigma,
  F_statistic = unname(summary(ols_model)$fstatistic[1]),
  Model_p_value = pf(
    summary(ols_model)$fstatistic[1],
    summary(ols_model)$fstatistic[2],
    summary(ols_model)$fstatistic[3],
    lower.tail = FALSE
  )
)

write.csv(
  model_fit,
  file.path(output_dir, "05_OLS_model_fit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------------------------
# 6. OLS 假设诊断
# ------------------------------------------------------------------------------
resid_ols <- residuals(ols_model)

ad_result <- nortest::ad.test(resid_ols)
bp_result <- lmtest::bptest(ols_model)

diagnostic_summary <- data.frame(
  Diagnostic = c(
    "Anderson-Darling normality",
    "Breusch-Pagan heteroskedasticity"
  ),
  Statistic = c(
    unname(ad_result$statistic),
    unname(bp_result$statistic)
  ),
  P_value = c(ad_result$p.value, bp_result$p.value)
)

write.csv(
  diagnostic_summary,
  file.path(output_dir, "06_model_diagnostics.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

vif_result <- car::vif(ols_model)

if (is.matrix(vif_result)) {
  vif_output <- data.frame(
    Variable = rownames(vif_result),
    vif_result,
    row.names = NULL,
    check.names = FALSE
  )
} else {
  vif_output <- data.frame(
    Variable = names(vif_result),
    VIF = as.numeric(vif_result)
  )
}

write.csv(
  vif_output,
  file.path(output_dir, "07_VIF.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

p_qq <- ggplot(
  data.frame(residual = resid_ols),
  aes(sample = residual)
) +
  stat_qq(alpha = 0.25) +
  stat_qq_line() +
  theme_minimal(base_size = 12) +
  labs(
    title = "OLS 残差 Q-Q 图",
    x = "理论分位数",
    y = "样本分位数"
  )

ggsave(
  file.path(output_dir, "02_OLS_residual_QQ.png"),
  p_qq, width = 7, height = 6, dpi = 180
)

p_resid <- ggplot(
  data.frame(
    fitted = fitted(ols_model),
    residual = resid_ols
  ),
  aes(x = fitted, y = residual)
) +
  geom_point(alpha = 0.15, size = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal(base_size = 12) +
  labs(
    title = "OLS 残差与拟合值",
    x = "拟合值",
    y = "残差"
  )

ggsave(
  file.path(output_dir, "03_OLS_residual_vs_fitted.png"),
  p_resid, width = 8, height = 6, dpi = 180
)

# ------------------------------------------------------------------------------
# 7. HC3 稳健标准误
# ------------------------------------------------------------------------------
coef_ols <- summary(ols_model)$coefficients
hc3_vcov <- sandwich::vcovHC(ols_model, type = "HC3")
hc3_test <- lmtest::coeftest(ols_model, vcov. = hc3_vcov)

hc3_se <- sqrt(diag(hc3_vcov))
hc3_ci <- cbind(
  Estimate = coef(ols_model),
  HC3_Lower = coef(ols_model) - 1.96 * hc3_se,
  HC3_Upper = coef(ols_model) + 1.96 * hc3_se
)

hc3_output <- data.frame(
  Term = names(coef(ols_model)),
  Estimate = coef(ols_model),
  OLS_SE = coef_ols[, "Std. Error"],
  OLS_P = coef_ols[, "Pr(>|t|)"],
  HC3_SE = hc3_test[, "Std. Error"],
  HC3_P = hc3_test[, "Pr(>|t|)"],
  HC3_Lower = hc3_ci[, "HC3_Lower"],
  HC3_Upper = hc3_ci[, "HC3_Upper"],
  row.names = NULL
)

write.csv(
  hc3_output,
  file.path(output_dir, "08_OLS_vs_HC3.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------------------------
# 8. Case Bootstrap
# ------------------------------------------------------------------------------
# 直接对观测行有放回重抽样，不要求所有残差共享同一方差结构。
# 该方法仍依赖观测之间近似独立，因此结果用于稳健性比较而非因果推断。
set.seed(SEED)

coef_names <- names(coef(ols_model))
boot_coef <- matrix(
  NA_real_,
  nrow = BOOT_B,
  ncol = length(coef_names),
  dimnames = list(NULL, coef_names)
)

n <- nrow(model_data)

for (b in seq_len(BOOT_B)) {
  idx <- sample.int(n, size = n, replace = TRUE)
  boot_data <- model_data[idx, , drop = FALSE]

  boot_fit <- try(
    lm(model_formula, data = boot_data),
    silent = TRUE
  )

  if (!inherits(boot_fit, "try-error")) {
    current_coef <- coef(boot_fit)
    boot_coef[b, names(current_coef)] <- current_coef
  }

  if (b %% 100 == 0) {
    cat("Bootstrap:", b, "/", BOOT_B, "\n")
  }
}

bootstrap_ci <- t(
  apply(
    boot_coef,
    2,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
)

ols_ci <- confint(ols_model, level = 0.95)

interval_compare <- data.frame(
  Term = coef_names,
  Estimate = coef(ols_model),
  OLS_Lower = ols_ci[, 1],
  OLS_Upper = ols_ci[, 2],
  HC3_Lower = hc3_ci[, "HC3_Lower"],
  HC3_Upper = hc3_ci[, "HC3_Upper"],
  Bootstrap_Lower = bootstrap_ci[, 1],
  Bootstrap_Upper = bootstrap_ci[, 2],
  Bootstrap_Valid_Runs = colSums(!is.na(boot_coef)),
  row.names = NULL
)

interval_compare <- interval_compare %>%
  mutate(
    OLS_Significant = OLS_Lower * OLS_Upper > 0,
    HC3_Significant = HC3_Lower * HC3_Upper > 0,
    Bootstrap_Significant = Bootstrap_Lower * Bootstrap_Upper > 0
  )

write.csv(
  interval_compare,
  file.path(output_dir, "09_interval_comparison.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------------------------
# 9. 置信区间比较图
# ------------------------------------------------------------------------------
plot_intervals <- interval_compare %>%
  filter(Term != "(Intercept)") %>%
  select(
    Term, Estimate,
    OLS_Lower, OLS_Upper,
    HC3_Lower, HC3_Upper,
    Bootstrap_Lower, Bootstrap_Upper
  ) %>%
  pivot_longer(
    cols = -c(Term, Estimate),
    names_to = c("Method", ".value"),
    names_pattern = "(OLS|HC3|Bootstrap)_(Lower|Upper)"
  )

p_ci <- ggplot(
  plot_intervals,
  aes(
    x = Estimate,
    y = reorder(Term, Estimate),
    xmin = Lower,
    xmax = Upper,
    linetype = Method
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(
    position = position_dodge(width = 0.6),
    height = 0
  ) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 1.5
  ) +
  theme_minimal(base_size = 10) +
  labs(
    title = "OLS、HC3 与 Case Bootstrap 95% 置信区间",
    x = "回归系数",
    y = NULL,
    linetype = "方法"
  )

ggsave(
  file.path(output_dir, "04_置信区间稳健性比较.png"),
  p_ci, width = 11, height = 9, dpi = 180
)

# ------------------------------------------------------------------------------
# 10. 稳健性汇总
# ------------------------------------------------------------------------------
robustness_summary <- interval_compare %>%
  filter(Term != "(Intercept)") %>%
  summarise(
    Total_terms = n(),
    OLS_significant = sum(OLS_Significant, na.rm = TRUE),
    HC3_significant = sum(HC3_Significant, na.rm = TRUE),
    Bootstrap_significant = sum(Bootstrap_Significant, na.rm = TRUE),
    All_three_agree = sum(
      OLS_Significant == HC3_Significant &
        HC3_Significant == Bootstrap_Significant,
      na.rm = TRUE
    )
  )

write.csv(
  robustness_summary,
  file.path(output_dir, "10_robustness_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\n============================================================\n")
cat("分析完成\n")
cat("============================================================\n")
cat(sprintf(
  "N = %d；Adjusted R² = %.4f\n",
  nobs(ols_model),
  summary(ols_model)$adj.r.squared
))
cat(sprintf(
  "Anderson-Darling p = %.4g；Breusch-Pagan p = %.4g\n",
  ad_result$p.value,
  bp_result$p.value
))
cat("结果目录：", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n")

sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

