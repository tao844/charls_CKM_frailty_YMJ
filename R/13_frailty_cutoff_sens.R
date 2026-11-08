## ============================================================
## 13_frailty_cutoff_sens.R —— 衰弱截断值敏感性（0.20 / 0.25 vs 主口径 0.21）
## 复核记录（2026-11-08）：六格模型同主规格（年龄时间轴、10 协变量、
## householdID 聚类），样本 = 12,182 模型样本。
##   0.20: 1.351 / 2.040 / 1.022 / 1.376 / 1.884（与主口径几乎逐位相同）
##   0.21: 1.349 / 2.063 / 1.023 / 1.384 / 1.882（主表）
##   0.25: 1.399 / 2.310 / 1.026 / 1.420 / 2.078（衰弱格 HR 略升，模式不变）
## 模型样本内衰弱占比：0.20 → 25.9%；0.21 → 25.2%；0.25 → 15.9%
## 产物：仅控制台汇总，无个体数据
## ============================================================
suppressMessages({library(dplyr); library(survival)})
s2 <- readRDS("output/rebuild/bl_cov_step2.rds")
s3 <- readRDS("output/rebuild/exposures_step3.rds")
dd <- s2 %>% inner_join(s3, by = "ID")
covs <- c("gender","edu_bl","rural_bl","marry_bl","smoken_bl","drinkev_bl",
          "sleep_bl","ln_hhcperc","ins_bl","pension_bl")
mok <- complete.cases(dd[, covs]) & !is.na(dd$FI41) & !is.na(dd$CKM_Stage)

run_cut <- function(cut) {
  fr <- cut(dd$FI41, breaks = c(-Inf, 0.10 - 1e-9, cut + 1e-9, Inf),
            labels = c("robust", "prefrail", "frail"))
  dd$j6 <- factor(paste(dd$CKM2, as.character(fr), sep = "_"),
                  levels = c("0-2_robust","0-2_prefrail","0-2_frail",
                             "3-4_robust","3-4_prefrail","3-4_frail"))
  dd$j6 <- relevel(dd$j6, ref = "0-2_robust")
  m <- coxph(Surv(age, age + dtime_yr, died) ~ j6 + gender + edu_bl + rural_bl +
             marry_bl + smoken_bl + drinkev_bl + sleep_bl + ln_hhcperc + ins_bl +
             pension_bl, data = dd[mok, ], cluster = householdID)
  cf <- summary(m)$conf.int
  cat("cutoff", cut, " 衰弱占比", round(mean(fr[mok] == "frail"), 3), "\n")
  print(round(cf[grep("j6", rownames(cf)), c("exp(coef)", "lower .95", "upper .95")], 3))
}

run_cut(0.20); run_cut(0.21); run_cut(0.25)
cat("13 完成\n")
