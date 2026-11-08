## ============================================================
## 12_fi_nockm_and_2wave_sens.R —— FI_noCKM 六格敏感性 + 双波连续确认敏感性
## （2026-11-08 方法学漏洞修复轮；复核头注见下）
##
## 复核记录：
##   FI_noCKM 由 00e_exposures_rebuild.R 直接产出（36 项、k36/n36、n36>=29），
##   本脚本仅拟合模型。主模型同规格：年龄时间轴 + 延迟进入、10 协变量、
##   householdID 聚类稳健方差。样本 = 协变量完整 & FI_noCKM/CKM 可判。
##   结果：n = 11,752，死亡 1,415；
##     0-2 prefrail 1.176 (0.932-1.484) | 0-2 frail 1.834 (1.439-2.337)
##     3-4 robust 1.086 (0.781-1.510)   | 3-4 prefrail 1.280 (1.007-1.627)
##     3-4 frail 1.710 (1.337-2.186)
##   双波确认：用 _cog_long_for_review.rds 的波次特异 cog（cog_h 为基线值），
##   事件 = 首次连续两波 cog<=thr；时间 = 确认波 t - 基线年；死亡竞争删失。
##   结果：n = 8,428（至少一次随访认知），确认事件 341，竞争死亡 894；
##     五格点估计 1.145-1.473，95%CI 全部跨 1（null 模式不变）。
## 产物：output/tables/_fi_nockm_sens_v2.rds（仅汇总，无个体数据）
## ============================================================
suppressMessages({library(dplyr); library(survival)})

s2 <- readRDS("output/rebuild/bl_cov_step2.rds")
s3 <- readRDS("output/rebuild/exposures_step3.rds")
dd <- s2 %>% inner_join(s3, by = "ID")

covs <- c("gender","edu_bl","rural_bl","marry_bl","smoken_bl","drinkev_bl",
          "sleep_bl","ln_hhcperc","ins_bl","pension_bl")  # age 为时间轴
cc <- complete.cases(dd[, covs]) & !is.na(dd$FI_noCKM) & !is.na(dd$CKM_Stage)
dd$j6nc <- factor(paste(dd$CKM2, dd$frail3_nc, sep = "_"),
                  levels = c("0-2_robust","0-2_prefrail","0-2_frail",
                             "3-4_robust","3-4_prefrail","3-4_frail"))
dd$j6nc <- relevel(dd$j6nc, ref = "0-2_robust")

m <- coxph(Surv(age, age + dtime_yr, died) ~ j6nc + gender + edu_bl + rural_bl +
           marry_bl + smoken_bl + drinkev_bl + sleep_bl + ln_hhcperc + ins_bl +
           pension_bl, data = dd[cc, ], cluster = householdID)
est <- summary(m)$conf.int
est <- est[grep("j6nc", rownames(est)), c("exp(coef)", "lower .95", "upper .95")]
cat("[FI_noCKM] n =", sum(cc), " deaths =", sum(dd$died[cc]), "\n")
print(round(est, 3))
saveRDS(list(n = sum(cc), deaths = sum(dd$died[cc]), est = est),
        "output/tables/_fi_nockm_sens_v2.rds")

## ---- 双波连续确认 ----
L <- readRDS("output/rebuild/_cog_long_for_review.rds")
b <- L$baseline
f <- L$followup
f <- f[!is.na(f$cog), ]                      # 波次特异分数；cog_h 为基线值勿用
f <- f[order(f$ID, f$t), ]
f$imp <- f$cog <= f$thr
f$prev_imp <- c(NA, head(f$imp, -1))
f$prev_imp[c(NA, head(f$ID, -1)) != f$ID] <- NA
conf <- f[f$imp & !is.na(f$prev_imp) & f$prev_imp, ]
conf <- conf[!duplicated(conf$ID), c("ID", "t", "bl_time.x")]
conf$ctime_c <- conf$t - conf$bl_time.x
lastfu <- f %>% group_by(ID) %>%
  summarise(lastfu = max(t - bl_time.x, na.rm = TRUE), .groups = "drop")
rs <- merge(b, conf[, c("ID", "ctime_c")], by = "ID", all.x = TRUE)
rs <- merge(rs, lastfu, by = "ID", all.x = TRUE)
rs$ci_evc <- as.integer(!is.na(rs$ctime_c))
rs$ci_tmc <- ifelse(rs$ci_evc == 1, rs$ctime_c,
             ifelse(rs$died == 1, pmin(rs$lastfu, rs$dtime_yr, na.rm = TRUE),
                    rs$lastfu))
mok <- !is.na(rs$ci_tmc) & rs$ci_tmc > 0
rs$j6 <- factor(rs$joint6, levels = c("0-2 | robust","0-2 | prefrail","0-2 | frail",
                                      "3-4 | robust","3-4 | prefrail","3-4 | frail"))
rs$j6 <- relevel(rs$j6, ref = "0-2 | robust")
m2 <- coxph(Surv(ci_tmc, ci_evc) ~ j6 + gender + edu_bl + rural_bl + marry_bl +
            smoken_bl + drinkev_bl + sleep_bl + ln_hhcperc + ins_bl + pension_bl,
            data = rs[mok, ], cluster = householdID)
est2 <- summary(m2)$conf.int
est2 <- est2[grep("j6", rownames(est2)), c("exp(coef)", "lower .95", "upper .95")]
cat("[two-wave confirmed] n =", sum(mok), " events =", sum(rs$ci_evc[mok]),
    " competing deaths =", sum(rs$ci_evc[mok] == 0 & rs$died[mok] == 1), "\n")
print(round(est2, 3))
cat("12 完成\n")
