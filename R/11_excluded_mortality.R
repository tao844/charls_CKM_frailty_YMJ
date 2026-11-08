## ============================================================
## 11_excluded_mortality.R —— 排除组成亡核算（Methods/Limitations 量化依据）
## 回答：3,472 名因基线认知缺失（代理访谈/题目缺失）被排除者的死亡率
## 口径：复现 00b 排除链前三步（年龄、卒中、认知缺失）→ ID →
##       匹配死亡库（CHARLS_exit_csv ... analysis_db_mortality.xlsx，
##       Harmonized CHARLS Version D 推导）
## 复核记录（2026-11-08）：可链接 2,165/3,472，死亡 468 = 21.6%；
##   其中代理访谈 171/855 = 20.0%，非代理缺失 285/1,261 = 22.6%；
##   对照 13,256 分析队列 1,546 = 11.7%（bl_cov_step2$died）
## 方法学落稿句："among the 3,472 so excluded, the 2,165 linkable to the
##   CHARLS mortality follow-up had a mortality of 21.6%, versus 11.7% in the
##   analytic cohort"
## ============================================================
suppressMessages({library(haven); library(dplyr); library(readxl)})

proj  <- "D:/00 陶万清的工作/07 working stata/01 CKM_frailty_revision"
pkg11 <- file.path(proj, "CHARLS_official/household_and_community_questionnaire_data")
id12 <- function(x) paste0(substr(x, 1, 9), "0", substr(x, 10, 11))

## ---- 1. 官方2011底表 + 认知计分（与 00b_rebuild_cohort.R 同口径）----
demo <- read_dta(file.path(pkg11, "demographic_background.dta")); demo$ID12 <- id12(demo$ID)
hs   <- read_dta(file.path(pkg11, "health_status_and_functioning.dta")); hs$ID12 <- id12(hs$ID)
d <- demo %>% select(ID = ID12, birth_yr = ba002_1) %>%
  left_join(hs %>% select(ID = ID12, proxy, da007_8_,
    dc001s1, dc001s2, dc001s3, dc002, dc003, dc019, dc020, dc021, dc022, dc023, dc025,
    all_of(paste0("dc006s", 1:11)), all_of(paste0("dc027s", 1:11))), by = "ID")
d$age <- 2011 - as.numeric(d$birth_yr)

imm_w <- rowSums(!is.na(d[paste0("dc006s", 1:10)])); del_w <- rowSums(!is.na(d[paste0("dc027s", 1:10)]))
s11i <- !is.na(d$dc006s11); s11d <- !is.na(d$dc027s11)
imm <- ifelse(imm_w > 0, imm_w, ifelse(s11i, 0, NA_real_))
del <- ifelse(del_w > 0, del_w, ifelse(s11d, 0, NA_real_))
memeory <- (imm + del) / 2
dow <- ifelse(is.na(d$dc002), NA, as.numeric(d$dc002 == 1))
season <- ifelse(is.na(d$dc003), NA, as.numeric(d$dc003 == 1))
s7a <- as.numeric(d$dc019); s7b <- as.numeric(d$dc020); s7c <- as.numeric(d$dc021)
s7d <- as.numeric(d$dc022); s7e <- as.numeric(d$dc023)
s7_item <- cbind(s7a == 93, s7b == 86, s7c == 79, s7d == 72, s7e == 65)
tics <- cbind(!is.na(d$dc001s1), !is.na(d$dc001s2), !is.na(d$dc001s3),
              dow, season, s7_item, d$dc025 == 1)
admin <- cbind(!is.na(d$dc001s1), !is.na(d$dc001s2), !is.na(d$dc001s3),
               !is.na(d$dc002), !is.na(d$dc003),
               !is.na(s7a), !is.na(s7b), !is.na(s7c), !is.na(s7d), !is.na(s7e),
               !is.na(d$dc025))
admin[is.na(admin)] <- FALSE; tics[is.na(tics)] <- FALSE
executive <- ifelse(rowSums(admin) > 0, rowSums(ifelse(admin, tics, 0)), NA_real_)
total_cog <- ifelse(is.na(memeory) | is.na(executive), NA, memeory + executive)

## ---- 2. 排除链（与 00b 一致：st1 年龄, st2 卒中, st3 认知缺失）----
st1 <- is.na(d$age) | d$age < 45
stroke <- ifelse(is.na(d$da007_8_), NA, as.numeric(d$da007_8_ == 1))
st2 <- !st1 & !is.na(stroke) & stroke == 1
st3 <- !st1 & !st2 & is.na(total_cog)
stopifnot(nrow(d) == 17705, sum(st3) == 3472, sum(d$proxy[st3] == 1, na.rm = TRUE) == 1235)

## ---- 3. 死亡库匹配 ----
excl <- d$ID[st3]
mort <- read_excel(file.path(proj, "CHARLS_exit_csv",
  "CHARLS — 死亡（Mortality） (dbcodebook.cn)", "analysis_db_mortality.xlsx"))
mm <- mort[match(excl, mort$id), ]
linked <- !is.na(mm$Die_status)
cat("认知排除 3,472 人：可链接", sum(linked), "；死亡", sum(mm$Die_status == 1, na.rm = TRUE),
    sprintf("(%.1f%%)", mean(mm$Die_status[linked] == 1) * 100), "\n")

proxy_sub <- !is.na(d$proxy[st3]) & d$proxy[st3] == 1
for (i in 1:2) {
  sel <- linked & (if (i == 1) proxy_sub else !proxy_sub & !is.na(d$proxy[st3]))
  sub <- mm[sel, ]
  cat(sprintf("  %s: n=%d 死亡=%d (%.1f%%) 死亡时间中位=%.1f年\n",
      ifelse(i == 1, "代理访谈", "非代理题目缺失"), sum(sel),
      sum(sub$Die_status == 1), mean(sub$Die_status == 1) * 100,
      median(sub$Die_time_months[sub$Die_status == 1], na.rm = TRUE) / 12))
}

## ---- 4. 对照：13,256 分析队列同口径死亡率 ----
s2 <- readRDS("output/rebuild/bl_cov_step2.rds")
cat("对照 13,256 分析队列: 死亡", sum(s2$died),
    sprintf("(%.1f%%)", mean(s2$died) * 100), "\n")
