# ============================================================
# Reviewer Comment 5 fix (root cause found):
# The survival analysis was broken because the LUAD code built OS from
# clinicalMatrix with  OS = ifelse(vital_status == "Dead", 1, 0),
# but TCGA uses "DECEASED" / "LIVING". Result: OS = 0 for everyone,
# no events -> degenerate KM curves -> "same p-value and at-risk counts".
#
# Fix: use the Xena survival files directly (they already have correct
# OS/OS.time) for BOTH cohorts.
# ============================================================
library(survival)
library(survminer)
library(pROC)

do_cohort <- function(expr_file, surv_file, id_col, tag) {
  # --- expression ---
  expr <- read.delim(expr_file, row.names = 1, check.names = FALSE)
  e <- as.data.frame(t(expr[c("PLA2G1B", "DSG2"), , drop = FALSE]))
  e$sampleID <- rownames(e)
  e$sampleID_short <- substr(e$sampleID, 1, 15)
  e$code <- substr(e$sampleID, 14, 15)
  e$group <- ifelse(e$code %in% c("01","02","03","04","05","06","07","08","09"), "Tumor",
              ifelse(e$code %in% c("10","11","12","13","14","15","16","17","18","19"), "Normal", "Other"))
  e <- e[e$group %in% c("Tumor", "Normal"), ]

  # --- survival (Xena file already has correct OS / OS.time) ---
  s <- read.delim(surv_file, check.names = FALSE)
  names(s)[names(s) == id_col] <- "sampleID"
  s$sampleID_short <- substr(s$sampleID, 1, 15)
  s$OS <- as.numeric(s$OS)
  s$OS.time <- as.numeric(s$OS.time)
  s <- s[!is.na(s$OS) & !is.na(s$OS.time) & s$OS.time > 0, ]

  m <- merge(e, s[, c("sampleID_short", "OS", "OS.time")], by = "sampleID_short")
  tumor <- m[m$group == "Tumor", ]
  cat("\n=====", tag, "tumor n =", nrow(tumor),
      " events =", sum(tumor$OS == 1), "=====\n")

  # --- Cox + risk split (identical to original intent) ---
  cox <- coxph(Surv(OS.time, OS) ~ PLA2G1B + DSG2, data = tumor)
  tumor$risk <- predict(cox, type = "risk")
  tumor$grp <- ifelse(tumor$risk > median(tumor$risk), "High", "Low")
  sd <- survdiff(Surv(OS.time, OS) ~ grp, data = tumor)
  pv <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  cat("log-rank p =", format(pv, digits = 6),
      "  High n =", sum(tumor$grp == "High"), " Low n =", sum(tumor$grp == "Low"), "\n")
  print(round(summary(cox)$conf.int, 3))

  fit <- survfit(Surv(OS.time, OS) ~ grp, data = tumor)
  g <- ggsurvplot(fit, data = tumor, pval = TRUE, risk.table = TRUE,
                  title = paste("TCGA", tag, "Overall Survival"))
  pdf(paste0(tag, "_KM_corrected.pdf"), width = 8, height = 7)
  print(g)
  dev.off()
  write.csv(tumor, paste0(tag, "_survival_corrected.csv"), row.names = FALSE)
}

# ---- LUAD ----
do_cohort("TCGA.LUAD.sampleMap_HiSeqV2.gz", "LUAD生存分析.txt", "sample", "LUAD")
# ---- LUSC ----
do_cohort("TCGA.LUSC.sampleMap_HiSeqV2.gz", "LUSC生存分析.txt", "xena_sample", "LUSC")
