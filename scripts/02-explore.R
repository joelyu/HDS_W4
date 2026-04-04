# =============================================================================
# 02-explore.R
# Exploratory analysis: cross-cohort immune infiltration & survival
#
# Focus:
#   1. Merge scores with clinical data, harmonise IDs and PAM50 labels
#   2. Univariate Cox screen (unadjusted) across 14 cell types — both cohorts
#   3. PAM50-adjusted Cox screen with FDR correction — identifies top predictors
#
# Expects: scores_metabric, scores_tcga from 01-immune-scoring.R
#          proc_dir defined by parent qmd
# =============================================================================

library(survival)

if (!exists("proc_dir")) proc_dir <- file.path("data", "processed")

# --- 1. Merge scores with clinical data --------------------------------------
clin_metabric <- read.csv(file.path(proc_dir, "metabric_clinical.csv"))
clin_tcga <- read.csv(file.path(proc_dir, "tcga_clinical.csv"))

# Standardise patient_id matching (TCGA sample IDs may have extra suffixes)
# METABRIC: patient_id matches directly
# TCGA: sampleId in expression data may differ from patientId in clinical data
# Try direct match first, then truncate to first 12 chars (TCGA barcode patient portion)
merge_scores_clinical <- function(scores, clinical, label) {
  # Direct merge
  merged <- inner_join(scores, clinical, by = c("patient_id", "cohort"))

  if (nrow(merged) == 0) {
    # TCGA barcodes: sample IDs are like TCGA-XX-XXXX-01, patient IDs TCGA-XX-XXXX
    scores$patient_id_short <- substr(scores$patient_id, 1, 12)
    clinical$patient_id_short <- substr(clinical$patient_id, 1, 12)
    merged <- inner_join(scores, clinical,
      by = c("patient_id_short", "cohort"),
      suffix = c("", ".clin")
    ) %>%
      select(-patient_id_short, -any_of("patient_id.clin"))
  }

  message(sprintf(
    "[%s] Merged: %d/%d patients matched to clinical",
    label, nrow(merged), nrow(scores)
  ))
  merged
}

surv_metabric <- merge_scores_clinical(scores_metabric, clin_metabric, "metabric")
surv_tcga <- merge_scores_clinical(scores_tcga, clin_tcga, "tcga")

# Parse OS event (cBioPortal format: "1:DECEASED" / "0:LIVING")
parse_os_event <- function(os_status) {
  as.integer(sub(":.*", "", os_status))
}

# Harmonise PAM50 labels (TCGA uses "BRCA_LumA", METABRIC uses "LumA")
clean_pam50 <- function(x) {
  x <- sub("^BRCA_", "", x)
  x <- sub("^Basal$", "Basal-like", x)
  x <- sub("^Her2$", "HER2-enriched", x)
  x <- sub("^LumA$", "Luminal A", x)
  x <- sub("^LumB$", "Luminal B", x)
  x <- sub("^Normal$", "Normal-like", x)
  x
}

surv_metabric <- surv_metabric %>%
  filter(!is.na(os_months) & !is.na(os_status) & os_status != "") %>%
  mutate(
    os_event = parse_os_event(os_status),
    pam50_subtype = clean_pam50(pam50_subtype)
  ) %>%
  filter(!is.na(pam50_subtype) & pam50_subtype != "" &
    pam50_subtype != "NC" & pam50_subtype != "claudin-low")

surv_tcga <- surv_tcga %>%
  filter(!is.na(os_months) & !is.na(os_status) & os_status != "") %>%
  mutate(
    os_event = parse_os_event(os_status),
    pam50_subtype = clean_pam50(pam50_subtype)
  ) %>%
  filter(!is.na(pam50_subtype) & pam50_subtype != "" &
    pam50_subtype != "NC" & pam50_subtype != "claudin-low")

message(sprintf(
  "\nSurvival cohorts: METABRIC n=%d (events=%d) | TCGA n=%d (events=%d)",
  nrow(surv_metabric), sum(surv_metabric$os_event),
  nrow(surv_tcga), sum(surv_tcga$os_event)
))

# --- 2. Univariate Cox: individual cell types (per cohort) ------------------
run_cox_screen <- function(surv_df, score_cols, label) {
  results <- lapply(score_cols, function(ct) {
    surv_df$score_z <- scale(surv_df[[ct]])[, 1]
    fit <- coxph(Surv(os_months, os_event) ~ score_z, data = surv_df)
    s <- summary(fit)
    data.frame(
      cohort = label,
      cell_type = ct,
      hr = s$conf.int[1, 1],
      hr_lower = s$conf.int[1, 3],
      hr_upper = s$conf.int[1, 4],
      p_value = s$coefficients[1, 5],
      concordance = concordance(fit)$concordance,
      n = s$n,
      events = s$nevent,
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    mutate(fdr = p.adjust(p_value, method = "BH"))
  # Note: score_z is created on a local copy of surv_df inside each lapply closure.
  # R's pass-by-value semantics mean the caller's data frame is never mutated.

  message(sprintf("\n[%s] Univariate Cox (unadjusted):", label))
  message(sprintf(
    "  %d/%d cell types FDR < 0.05",
    sum(results$fdr < 0.05), nrow(results)
  ))

  # Print significant results
  sig <- results %>%
    filter(fdr < 0.05) %>%
    arrange(p_value)
  if (nrow(sig) > 0) {
    for (i in seq_len(nrow(sig))) {
      message(sprintf(
        "  %s: HR=%.3f [%.3f-%.3f] p=%.2e",
        sig$cell_type[i], sig$hr[i], sig$hr_lower[i],
        sig$hr_upper[i], sig$p_value[i]
      ))
    }
  }

  results
}

# All 14 cell types (excluding T_cells aggregate which is sum of subsets)
score_cols <- setdiff(
  colnames(scores_metabric),
  c("patient_id", "cohort", "T_cells")
)

stopifnot(all(score_cols %in% colnames(surv_metabric)))
stopifnot(all(score_cols %in% colnames(surv_tcga)))

cox_metabric <- run_cox_screen(surv_metabric, score_cols, "metabric")
cox_tcga <- run_cox_screen(surv_tcga, score_cols, "tcga")

cox_both <- bind_rows(cox_metabric, cox_tcga)

# --- 3. PAM50-adjusted Cox (per cohort) --------------------------------------
run_cox_adjusted <- function(surv_df, score_cols, label) {
  results <- lapply(score_cols, function(ct) {
    surv_df$score_z <- scale(surv_df[[ct]])[, 1]
    fit <- coxph(Surv(os_months, os_event) ~ score_z + pam50_subtype, data = surv_df)
    s <- summary(fit)
    # Extract the score_z coefficient (first row)
    data.frame(
      cohort = label,
      cell_type = ct,
      hr = s$conf.int[1, 1],
      hr_lower = s$conf.int[1, 3],
      hr_upper = s$conf.int[1, 4],
      p_value = s$coefficients[1, 5],
      concordance = concordance(fit)$concordance,
      n = s$n,
      events = s$nevent,
      model = "PAM50-adjusted",
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    mutate(fdr = p.adjust(p_value, method = "BH"))

  message(sprintf("\n[%s] PAM50-adjusted Cox:", label))
  message(sprintf(
    "  %d/%d cell types FDR < 0.05",
    sum(results$fdr < 0.05), nrow(results)
  ))

  sig <- results %>%
    filter(fdr < 0.05) %>%
    arrange(p_value)
  if (nrow(sig) > 0) {
    for (i in seq_len(nrow(sig))) {
      message(sprintf(
        "  %s: HR=%.3f [%.3f-%.3f] p=%.2e",
        sig$cell_type[i], sig$hr[i], sig$hr_lower[i],
        sig$hr_upper[i], sig$p_value[i]
      ))
    }
  }

  results
}

adj_metabric <- run_cox_adjusted(surv_metabric, score_cols, "metabric")
adj_tcga <- run_cox_adjusted(surv_tcga, score_cols, "tcga")

adj_both <- bind_rows(adj_metabric, adj_tcga)

# --- 4. Summary table: side-by-side HRs -------------------------------------
comparison <- cox_both %>%
  select(cohort, cell_type, hr, hr_lower, hr_upper, p_value, fdr) %>%
  mutate(model = "unadjusted") %>%
  bind_rows(
    adj_both %>%
      select(cohort, cell_type, hr, hr_lower, hr_upper, p_value, fdr) %>%
      mutate(model = "PAM50-adjusted")
  )

# Save for later plotting
write.csv(comparison, file.path(proc_dir, "cox_comparison.csv"), row.names = FALSE)
write.csv(surv_metabric, file.path(proc_dir, "metabric_surv_ready.csv"), row.names = FALSE)
write.csv(surv_tcga, file.path(proc_dir, "tcga_surv_ready.csv"), row.names = FALSE)

message("\n=== Exploration complete ===")
message(sprintf("Output files in %s:", proc_dir))
for (f in list.files(proc_dir, pattern = "\\.csv$")) {
  message(sprintf("  %s", f))
}
