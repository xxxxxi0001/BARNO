library(survival)
n<-5
tcga_data<-readRDS("TCGA_SKCM.rds")
top_tf_Treg<-readRDS("Top_TF_of_Treg_punished.rds")
top_tf_Th1<-readRDS("Top_TF_of_Th1_punished.rds")
top_tf_Tfh<-readRDS("Top_TF_of_Tfh_punished.rds")
top_tf_TEMemory<-readRDS("Top_TF_of_TEMemory_punished.rds")
Treg_pos<-as.character(head(top_tf_Treg$top_positive_summary,n)$TF)
Treg_neg<-as.character(head(top_tf_Treg$top_negative_summary,n)$TF)
Th1_pos<-as.character(head(top_tf_Th1$top_positive_summary,n)$TF)
Th1_neg<-as.character(head(top_tf_Th1$top_negative_summary,n)$TF)
Tfh_pos<-as.character(head(top_tf_Tfh$top_positive_summary,n)$TF)
Tfh_neg<-as.character(head(top_tf_Tfh$top_negative_summary,n)$TF)
TEMemory_pos<-as.character(head(top_tf_TEMemory$top_positive_summary,n)$TF)
TEMemory_neg<-as.character(head(top_tf_TEMemory$top_negative_summary,n)$TF)
top_tf_cytotoxic<-readRDS("Top_TF_of_Cytotoxic_T_punished.rds")
top_tf_exhausted<-readRDS("Top_TF_of_Exhausted_T_punished.rds")
cytotoxic_pos<-as.character(head(top_tf_cytotoxic$top_positive_summary,n)$TF)
cytotoxic_neg<-as.character(head(top_tf_cytotoxic$top_negative_summary,n)$TF)
exhausted_pos<-as.character(head(top_tf_exhausted$top_positive_summary,n)$TF)
exhausted_neg<-as.character(head(top_tf_exhausted$top_negative_summary,n)$TF)
top_tf_Melanocytic<-readRDS("Top_TF_of_Melanocytic_punished.rds")
top_tf_stress_melanoma<-readRDS("Top_TF_of_stress_melanoma_punished.rds")
top_tf_proliferate_melanoma<-readRDS("Top_TF_of_proliferate_melanoma_punished.rds")
top_tf_invasive_melanoma<-readRDS("Top_TF_of_invasive_melanoma_punished.rds")
melanocytic_pos<-as.character(head(top_tf_Melanocytic$top_positive_summary,n)$TF)
melanocytic_neg<-as.character(head(top_tf_Melanocytic$top_negative_summary,n)$TF)
stress_melanoma_pos<-as.character(head(top_tf_stress_melanoma$top_positive_summary,n)$TF)
stress_melanoma_neg<-as.character(head(top_tf_stress_melanoma$top_negative_summary,n)$TF)
proliferate_melanoma_pos<-as.character(head(top_tf_proliferate_melanoma$top_positive_summary,n)$TF)
proliferate_melanoma_neg<-as.character(head(top_tf_proliferate_melanoma$top_negative_summary,n)$TF)
invasive_melanoma_pos<-as.character(head(top_tf_invasive_melanoma$top_positive_summary,n)$TF)
invasive_melanoma_neg<-as.character(head(top_tf_invasive_melanoma$top_negative_summary,n)$TF)
module_all<-list(
  melanocytic=list(pos=melanocytic_pos,neg=melanocytic_neg),
  stress_melanoma=list(pos=stress_melanoma_pos,neg=stress_melanoma_neg),
  proliferate_melanoma=list(pos=proliferate_melanoma_pos,neg=proliferate_melanoma_neg),
  invasive_melanoma=list(pos=invasive_melanoma_pos,neg=invasive_melanoma_neg),
  cytotoxic=list(pos=cytotoxic_pos,neg=cytotoxic_neg),
  exhausted=list(pos=exhausted_pos,neg=exhausted_neg),
  Treg=list(pos=Treg_pos,neg=Treg_neg),
  Th1=list(pos=Th1_pos,neg=Th1_neg),
  Tfh=list(pos=Tfh_pos,neg=Tfh_neg),
  TEMemory=list(pos=TEMemory_pos,neg=TEMemory_neg)
)
df<-readRDS("TCGA_SKCM.rds")
df_tpm<-df$tpm
calculate_module_score<-function(genes, df_tpm) {
  genes_use <- intersect(unique(genes), rownames(df_tpm))
  if (length(genes_use) == 0) {
    return(rep(NA, ncol(df_tpm)))
  }
  colMeans(df_tpm[genes_use, , drop = FALSE], na.rm = TRUE)
}

# initialization
module_scores <- list()

# for each module, calculate module score by mean
for (module_name in names(module_all)) {
  for (direction in names(module_all[[module_name]])) {
    genes <- module_all[[module_name]][[direction]]
    score_name <- paste(module_name, direction, sep = "_")
    module_scores[[score_name]] <- calculate_module_score(genes, df_tpm)
  }
}

# into data frame
module_score_df <- as.data.frame(module_scores)

# and save survival data
module_score_df$survival <- tcga_data$survival

# for each module run cox
module_names <- names(module_scores)
module_cox_results <- lapply(module_names, function(score_name) {
  fit <- coxph(
    as.formula(paste0("survival ~ ", score_name)),
    data = module_score_df
  )
  
  cox_summary <- summary(fit)
  # and store result
  data.frame(
    module = score_name,
    coefficient = cox_summary$coefficients[1, "coef"],
    HR = cox_summary$coefficients[1, "exp(coef)"],
    z_score = cox_summary$coefficients[1, "z"],
    p.value = cox_summary$coefficients[1, "Pr(>|z|)"],
    lower95 = cox_summary$conf.int[1, "lower .95"],
    upper95 = cox_summary$conf.int[1, "upper .95"]
  )
})

# bind result 
module_cox_results <- do.call(rbind, module_cox_results)
# adjusted p
module_cox_results$FDR <- p.adjust(module_cox_results$p.value, method = "fdr")
# get additional information
module_cox_results$direction <- ifelse(module_cox_results$HR > 1, "pro-death", "pro-survival")
# order based on p
module_cox_results <- module_cox_results[order(module_cox_results$p.value), ]
# save
write.csv(module_cox_results,"single_module_cox_5.csv")
