set.seed(123)
melanocytic_genes <- c("TYR", "DCT", "MLANA", "PMEL", "MITF", "SOX10", "MIA", "GPNMB", "S100B", "ERBB3")
weight_matrix_melanocytic<-TR_weight("Mel.malignant.rds",melanocytic_genes,genie_file = "genie_matrix.rds")

saveRDS(weight_matrix_melanocytic,file="TR_Weight_Melanocytic_punished.rds")