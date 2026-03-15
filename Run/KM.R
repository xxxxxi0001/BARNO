weight_matrix_melanocytic<-readRDS("TR_Weight_Melanocytic.rds")
signature_genes<-c("APEX1", "DLX5", "ENO1", "ILF2", "ILF3", "NME1", "NONO", "PARP1", "RAN", "ZNF207")
proliferate_melanoma_APEX1_fit<-km_curve(signature_genes,"TCGA_SKCM.rds","Top 10 TF of Proliferate Melanoma Negative")