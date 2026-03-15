weight_matrix_melanocytic<-readRDS("TR_Weight_mel_proliferate_punished.rds")
mst_proliferate_melanoma_RAN<-layered_graph("RAN",weight_matrix_melanocytic)
directed_mst_visualization(mst_proliferate_melanoma_RAN,"RAN",aspect_ratio = 0.5,text_size=0.8,text_centrality=0,arrow_size=1)
