#' Construct GRN for selected tf
#'
#' @param TF Selected TF for GRN construction
#' @param weight_matrix Weight Matrix's name
#' @param layer1_top How many TF you want to see in first layer
#' @param layer2_top How many TF you want to see in second layer
#' 
#' @return Layered graph
#' 
#' @export
layered_graph<-function(TF, 
                        weight_matrix,
                        layer1_top=10,
                        layer2_top=5){

  distance<-weight_matrix$distance_matrix
  distance_df<-as.data.frame(as.table(distance))
  colnames(distance_df) <- c("TF_name", "Target", "TR_Distance")

  weight<-weight_matrix$TR_weight_matrix
  weight_df<-as.data.frame(as.table(weight))
  colnames(weight_df) <- c("TF_name", "Target", "TR_Weight")

  tr_df<-merge(weight_df, distance_df, by = c("TF_name", "Target"))
  
  layer1_all<-tr_df[tr_df$TF_name==TF,]
  layer1_sort<-layer1_all[order(layer1_all$TR_Distance),]
  layer1<-head(layer1_sort,layer1_top)

  all_tf<-unique(as.character(tr_df$TF_name))
  intermediate_tf<-intersect(as.character(layer1$Target), all_tf)

  if (length(intermediate_tf)>0) {
    layer2_list<-list()
    for (i in seq_along(intermediate_tf)) {
      intermediate<-intermediate_tf[[i]]
      df_intermediate_tf<-tr_df[tr_df$TF_name==intermediate,]
      df_intermediate_tf_sorted<- df_intermediate_tf[order(df_intermediate_tf$TR_Distance),]
      layer2_list[[i]]<-head(df_intermediate_tf_sorted,layer2_top)
    }
    layer2<-do.call(rbind,layer2_list)
  }
  if (exists("layer2")){
    df_all_layers<-rbind(layer1,layer2)
    df_all_layers<-unique(df_all_layers)
  }
  else(df_all_layers<-layer1)
  
  graph <- igraph::graph_from_data_frame(df_all_layers, directed = TRUE)
  return(graph)
}