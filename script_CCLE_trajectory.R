library(tidyverse)

#Trajectory analysis
#download matrixes 
viper_es <- readRDS("./viper_es.rds")
ccle_mat_hub_ahr <- readRDS("./ccle_mat_hub_ahr.rds")
ccle_mat_hub_ahr <- ccle_mat_hub_ahr[!duplicated(rownames(ccle_mat_hub_ahr)), ]
ccle_mat_hub_ahr_score <- readRDS("./ccle_mat_hub_ahr_score.rds")
ccle_mat_hub_ahr_score <- ccle_mat_hub_ahr_score[!duplicated(rownames(ccle_mat_hub_ahr_score)), ]
ccle_mat_hub_ahr_common <- readRDS("./ccle_mat_hub_ahr_common.rds")
ccle_mat_hub_ahr_common <- ccle_mat_hub_ahr_common[!duplicated(rownames(ccle_mat_hub_ahr_common)), ]
ann <- readRDS("./ann.rds")
ann$

#wrap expression data
viper_es_w <- wrap_expression(
  expression = t(viper_es),
  counts = t(viper_es)
)
#add prior information
viper_es_w <- add_prior_information(
  viper_es_w,
  groups_id = data.frame(cell_id= colnames(viper_es), group_id = ann$ahr_go)
)
#wrap expression data
ccle_mat_hub_ahr_w <- wrap_expression(
  expression = t(ccle_mat_hub_ahr),
  counts = t(ccle_mat_hub_ahr)
)
#add prior information
ccle_mat_hub_ahr_w <- add_prior_information(
  ccle_mat_hub_ahr_w,
  groups_id = data.frame(cell_id= colnames(ccle_mat_hub_ahr), group_id = ann$ahr_go)
)
#wrap expression data
ccle_mat_hub_ahr_score_w <- wrap_expression(
  expression = t(ccle_mat_hub_ahr_score),
  counts = t(ccle_mat_hub_ahr_score)
)
#add prior information
ccle_mat_hub_ahr_score_w <- add_prior_information(
  ccle_mat_hub_ahr_score_w,
  groups_id = data.frame(cell_id= colnames(ccle_mat_hub_ahr_score), group_id = ann$ahr_go)
)
#wrap expression data
ccle_mat_hub_ahr_common_w <- wrap_expression(
  expression = t(ccle_mat_hub_ahr_common),
  counts = t(ccle_mat_hub_ahr_common)
)
#add prior information
ccle_mat_hub_ahr_common_w <- add_prior_information(
  ccle_mat_hub_ahr_common_w,
  groups_id = data.frame(cell_id= colnames(ccle_mat_hub_ahr_common), group_id = ann$ahr_go)
)

#run a method slingshot
model_viper_es <- infer_trajectory(viper_es_w, "slingshot")
saveRDS(model, "./model_viper_es_slingshot.rds")
model_viper_es_umap <- model_viper_es %>% add_dimred(dyndimred::dimred_umap, expression_source = viper_es_w$expression)

png("./trajectory_slingshot_ahr_go_gr.png", width = 900, height = 700)
plot_dimred(
  model_viper_es_umap, 
  expression_source = viper_es_w$expression, 
  grouping = ann$ahr_go
)
dev.off()

png("./trajectory_slingshot_ahr_go_tissue.png", width = 900, height = 700)
plot_dimred(
  model_viper_es_umap, 
  expression_source = viper_es_w$expression, 
  grouping = ann$tissue
)
dev.off()


#run a method slingshot
model_ccle_mat_hub_ahr <- infer_trajectory(ccle_mat_hub_ahr_w, "slingshot")
model_ccle_mat_hub_ahr_umap <- model_ccle_mat_hub_ahr %>% add_dimred(dyndimred::dimred_umap, expression_source = ccle_mat_hub_ahr_w$expression)
saveRDS(model, "./model_ccle_mat_hub_ahr_slingshot.rds")

png("./trajectory_slingshot_ahr_hub_ahr_gr.png", width = 900, height = 700)
plot_dimred(
  model_ccle_mat_hub_ahr_umap, 
  expression_source = ccle_mat_hub_ahr_w$expression, 
  grouping = ann$ahr_hub_ahr
)
dev.off()

png("./trajectory_slingshot_ahr_hub_ahr_tissue.png", width = 900, height = 700)
plot_dimred(
  model_ccle_mat_hub_ahr_umap, 
  expression_source = ccle_mat_hub_ahr_w$expression, 
  grouping = ann$tissue
)
dev.off()

#run a method slingshot
model_ccle_mat_hub_ahr_score <- infer_trajectory(ccle_mat_hub_ahr_score_w, "slingshot")
model_ccle_mat_hub_ahr_score_umap <- model_ccle_mat_hub_ahr_score %>% add_dimred(dyndimred::dimred_umap, expression_source = ccle_mat_hub_ahr_score_w$expression)
saveRDS(model, "./model_ccle_mat_hub_ahr_score_slingshot.rds")

png("./trajectory_slingshot_ahr_hub_ahr_score_gr.png", width = 900, height = 700)
plot_dimred(
  model_ccle_mat_hub_ahr_score_umap, 
  expression_source = ccle_mat_hub_ahr_score_w$expression, 
  grouping = ann$ahr_hub_ahr_score
)
dev.off()

png("./trajectory_slingshot_ahr_hub_ahr_score_tissue.png", width = 900, height = 700)
plot_dimred(
  model_ccle_mat_hub_ahr_score_umap, 
  expression_source = ccle_mat_hub_ahr_score_w$expression, 
  grouping = ann$tissue
)
dev.off()

#run a method slingshot
model_ccle_mat_hub_ahr_common <- infer_trajectory(ccle_mat_hub_ahr_common_w, "slingshot")
model_ccle_mat_hub_ahr_common_umap <- model_ccle_mat_hub_ahr_common %>% add_dimred(dyndimred::dimred_umap, expression_source = ccle_mat_hub_ahr_common_w$expression)
saveRDS(model, "./model_ccle_mat_hub_ahr_common_slingshot.rds")

png("./trajectory_slingshot_ahr_hub_ahr_common_gr.png", width = 900, height = 700)
plot_dimred(
  model_ccle_mat_hub_ahr_common_umap, 
  expression_source = ccle_mat_hub_ahr_common_w$expression, 
  grouping = ann$ahr_hub_ahr_common
)
dev.off()

png("./trajectory_slingshot_ahr_hub_ahr_common_tissue.png", width = 900, height = 700)
plot_dimred(
  model_ccle_mat_hub_ahr_common_umap, 
  expression_source = ccle_mat_hub_ahr_common_w$expression, 
  grouping = ann$tissue
)
dev.off()