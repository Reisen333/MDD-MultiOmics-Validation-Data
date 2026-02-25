#figure 2A bubble plot
library(dplyr)
library(ggplot2)
library(forcats)
library(readr)
library(scales)
allDEGresult <- read_csv("6个数据的DEG并集通路富集结果.csv") %>% mutate(method = "AllDEG")
RRAresult <- read_csv("RRA方法DEG通路富集结果.csv") %>% mutate(method = "RRA method")
geneintresult <- read_csv("基因交集法DEG通路富集结果.csv") %>% mutate(method = "Overlap method")
combined_results <- bind_rows(allDEGresult, RRAresult, geneintresult)
combined_results <- combined_results %>%
  mutate(database = if_else(source %in% c("GO:BP", "GO:MF", "GO:CC"), "GO", source))
significant_results <- combined_results %>%
  filter(adjusted_p_value < 0.05)
pathway_counts <- significant_results %>%
  group_by(term_name) %>%
  summarise(method_count = n_distinct(method), .groups = 'drop')
intersection_pathways <- pathway_counts %>%
  filter(method_count == 3) %>%
  pull(term_name)
top_intersection_data <- significant_results %>%
  filter(term_name %in% intersection_pathways)
pathway_avg_significance <- top_intersection_data %>%
  group_by(term_name, database) %>%
  summarise(avg_neg_log_p = mean(negative_log10_of_adjusted_p_value), .groups = 'drop')
top_per_database <- pathway_avg_significance %>%
  group_by(database) %>%
  top_n(20, wt = avg_neg_log_p) %>%
  ungroup()
final_pathway_list <- unique(top_per_database$term_name)
plot_data <- significant_results %>%
  filter(term_name %in% final_pathway_list) %>%
  mutate(method = factor(method, levels = c("Overlap method", "RRA method", "AllDEG")),
         database = factor(database, levels = c("GO", "KEGG", "REAC"))) %>%
  arrange(database, term_name) %>%
  mutate(term_name = fct_inorder(term_name))
method_gradients <- list(
  "Overlap method" = seq_gradient_pal("#f287a5", "#e50f4b"),  
  "RRA method"     = seq_gradient_pal("#a4bee1", "#497cc2"),  
  "AllDEG"         = seq_gradient_pal("#afa3cc", "#5f4798"))  
plot_data <- plot_data %>%
  mutate(color_value = negative_log10_of_adjusted_p_value)
plot_data <- plot_data %>%
  group_by(method) %>%
  mutate(normalized = (color_value - min(color_value)) / (max(color_value) - min(color_value)),
         color_hex = case_when(
           method == "Overlap method" ~ method_gradients[["Overlap method"]](normalized),
           method == "RRA method"     ~ method_gradients[["RRA method"]](normalized),
           method == "AllDEG"         ~ method_gradients[["AllDEG"]](normalized))) %>% ungroup()
bubble_plot <- ggplot(plot_data,aes(x = term_name,
                                    y = method,
                                    size = negative_log10_of_adjusted_p_value,
                                    color = color_hex)) + geom_point(alpha = 1) +
  facet_grid(~ database, scales = "free_x", space = "free_x", switch = "x") +
  scale_size_continuous(range = c(3, 10), name = "-log10(P.adj)") +
  scale_color_identity() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.placement = "outside",
        strip.background = element_rect(fill = "gray90", colour = "gray50"),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_line(colour = "grey80", linetype = "dashed"),
        panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "grey")) +
  labs(x = "Pathway", y = "", title = "Top 20 pathway results in different methods")
print(bubble_plot)


#preprocessing of sc-RNA data (here is an example of GSE201687, same procedure is also applied to GSE144136.)
library(Seurat)
library(dplyr)
library(tidyverse)
library(harmony)
folders <- list.dirs(full.names = FALSE, recursive = FALSE)
seurat_objects <- list()
for (folder in folders) {
  data_path <- file.path("C:/R/Rworking/抑郁症数据/单细胞 雌性猕猴的小胶质细胞特异性反应/单细胞", folder)
  data <- Read10X(data.dir = data_path)
  seurat_obj <- CreateSeuratObject(counts = data, project = folder)
  seurat_objects[[folder]] <- seurat_obj}
saveRDS(seurat_objects, file = "数据集体读取.rds")
sce.all <- merge(seurat_objects[[1]], seurat_objects[-1])
rm(depressed_list)
rm(seurat_objects)
gc()
sce.all <- JoinLayers(sce.all)
sce.all <- JoinLayers(sce.all, layer = "data")
min.cells <- 3
min.features <- 100
counts <- GetAssayData(sce.all,assay = "RNA", slot = "counts")
gene_counts <- colSums(counts > 0) 
filtered_genes <- names(gene_counts[gene_counts >= min.cells])  
sce.all <- sce.all[, filtered_genes] 
cell_counts <- rowSums(counts > 0) 
filtered_cells <- names(cell_counts[cell_counts >= min.features])  
sce.all <- sce.all[filtered_cells, ] 
rm(counts)
gc()
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^mt-")
pctMT <- 5
sce.all <- subset(sce.all, subset = percent.mt < pctMT)
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst", nfeatures = 2000)
sce.all <- NormalizeData(sce.all)
sce.all <- ScaleData(sce.all, features = rownames(sce.all))
sce.all <- RunPCA(sce.all, pc.genes = VariableFeatures(sce.all),seed.use=3)
sce.all <- RunHarmony(sce.all, group.by = "orig.ident")
sce.all <- RunUMAP(sce.all, dims = 1:10, reduction = "harmony")
sce.all <- FindNeighbors(sce.all, dims = 1:30, reduction = "harmony")
sce.all <- FindClusters(sce.all, resolution = 0.8)
saveRDS(sce.all, file = "抑郁组和健康组整合.rds")#The data of GSE144136 saved as"同猕猴预处理.rds"


#figure 2B UMAP plot
MDD <- readRDS("C:/R/Rworking/抑郁症数据/单细胞 重度抑郁症前额叶皮层的单核转录组学涉及少突胶质细胞前体细胞和兴奋性神经元/同猕猴预处理.rds")
yiyv <- readRDS("C:/R/Rworking/抑郁症数据/单细胞 雌性猕猴的小胶质细胞特异性反应/单细胞/抑郁组和健康组整合.rds")
#unifying the cell type names of the two datasets
mdd_meta <- MDD@meta.data
mdd_meta <- mdd_meta %>%
  mutate(
    cell_type_simplified = case_when(
      str_starts(cell_type, "Astros_") ~ "AST",
      str_starts(cell_type, "Ex_") ~ "EXN",
      str_starts(cell_type, "Inhib_") ~ "INH",
      str_starts(cell_type, "Oligos_") ~ "OLI",
      str_starts(cell_type, "OPCs_") ~ "OPC",
      cell_type == "Micro/Macro" ~ "MIC",
      TRUE ~ "OTHER"),
    condition = factor(ifelse(group == "Control", "Control", "Disease"), levels = c("Control", "Disease")))
MDD <- AddMetaData(MDD, metadata = mdd_meta)
DimPlot(yiyv, reduction = "umap", group.by = "cell_type",label = T)
DimPlot(MDD, reduction = "umap", group.by = "cell_type_simplified",label = T)


#figure 2C stacked bar plot
library(dplyr)
library(ggplot2)
library(patchwork) 
library(scales)
my_colors <- c(
  "AST" = "#9999FF", 
  "EXN" = "#39c5bb", 
  "INH" = "#2A9D8E", 
  "INT" = "#2A9D8E",
  "MIC" = "#497cc2", 
  "OLI" = "#e50f4b",
  "OPC" = "#ff8800")
#create general functions for proportion calculation and drawing
create_proportion_plot <- function(metadata, condition_col, cell_type_col, plot_title) {
  proportion_data <- metadata %>%
    select(
      condition = all_of(condition_col),
      cell_type = all_of(cell_type_col)) %>%
    #standardize cell type names
    mutate(cell_type = case_when(
      cell_type == "INT" ~ "INH",
      cell_type == "Inhib-VIP" ~ "INH", 
      TRUE ~ cell_type )) %>%
    group_by(condition, cell_type) %>%
    summarise(count = n(), .groups = 'drop_last') %>%
    mutate(proportion = count / sum(count)) %>%
    mutate(cell_type = factor(cell_type, levels = names(my_colors))) %>%
    arrange(desc(cell_type)) %>%
    group_by(condition) %>%
    mutate(label_y_pos = cumsum(proportion) - 0.5 * proportion) %>%
    ungroup()
  p <- ggplot( proportion_data, aes(x = condition, y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack", color = "white") +
    geom_text(
      data = . %>% filter(proportion > 0.015), 
      aes(y = label_y_pos, label = percent(proportion, accuracy = 0.1)),
      color = "black",
      size = 3.5) +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = my_colors) +
    labs( title = plot_title,
          x = "group",
          y = "cell types proportion",
          fill = "cell type" ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_blank(), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
  return(p)}
#GSE201687
p_yiyv <- create_proportion_plot(
  metadata = yiyv@meta.data,
  condition_col = "group",
  cell_type_col = "cell_type",
  plot_title = "Monkey data")
#GSE144136
p_mdd <- create_proportion_plot(
  metadata = MDD@meta.data,
  condition_col = "group",
  cell_type_col = "cell_type_simplified",
  plot_title = "Human data")
final_plot <- p_mdd + p_yiyv + 
  plot_layout(guides = 'collect') & 
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)))
print(final_plot)


#figure 2D violin plot
library(viridis)
library(tidyr)
library(tibble)
library(scales)
#need to run the "gene intersection analysis" part of the code in "figure 1 code" first
#gene enrichment in each cell type of GSE144136
#check how many of these genes are present in the expression matrix of the new dataset
up_genes_in_MDD <- intersect(intersect_up_degs, rownames(MDD))
down_genes_in_MDD <- intersect(intersect_down_degs, rownames(MDD))
MDD <- AddModuleScore(
  object = MDD,
  features = list(up_genes_in_MDD),
  name = "Bulk_UP_Score_on_MDD")
MDD <- AddModuleScore(
  object = MDD,
  features = list(down_genes_in_MDD),
  name = "Bulk_DOWN_Score_on_MDD")
up_medians <- tapply(MDD$Bulk_UP_Score_on_MDD1, 
                     MDD$cell_type_simplified,median)
up_sorted_levels <- names(sort(up_medians, decreasing = TRUE))
MDD$cell_type_sorted_up <- factor(MDD$cell_type_simplified, levels = up_sorted_levels)
down_medians <- tapply(MDD$Bulk_DOWN_Score_on_MDD1, MDD$cell_type_simplified, median)
down_sorted_levels <- names(sort(down_medians, decreasing = TRUE))
MDD$cell_type_sorted_down <- factor(MDD$cell_type_simplified, levels = down_sorted_levels)
p_vln_up_MDD_sorted <- VlnPlot(
  object = MDD,
  features = "Bulk_UP_Score_on_MDD1",
  group.by = "cell_type_sorted_up", 
  cols = my_colors,pt.size = 0) +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
  labs(title = "Common Up-regulated DEGs Score in Human data", 
       y = "Enrichment Score", 
       x = "Identity (Sorted by Score)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")
p_vln_down_MDD_sorted <- VlnPlot(
  object = MDD,
  features = "Bulk_DOWN_Score_on_MDD1",
  group.by = "cell_type_sorted_down",
  cols = my_colors,pt.size = 0) + 
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
  labs(title = "Common Down-regulated DEGs Score in Human data", 
       y = "Enrichment Score", 
       x = "Identity (Sorted by Score)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")
print(p_vln_up_MDD_sorted | p_vln_down_MDD_sorted)
#gene enrichment in each cell type of GSE201687
up_genes_in_yiyv <- intersect(intersect_up_degs, rownames(yiyv))
down_genes_in_yiyv <- intersect(intersect_down_degs, rownames(yiyv))
yiyv <- AddModuleScore(
  object = yiyv,
  features = list(up_genes_in_yiyv),
  name = "Bulk_UP_Score_on_yiyv")
yiyv <- AddModuleScore(
  object = yiyv,
  features = list(down_genes_in_yiyv),
  name = "Bulk_DOWN_Score_on_yiyv")
up_medians <- tapply(yiyv$Bulk_UP_Score_on_yiyv1, 
                     yiyv$cell_type,median)
up_sorted_levels <- names(sort(up_medians, decreasing = TRUE))
yiyv$cell_type_sorted_up <- factor(yiyv$cell_type, levels = up_sorted_levels)
down_medians <- tapply(yiyv$Bulk_DOWN_Score_on_yiyv1, yiyv$cell_type, median)
down_sorted_levels <- names(sort(down_medians, decreasing = TRUE))
yiyv$cell_type_sorted_down <- factor(yiyv$cell_type, levels = down_sorted_levels)
p_vln_up_yiyv_sorted <- VlnPlot(
  object = yiyv,
  features = "Bulk_UP_Score_on_yiyv1",
  group.by = "cell_type_sorted_up", 
  cols = my_colors,pt.size = 0) +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
  labs(title = "Common Up-regulated DEGs Score in Monkey data", 
       y = "Enrichment Score", 
       x = "Identity (Sorted by Score)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")
p_vln_down_yiyv_sorted <- VlnPlot(
  object = yiyv,
  features = "Bulk_DOWN_Score_on_yiyv1",
  group.by = "cell_type_sorted_down",
  cols = my_colors,pt.size = 0) + 
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
  labs(title = "Common Down-regulated DEGs Score in Monkey data", 
       y = "Enrichment Score", 
       x = "Identity (Sorted by Score)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")
print(p_vln_up_yiyv_sorted | p_vln_down_yiyv_sorted)


#bayesian deconvolution analysis
#need to run the "preprocessing of sc-RNA data" and "figure 2B UMAP plot" part of the code in "figure 2 code" first
#prepare sc-RNA input data
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(tidytext)
library(tidyr)
library(BayesPrism)
library(reshape2)
library(future)
sc_counts <- GetAssayData(MDD, assay = "RNA", layer = "counts")
cell_types <- MDD$cell_type_simplified
names(cell_types) <- colnames(sc_counts)
#prepare bulk RNA input data (taking GSE98793 as example)
bulk_matrix_log2 <- as.matrix(GSE98793)
#if the data is raw count data, then skip the following "bulk_matrix_linear <- 2^bulk_matrix_log2"
bulk_matrix_linear <- 2^bulk_matrix_log2
common_genes <- intersect(rownames(sc_counts), rownames(bulk_matrix_linear))
sc_counts_filtered <- sc_counts[common_genes, ]
bulk_matrix_filtered <- bulk_matrix_linear[common_genes, ]
sc_counts_filtered <- sc_counts_filtered[rowSums(sc_counts_filtered) > 0, ]
common_genes_final <- rownames(sc_counts_filtered)
bulk_matrix_filtered <- bulk_matrix_filtered[common_genes_final, ]
n_cores_to_use <- availableCores() - 1
plan("multisession", workers = n_cores_to_use)
#create bayesPrism object
sc_counts_transposed <- t(sc_counts_filtered)
bulk_matrix_transposed <- t(bulk_matrix_filtered)
final_gene_order <- intersect(colnames(sc_counts_transposed), colnames(bulk_matrix_transposed))
sc_final <- sc_counts_transposed[, final_gene_order]
bulk_final <- bulk_matrix_transposed[, final_gene_order]
prism_obj <- new.prism(
  reference = sc_final,
  mixture = bulk_final,
  cell.type.labels = cell_types,
  cell.state.labels = cell_types,
  key = NULL,
  outlier.cut = 0.01,
  outlier.fraction = 0.1)
#run bayesprism
start_time <- Sys.time()
bp_results <- run.prism(prism = prism_obj, n.cores = n_cores_to_use)
end_time <- Sys.time()
plan("sequential")
saveRDS(bp_results, file = "MDD对GSE98793贝叶斯结果.rds")#save each bulk data


#figure 2E box plot (taking bayesian result of GSE76826 as example)
library(BayesPrism)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
bp.result <- MDD对GSE76826贝叶斯结果
theta <- bp.result@posterior.theta_f@theta
theta_df <- as.data.frame(theta) %>%
  tibble::rownames_to_column("sample") %>%
  mutate(group = ifelse(grepl("control", sample), "control", "MDD"))
theta_long_df <- theta_df %>%
  pivot_longer(
    cols = -c(sample, group),
    names_to = "cell_type",
    values_to = "proportion")
head(theta_long_df)
cell_to_plot <- "OLI"    #here can changed cell type
plot_data <- theta_long_df %>%
  filter(cell_type == cell_to_plot)
ggplot(plot_data, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA,alpha = 1) +
  geom_jitter(width = 0.2, size = 2,alpha = 0.3) +
  scale_fill_manual(values = c("control" = "#497cc2", "MDD" = "#e50f4b")) +
  stat_compare_means(method = "wilcox.test", 
                     label.x.npc = "center",
                     label.y.npc = 0.9, 
                     aes(label = ..p.format..)) +
  labs(title = paste("Proportion of", cell_to_plot),
       x = "Group",y = "Cell Type Proportion") +
  theme_classic() +
  theme(legend.position = "none")


#figure 2F line plot
library(BayesPrism)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
bp.result <- MDD对益生菌bulk贝叶斯结果
theta <- bp.result@posterior.theta_f@theta
theta_df <- as.data.frame(theta) %>%
  tibble::rownames_to_column("sample") %>%
  separate(sample, into = c("Rep", "Region", "Group"), sep = "-") %>%
  pivot_longer(
    cols = -c(Rep, Region, Group),
    names_to = "cell_type",
    values_to = "proportion") %>%
  mutate( Rep = str_remove(Rep, "Rep"),
          Group = case_when(
            Group == "FT" ~ "Fluoxetine",
            Group == "Bupropion" ~ "Bupropion",
            Group == "Desipramine" ~ "Desipramine",
            Group == "Probiotics" ~ "Probiotics",
            Group == "Maltodextrin" ~ "Maltodextrin",
            Group == "Control" ~ "Control",
            TRUE ~ Group)) %>%
  mutate(Rep = factor(Rep, levels = as.character(1:5)),
         Region = factor(Region),
         Group = factor(Group, levels = c("Control", "Maltodextrin", "Probiotics", "Fluoxetine", "Bupropion", "Desipramine")))
cell_type_summary <- theta_df %>% 
  group_by(Group, cell_type) %>% 
  summarise(mean_proportion = mean(proportion),
            se = sd(proportion) / sqrt(n()),.groups = "drop")
cell_type_summary$Group <- factor(cell_type_summary$Group,
                                  levels = c("Control","Fluoxetine","Maltodextrin","Probiotics","Bupropion","Desipramine"))
ggplot(cell_type_summary, aes(x = Group, y = mean_proportion, group = cell_type, color = cell_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_proportion - se, ymax = mean_proportion + se),
                width = 0.1, linewidth = 0.6) +
  scale_color_manual(values = my_colors) +
  labs( title = "Cell Type Proportions Across Treatment Groups",
        x = "Treatment Group",
        y = "Mean Proportion",
        color = "Cell Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16))


#figure 2F stacked bar plot
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
#create a general function
create_bayes_plot <- function(bp_results, group_expression, dataset_name, color_palette) {
  cell_fractions <- get.fraction(bp = bp_results, which.theta = "final", state.or.type = "type")
  cell_fractions_df <- as.data.frame(cell_fractions)
  cell_fractions_df$Sample <- rownames(cell_fractions_df)
  cell_fractions_df <- cell_fractions_df %>%
    mutate(Group = !!group_expression)
  print(table(cell_fractions_df$Group))
  aggregated_fractions <- cell_fractions_df %>%
    dplyr::select(-Sample) %>%
    group_by(Group) %>%
    summarise(across(everything(), mean), .groups = 'drop')
  plot_data_aggregated <- aggregated_fractions %>%
    pivot_longer(cols = -Group,
                 names_to = "CellType",
                 values_to = "MeanFraction")
  p <- ggplot(plot_data_aggregated, aes(x = Group, y = MeanFraction, fill = CellType)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", MeanFraction * 100)),
              position = position_stack(vjust = 0.5),
              color = "black", size = 4) +
    theme_classic() +
    labs(title = paste("BayesPrism result:", dataset_name),
         x = "group",y = "cell type proportion",
         fill = "cell type") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_palette) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  return(p)}
#use the drawing function
fluoxetine_groups <- quote(case_when(
  grepl("Sham", Sample, ignore.case = TRUE) ~ "Sham",
  grepl("Fluoxetine", Sample, ignore.case = TRUE) ~ "Fluoxetine"))
plot_fluoxetine <- create_bayes_plot(MDD对氟西汀bulk数据贝叶斯结果, fluoxetine_groups, "GSE194289", my_colors)
print(plot_fluoxetine)


#cellchat analysis (here is an example of GSE201687, same procedure is also applied to GSE144136.)
library(CellChat)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(tidytext)
library(tidyr)
yiyv <- readRDS("C:/R/Rworking/抑郁症数据/单细胞 雌性猕猴的小胶质细胞特异性反应/单细胞/抑郁组和健康组整合.rds")
seurat_depression <- subset(yiyv, subset = group == "depression")
seurat_healthy <- subset(yiyv, subset = group == "healthy")
rm(yiyv)
gc()
#depression group analysis
data_input_dep <- GetAssayData(seurat_depression, assay = "RNA", layer = "data") # 或者 slot = "data"
meta_dep <- seurat_depression@meta.data
cellchat_dep <- createCellChat(object = data_input_dep, meta = meta_dep, group.by = "cell_type")
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat_dep@DB <- CellChatDB.use
cellchat_dep <- subsetData(cellchat_dep)
cellchat_dep <- identifyOverExpressedGenes(cellchat_dep)
cellchat_dep <- identifyOverExpressedInteractions(cellchat_dep)
#healthy group analysis
data_input_healthy <- GetAssayData(seurat_healthy, assay = "RNA", layer = "data")
meta_healthy <- seurat_healthy@meta.data
cellchat_healthy <- createCellChat(object = data_input_healthy, meta = meta_healthy, group.by = "cell_type")
cellchat_healthy@DB <- CellChatDB.use
cellchat_healthy <- subsetData(cellchat_healthy)
cellchat_healthy <- identifyOverExpressedGenes(cellchat_healthy)
cellchat_healthy <- identifyOverExpressedInteractions(cellchat_healthy)
cellchat_dep <- computeCommunProb(cellchat_dep)
cellchat_dep <- filterCommunication(cellchat_dep, min.cells = 10)
cellchat_dep <- computeCommunProbPathway(cellchat_dep)
cellchat_dep <- aggregateNet(cellchat_dep)
cellchat_healthy <- computeCommunProb(cellchat_healthy)
cellchat_healthy <- filterCommunication(cellchat_healthy, min.cells = 10)
cellchat_healthy <- computeCommunProbPathway(cellchat_healthy)
cellchat_healthy <- aggregateNet(cellchat_healthy)


#figure 2G circle plot
cell_types_healthy <- levels(cellchat_healthy@idents)
colors_for_healthy <- my_colors[cell_types_healthy]
cell_types_dep <- levels(cellchat_dep@idents)
colors_for_dep <- my_colors[cell_types_dep]
groupSize_healthy <- as.numeric(table(cellchat_healthy@idents))
groupSize_dep <- as.numeric(table(cellchat_dep@idents))
par(mfrow = c(1, 2), xpd=TRUE)
netVisual_circle(
  cellchat_healthy@net$count,
  vertex.weight = groupSize_healthy,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = colors_for_healthy,
  title.name = "Healthy: Number of Interactions")
par(mfrow = c(1, 2), xpd=TRUE)
netVisual_circle(
  cellchat_dep@net$count,
  vertex.weight = groupSize_dep,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = colors_for_dep,
  title.name = "Depression: Number of Interactions")


