#validation cohort data preprocessing
data <- read.delim("C:/R/Rworking/抑郁症数据/TR202509281050J5EB/input/4_expression/4.1_expression/Expression_with_annotation_fpkm.xls", row.names = 1, sep = '\t', check.names = FALSE)
library(dplyr)
library(tibble)
processed_data_dplyr <- data %>%
  rownames_to_column("EnsemblID") %>%
  select(EnsemblID, Name, contains(":read count")) %>%
  mutate(Gene = ifelse(is.na(Name) | Name == "" | Name == "-", EnsemblID, Name)) %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), mean), .groups = 'drop') %>%
  mutate(across(where(is.numeric), round)) %>%
  column_to_rownames("Gene") %>%
  rename_with(~gsub(":read count", "", .x), contains(":read count"))
head(processed_data_dplyr)
#DEG analysis
library(DESeq2)
count_data <- as.matrix(round(processed_data_dplyr))
head(count_data)
sample_names <- colnames(count_data)
sample_groups <- gsub("_.*", "", sample_names)
colData <- data.frame(
  row.names = sample_names,
  condition = factor(sample_groups))
print(colData)
all(colnames(count_data) == rownames(colData))
dds <- DESeqDataSetFromMatrix(
  countData = count_data,colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "MDD", "HC"))
head(res)
res_df <- as.data.frame(res)
pvalue_threshold <- 0.05  
log2fc_threshold <- log2(1.1)
significant_genes_pvalue <- subset(
  res_df, !is.na(pvalue) & pvalue < pvalue_threshold & abs(log2FoldChange) > log2fc_threshold)
upregulated_genes <- subset(significant_genes_pvalue, log2FoldChange > 0)
downregulated_genes <- subset(significant_genes_pvalue, log2FoldChange < 0)


#figure 5B violin plot
#need to run the "preprocessing of sc-RNA data" part of the code in "figure 2 code"
my_colors <- c(
  "AST" = "#9999FF", 
  "EXN" = "#39c5bb", 
  "INH" = "#2A9D8E", 
  "INT" = "#2A9D8E",
  "MIC" = "#497cc2", 
  "OLI" = "#e50f4b",
  "OPC" = "#ff8800")
bulk_up_genes <- rownames(upregulated_genes)
bulk_down_genes <- rownames(downregulated_genes)
#create enrichment scoring function
run_seurat_enrichment <- function(seurat_obj, up_genes, down_genes, cell_type_col, dataset_name, color_palette) {
  genes_in_obj <- rownames(seurat_obj)
  up_genes_in_obj <- intersect(up_genes, genes_in_obj)
  down_genes_in_obj <- intersect(down_genes, genes_in_obj)
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = list(up_genes_in_obj),
    name = "UP_Score")
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = list(down_genes_in_obj),
    name = "DOWN_Score")
  up_medians <- tapply(seurat_obj[["UP_Score1"]][,1], seurat_obj[[cell_type_col]], median)
  up_sorted_levels <- names(sort(up_medians, decreasing = TRUE))
  seurat_obj$temp_sorted_up <- factor(seurat_obj[[cell_type_col]][,1], levels = up_sorted_levels)
  down_medians <- tapply(seurat_obj[["DOWN_Score1"]][,1], seurat_obj[[cell_type_col]], median)
  down_sorted_levels <- names(sort(down_medians, decreasing = TRUE))
  seurat_obj$temp_sorted_down <- factor(seurat_obj[[cell_type_col]][,1], levels = down_sorted_levels)
  p_up <- VlnPlot(
    object = seurat_obj,
    features = "UP_Score1",
    group.by = "temp_sorted_up",
    cols = color_palette, 
    pt.size = 0) +
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
    labs(title = paste("Up-regulated Genes Score in", dataset_name),
         y = "Enrichment Score", x = "Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  p_down <- VlnPlot(
    object = seurat_obj,
    features = "DOWN_Score1",
    group.by = "temp_sorted_down",
    cols = my_colors,
    pt.size = 0) + 
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
    labs(title = paste("Down-regulated Genes Score in", dataset_name),
         y = "Enrichment Score", 
         x = "Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  return(list(plot_up = p_up, plot_down = p_down))}
#execute function
mdd_plots <- run_seurat_enrichment( 
  seurat_obj = MDD,
  up_genes = bulk_up_genes,
  down_genes = bulk_down_genes,
  cell_type_col = "cell_type_simplified",
  dataset_name = "Human Data",
  color_palette = my_colors) 
yiyv_plots <- run_seurat_enrichment(
  seurat_obj = yiyv,
  up_genes = bulk_up_genes,
  down_genes = bulk_down_genes,
  cell_type_col = "cell_type",
  dataset_name = "Monkey Data",
  color_palette = my_colors)
print(mdd_plots$plot_up | mdd_plots$plot_down)
print(yiyv_plots$plot_up | yiyv_plots$plot_down)


#figure 5C bubble plot
library(dplyr)
library(ggplot2)
library(forcats)
library(readr)
library(patchwork)
pathway <- read_csv("验证实验通路富集结果.csv") %>% mutate(method = "验证实验")
pathway <- pathway %>%
  mutate(database = if_else(source %in% c("GO:BP", "GO:MF", "GO:CC"), "GO", source))
significant_results <- pathway %>%
  filter(adjusted_p_value < 0.05)
top_n_value <- 20#top20
top_pathways_data <- significant_results %>%
  group_by(database) %>%
  top_n(top_n_value, wt = negative_log10_of_adjusted_p_value) %>%
  ungroup()
db_levels <- c("GO", "KEGG", "REAC")
plot_list <- list()
for (db in db_levels) {
  local_data <- top_pathways_data %>%
    filter(database == db)
  if(nrow(local_data) == 0) next
  local_data <- local_data %>%
    arrange(intersection_size) %>%
    mutate(term_name = fct_inorder(term_name))
  p <- ggplot(local_data, aes(x = intersection_size, y = term_name, 
                              color = negative_log10_of_adjusted_p_value,
                              size = intersection_size)) +
    geom_point(alpha = 0.8) +
    scale_color_gradient(
      low = "#497cc2", high = "#e50f4b", 
      name = "-log10(P.adj)") +
    scale_size_continuous(range = c(3, 10), name = "Gene Count") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_blank(), 
          panel.grid.major.y = element_line(colour = "grey90", linetype = "dashed"),
          panel.grid.major.x = element_line(colour = "grey95"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 1)) +
    labs(x = NULL, title = db)
  plot_list[[db]] <- p}
final_plot <- wrap_plots(plot_list, ncol = 1) + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = paste0("Top ", top_n_value, " Enriched Pathways"),
                  theme = theme(plot.title = element_text(size = 14, hjust = 0.5))) &
  theme(plot.margin = margin(5, 5, 5, 5))
last_db <- names(plot_list)[length(plot_list)]
plot_list[[last_db]] <- plot_list[[last_db]] + labs(x = "Gene Count (Intersection Size)")
final_plot <- wrap_plots(plot_list, ncol = 1) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = paste0("Top ", top_n_value, " Enriched Pathways"))
print(final_plot)


#figure 5D bar plot
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(ggpubr)
mdd_samples_to_keep <- c("MDD_1", "MDD_2", "MDD_4", "MDD_6")
hc_samples_to_keep <- c("HC_1", "HC_2", "HC_3")
all_samples_to_keep <- c(mdd_samples_to_keep, hc_samples_to_keep)
normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts_subset <- normalized_counts[, all_samples_to_keep]
target_genes <- c("ANO10", "HK2", "CCDC50")
plot_data_subset <- as.data.frame(normalized_counts_subset) %>%
  rownames_to_column("Gene") %>%
  filter(Gene %in% target_genes) %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  mutate(condition = factor(ifelse(grepl("MDD", Sample), "MDD", "HC"), levels = c("HC", "MDD")))
head(plot_data_subset)
ggplot(plot_data_subset, aes(x = condition, y = Expression, fill = condition)) +
  geom_rect(data = . %>% distinct(condition),
            aes(xmin = as.numeric(condition) - 0.5,
                xmax = as.numeric(condition) + 0.5,
                ymin = -Inf,ymax = Inf,fill = condition),
            alpha = 0.2, 
            inherit.aes = FALSE) + 
  geom_violin(trim = FALSE, alpha = 0.7) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(title = "Gene Expression in MDD vs HC Subgroups",
       x = "Group",
       y = "Normalized Expression Count") +
  scale_fill_manual(values = c("MDD" = "#e50f4b", "HC" = "#497cc2")) + 
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_blank(), 
        strip.text = element_text(face = "bold", size = 12))


#figure 5E bar plot
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
samples_to_include <- c(
  "MDD_1", "MDD_2", "MDD_4", "MDD_6",  
  "HC_1", "HC_2", "HC_3" )              
calibration_method <- "mean" 
control_group_name <- "HC"
calibrator_sample_id <- "HC_1"
data_dct <- read_csv("qPCR_delta_ct.csv")
data_filtered <- data_dct %>%
  filter(Group1 %in% samples_to_include)
if (calibration_method == "mean") {
  calibrator_dcts <- data_filtered %>%
    filter(Group == control_group_name) %>%
    group_by(Gene) %>%
    summarise(calibrator_val = mean(delta_Ct))
  plot_title <- paste(control_group_name)} else if (calibration_method == "single") {
    calibrator_dcts <- data_filtered %>%
      filter(Group1 == calibrator_sample_id) %>%
      select(Gene, calibrator_val = delta_Ct)
    plot_title <- paste(calibrator_sample_id)} else {
      stop("WARNING")}
print(calibrator_dcts)
#calculate ΔΔCt and relative expression levels
data_recalculated <- data_filtered %>%
  left_join(calibrator_dcts, by = "Gene") %>%
  mutate(delta_delta_ct = delta_Ct - calibrator_val,
         relative_expression = 2^(-delta_delta_ct))
summary_data <- data_recalculated %>%
  group_by(Group, Gene) %>%
  summarise(N = n(),
            mean_expression = mean(relative_expression),
            sem = sd(relative_expression) / sqrt(N),.groups = 'drop')
summary_data$Group <- factor(summary_data$Group, levels = c("HC", "MDD"))
data_recalculated$Group <- factor(data_recalculated$Group, levels = c("HC", "MDD"))
group_colors <- c("HC" = "#497cc2", "MDD" = "#e50f4b")
final_plot <- ggplot(data = summary_data, 
                     aes(x = Gene, y = mean_expression, fill = Group)) +
  geom_bar(stat = "identity",
           position = position_dodge(0.8),
           width = 0.7, color = "black", 
           linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_expression - sem, ymax = mean_expression + sem),
                position = position_dodge(0.8),
                width = 0.25,linewidth = 0.8) +
  stat_compare_means(data = data_recalculated,
                     aes(y = relative_expression, group = Group),
                     label = "p.signif",method = "t.test",
                     hide.ns = FALSE) +
  scale_fill_manual(values = group_colors) +
  labs(title = plot_title,x = "Gene",
       y = expression("Relative Expression (2"^-Delta*Delta*C[t]*")"),fill = "Group") +
  theme_classic(base_size = 14) +
  theme( plot.title = element_text(hjust = 0.5, face = "bold"),
         axis.text = element_text(color = "black", size = 12),
         axis.title = element_text(face = "bold"),
         legend.position = "top")
print(final_plot)
