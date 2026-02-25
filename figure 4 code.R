#figure 4A volcano plot
#need to run the "RRA analysis" part of the code in "figure 1 code"
library(ggplot2)
library(ggrepel)
library(dplyr)
full_rra_results_raw <- rbind(rra_up_results, rra_down_results)
full_rra_results_dedup <- full_rra_results_raw %>%
  group_by(Name) %>%
  slice_min(order_by = Score, n = 1, with_ties = FALSE) %>%
  ungroup()
avg_logfc_all <- sapply(full_rra_results_dedup$Name, calculate_average_logfc, results_list = all_limma_results)
full_rra_results_dedup$Avg_logFC <- avg_logfc_all
p_threshold_rra <- 0.05
logfc_threshold_rra <- log2(1.1)
plot_data_rra <- full_rra_results_dedup %>%
  filter(!is.na(Avg_logFC) & !is.na(Score)) %>%
  mutate(neg_log10_Score = -log10(Score),
         Status = case_when(
           Avg_logFC > logfc_threshold_rra & Score < p_threshold_rra  ~ "Upregulated",
           Avg_logFC < -logfc_threshold_rra & Score < p_threshold_rra ~ "Downregulated",
           TRUE                                                       ~ "Not Significant"))
genes_to_highlight <- c("ANO10", "HK2", "CCDC50")
highlight_subset_rra <- subset(plot_data_rra, Name %in% genes_to_highlight)
ano10_rank <- which(core_up_genes_df$Name == "ANO10")
hk2_rank <- which(core_up_genes_df$Name == "HK2")
ccdc50_rank <- which(core_down_genes_df$Name == "CCDC50")
total_up <- nrow(core_up_genes_df)
total_down <- nrow(core_down_genes_df)
ano10_stats <- subset(core_up_genes_df, Name == "ANO10")
hk2_stats <- subset(core_up_genes_df, Name == "HK2")
ccdc50_stats <- subset(core_down_genes_df, Name == "CCDC50")
summary_text <- paste(
  "Key up gene (total: ", total_up, "):\n",
  "  - ANO10: Rank ", ano10_rank, "\n (AvgLogFC=", round(ano10_stats$Avg_logFC, 2), ", Score=", formatC(ano10_stats$Score, format = "e", digits = 2), ")\n",
  "  - HK2:   Rank ", hk2_rank, "\n (AvgLogFC=", round(hk2_stats$Avg_logFC, 2), ", Score=", formatC(hk2_stats$Score, format = "e", digits = 2), ")\n\n",
  "Key down gene (total: ", total_down, "):\n",
  "  - CCDC50: Rank ", ccdc50_rank, "\n (AvgLogFC=", round(ccdc50_stats$Avg_logFC, 2), ", Score=", formatC(ccdc50_stats$Score, format = "e", digits = 2), ")",
  sep = "")
meta_volcano_plot <- ggplot(
  plot_data_rra, 
  aes(x = Avg_logFC, y = neg_log10_Score)) +
  geom_point(aes(color = Status), alpha = 0.5, size = 2) +
  geom_point(data = highlight_subset_rra, aes(fill = Status), size = 5, shape = 21, color = "black", stroke = 1) +
  geom_text_repel(data = highlight_subset_rra, aes(label = Name), box.padding = 1.2, point.padding = 0.6, 
                  segment.color = 'black', fontface = "bold", size = 5) +
  annotate(geom = "label", x = -Inf, y = Inf,label = summary_text, hjust = 0, vjust = 1,
           size = 3.5, fill = "white", alpha = 0.8, label.padding = unit(0.5, "lines")) +
  scale_color_manual(values = c("Upregulated" = "#e50f4b", "Downregulated" = "#497cc2", "Not Significant" = "grey80"), 
                     name = "Gene Status") +
  scale_fill_manual(values = c("Upregulated" = "#e50f4b", "Downregulated" = "#497cc2"), guide = "none") +
  geom_hline(yintercept = -log10(p_threshold_rra), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-logfc_threshold_rra, logfc_threshold_rra), linetype = "dashed", color = "grey50") +
  labs(title = "Meta-analysis Volcano Plot: Focus on ANO10, HK2, CCDC50", x = "Average Log2(Fold Change)", 
       y = "-Log10(RRA Score)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold")) + coord_cartesian(xlim = c(-5000, 5000))
print(meta_volcano_plot)


#figure 4B heatmap plot
#need to run the "figure 1D volcano plot" part of the code in "figure 1 code"
top_n <- 20
genes_to_force_include <- c("ANO10", "HK2", "CCDC50")
top_n_genes_df <- all_core_genes_df %>%
  filter(Score < 0.05) %>%
  mutate(Status = ifelse(Avg_logFC > 0, "Upregulated", "Downregulated")) %>%
  group_by(Status) %>%
  arrange(Score, .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()
forced_genes_df <- all_core_genes_df %>%
  filter(Name %in% genes_to_force_include) %>%
  mutate(Status = ifelse(Avg_logFC > 0, "Upregulated", "Downregulated"))
top_heatmap_genes_df <- bind_rows(top_n_genes_df, forced_genes_df) %>%
  distinct(Name, .keep_all = TRUE) %>%
  arrange(desc(Status), Score) 
heatmap_matrix <- sapply(all_limma_results, function(study_result) {
  full_table <- study_result$full_table
  full_table[match(top_heatmap_genes_df$Name, rownames(full_table)), "logFC"]})
rownames(heatmap_matrix) <- top_heatmap_genes_df$Name
annotation_row_df <- data.frame(
  Avg_logFC = top_heatmap_genes_df$Avg_logFC,
  `-log10(Score)` = -log10(top_heatmap_genes_df$Score),
  row.names = top_heatmap_genes_df$Name)
annotation_row_df$Highlight <- factor(
  ifelse(rownames(annotation_row_df) %in% genes_to_force_include, "Core Gene", "Other"),
  levels = c("Core Gene", "Other"))
main_heatmap_colors <- colorRampPalette(c("#008000", "#f7f7f7", "#EE4B2B"))(100)
annotation_colors_list <- list(
  Avg_logFC = colorRampPalette(c("#4393c3", "white", "#d6604d"))(100),
  `-log10(Score)` = colorRampPalette(c("white", "#4393c3"))(100),
  Highlight = c(`Core Gene` = "gold", Other = "white"))
gap_position <- sum(top_heatmap_genes_df$Status == "Upregulated")
heatmap_plot_final <- pheatmap(
  mat = heatmap_matrix,
  na_col = "grey85",
  color = main_heatmap_colors, 
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",            
  fontsize_number = 8, 
  cluster_rows = FALSE,              
  cluster_cols = FALSE,              
  gaps_row = gap_position,
  annotation_row = annotation_row_df,
  annotation_colors = annotation_colors_list, 
  angle_col = "45",                  
  border_color = "white",            
  main = "Top DEGs with Core Genes Forced Included",
  fontsize_row = 10,
  fontsize_col = 10)
print(heatmap_plot_final)


#figure 4C UMAP plot
#need to run the "preprocessing of sc-RNA data" part of the code in "figure 2 code"
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
#define drawing function
plot_enhanced_feature <- function(seurat_object, gene, dataset_name, pt_size = 0.1) {
  if (!gene %in% rownames(seurat_object)) {
    return(ggplot() + labs(title = paste(dataset_name, "-", gene), subtitle = "Gene not found") + 
             theme_void())}
  p <- FeaturePlot(seurat_object,
                   features = gene,
                   order = TRUE, 
                   pt.size = pt_size,
                   cols = c("grey90", "#FFFF00", "#D93A49"), 
                   min.cutoff = 'q1', 
                   max.cutoff = 'q99') +
    labs(title = paste(dataset_name, "-", gene)) +
    theme(
      plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      panel.border = element_rect(colour = "black", fill=NA, size=1))
  return(p)}
genes_of_interest <- c("ANO10", "HK2", "CCDC50")
#execute function
mdd_expression_plots <- lapply(genes_of_interest, function(gene) {
  plot_enhanced_feature(seurat_object = MDD, gene = gene, dataset_name = "MDD")})
yiyv_expression_plots <- lapply(genes_of_interest, function(gene) {
  plot_enhanced_feature(seurat_object = yiyv, gene = gene, dataset_name = "yiyv")})
mdd_row <- wrap_plots(mdd_expression_plots, nrow = 1)
yiyv_row <- wrap_plots(yiyv_expression_plots, nrow = 1)
final_plot <- mdd_row / yiyv_row +
  plot_annotation(title = "Core Gene Expression in Single-Cell Datasets",
                  theme = theme(plot.title    = element_text(hjust = 0.5, size = 20, face = "bold"),
                                plot.subtitle = element_text(hjust = 0.5, size = 16)))
print(final_plot)


#figure 4D violin plot
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)
library(patchwork)
library(stringr)
genes_gse <- c("ANO10", "HK2", "CCDC50")
genes_rat <- c("Ano10", "Hk2", "Ccdc50")
my_color_palette <- c("Control"="#497cc2",
                      "Maltodextrin"="#39c5bb",
                      "Probiotics"="#ffb6c1",
                      "Desipramine"="#2A9D8E",
                      "Bupropion"="#9999FF",
                      "FT"="#ff8800")
fluoxetine_color_palette <- c("Sham"= "#497cc2","Fluoxetine" = "#ff8800")
#create plotting functions for GSE98793 and GSE76826
plot_gse_expression <- function(dataset, dataset_name, genes_to_plot) {
  genes_found <- genes_to_plot[genes_to_plot %in% rownames(dataset)]
  if (length(genes_found) == 0) {
    return(ggplot() + labs(title = paste(dataset_name, "- No Genes Found")) + theme_void())}
  sample_info <- data.frame(Sample = colnames(dataset))
  sample_info <- sample_info %>%
    mutate(Group = case_when(
      grepl("control", Sample, ignore.case = TRUE) ~ "Control",
      grepl("MDD", Sample, ignore.case = TRUE) ~ "MDD",
      TRUE ~ "Other")) %>%
    filter(Group %in% c("Control", "MDD")) %>%
    mutate(Group = factor(Group, levels = c("Control", "MDD")))
  plot_df <- dataset[genes_found, sample_info$Sample, drop = FALSE] %>%
    rownames_to_column("Gene") %>%
    pivot_longer( cols = -Gene,names_to = "Sample",
                  values_to = "Expression") %>%
    left_join(sample_info, by = "Sample")
  p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
    geom_rect(data = . %>% distinct(Group),
              aes(xmin = as.numeric(Group) - 0.5,xmax = as.numeric(Group) + 0.5,
                  ymin = -Inf,ymax = Inf,fill = Group),
              alpha = 0.2,inherit.aes = FALSE) +
    geom_violin(trim = FALSE, alpha = 1) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    facet_wrap(~ Gene, scales = "free_y", ncol = length(genes_found)) +
    stat_compare_means(comparisons = list(c("Control", "MDD")),
                       method = "wilcox.test",label = "p.signif",
                       label.y.npc = 0.9) +
    scale_fill_manual(values = c("Control" = "#497cc2", "MDD" = "#e50f4b")) +
    labs( title = paste("Gene Expression in", dataset_name),
          x = "Group",y = "Normalized Expression") +
    theme_bw(base_size = 12) +
    theme( plot.title = element_text(hjust = 0.5, face = "bold"),
           strip.background = element_blank(),
           strip.text = element_text(face = "bold.italic"),
           legend.position = "none",
           axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)}
#create plotting functions for GSE194289
plot_fluoxetine_expression <- function(dataset, dataset_name, genes_to_plot) {
  genes_found <- genes_to_plot[genes_to_plot %in% rownames(dataset)]
  if (length(genes_found) == 0) {
    return(ggplot() + labs(title = paste(dataset_name, "- No Genes Found")) + theme_void())}
  sample_info <- data.frame(Sample = colnames(dataset))
  sample_info <- sample_info %>%
    mutate(Group = case_when(
      grepl("Sham", Sample, ignore.case = TRUE) ~ "Sham",
      grepl("Fluoxetine", Sample, ignore.case = TRUE) ~ "Fluoxetine",
      TRUE ~ "Other")) %>%
    filter(Group %in% c("Sham", "Fluoxetine")) %>%
    mutate(Group = factor(Group, levels = c("Sham", "Fluoxetine")))
  plot_df <- dataset[genes_found, sample_info$Sample, drop = FALSE] %>%
    rownames_to_column("Gene") %>%
    pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
    left_join(sample_info, by = "Sample")
  p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
    geom_rect(data = . %>% distinct(Group),
              aes(xmin = as.numeric(Group) - 0.5,xmax = as.numeric(Group) + 0.5,
                  ymin = -Inf,ymax = Inf, fill = Group),
              alpha = 0.2, inherit.aes = FALSE) +
    geom_violin(trim = FALSE, alpha = 1) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    facet_wrap(~ Gene, scales = "free_y", ncol = length(genes_found)) +
    stat_compare_means(comparisons = list(c("Sham", "Fluoxetine")),
                       method = "wilcox.test",label = "p.signif",
                       label.y.npc = 0.9) +
    scale_fill_manual(values = fluoxetine_color_palette) +
    labs(title = paste("Gene Expression in", dataset_name),
         x = "Group", y = "Expression Count") +
    theme_bw(base_size = 12) +
    theme( plot.title = element_text(hjust = 0.5, face = "bold"),
           strip.background = element_blank(),
           strip.text = element_text(face = "bold.italic"),
           legend.position = "none",
           axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)}
#create plotting functions for GSE222756
plot_probiotics_expression <- function(dataset, dataset_name, genes_to_plot) {
  genes_found <- genes_to_plot[genes_to_plot %in% rownames(dataset)]
  if (length(genes_found) == 0) {
    return(ggplot() + labs(title = paste(dataset_name, "- No Genes Found")) + theme_void())}
  sample_info <- data.frame(Sample = colnames(dataset))
  sample_info <- sample_info %>%
    mutate(Group = str_extract(Sample, "(FT|Bupropion|Desipramine|Probiotics|Maltodextrin|Control)$"))
  group_levels <- c("Control", "Maltodextrin", "Probiotics", "Desipramine", "Bupropion", "FT")
  sample_info$Group <- factor(sample_info$Group, levels = group_levels)
  plot_df <- dataset[genes_found, sample_info$Sample, drop = FALSE] %>%
    rownames_to_column("Gene") %>%
    pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
    left_join(sample_info, by = "Sample")
  p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
    geom_rect(data = . %>% distinct(Group),
              aes(xmin = as.numeric(Group) - 0.5,xmax = as.numeric(Group) + 0.5,
                  ymin = -Inf,ymax = Inf,fill = Group),
              alpha = 0.2,inherit.aes = FALSE) +
    geom_violin(trim = FALSE, alpha = 1) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    facet_wrap(~ Gene, scales = "free_y", ncol = length(genes_found)) +
    stat_compare_means(method = "wilcox.test",ref.group = "Control",
                       label = "p.signif",hide.ns = FALSE,vjust = -0.2) +
    scale_fill_manual(values = my_color_palette) +
    labs(title = paste("Gene Expression in", dataset_name),
         x = "Group",y = "Expression Count") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold.italic"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)}
#execute all function
p_gse98793 <- plot_gse_expression(GSE98793, "GSE98793", genes_gse)
p_gse76826 <- plot_gse_expression(GSE76826, "GSE76826", genes_gse)
p_fluoxetine <- plot_fluoxetine_expression(`GSE194289`, "GSE194289", genes_rat)
p_probiotics <- plot_probiotics_expression(`GSE222756`, "GSE222756", genes_rat)
final_plot <- (p_gse98793 + p_gse76826) / (p_fluoxetine + p_probiotics)
final_plot <- final_plot + plot_annotation(
  title = 'Expression of ANO10,HK2,CCDC50 Across Datasets',
  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face="bold")))
print(final_plot)


#figure 4F sankey bubble plot
library(readr)
library(tidyverse)
library(ggnewscale) 
library(ggsankey) 
library(cowplot)
RRAresult <- read_csv("RRA方法DEG通路富集结果.csv") %>% mutate(method = "RRA method")
geneintresult <- read_csv("基因交集法DEG通路富集结果.csv") %>% mutate(method = "Overlap method")
combined_data <- bind_rows(RRAresult, geneintresult)
#select some pathways from the pathway results
selected_pathway_names <- c(
  # GO
  "cell death",
  "cellular component organization or biogenesis",
  "cellular response to cytokine stimulus",
  "multicellular organism development",
  "nervous system process",
  "protein localization",
  "response to cytokine",
  "response to stimulus",
  "response to stress",
  # KEGG
  "HIF-1 signaling pathway",
  "Insulin signaling pathway",
  "Metabolic pathways",
  # REAC
  "Disease",
  "Glucose metabolism",
  "Metabolism of carbohydrates and carbohydrate derivatives",
  "Programmed Cell Death")
target_genes_sankey <- c("ANO10", "HK2", "CCDC50", "CD163", "GBE1", "TIMP2")
clean_term_name <- function(name) {
  name %>% str_remove("^[A-Z]{2,4}:[A-Z0-9]+_") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()}
top_pathways <- combined_data %>%
  filter(term_name %in% selected_pathway_names) %>%
  group_by(term_name) %>%
  summarise(min_p_value = min(adjusted_p_value), .groups = "drop") %>%
  filter(term_name %in% selected_pathway_names)
pathway_counts <- combined_data %>%
  filter(term_name %in% top_pathways$term_name) %>%
  group_by(term_name) %>%
  slice_min(order_by = adjusted_p_value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(term_name, count = intersection_size)
pathway_y_positions <- top_pathways %>%
  left_join(pathway_counts, by = "term_name") %>%
  mutate(term_name_clean = clean_term_name(term_name)) %>%
  arrange(min_p_value) %>%
  slice(n():1) %>%
  mutate(ymax = cumsum(count),
         ymin = ymax - count,
         label_pos = (ymin + ymax) / 2)
sankey_links_data <- combined_data %>%
  filter(term_name %in% top_pathways$term_name) %>%
  separate_rows(intersections, sep = ",") %>%
  filter(intersections %in% target_genes_sankey) %>%
  mutate(term_name_clean = clean_term_name(term_name)) %>%
  select(gene = intersections, pathway = term_name_clean) %>%
  distinct() %>%
  bind_rows(tibble(
    gene    = NA_character_,
    pathway = setdiff(top_pathways$term_name %>% clean_term_name(), unique(.$pathway))))
plot_data_bubble <- combined_data %>%
  filter(term_name %in% top_pathways$term_name) %>%
  mutate(term_name_clean = clean_term_name(term_name)) %>%
  left_join(select(pathway_y_positions, term_name_clean, label_pos), by = "term_name_clean") %>%
  mutate(plot_value = if_else(method == "Overlap method",
                              -negative_log10_of_adjusted_p_value,
                              negative_log10_of_adjusted_p_value))
#draw sankey plot
if (any(!is.na(sankey_links_data$gene))) {
  sankey_df_long <- make_long(sankey_links_data, gene, pathway)
  pathway_order <- pathway_y_positions$term_name_clean
  gene_order <- sort(unique(sankey_links_data$gene), na.last = TRUE)
  sankey_df_long$node <- factor(sankey_df_long$node, levels = c(gene_order, pathway_order))
  sankey_df_long$next_node <- factor(sankey_df_long$next_node, levels = c(gene_order, pathway_order))
  p_sankey <- ggplot(sankey_df_long, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = "black", show.legend = FALSE) +
    geom_sankey_text(size = 3.5, color = "black") +
    theme_void() +
    theme(plot.margin = unit(c(1, 20, 1, 1), "cm"))
} else {
  p_sankey <- ggplot() + 
    annotate("text", x=1, y=1) +
    theme_void() +
    theme(plot.margin = unit(c(1, 20, 1, 1), "cm"))}
#draw bubble plot
p_bubble <- ggplot(plot_data_bubble, aes(x = plot_value, y = label_pos)) +
  geom_point(data = . %>% filter(method == "Overlap method"), aes(size = intersection_size, color = adjusted_p_value)) +
  scale_color_gradient(name = "Overlap method Adj. P-value", low = "#0e3bf0", high = "#cdc1f2", guide = guide_colorbar(order = 1)) +
  new_scale_color() +
  geom_point(data = . %>% filter(method == "RRA method"), aes(size = intersection_size, color = adjusted_p_value)) +
  scale_color_gradient(name = "RRA method Adj. P-value", low = "#de1a85", high = "#f2c4d5", guide = guide_colorbar(order = 2)) +
  scale_size_continuous(name = "Gene Count") +
  scale_x_continuous(name = "-log10(Adjusted P-value)", labels = abs, expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  labs(y = NULL) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank())
final_plot <- ggdraw() +
  draw_plot(p_sankey) +
  draw_plot(p_bubble, x = 0.45, y = 0.025, width = 0.55, height = 0.95)
print(final_plot)



