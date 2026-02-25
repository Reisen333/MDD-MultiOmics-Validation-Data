#Supplementary Fig.1 PCA plot
#need to run the "DEG analysis" part of the code in "figure 1 code" first
#create PCA plotting function
get_sex <- function(sample_names) {
  sex <- sub(".*_(M|F)_.*", "\\1", sample_names)
  sex[!sex %in% c("M", "F")] <- NA 
  return(sex)}
plot_pca_on_degs <- function(expression_data, deg_list, group_info, plot_title, sex_info = NULL) {
  if (length(deg_list) < 2) {
    return(ggplot() + labs(title = plot_title, subtitle = "DEG deficiency") + theme_void())}
  pca_input <- t(expression_data[deg_list, ])
  pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)
  pca_plot_data <- data.frame(pca_result$x)
  percent_variance <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
  pca_plot_data$group <- group_info
  if (!is.null(sex_info) && length(unique(na.omit(sex_info))) > 1) {
    pca_plot_data$sex <- factor(sex_info)
    p <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = group, shape = sex)) +
      geom_point(size = 4, alpha = 0.8) +
      scale_color_manual(values = c("Control" = "#497cc2", "MDD" = "#e50f4b"), name="Group") +
      scale_shape_manual(values = c("M" = 15, "F" = 17), name="Sex", na.translate = FALSE) + 
      guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))
  } else {
    p <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = group, shape = group)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(values = c("Control" = "#497cc2", "MDD" = "#e50f4b"), name="Group") +
      scale_shape_manual(values = c("Control" = 16, "MDD" = 17), name="Group")}
  p <- p + labs(
    title = plot_title,
    x = paste0("PC1 (", percent_variance[1], "%)"),
    y = paste0("PC2 (", percent_variance[2], "%)")
  ) +  theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom")
  return(p)}
pca_plots_list <- lapply(names(all_limma_results), function(name) {
  result <- all_limma_results[[name]]
  subset_data <- result$subset_data
  deg_list <- result$gene_lists$all
  samples <- colnames(subset_data)
  groups <- factor(ifelse(grepl(dataset_info[[which(sapply(dataset_info, `[[`, "name") == name)]]$control_pattern, samples), "Control", "MDD"))
  sex <- get_sex(samples)
  plot_pca_on_degs(
    expression_data = subset_data, 
    deg_list = deg_list, 
    group_info = groups, 
    plot_title = name, 
    sex_info = sex)})
final_pca_plot <- wrap_plots(pca_plots_list, ncol = 2) +
  plot_annotation(
    title = "PCA of Samples using DEGs (Color by Group, Shape by Sex)",
    tag_levels = 'A',
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
print(final_pca_plot)


#Supplementary Fig.2 volcano plot
#create volcano plotting function
plot_volcano <- function(limma_res, p_cutoff = 0.05, fc_cutoff = log2(1.1)) {
  df <- limma_res$full_table
  df$change <- "Stable"
  df$change[df$P.Value < p_cutoff & df$logFC > fc_cutoff] <- "Up-regulated"
  df$change[df$P.Value < p_cutoff & df$logFC < -fc_cutoff] <- "Down-regulated"
  df$change <- factor(df$change, levels = c("Up-regulated", "Down-regulated", "Stable"))
  my_colors <- c("Up-regulated" = "#e50f4b", 
                 "Down-regulated" = "#497cc2", 
                 "Stable" = "grey80")
  p <- ggplot(df, aes(x = logFC, y = -log10(P.Value), color = change)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = my_colors) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black", linewidth = 0.4) +
    labs(title = limma_res$dataset_name,
         x = expression(log[2]("Fold Change")),
         y = expression(-log[10]("P Value")),
         color = "Status") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  return(p)}
#execute function
volcano_plots_list <- lapply(all_limma_results, function(res) {
  plot_volcano(res, p_cutoff = 0.05, fc_cutoff = log2(1.1))})
combined_volcano <- wrap_plots(volcano_plots_list) + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "Volcano Plots of Differential Expression",
    subtitle = "Threshold: p < 0.05, |log2FC| > log2(1.1) (Red: Up, Blue: Down)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 12)))
print(combined_volcano)


#Supplementary Fig.3A box plot
#need to run the "bayesian deconvolution analysis" part of the code in "figure 2 code" first
library(BayesPrism)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
bp.result <- MDD对GSE98793贝叶斯结果
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
cell_to_plot <- "OLI"   #here can changed cell type
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


#Supplementary Fig.3B stacked bar plot
#need to run the "bayesian deconvolution analysis" part of the code in "figure 2 code" first
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
gse98793_groups <- quote(ifelse(grepl("control", Sample, ignore.case = TRUE), "Control", "MDD"))
plot_gse98793 <- create_bayes_plot(MDD对GSE98793贝叶斯结果, gse98793_groups, "GSE98793", my_colors)
print(plot_gse98793)
gse76826_groups <- quote(ifelse(grepl("control", Sample, ignore.case = TRUE), "Control", "MDD"))
plot_gse76826 <- create_bayes_plot(MDD对GSE76826贝叶斯结果, gse76826_groups, "GSE76826", my_colors)
print(plot_gse76826)


#Supplementary Fig.3C box plot
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
cell_to_plot <- "MIC"    #here can changed cell type
plot_data <- theta_df %>%
  filter(cell_type == cell_to_plot)
ggplot(plot_data, aes(x = Group, y = proportion, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 1) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.2) +
  scale_fill_manual(values = c("Control"="#497cc2", "Maltodextrin"="#39c5bb", "Probiotics"="#ffb6c1", 
                               "Desipramine"="#2A9D8E", "Bupropion"="#9999FF", "FT"="#ff8800", "Fluoxetine" = "#ff8800")) +
  stat_compare_means(method = "kruskal.test", label.y.npc = 0.9) + # 多组比较用Kruskal-Wallis检验
  labs(title = paste("Proportion of", cell_to_plot, "across treatments"),
       x = "Treatment Group",
       y = "Cell Type Proportion") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")


#Supplementary Fig.4 circle plot
#need to run the "cellchat analysis" part of the code in "figure 2 code" first
cell_types_healthy <- levels(cellchat_healthy@idents)
colors_for_healthy <- my_colors[cell_types_healthy]
cell_types_dep <- levels(cellchat_dep@idents)
colors_for_dep <- my_colors[cell_types_dep]
groupSize_healthy <- as.numeric(table(cellchat_healthy@idents))
groupSize_dep <- as.numeric(table(cellchat_dep@idents))
par(mfrow = c(1, 2), xpd=TRUE)
netVisual_circle(
  cellchat_healthy@net$weight,
  vertex.weight = groupSize_healthy,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = colors_for_healthy,
  title.name = "Healthy: Interaction Strength")
par(mfrow = c(1, 2), xpd=TRUE) 
netVisual_circle(
  cellchat_dep@net$weight,
  vertex.weight = groupSize_dep,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = colors_for_dep,
  title.name = "Depression: Interaction Strength")


#Supplementary Fig.5 and Fig.6 violin plot
#need to run the "DEG analysis" part of the code in "figure 1 code" first
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(ggpubr)
oli_genes_of_interest <- c("FUT8", "PTPRK", "ZNF532", "CHIC1", "CCDC50")
mic_genes_of_interest <- c("ANO10", "ASPH", "ATP8B4", "CD163", "FGD4", "JDP2", "IRAK3", 
                           "IL13RA1", "HK2", "MTMR3", "MPP1", "PLXDC2", "TFEB", 
                           "ST3GAL6", "SLC25A40", "GBE1", "MAPK14", "TIMP2")
all_genes_of_interest <- unique(c(oli_genes_of_interest, mic_genes_of_interest))
gene_bg_info <- data.frame(Gene = all_genes_of_interest) %>%
  mutate(Type = ifelse(Gene %in% oli_genes_of_interest, "OLI", "MIC")) %>%
  mutate(BgColor = ifelse(Type == "OLI", "#e50f4b", "#497cc2"))
ordered_genes <- c(sort(oli_genes_of_interest), sort(mic_genes_of_interest))
#create a drawing function
plot_gene_expression_combined <- function(dataset_name, limma_result, genes_to_plot, gene_order, bg_map) {
  expression_data <- limma_result$subset_data
  sample_names <- colnames(expression_data)
  genes_found <- genes_to_plot[genes_to_plot %in% rownames(expression_data)]
  if (length(genes_found) == 0) {
    return(ggplot() + labs(title = dataset_name, subtitle = "No target genes found") + theme_void())}
  info <- dataset_info[[which(sapply(dataset_info, `[[`, "name") == dataset_name)]]
  groups <- factor(
    ifelse(grepl(info$control_pattern, sample_names), "Control", "MDD"),
    levels = c("Control", "MDD"))
  plot_df <- expression_data[genes_found, , drop = FALSE] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    mutate(Group = groups) %>%
    pivot_longer(cols = all_of(genes_found), names_to = "Gene", values_to = "Expression")
  plot_df$Gene <- factor(plot_df$Gene, levels = gene_order)
  current_bg_data <- bg_map %>% filter(Gene %in% levels(plot_df$Gene))
  current_bg_data$Gene <- factor(current_bg_data$Gene, levels = levels(plot_df$Gene))
  p <- ggplot(plot_df, aes(x = Group, y = Expression)) +
    geom_rect(data = current_bg_data,
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = I(BgColor)),
              alpha = 0.15, inherit.aes = FALSE) +
    geom_violin(aes(fill = Group), trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    stat_compare_means(comparisons = list(c("Control", "MDD")), 
                       method = "wilcox.test", 
                       label = "p.signif") +
    facet_wrap(~ Gene, scales = "free_y", ncol = 7) + 
    scale_fill_manual(values = c("Control" = "#497cc2", "MDD" = "#e50f4b")) +
    labs(title = paste("Gene Expression in", dataset_name, "(Oli: Red BG, Mic: Blue BG)"),
         x = NULL,
         y = "Normalized Expression") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.background = element_rect(fill = "grey95", color = "black"),
          strip.text = element_text(face = "bold.italic"),
          legend.position = "none")
  return(p)}
#use the drawing function
p_gse98793_combined <- plot_gene_expression_combined(
  dataset_name = "GSE98793",
  limma_result = all_limma_results$GSE98793,
  genes_to_plot = all_genes_of_interest,
  gene_order = ordered_genes,           
  bg_map = gene_bg_info)
p_gse76826_combined <- plot_gene_expression_combined(
  dataset_name = "GSE76826",
  limma_result = all_limma_results$GSE76826,
  genes_to_plot = all_genes_of_interest,
  gene_order = ordered_genes,            
  bg_map = gene_bg_info)
final_combined_plot <- p_gse98793_combined / p_gse76826_combined +
  plot_annotation(
    title = 'Expression of OLI & MIC Genes in Human Blood Samples',
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)))
print(final_combined_plot)


#Supplementary Fig.7 and Fig.8 violin plot
library(tidyverse) 
library(patchwork) 
library(ggpubr)    
library(stringr)
GSE222756 <- read.delim("C:/R/Rworking/抑郁症数据/单细胞 抗抑郁药物和益生菌作用的共同及独特转录组特征/益生菌bulk数据修改.txt", row.names = 1, sep = '\t', check.names = FALSE)
GSE194289 <- readRDS("C:/R/Rworking/抑郁症数据/单细胞 氟西汀作用在 27 个大脑区域的综合多组学景观揭示了能量代谢和区域特异性染色质重塑的整体增加/GSE194289基因名称修改.rds")
rat_oli_genes <- c("Fut8", "Ptprk", "Zfp532", "Tsx", "Ccdc50")
rat_mic_genes <- c("Ano10", "Asph", "Atp8b4", "Cd163", "Fgd4", "Jdp2", "Irak3", 
                   "Il13ra1-ps1", "Il13ra1", "Hk2", "Mtmr3", "Plxdc2", "Tfeb", 
                   "St3gal6", "Slc25a40", "Gbe1", "Mapk14")
all_target_genes <- unique(c(rat_oli_genes, rat_mic_genes))
gene_bg_info <- data.frame(Gene = all_target_genes) %>%
  mutate(Type = ifelse(Gene %in% rat_oli_genes, "OLI", "MIC")) %>%
  mutate(BgColor = ifelse(Type == "OLI", "#e50f4b", "#497cc2"))
#create a drawing function
plot_bulk_expression_with_bg <- function(plot_df, title, stat_comparisons = NULL, color_palette = NULL, bg_map = NULL) {
  genes_found <- unique(plot_df$Gene)
  if (length(genes_found) == 0) {
    return(ggplot() + labs(title = title, subtitle = "No target genes found") + theme_void())}
  current_bg_data <- bg_map %>% filter(Gene %in% genes_found)
  if (is.factor(plot_df$Gene)) {
    current_bg_data$Gene <- factor(current_bg_data$Gene, levels = levels(plot_df$Gene))}
  p <- ggplot(plot_df, aes(x = Group, y = Expression)) +
    geom_rect(data = current_bg_data, 
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = I(BgColor)), 
              alpha = 0.15, 
              inherit.aes = FALSE) +
    geom_violin(aes(fill = Group), trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    facet_wrap(~ Gene, scales = "free_y", ncol = 6) +
    labs(title = title, x = NULL, y = "Log-Normalized Expression (log1p)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold.italic"),
          legend.position = "bottom",
          legend.title = element_blank())
  if (!is.null(stat_comparisons)) {
    p <- p + stat_compare_means(
      comparisons = stat_comparisons,
      method = "wilcox.test",
      label = "p.signif",
      vjust = 0.5, hide.ns = FALSE)}
  if (!is.null(color_palette)) {
    p <- p + scale_fill_manual(values = color_palette)}
  return(p)}
#use the drawing function
ordered_genes <- c(sort(rat_oli_genes), sort(rat_mic_genes))
fluoxetine_long_df <- GSE194289 %>%
  rownames_to_column("Gene") %>%
  filter(Gene %in% all_target_genes) %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  mutate(Group = factor(ifelse(str_detect(Sample, "Sham"), "Sham", "Fluoxetine"), levels = c("Sham", "Fluoxetine")),
         Expression = log1p(Expression))
fluoxetine_long_df$Gene <- factor(fluoxetine_long_df$Gene, levels = ordered_genes)
plot_fluoxetine <- plot_bulk_expression_with_bg(
  plot_df = fluoxetine_long_df,
  title = "Expression of Oli (Red) & Mic (Blue) Genes (GSE194289)",
  stat_comparisons = list(c("Sham", "Fluoxetine")),
  color_palette = fluoxetine_colors,
  bg_map = gene_bg_info)
probiotics_long_df <- GSE222756 %>%
  rownames_to_column("Gene") %>%
  filter(Gene %in% all_target_genes) %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  mutate(Group = case_when(
    str_detect(Sample, "Control") ~ "Control", str_detect(Sample, "Maltodextrin") ~ "Maltodextrin",
    str_detect(Sample, "Probiotics") ~ "Probiotics", str_detect(Sample, "Desipramine") ~ "Desipramine",
    str_detect(Sample, "Bupropion") ~ "Bupropion", str_detect(Sample, "FT") ~ "FT",
    TRUE ~ "Unknown"), Expression = log1p(Expression))
group_order <- c("Control", "Maltodextrin", "Probiotics", "Desipramine", "Bupropion", "FT")
probiotics_long_df$Group <- factor(probiotics_long_df$Group, levels = group_order)
probiotics_long_df$Gene <- factor(probiotics_long_df$Gene, levels = ordered_genes)
comparison_groups <- group_order[group_order != "Control"]
probiotics_comparisons <- lapply(comparison_groups, function(x) c("Control", x))
plot_probiotics <- plot_bulk_expression_with_bg(
  plot_df = probiotics_long_df,
  title = "Expression of Oli (Red) & Mic (Blue) Genes (GSE222756)",
  stat_comparisons = probiotics_comparisons,
  color_palette = probiotics_colors,
  bg_map = gene_bg_info)
final_combined_plot <- plot_fluoxetine / plot_probiotics +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(
    title = "Comparative Gene Expression with Cell-Type Backgrounds",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)))
print(final_combined_plot)


