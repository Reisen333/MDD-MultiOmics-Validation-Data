#figure 3A tile plot
#need to run the "preprocessing of sc-RNA data" and "figure 2B UMAP plot" part of the code in "figure 2 code" ,and "gene intersection analysis" part of the code in "figure 1 code"
Idents(MDD) <- "cell_type_simplified" 
mdd_markers <- FindAllMarkers(MDD, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
Idents(yiyv) <- "cell_type"
yiyv_markers <- FindAllMarkers(yiyv, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
get_overlap_for_plotting <- function(marker_df, deg_list, n_top = 500) {
  overlap_df <- marker_df %>%
    dplyr::filter(p_val_adj < 0.05 & gene %in% deg_list) %>% 
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n = n_top, order_by = avg_log2FC) %>%
    dplyr::ungroup() %>% dplyr::select(gene, cluster) %>% dplyr::distinct()
  return(overlap_df)}
#create drawing functions
create_overlap_plot <- function(overlap_df, title) {
  if (nrow(overlap_df) == 0) {
    return(ggplot() + labs(title = title, subtitle = "No overlapping genes found") + theme_void())}
  count_df <- overlap_df %>% dplyr::count(cluster, name = "count")
  plot_df <- left_join(overlap_df, count_df, by = "cluster")
  p <- ggplot(plot_df, aes(x = cluster, y = gene, fill = count)) +
    geom_tile(color = "white", lwd = 0.5) +
    scale_fill_viridis(name = "Overlap Count") +
    labs(
      title = title,
      x = "Cell Type",
      y = "Bulk DEG Gene") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1,size = 11),
      axis.text.y = element_text(size = 11))
  return(p)}
mdd_up_overlap <- get_overlap_for_plotting(mdd_markers, intersect_up_degs)
mdd_down_overlap <- get_overlap_for_plotting(mdd_markers, intersect_down_degs)
yiyv_up_overlap <- get_overlap_for_plotting(yiyv_markers, intersect_up_degs)
yiyv_down_overlap <- get_overlap_for_plotting(yiyv_markers, intersect_down_degs)
p_mdd_up <- create_overlap_plot(mdd_up_overlap, "MDD: Upregulated DEGs vs Cell Markers")
p_mdd_down <- create_overlap_plot(mdd_down_overlap, "MDD: Downregulated DEGs vs Cell Markers")
p_yiyv_up <- create_overlap_plot(yiyv_up_overlap, "yiyv: Upregulated DEGs vs Cell Markers")
p_yiyv_down <- create_overlap_plot(yiyv_down_overlap, "yiyv: Downregulated DEGs vs Cell Markers")
(p_mdd_up | p_yiyv_up) / (p_mdd_down | p_yiyv_down)


#figure 3B heatmap plot
#need to run the "RRA analysis" part of the code in "figure 1 code"
library(dplyr)
library(pheatmap)
oli_genes <- c("FUT8", "PTPRK", "ZNF532","CHIC1","CCDC50")
mic_genes <- c("ANO10", "ASPH", "ATP8B4", "CD163", "FGD4", "JDP2", "IRAK3",
               "IL13RA1", "HK2", "MTMR3", "MPP1", "PLXDC2", "TFEB",
               "ST3GAL6", "SLC25A40", "GBE1", "MAPK14", "TIMP2")
mic_genes_df <- subset(final_up_summary, Name %in% mic_genes) %>%
  mutate(Rank = match(Name, rra_up_results$Name),
         Status = "Upregulated")                  
oli_genes_df <- subset(final_down_summary, Name %in% oli_genes) %>%
  mutate(Rank = match(Name, rra_down_results$Name),
         Status = "Downregulated"   ) 
custom_heatmap_genes_df <- rbind(mic_genes_df, oli_genes_df) %>%
  arrange(desc(Status)) 
heatmap_matrix_custom <- sapply(all_limma_results, function(study_result) {
  full_table <- study_result$full_table
  full_table[match(custom_heatmap_genes_df$Name, rownames(full_table)), "logFC"]})
rownames(heatmap_matrix_custom) <- custom_heatmap_genes_df$Name
annotation_row_df_custom <- data.frame(
  `Avg_logFC` = custom_heatmap_genes_df$Avg_logFC,
  `-log10(Score)` = -log10(custom_heatmap_genes_df$Score),
  row.names = custom_heatmap_genes_df$Name)
main_heatmap_colors <- colorRampPalette(c("#008000", "#f7f7f7", "#EE4B2B"))(100)
annotation_colors_list_custom <- list(
  `Avg_logFC` = colorRampPalette(c("#4393c3", "white", "#d6604d"))(100),
  `-log10(Score)` = colorRampPalette(c("white", "#4393c3"))(100))
custom_row_labels <- paste0(custom_heatmap_genes_df$Name, " (Rank: ", custom_heatmap_genes_df$Rank, ")")
gap_position_custom <- length(mic_genes)
heatmap_plot_custom <- pheatmap(
  mat = heatmap_matrix_custom,
  na_col = "grey85",
  color = main_heatmap_colors,     
  cluster_rows = FALSE,                
  cluster_cols = TRUE,                 
  gaps_row = gap_position_custom,     
  annotation_row = annotation_row_df_custom,
  annotation_colors = annotation_colors_list_custom, 
  labels_row = custom_row_labels,    
  angle_col = "45",
  border_color = "white",
  main = "Heatmap of Specific Genes with RRA Ranks",
  fontsize_row = 10,
  fontsize_col = 10)
print(heatmap_plot_custom)


#figure 3C line plot
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(tibble)
library(stringr)
library(rstatix)
GSE98793 <- readRDS("C:/R/Rworking/抑郁症数据/血液样本/GSE98793.rds")
GSE76826 <- readRDS("C:/R/Rworking/抑郁症数据/血液样本/GSE76826.rds")
GSE222756 <- read.delim("C:/R/Rworking/抑郁症数据/单细胞 抗抑郁药物和益生菌作用的共同及独特转录组特征/益生菌bulk数据修改.txt", row.names = 1, sep = '\t', check.names = FALSE)
GSE194289 <- readRDS("C:/R/Rworking/抑郁症数据/单细胞 氟西汀作用在 27 个大脑区域的综合多组学景观揭示了能量代谢和区域特异性染色质重塑的整体增加/GSE194289基因名称修改.rds")
mymicgenes <- c("FUT8", "CCDC50", "ZNF532")
myoligenes <- c("ANO10", "CD163", "IRAK3", "IL13RA1", "HK2", "MTMR3", "PLXDC2", "TFEB", "ST3GAL6", "SLC25A40", "MAPK14")
mygenes <- c(mymicgenes, myoligenes)
myratmicegenes <- c("Fut8", "Ccdc50", "Zfp532")
myratoligenes <- c("Ano10", "Cd163", "Irak3", "Il13ra1", "Hk2", "Mtmr3", "Plxdc2", "Tfeb", "St3gal6", "Slc25a40", "Mapk14")
myratgenes <- c(myratmicegenes, myratoligenes)
human_colors <- c("Control" = "#497cc2", "MDD" = "#e50f4b")
fluoxetine_colors <- c("Sham" = "#497cc2", "Fluoxetine" = "#ff8800")
probiotics_colors <- c("Control"="#497cc2", "Maltodextrin"="#39c5bb", "Probiotics"="#ffb6c1", 
                       "Desipramine"="#2A9D8E", "Bupropion"="#9999FF", "FT"="#ff8800")
#creat a plotting function for GSE98793,GSE76826,GSE194289
plot_line_expression_simple <- function(plot_df, title, genes_cat1, genes_cat2, ref_group, color_palette = NULL) {
  gene_order <- c(genes_cat1, genes_cat2)
  plot_df_ordered <- plot_df %>%
    mutate(Gene = factor(Gene, levels = gene_order)) %>% filter(!is.na(Gene))
  summary_df <- plot_df_ordered %>%
    group_by(Gene, Group) %>%
    summarise(mean_expression = mean(Expression, na.rm = TRUE),
              se_expression = sd(Expression, na.rm = TRUE) / sqrt(n()), .groups = 'drop')
  stat.test <- plot_df_ordered %>%
    group_by(Gene) %>%
    wilcox_test(Expression ~ Group, ref.group = ref_group) %>%
    add_significance("p") %>%
    left_join(summary_df %>% group_by(Gene) %>% summarise(y.position = max(mean_expression + se_expression, na.rm = TRUE)), by = "Gene")
  p <- ggplot(data = summary_df, aes(x = Gene, y = mean_expression, group = Group, color = Group)) +
    geom_rect(data = data.frame(xmin = c(0.5, length(genes_cat1) + 0.5), xmax = c(length(genes_cat1) + 0.5, Inf), category = c("cat1", "cat2")),
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = category), inherit.aes = FALSE, alpha = 0.2) +
    scale_fill_manual(values = c("cat1" = "#e50f4b", "cat2" = "#497cc2"), guide = "none") +
    geom_errorbar(aes(ymin = mean_expression - se_expression, ymax = mean_expression + se_expression), width = 0.5, linewidth = 1) +
    geom_line(linewidth = 1) +
    geom_point(size = 3, shape = 21, fill = "white", stroke = 1.5) +
    geom_text(data = stat.test, aes(x = Gene, y = y.position, label = p.signif), 
              inherit.aes = FALSE, size = 5, vjust = -0.5) +
    labs(title = title, x = "Gene", y = "Mean Normalized Expression", color = "Group") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "grey85", linetype = "dashed"),#dashed虚线，solid实线
          panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"))
  if (!is.null(color_palette)) { p <- p + scale_color_manual(values = color_palette) }
  return(p)}
#creat a plotting function for GSE222756
plot_line_expression_probiotics <- function(plot_df, title, genes_cat1, genes_cat2, ref_group, color_palette = NULL, label_spacing_factor = 0.05, hide_ns = TRUE) {
  gene_order <- c(genes_cat1, genes_cat2)
  plot_df_ordered <- plot_df %>%
    mutate(Gene = factor(Gene, levels = gene_order)) %>% 
    filter(!is.na(Gene))
  summary_df <- plot_df_ordered %>%
    group_by(Gene, Group) %>%
    summarise(mean_expression = mean(Expression, na.rm = TRUE),
              se_expression = sd(Expression, na.rm = TRUE) / sqrt(n()), .groups = 'drop')
  stat.test <- plot_df_ordered %>%
    group_by(Gene) %>%
    wilcox_test(Expression ~ Group, ref.group = ref_group) %>%
    add_significance("p")
  if (hide_ns) {stat.test <- stat.test %>% filter(p.signif != "ns")}
  y_positions <- summary_df %>%
    group_by(Gene) %>%
    summarise(y_ceiling = max(mean_expression + se_expression, na.rm = TRUE),
              y_range = max(mean_expression) - min(mean_expression)) %>%
    mutate(y_range = ifelse(y_range == 0, 1, y_range))
  stat.test <- stat.test %>%
    left_join(y_positions, by = "Gene") %>%
    # 计算固定的间隔
    mutate(label_interval = y_range * label_spacing_factor) %>%
    arrange(Gene, group2) %>%
    group_by(Gene) %>%
    mutate(y.position = y_ceiling + (row_number() * label_interval)) %>%
    ungroup()
  p <- ggplot(data = summary_df, aes(x = Gene, y = mean_expression, group = Group, color = Group)) +
    geom_rect(data = data.frame(xmin = c(0.5, length(genes_cat1) + 0.5), xmax = c(length(genes_cat1) + 0.5, Inf), category = c("cat1", "cat2")),
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = category), inherit.aes = FALSE, alpha = 0.2) +
    scale_fill_manual(values = c("cat1" = "#e50f4b", "cat2" = "#497cc2"), guide = "none") +
    geom_errorbar(aes(ymin = mean_expression - se_expression, ymax = mean_expression + se_expression), width = 0.5, linewidth = 1) +
    geom_line(linewidth = 1) +
    geom_point(size = 3, shape = 21, fill = "white", stroke = 1.5) +
    # 只有当 stat.test 不为空时才绘制文本
    {if(nrow(stat.test) > 0) 
      geom_text(data = stat.test, aes(x = Gene, y = y.position, label = p.signif, color = group2),
                inherit.aes = FALSE, size = 4, vjust = 0)} + # vjust设为0防止压住线
    labs(title = title, x = "Gene", y = "Mean Normalized Expression", color = "Group") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "grey85", linetype = "dashed"),
          panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"))
  if (!is.null(color_palette)) { p <- p + scale_color_manual(values = color_palette, name = "Group") }
  return(p)}
# GSE98793
plot_df_gse98793 <- GSE98793[mygenes, ] %>% t() %>% as.data.frame() %>% mutate(Group = factor(ifelse(grepl("control", rownames(.), ignore.case = TRUE), "Control", "MDD"))) %>% pivot_longer(cols = -Group, names_to = "Gene", values_to = "Expression")
line_gse98793 <- plot_line_expression_simple(plot_df = plot_df_gse98793, title = "Human Blood: GSE98793", genes_cat1 = mymicgenes, genes_cat2 = myoligenes, ref_group = "Control", color_palette = human_colors)
# GSE76826
plot_df_gse76826 <- GSE76826[mygenes, ] %>% t() %>% as.data.frame() %>% mutate(Group = factor(ifelse(grepl("control", rownames(.), ignore.case = TRUE), "Control", "MDD"))) %>% pivot_longer(cols = -Group, names_to = "Gene", values_to = "Expression")
line_gse76826 <- plot_line_expression_simple(plot_df = plot_df_gse76826, title = "Human Blood: GSE76826", genes_cat1 = mymicgenes, genes_cat2 = myoligenes, ref_group = "Control", color_palette = human_colors)
# GSE194289
plot_df_fluoxetine <- GSE194289 %>% rownames_to_column("Gene") %>% filter(Gene %in% myratgenes) %>% pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>% mutate(Group = factor(ifelse(str_detect(Sample, "Sham"), "Sham", "Fluoxetine")))
line_fluoxetine <- plot_line_expression_simple(plot_df = plot_df_fluoxetine, title = "Rat Drug Treatment: GSE194289", genes_cat1 = myratmicegenes, genes_cat2 = myratoligenes, ref_group = "Sham", color_palette = fluoxetine_colors)
# GSE222756
group_order <- c("Control", "Maltodextrin", "Probiotics", "Desipramine", "Bupropion", "FT")
plot_df_probiotics <- GSE222756 %>% 
  rownames_to_column("Gene") %>% 
  filter(Gene %in% myratgenes) %>% 
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>% 
  mutate(Group = factor(str_extract(Sample, paste(group_order, collapse = "|")), levels = group_order))
line_probiotics <- plot_line_expression_probiotics(
  plot_df = plot_df_probiotics, 
  title = "Rat Drug Treatment: GSE222756",
  genes_cat1 = myratmicegenes, 
  genes_cat2 = myratoligenes,
  ref_group = "Control",
  color_palette = probiotics_colors,
  label_spacing_factor = 0.08,hide_ns = F)
final_line_plot <- (line_gse98793 | line_gse76826) / (line_fluoxetine | line_probiotics) +
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)))
print(final_line_plot)



