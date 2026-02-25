options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 32000 * 1024^2)
set.seed(12345)
library(limma)
library(ggplot2)
library(patchwork)
library(ggvenn)
library(clusterProfiler)
library(enrichplot)
dataset_info <- list(
  list(name = "GSE98793", data = GSE98793, control_pattern = "^control", case_pattern = "^MDD"),
  list(name = "GSE76826", data = GSE76826, control_pattern = "^control", case_pattern = "^MDD"),
  list(name = "GSE290797", data = GSE290797, control_pattern = "^control", case_pattern = "^MDD"),
  list(name = "GSE247998", data = GSE247998, control_pattern = "^control", case_pattern = "^MDD"),
  list(name = "GSE99725", data = GSE99725, control_pattern = "controlbefore", case_pattern = "MDDbefore"),
  list(name = "GSE185855", data = GSE185855, control_pattern = "^control", case_pattern = "MDDPRE"))

#DEG analysis
#create limma analysis function
run_unified_limma_analysis <- function(dataset_name, expression_data, control_pattern, case_pattern, 
                                       p_val_cutoff = 0.05, log_fc_cutoff = log2(1.1)) {# DEG threshold p<0.05,log2|FC|>log2(1.1)
  #data preprocessing
  all_samples <- colnames(expression_data)
  relevant_samples_mask <- grepl(paste(control_pattern, case_pattern, sep = "|"), all_samples)
  subset_data <- expression_data[, relevant_samples_mask]
  #group
  groups <- factor(
    ifelse(grepl(control_pattern, colnames(subset_data)), "control", "case"),
    levels = c("control", "case"))
  #limma analysis
  design <- model.matrix(~0 + groups)
  colnames(design) <- levels(groups)
  fit <- lmFit(subset_data, design)
  contrast.matrix <- makeContrasts(case_vs_control = case - control, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  #differential genes result
  full_results_table <- topTable(fit2, coef = "case_vs_control", number = Inf)
  #exceeding the threshold DEGs
  deg_significant <- subset(full_results_table, P.Value < p_val_cutoff & abs(logFC) > log_fc_cutoff)
  #organize into a list
  results <- list(
    dataset_name = dataset_name,
    full_table = full_results_table,
    significant_degs = deg_significant,
    gene_lists = list(
      all = rownames(deg_significant),
      up = rownames(subset(deg_significant, logFC > 0)),
      down = rownames(subset(deg_significant, logFC < 0))),
    deg_counts = c(
      total = nrow(deg_significant),
      up = sum(deg_significant$logFC > 0),
      down = sum(deg_significant$logFC < 0)
    ),subset_data = subset_data)
  cat(paste(dataset_name, "found", results$deg_counts["total"], "DEGs。up-regulation:", results$deg_counts["up"], "down-regulation:", results$deg_counts["down"], "\n"))
  return(results)}
#execute function
all_limma_results <- lapply(dataset_info, function(info) {
  run_unified_limma_analysis(
    dataset_name = info$name,
    expression_data = info$data,
    control_pattern = info$control_pattern,
    case_pattern = info$case_pattern)})
#name DEG list
names(all_limma_results) <- sapply(dataset_info, `[[`, "name")


#figure 1B bar plot
#drawing data preparation
plot_data_list <- lapply(names(all_limma_results), function(name) {
  counts <- all_limma_results[[name]]$deg_counts
  data.frame(
    Dataset = name,
    Regulation = c("Up-regulated", "Down-regulated"),
    Count = c(counts["up"], -counts["down"]))})
plot_data <- do.call(rbind, plot_data_list)
#bar plot
deg_count_plot <- ggplot(plot_data, aes(x = Dataset, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_y_continuous(labels = abs, name = "Number of DEGs") +
  geom_hline(yintercept = 0, colour = "black") +
  labs(
    title = "Number of Differentially Expressed Genes (DEGs)",
    subtitle = "Threshold: p < 0.05, |log2FC| > log2(1.1)",
    fill = "Regulation Status") +
  scale_fill_manual(values = c("Up-regulated" = "#e50f4b", "Down-regulated" = "#497cc2")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5))
print(deg_count_plot)


#gene intersection analysis
all_deg_list <- lapply(all_limma_results, function(res) res$gene_lists$all)
up_deg_list <- lapply(all_limma_results, function(res) res$gene_lists$up)
down_deg_list <- lapply(all_limma_results, function(res) res$gene_lists$down)
#extract upregulated and downregulated gene lists from GSE98793 and GSE76826
up_degs_gse98793 <- all_limma_results$GSE98793$gene_lists$up
up_degs_gse76826 <- all_limma_results$GSE76826$gene_lists$up
all_degs_gse98793 <- all_limma_results$GSE98793$gene_lists$all
down_degs_gse98793 <- all_limma_results$GSE98793$gene_lists$down
down_degs_gse76826 <- all_limma_results$GSE76826$gene_lists$down
all_degs_gse76826 <- all_limma_results$GSE76826$gene_lists$all
#identify the intersection of up-regulated genes
intersect_up_degs <- intersect(up_degs_gse98793, up_degs_gse76826)
cat(paste("both up-regulated genes number:", length(intersect_up_degs), "\n"))
print(intersect_up_degs)
#identify the intersection of down-regulated genes
intersect_down_degs <- intersect(down_degs_gse98793, down_degs_gse76826)
cat(paste("both down-regulated genes number:", length(intersect_down_degs), "\n"))
print(intersect_down_degs)
#find the common DEG intersection of these two sets of data
intersect_all_degs <- intersect(all_degs_gse98793, all_degs_gse76826)
cat(paste("total DEG number:", length(intersect_all_degs), "\n"))
print(intersect_all_degs)


#figure 1C upset plot
library(UpSetR)
upset(fromList(all_deg_list), nsets = 6, order.by = "freq", text.scale = 2)
#figure 1C venn plot
#prepare data
venn_up_list <- list(GSE98793 = up_degs_gse98793,GSE76826 = up_degs_gse76826)
venn_down_list <- list(GSE98793 = down_degs_gse98793,GSE76826 = down_degs_gse76826)
#up-regulated gene venn plot
venn_up_plot <- ggvenn(
  venn_up_list,
  fill_color = c("#e50f4b", "#f287a5"),
  stroke_size = 0.5, set_name_size = 5
) + labs(title = "Common Up-regulated Genes") + theme(plot.title = element_text(hjust = 0.5))
#up-regulated gene venn plot
venn_down_plot <- ggvenn(
  venn_down_list,
  fill_color = c("#497cc2", "#a4bee1"),
  stroke_size = 0.5, set_name_size = 5
) + labs(title = "Common Down-regulated Genes") + theme(plot.title = element_text(hjust = 0.5))
# final venn plot
final_plot <- volcano_plot / (venn_up_plot + venn_down_plot) +
  plot_annotation(
    title = 'Combined DEG Analysis of MDD Datasets GSE98793 & GSE76826',
    caption = 'Figure generated using ggplot2, ggvenn, and patchwork') & 
  theme(plot.title = element_text(hjust = 0.5))
print(final_plot)


#RRA analysis
library(RobustRankAggreg)
#prepare data(need "all_limma_results")
#create a list of up-regulated and down-regulated DEG results from "all_limma_results"
rra_input_up_lists <- list()
rra_input_down_lists <- list()
for (dataset_name in names(all_limma_results)) {
  full_table <- all_limma_results[[dataset_name]]$full_table
  #select and sort the list of up-regulated genes
  up_genes <- subset(full_table, logFC > 0)
  up_genes_sorted <- up_genes[order(up_genes$P.Value), ]
  rra_input_up_lists[[dataset_name]] <- rownames(up_genes_sorted)
  #select and sort the list of up-regulated genes
  down_genes <- subset(full_table, logFC < 0)
  down_genes_sorted <- down_genes[order(down_genes$P.Value), ]
  rra_input_down_lists[[dataset_name]] <- rownames(down_genes_sorted)}
#cleaning up-regulated genes
cleaned_up_lists <- list()
for (dataset_name in names(rra_input_up_lists)) {
  original_list <- rra_input_up_lists[[dataset_name]]
  cat(paste("Original gene number:", length(original_list), "\n"))
  #check and remove NA and empty strings
  valid_list <- original_list[!is.na(original_list) & original_list != ""]
  #check and remove duplicates
  final_list <- unique(valid_list)
  cat(paste("clean gene number:", length(final_list), "\n"))
  cleaned_up_lists[[dataset_name]] <- final_list}
#cleaning down-regulated genes
cleaned_down_lists <- list()
for (dataset_name in names(rra_input_down_lists)) {
  original_list <- rra_input_down_lists[[dataset_name]]
  cat(paste("Original gene number:", length(original_list), "\n"))
  #check and remove NA and empty strings
  valid_list <- original_list[!is.na(original_list) & original_list != ""]
  #check and remove duplicates
  final_list <- unique(valid_list)
  cat(paste("clean gene number:", length(final_list), "\n"))
  cleaned_down_lists[[dataset_name]] <- final_list}
#execute RRA analysis
#up-regulated gene RRA
rra_up_results <- aggregateRanks(glist = cleaned_up_lists)
#down-regulated gene RRA
rra_down_results <- aggregateRanks(glist = cleaned_down_lists)
cat("\n\nRRA analysis complete\n")
#selected the final list of genes based on the integrated p-values (Score)
rra_final_up_genes <- subset(rra_up_results, Score < 0.05)$Name
rra_final_down_genes <- subset(rra_down_results, Score < 0.05)$Name
cat(paste("\np < 0.05，found", length(rra_final_up_genes), "up-regulated genes\n"))
cat(paste("p < 0.05，found", length(rra_final_down_genes), "down-regulated genes\n"))
#create a function to calculate the average logFC of the specified gene.
calculate_average_logfc <- function(gene_name, results_list) {
  logfc_values <- c()
  for (dataset_name in names(results_list)) {
    full_table <- results_list[[dataset_name]]$full_table
    if (gene_name %in% rownames(full_table)) {
      logfc_values <- c(logfc_values, full_table[gene_name, "logFC"])}}
  return(mean(logfc_values, na.rm = TRUE))}
#calculate the average logFC for all the genes with significant RRA results
#up-regulated RRA gene
avg_logfc_up <- sapply(rra_final_up_genes, calculate_average_logfc, results_list = all_limma_results)
#down-regulated RRA gene
avg_logfc_down <- sapply(rra_final_down_genes, calculate_average_logfc, results_list = all_limma_results)
#create a summary table of RRA results for up-regulated genes
final_up_summary <- data.frame(
  Name = rra_final_up_genes,
  Avg_logFC = avg_logfc_up)
final_up_summary <- merge(final_up_summary, rra_up_results, by = "Name")
final_up_summary <- final_up_summary[order(final_up_summary$Score), c("Name", "Score", "Avg_logFC")]
#create a summary table of RRA results for up-regulated genes
final_down_summary <- data.frame(
  Name = rra_final_down_genes,
  Avg_logFC = avg_logfc_down)
final_down_summary <- merge(final_down_summary, rra_down_results, by = "Name")
final_down_summary <- final_down_summary[order(final_down_summary$Score), c("Name", "Score", "Avg_logFC")]
print(head(final_up_summary))
print(head(final_down_summary))
#screen the genes corresponding to the RRA results (Score < 0.05,Avg_logFC > log2(1.1))
logfc_threshold <- log2(1.1)
#screening up-regulated RRA genes
core_up_genes_df <- subset(final_up_summary, Score < 0.05 & Avg_logFC > logfc_threshold)
#screening down-regulated RRA genes
core_down_genes_df <- subset(final_down_summary, Score < 0.05 & Avg_logFC < -logfc_threshold)


#figure 1D volcano plot
library(ggplot2)
library(dplyr)
library(ggrepel)
core_up_genes_df$Status <- "Upregulated"
core_down_genes_df$Status <- "Downregulated"
all_core_genes_df <- rbind(core_up_genes_df, core_down_genes_df)
all_core_genes_df$log10Score <- -log10(all_core_genes_df$Score)
all_core_genes_df$label <- ifelse(all_core_genes_df$log10Score > 5 & abs(all_core_genes_df$Avg_logFC) > 50, all_core_genes_df$Name, "")
plot_data <- all_core_genes_df %>%
  mutate(neg_log10_Score = -log10(Score),Status = case_when(
    Score < 0.01 & Avg_logFC > log2(1.1)   ~ "Upregulated",
    Score < 0.01 & Avg_logFC < -log2(1.1)  ~ "Downregulated",
    TRUE                           ~ "Not Significant"))
genes_to_label <- plot_data %>%
  filter(Status != "Not Significant") %>%
  arrange(Score) %>%
  group_by(Status) %>%
  slice_head(n = 5) %>%
  ungroup()
up_count <- sum(plot_data$Status == "Upregulated")
down_count <- sum(plot_data$Status == "Downregulated")
volcano_plot_optimized <- ggplot(plot_data, aes(x = Avg_logFC, y = neg_log10_Score)) +
  geom_point(aes(color = Status), alpha = 0.5, size = 2) +
  scale_color_manual(
    values = c("Upregulated" = "#e50f4b", "Downregulated" = "#497cc2", "Not Significant" = "grey"),
    name = "Gene Status") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_text_repel(data = genes_to_label,aes(label = Name),
                  box.padding = 0.5,point.padding = 0.2,
                  segment.color = 'grey50',max.overlaps = Inf) +# 尽量显示所有标签# 设置主题和标签
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(linetype = "dashed"),
        panel.grid.minor = element_blank()) +
  labs( title = "Meta-Analysis of Differentiatedly Expressed Genes",
        subtitle = paste(length(all_limma_results), "datasets integrated via RRA"),
        x = "Average log2(Fold Change)", y = "-log10(RRA Score)") +
  coord_cartesian(xlim = c(-5000, 5000))
print(volcano_plot_optimized)



