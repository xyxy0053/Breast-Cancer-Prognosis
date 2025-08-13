# 加载必要的 R包
library(ggplot2)
library(ggsci) # 用于更美观的颜色方案（可选）

# 读取输入文件（假设计价为你是的是制表符分隔文本文件）
input_file <- "volcano.txt"
# 使用 read.table，假设有列名
data <- read.table(input_file, header = TRUE)

# 设置显著性阈值
log2fc_threshold <- 1 # 对数倍数变化的阈值
pvalue_threshold <- 0.05 # p 值的阈值

# 添加基因分组信息（用于上色）
data$group <- "Not significant" # 初始化分组为不显著
data$group[data$padj < pvalue_threshold & abs(data$log2FoldChange) > log2fc_threshold] <- "Significant"
data$group[data$log2FoldChange > log2fc_threshold & data$padj < pvalue_threshold] <- "Up-regulated"
data$group[data$log2FoldChange < -log2fc_threshold & data$padj < pvalue_threshold] <- "Down-regulated"

# 绘制火山图
volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(alpha = 0.3, size = 2.5) + # 绘制散点图
  scale_color_manual(values = c("Up-regulated" = "#D9352A", # 设置颜色
                                "Down-regulated" = "#4575B4",
                                "Significant" = "#000000",
                                "Not significant" = "#999999")) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "red") + # 添加 p 值阈值的水平虚线
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "red") + # 添加倍数变化阈值的垂直虚线
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10(Adjusted P-value)",
       color = "Gene Expression") +
  theme_minimal(base_size = 12) + # 设置主题风格
  theme(panel.grid.major = element_line(color = "lightgray", linetype = "dotted"), # 添加网格线
        legend.position = "bottom") # 将图例放置在底部

# 显示火山图
print(volcano_plot)
