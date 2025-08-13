require(survival)
require(glmnet)
require(survminer)
require(ggplot2)
require(dplyr)
require(litedown)
require(markdown)


# 加载必要的R包
library(survival)
library(glmnet)
library(survminer)
library(ggplot2)
library(dplyr)
library(litedown)
library(markdown)

# 1. 数据准备
# 读取表达矩阵数据
expression_data <- read.delim("immune_DEGs.fpkm.txt", header = TRUE, row.names = 1, check.names = FALSE)

# 提取肿瘤样本（1102个肿瘤样本）
tumor_samples <- grep("^Tumor_", colnames(expression_data), value = TRUE)
tumor_expression <- t(expression_data[, tumor_samples])

# 读取临床数据
clinical_data <- read.delim("clinical_data.deal.txt", header = TRUE)

# 确保样本ID匹配
# 将临床数据中的ID格式化为与表达矩阵相同的格式
clinical_data$ID <- gsub("-", "_", clinical_data$ID)  # 移除可能的连字符
rownames(clinical_data) <- clinical_data$ID

# 只保留有临床数据的肿瘤样本
common_samples <- intersect(rownames(tumor_expression), clinical_data$ID)
tumor_expression <- tumor_expression[common_samples, ]
clinical_data <- clinical_data[common_samples, ]

# 准备生存数据
surv_data <- data.frame(
  time = clinical_data$OS_time,
  status = clinical_data$OS_status
)
rownames(surv_data) <- common_samples

# 2. LASSO Cox回归分析
# 准备数据
x <- as.matrix(tumor_expression)
y <- Surv(surv_data$time, surv_data$status)

# 执行交叉验证的LASSO Cox回归
set.seed(123)  # 确保结果可重复
cv_fit <- cv.glmnet(
  x = x,
  y = y,
  family = "cox",
  type.measure = "C",
  alpha = 1,  # LASSO回归
  nfolds = 10  # 10折交叉验证
)

# 3. 结果可视化
# 绘制交叉验证曲线
plot(cv_fit, main = "Cross-Validation for Lambda Selection")

# 绘制系数路径图
plot(cv_fit$glmnet.fit, xvar = "lambda", 
     main = "Coefficient Paths in LASSO Cox Regression")

# 4. 提取重要基因
# 选择最优lambda值 (lambda.min或lambda.1se)
#optimal_lambda <- cv_fit$lambda.min  # 更宽松的选择
 optimal_lambda <- cv_fit$lambda.1se  # 更严格的选择

# 提取系数不为零的基因
coefs <- coef(cv_fit, s = optimal_lambda)
selected_genes <- rownames(coefs)[which(coefs != 0)]
gene_coefs <- coefs[which(coefs != 0)]

selected_genes
optimal_lambda

# 创建基因特征
gene_signature <- data.frame(
  Gene = selected_genes,
  Coefficient = as.numeric(gene_coefs))
rownames(gene_signature) <- NULL

# 5. 构建风险评分模型
# 计算每个样本的风险评分
risk_scores <- x[, selected_genes] %*% gene_coefs
colnames(risk_scores) <- "RiskScore"

# 将风险评分与生存数据合并
final_data <- cbind(surv_data, RiskScore = as.numeric(risk_scores))

# 6. 评估基因特征
# 按中位数分组
risk_group <- ifelse(final_data$RiskScore > median(final_data$RiskScore), 
                   "High Risk", "Low Risk")

# 绘制Kaplan-Meier曲线
fit <- survfit(Surv(time, status) ~ risk_group, data = final_data)
ggsurvplot(fit, 
           data = final_data,
           pval = TRUE,
           risk.table = TRUE,
           title = "Survival Analysis by Risk Group",
           xlab = "Time (Months)")

# 7. 保存结果
# 保存基因特征
write.csv(gene_signature, "immune_gene_signature.csv", row.names = FALSE)

# 保存风险评分
write.csv(final_data, "patient_risk_scores.csv", row.names = TRUE)

# 8. 输出重要信息
cat("Optimal lambda:", optimal_lambda, "\n")
cat("Number of selected genes:", length(selected_genes), "\n")
cat("Selected genes:\n")
print(gene_signature)



# 提取系数不为零的基因及其系数
coefs <- coef(cv_fit, s = optimal_lambda)
selected_genes <- coefs@Dimnames[[1]][coefs@i + 1]  # 获取基因名称
gene_coefs <- coefs@x  # 获取系数值

# 创建基因特征数据框
gene_signature <- data.frame(
  Gene = selected_genes,
  Coefficient = gene_coefs
)

# 打印基因特征
print(gene_signature)



# 生成公式字符串
formula_parts <- paste0(
  sprintf("%.2f", gene_signature$Coefficient),
  " × ",
  gene_signature$Gene
)

# 组合成完整公式
risk_formula <- paste("RiskScore =", paste(formula_parts, collapse = " + "))

# 打印风险评分公式
cat(risk_formula, "\n")
RiskScore = 0.10 × ENSG00000205186 + -0.00 × ENSG00000138755 + -0.00 × ENSG00000173432 + 0.01 × ENSG00000143171 + -0.00 × ENSG00000170345 + -0.33 × ENSG00000145777 + 0.02 × ENSG00000113389 + -0.02 × ENSG00000136011 + 0.01 × ENSG00000139874 + 0.05 × ENSG00000106178 + -0.00 × ENSG00000170893 


# 确保表达矩阵只包含signature genes
signature_expression <- tumor_expression[, gene_signature$Gene]

# 计算风险评分
risk_scores <- apply(signature_expression, 1, function(x) {
  sum(x * gene_signature$Coefficient)
})

# 添加到临床数据
clinical_data$RiskScore <- risk_scores[rownames(clinical_data)]



library(ggplot2)
library(tidyr)

# 选择影响最大的5个基因（按绝对值）
top_genes <- gene_signature[order(abs(gene_signature$Coefficient), decreasing = TRUE), ]
top_genes <- head(top_genes, 5)

# 准备数据
plot_data <- tumor_expression[, top_genes$Gene] %>% 
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Gene", values_to = "Expression")

# 添加系数信息
plot_data <- merge(plot_data, top_genes, by = "Gene")

# 绘制小提琴图
ggplot(plot_data, aes(x = Gene, y = log2(Expression + 1), fill = Gene)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, fill = "white") +
  facet_wrap(~ Gene, scales = "free_x", nrow = 1) +
  labs(title = "Expression Distribution of Key Signature Genes",
       subtitle = paste("Positive coefficients increase risk, negative decrease risk"),
       y = "log2(FPKM + 1)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
		
		
		
		
forest_plot <- ggplot(forest_data, aes(x = HR, y = reorder(gene, HR))) +
  geom_point(aes(size = -log10(`Pr(>|z|)`)), 
             color = ifelse(forest_data$HR > 1, "red", "blue")) +
  geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), 
                 height = 0.2, color = "grey50") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_x_log10() +
  labs(
    title = "Multivariate Cox Regression Results",
    x = "Hazard Ratio (log scale)",
    y = "Gene Symbol"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )


		