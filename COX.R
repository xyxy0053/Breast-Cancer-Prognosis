setwd("D:\\test")

AssayData <- read.table("immune_DEGs.fpkm.filter.txt",sep="\t",header=T,row.names=1)
AssayData <- AssayData[order(colnames(AssayData))]
clinical <- read.table("clinical_data.all.filter.txt",head=T,sep="\t")
AssayData <- AssayData[,colnames(AssayData)%in%t(clinical[,1])]
phenoData <- clinical[clinical[,1]%in%colnames(AssayData),]
dim(AssayData)
dim(phenoData)

##筛选表达在一半以上样本低表达的基因
keep.exprs <- rowSums(AssayData>0)>=0.5*ncol(AssayData)
AssayData <- AssayData[keep.exprs,]
dim(AssayData)

t_exp <- AssayData
pheno <- phenoData

library(survival)
## Eliminate more than half of the gene data not expressed
t_exp <- t_exp[apply(t_exp, 1, function(x) sum(x == 0))<(ncol(t_exp)*0.49),]
dim(t_exp)
group_data <- apply(t_exp , 1 , function(gene){
  name <- rownames(gene)
  gene <- unlist(gene)
  group <- ifelse(gene >= median(gene), '2', '1')
  names(group) <- name
  return(group)
})
##single
group_data <- as.data.frame(group_data, stringsAsFactors = F)
survival_dat <- data.frame(
                           age = pheno$age,
                           gender = pheno$gender,
                           stage = pheno$stage,
                           status = pheno$status,
                           T = pheno$T,
                           N = pheno$N,
                           M = pheno$M,
                           time = pheno$time,
                           stringsAsFactors = F)
survival_dat <- cbind(group_data, survival_dat)
colnames(survival_dat) <- sub("\\-", "", colnames(survival_dat))
covariates <- as.character(colnames(survival_dat))
univ_formulas <- sapply(covariates,
                        function(x){
                          ##print(x)
                          as.formula(paste('Surv(time, status)~', x))
                        })
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = survival_dat)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value <- signif(x$wald["pvalue"], digits = 2)
                         beta <- signif(x$coef[1], digits = 2)
                         HR <- signif(x$coef[2], digits = 2)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res <- c(beta, HR, p.value)
                         names(res) <- c("coef", "HR (95% CI for HR)", "p.value")
                         return(res)
                       })
res_single <- as.data.frame(t(do.call(cbind, univ_results)))
write.table(res_single,"COX-res-single.txt",sep="\t",row.names=T,col.names=T)

##保存单因素结果
table(as.matrix(res_single$p.value) <= 0.05)
res_single <- res_single[as.matrix(res_single$p.value) <= 0.05, ]
res_single <- res_single[order(res_single$p.value), ]
single_pick <- rownames(res_single)

## multi
multi_results <- apply(t_exp , 1 , function(gene){
  ## gene <- t_exp[1, ]
  gene <- unlist(gene)
  group <- ifelse(gene >= median(gene), '2', '1')
  ##这里选取在单因素里p < 0.05的进行多因素回归（例如age,stage等）
  survival_dat <- data.frame(group = group, 
                             age = pheno$age,
                             stage = pheno$stage,
                             T = pheno$T,
                             N = pheno$N,
                             M = pheno$M,
                             status = pheno$status,
                             time = pheno$time,
                             stringsAsFactors = F)
  res.cox <- coxph(Surv(time, status) ~ age + stage + T + N + M + group, 
                   data =  survival_dat)
  ## summary(res.cox)
  beta <- coef(res.cox)
  se <- sqrt(diag(vcov(res.cox)))
  HR <- exp(beta)
  HRse <- HR * se
  #summary(m)
  res <- as.data.frame(round(cbind(coef = beta,
                     se = se,
                     z = beta/se,
                     p.value = 1 - pchisq((beta/se)^2, 1),
                     HR = HR,
                     HRse = HRse,
                     HRz = (HR - 1) / HRse,
                     HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 5))
  return(res['group2',])
})
multi_results <- do.call(rbind, multi_results)
table(multi_results$p.value <= 0.05)
res_multi <- multi_results[multi_results$p.value <= 0.05, ]
res_multi <- res_multi[order(res_multi$p.value), ]
write.table(res_multi,"COX-res-multiple.txt",sep="\t",row.names=T,col.names=T)

multi_pick <- rownames(res_multi)
overgene <- intersect(multi_pick, single_pick)

model_exp <- group_data[,overgene] 
colnames(model_exp) <- overgene
dat <- cbind(pheno, model_exp)

library("survminer")
colnames(dat)
res_multi1 <- res_multi[rownames(res_multi)%in%overgene,]
##res_multi1 <- res_multi1[res_multi1$HR > 1, ]
res_multi1 <- res_multi1[order(res_multi1$HR),]

dim(res_multi1)
pick_genes <- rownames(res_multi1)
write.table(res_multi1,"COX-res-multi.txt",sep="\t",row.names=T,col.names=T)

sample1 <- read.table("COX-res-multi.txt",header=T,sep="\t")
sample1<-data.frame(sample1,stringsAsFactors=FALSE)

library(ggplot2)
pdf("多因素森林图.pdf",width=7,height=12)
p <- ggplot(sample1,aes(x=HR,y=ID,color=p.value))+
  geom_errorbarh(aes(xmax=HRCILL,xmin=HRCIUL),color="black",height=0,size=0.8)+#加线
  geom_point(aes(x=HR,y=ID),size=4,shape=18)+###点
  geom_vline(xintercept=1,linetype="dashed",size=1.2)+##垂直线
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6))+
  ##scale_y_discrete(labels=c(as.character(sample[,1])))+
  ##coord_trans(x="log2")+
  scale_color_continuous(low="#EDE972",high="#217021")+
  xlab("Hazard ratios")+
  labs(color="P value",title="")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+  
  theme(axis.text=element_text(colour="black",size=12))+  
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=12))+
  theme(axis.title=element_text(colour="black",size=12))
p
dev.off()

library(survminer)
model_exp <- t_exp[overgene, ]
for (i in 1:nrow(model_exp)) {
  gene <- model_exp[i, ]
  name <- rownames(gene)
  gene <- unlist(gene)
  group <- ifelse(gene >= median(gene), 'high', 'low')
  survival_dat <- data.frame(group = group,
                             status = pheno$status,
                             time = pheno$time,
                             stringsAsFactors = F)
  fit <- survfit(Surv(time, status) ~ group, data = survival_dat)
  ggsurvplot(fit, data = survival_dat,
             surv.median.line = "hv",
             legend.title = "Group",
             legend.labs = c("High", "Low"),
             pval = TRUE,
             conf.int = TRUE,
             palette = "jco",
             ggtheme = theme_bw()
  )
  ggsave(filename = paste(name, '.pdf', sep = ''))
}

























