library(DESeq2)
library(readxl)
library(tidyverse)
library(cowplot)
library(PCAtools)
library(ggfortify)
library(ggsci)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggVennDiagram)
library(tidyHeatmap)
count_data <- read.table("Vincent_gene_counts.txt",header = TRUE)
rownames(count_data) <- count_data$gene_id
count_data <- count_data[-1]

sample_metadata <- read.delim2('sampleinfo.txt') %>% as.data.frame() 
table(colnames(count_data)==sample_metadata$Sample)

sample_metadata$Group <- factor(sample_metadata$Group,levels = c('Ctrl', 'NaVPA', 'AP', 'NaVPA_AP'))
design_matrix <- model.matrix(~Group + Batch, data=sample_metadata)
colnames(design_matrix) <- gsub("Group", "", colnames(design_matrix))

sample_metadata$Group <- factor(sample_metadata$Group)
sample_metadata$Batch <- factor(sample_metadata$Batch)
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_metadata,
                              design = ~Batch + Group)
##初步过滤掉杂音
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)
# 可视化
vsd <- vst(dds, blind=FALSE)

mat <- assay(vsd)
mm <- model.matrix(~Group, colData(vsd))
mat <-limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
assay(vsd) <- mat
pca_removebatch  <- plotPCA(vsd,intgroup = "Group", ntop = 500)+
  geom_point(size = 1) +
  geom_text_repel(aes(color = Group, label = colnames(vsd)), max.overlaps = 20) +
  scale_color_bmj() +
  theme_cowplot()+
  theme(plot.title = element_text(hjust = 0.5)) +  
  guides(color = guide_legend(title = NULL),
         shape = guide_legend(title = NULL))

ggsave('pca_removebatch.pdf',pca_removebatch,limitsize = T,width = 5, height = 8,units = 'in')


# 差异分析
gene_id <- bitr(
  rownames(dds),
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = org.Mm.eg.db
)
dds <- dds[gene_id$ENSEMBL]
table(rownames(dds) == gene_id$ENSEMBL)
rownames(dds) <- gene_id$SYMBOL
resultsNames(dds)
dds <- DESeq(dds)

res_AP <- results(dds, contrast=c("Group","NaVPA_AP","AP"))
res_ctrl <- results(dds, contrast=c("Group","NaVPA","Ctrl"))

res_AP <- res_AP[order(res_AP$log2FoldChange,decreasing = c(TRUE,FALSE )), ]
res_AP <- res_AP[!is.na(res_AP$pvalue),]

res_ctrl  <- res_ctrl[order(res_ctrl$log2FoldChange,decreasing = c(TRUE,FALSE )), ]
res_ctrl  <- res_ctrl[!is.na(res_ctrl$pvalue),]

sum(res_AP$pvalue < 0.05 & res_AP$log2FoldChange >0, na.rm=TRUE)
sum(res_AP$pvalue < 0.05 & res_AP$log2FoldChange <0, na.rm=TRUE)

sum(res_ctrl$pvalue < 0.05 & res_AP$log2FoldChange >0, na.rm=TRUE)
sum(res_ctrl$pvalue < 0.05 & res_AP$log2FoldChange <0, na.rm=TRUE)



res_AP$label <- ifelse(rownames(res_AP) %in% c("Cbs", "Sardh", "Mthfr", "Mat1a", "Aldh1l2"), 
                       rownames(res_AP), NA)

vol_AP <- ggplot(res_AP, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = case_when(
    pvalue < 0.05 & log2FoldChange > 0 ~ "Upregulated",
    pvalue < 0.05 & log2FoldChange < 0 ~ "Downregulated",
    TRUE ~ "Not Significant"
  )), size = 1) +  # 设置点的大小
  scale_color_manual(values = c("Upregulated" = "#ED0000FF",
                                "Downregulated" = "#00468BFF",
                                "Not Significant" = "#ADB6B6FF")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # 添加垂直虚线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # 添加水平虚线
  labs(title = "Volcano Plot of NaVPA_AP vs AP",
       x = "Log2 Fold Change (NaVPA_AP vs AP)",
       y = "-log10(pvalue)",
       color = "Significance") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_blank()  # 移除legend标题
  )

ggsave('vol_AP.pdf',vol_AP,limitsize = T,width = 7,height = 5,units = 'in')

vol_ctrl <- ggplot(res_ctrl , aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = case_when(
    pvalue < 0.05 & log2FoldChange > 0 ~ "Upregulated",
    pvalue < 0.05 & log2FoldChange < 0 ~ "Downregulated",
    TRUE ~ "Not Significant"
  )), size = 1) +  # 设置点的大小
  scale_color_manual(values = c("Upregulated" = "#ED0000FF",
                                "Downregulated" = "#00468BFF",
                                "Not Significant" = "#ADB6B6FF")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # 添加垂直虚线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # 添加水平虚线
  labs(title = "Volcano Plot of NaVPA vs ctrl",
       x = "Log2 Fold Change (NaVPA vs ctrl)",
       y = "-log10(pvalue)",
       color = "Significance") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_blank()  # 移除legend标题
  )
ggsave('vol_ctrl.pdf',vol_ctrl,limitsize = T,width = 7,height = 5,units = 'in')


res_ctrl_id <-  bitr(
  rownames(res_ctrl),         # 输入基因符号（行名）
  fromType = "SYMBOL",        # 输入 ID 类型
  toType ="ENTREZID",        # 目标 ID 类型
  OrgDb = org.Mm.eg.db        # 小鼠注释数据库
)
res_ctrl <- res_ctrl %>% as.data.frame() 
res_ctrl <- res_ctrl[res_ctrl_id$SYMBOL,]
res_ctrl <- res_ctrl %>%
  mutate(ENTREZID= res_ctrl_id$ENTREZID)

res_AP_id <-  bitr(
  rownames(res_AP),         # 输入基因符号（行名）
  fromType = "SYMBOL",        # 输入 ID 类型
  toType ="ENTREZID",        # 目标 ID 类型
  OrgDb = org.Mm.eg.db        # 小鼠注释数据库
)
res_AP <- res_AP %>% as.data.frame() 
res_AP <- res_AP[res_AP_id$SYMBOL,]
res_AP <- res_AP %>%
  mutate(ENTREZID= res_AP_id$ENTREZID)
write.csv(res_AP,'res_AP.csv')
write.csv(res_ctrl,'res_ctrl.csv')

res_ctrl_sig <- res_ctrl[res_ctrl$pvalue < 0.05, ]
sum(res_ctrl_sig$log2FoldChange>0)

res_AP_sig <- res_AP[res_AP$pvalue < 0.05, ]
sum(res_AP_sig$log2FoldChange>0)
# 获取基因名称
genes_ctrl <- rownames(res_ctrl_sig)
genes_AP <- rownames(res_AP_sig)
genes_ctrl <- unique(genes_ctrl)
genes_AP <- unique(genes_AP)
x<-list(genes_AP,genes_ctrl)
ggVennDiagram(list(NaVPA_AP_vs_AP = genes_AP, NaVPA_vs_Ctrl = genes_ctrl),
              category.names = c("NaVPA_AP_vs_AP","NaVPA_vs_Ctrl"),
              label = "count", 
              label_color = "black",
              label_alpha = 0,
              edge_lty = "solid",  # 改为实线边框
              edge_size = 1) +
  scale_fill_gradient(low="#66C2A5",high = "#FC8D62",name = "gene count") +
  coord_flip()  # 添加横向显示

genes_ctrl_up <- rownames(res_ctrl_sig[res_ctrl_sig$log2FoldChange>0,])
genes_ctrl_down <- rownames(res_ctrl_sig[res_ctrl_sig$log2FoldChange<0,])
genes_AP_up <- rownames(res_AP_sig[res_AP_sig$log2FoldChange>0,])
genes_AP_down <- rownames(res_AP_sig[res_AP_sig$log2FoldChange<0,])


ggVennDiagram(list(NaVPA_AP_vs_AP_up = genes_AP_up, NaVPA_vs_Ctrl_up = genes_ctrl_up),
              category.names = c("NaVPA_AP_vs_AP_up","NaVPA_vs_Ctrl_up"),
              label = "count", 
              label_color = "black",
              label_alpha = 0,
              edge_lty = "solid",  # 改为实线边框
              edge_size = 1) +
  scale_fill_gradient(low="#66C2A5",high = "#FC8D62",name = "gene count") +
  coord_flip()  # 添加横向显示

ggVennDiagram(list(NaVPA_AP_vs_AP_down = genes_AP_down, NaVPA_vs_Ctrl_down = genes_ctrl_down),
              category.names = c("NaVPA_AP_vs_AP_down","NaVPA_vs_Ctrl_down"),
              label = "count", 
              label_color = "black",
              label_alpha = 0,
              edge_lty = "solid",  # 改为实线边框
              edge_size = 1) +
  scale_fill_gradient(low="#66C2A5",high = "#FC8D62",name = "gene count") +
  coord_flip()  # 添加横向显示



co_DEG <- intersect(genes_AP,genes_ctrl)
co_DEG <- unique(c(genes_AP,genes_ctrl)) 


#### 热图绘制
mat
mat_id <- bitr(
  rownames(mat),         # 输入基因符号（行名）
  fromType = "ENSEMBL",        # 输入 ID 类型
  toType = "SYMBOL",        # 目标 ID 类型
  OrgDb = org.Mm.eg.db        # 小鼠注释数据库
)
mat <- mat[mat_id$ENSEMBL,]
rownames(mat) <- mat_id$SYMBOL

mat
heatmatrix_long_ALL <- mat %>%as.data.frame() %>% 
  rownames_to_column(var = "gene") %>%  
  pivot_longer(
    cols = -gene,                             
    names_to = "sampleid",                    
    values_to = "expression"             
  )

heatmatrix_long_ALL_plot <- heatmatrix_long_ALL %>%as.data.frame() %>% 
  mutate(Group = case_when(
    str_detect(sampleid, "^AP\\d+$") ~ "AP",
    str_detect(sampleid, "^NaVPA_AP\\d+$") ~ "NaVPA_AP",
    str_detect(sampleid, "^NaVPA\\d+$") ~ "NaVPA",
    str_detect(sampleid, "^Ctrl\\d+$") ~ "Ctrl",
    TRUE ~ "Other"  # 处理不符合上述模式的样本
  )) 

heatmatrix_long_ALL_plot$Group <- factor(heatmatrix_long_ALL_plot$Group, levels = c("Ctrl","NaVPA","AP","NaVPA_AP"))

heatmatrix_long_ALL_plot$sampleid <- factor(heatmatrix_long_ALL_plot$sampleid,
                                        levels = c("Ctrl1","Ctrl2","Ctrl3","AP1","AP2","AP3",
                                                   "NaVPA1","NaVPA2","NaVPA3",
                                                   "NaVPA_AP1","NaVPA_AP2","NaVPA_AP3"))


heatmatrix_long_ALL_plot|> 
  group_by(Group) |>
  heatmap(gene,sampleid, expression,na_col = "#F5F4F0",
          scale='both',
          palette_value = circlize::colorRamp2(c(-2,0,2), c("#00468BFF", "#F5F4F0","#ED0000FF")),
          cluster_col = FALSE,cluster_rows = T,
          palette_grouping = list(c('#ECECEC','#244B84', '#DC716E', '#850E00'))
  ) |> 
  save_pdf("rnaseq_ALL_heatmap.pdf",width = 10,height = 25,units = c("in"))


target_gene <- c('Dhfr','Mthfr','Ggh','Shmt1','Shmt2','Tyms','Mthfd1','Mthfd2',#Folate cycle
                 'Mthfd1l','Mthfd2l','Gart','Atic',#Folate cycle
                 'Mtr','Mtrr','Mat1a','Mat2a','Mat2b','Gnmt','Sardh','Chdh',#Methionine cycle
                 'Ahcy','Amd1','Odc1','Sms','Srm','Mtap',#Methionine cycle
                 'Cbs','Cth','Gss','Gclc','Gclm','Gsr','Gpx1','Gpx3','Gpx4'#Transsulfurantion and GSH metabolism
                 )
target_gene <- target_gene[target_gene %in% rownames(mat)]
heatmatrix_long <- mat[target_gene,]
rownames(heatmatrix_long) <- target_gene
heatmatrix_long <- heatmatrix_long %>%as.data.frame() %>% 
  rownames_to_column(var = "gene") %>%  
  pivot_longer(
    cols = -gene,                             
    names_to = "sampleid",                    
    values_to = "expression"             
  )


heatmatrix_long_plot <- heatmatrix_long %>%as.data.frame() %>% 
  mutate(Group = case_when(
    str_detect(sampleid, "^AP\\d+$") ~ "AP",
    str_detect(sampleid, "^NaVPA_AP\\d+$") ~ "NaVPA_AP",
    str_detect(sampleid, "^NaVPA\\d+$") ~ "NaVPA",
    str_detect(sampleid, "^Ctrl\\d+$") ~ "Ctrl",
    TRUE ~ "Other"  # 处理不符合上述模式的样本
  )) %>% 
  mutate(Gene_Label = case_when(
      gene %in% c('Dhfr','Mthfr','Ggh','Shmt1','Shmt2','Tyms','Mthfd1','Mthfd2',
                  'Mthfd1l','Mthfd2l','Gart','Atic') ~ "Folate cycle",
      gene %in% c('Mtr','Mtrr','Mat1a','Mat2a','Mat2b','Gnmt','Sardh','Chdh',
                  'Ahcy','Amd1','Odc1','Sms','Srm','Mtap') ~ "Methionine cycle",
      gene %in% c('Cbs','Cth','Gss','Gclc','Gclm','Gsr','Gpx1','Gpx3','Gpx4') ~ 'Transsulfurantion and GSH metabolism',
    ))

heatmatrix_long_plot$Group <- factor(heatmatrix_long_plot$Group, levels = c("Ctrl","NaVPA","AP","NaVPA_AP"))

heatmatrix_long_plot$Gene_Label <- factor(heatmatrix_long_plot$Gene_Label,
                                                levels = c("Folate cycle","Methionine cycle","Transsulfurantion and GSH metabolism"))

heatmatrix_long_plot$sampleid <- factor(heatmatrix_long_plot$sampleid,
                                        levels = c("Ctrl1","Ctrl2","Ctrl3","AP1","AP2","AP3",
                                                   "NaVPA1","NaVPA2","NaVPA3",
                                                   "NaVPA_AP1","NaVPA_AP2","NaVPA_AP3"))


heatmatrix_long_plot|> 
  group_by(Gene_Label,Group) |>
  heatmap(gene,sampleid, expression,na_col = "#F5F4F0",
          scale='row',
          palette_value = circlize::colorRamp2(c(-2,0,2), c("#00468BFF", "#F5F4F0","#ED0000FF")),
          cluster_col = FALSE,cluster_rows = F,
          palette_grouping = list(
            # For first grouping (vs)
            c("#66C2A5", "#b58b4c","#FC8D62"), 
            # For second grouping (property_group)
            c('#ECECEC','#244B84', '#DC716E', '#850E00')
          )
  ) |> 
  save_pdf("rnaseq_heatmap.pdf",width = 6,height = 10,units = c("in"))




