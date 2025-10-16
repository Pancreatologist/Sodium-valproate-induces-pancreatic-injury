library(data.table)
library(ggplot2)
library(cowplot)
library(plotROC)
library(tidyverse)
library(limma)
library(ggrepel)
#加载R包
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(dplyr) # A Grammar of Data Manipulation
library(RColorBrewer) # ColorBrewer Palettes
library(grid) # The Grid Graphics Package
library(scales) # Scale Functions for Visualization

df_pancreas <- readxl::read_xlsx('Vicent代谢组学.xlsx',sheet = 1) %>% t() %>% as.data.frame() 
colnames(df_pancreas) <- df_pancreas[1,]
df_pancreas <- df_pancreas[-1,]
df_pancreas <- df_pancreas %>%
  mutate(across(everything(), as.numeric, na.rm = TRUE)) %>%
  mutate(across(everything(), log))
sampleinfo_pancreas <- readxl::read_xlsx('Vicent代谢组学.xlsx',sheet = 2)
metadata_pancreas <- readxl::read_xlsx('Vicent代谢组学.xlsx',sheet = 3)
sampleinfo_pancreas$group <- factor(sampleinfo_pancreas$group, 
                                    levels = c("ctrl", "Na.VPA", "CER.AP", "Na.VPA.CER.AP"))

##PCR
pca_data <- t(df_pancreas)
pca_data  <- prcomp(pca_data )
summary(pca_data)
# 使用ggplot2包绘制PCA图
rownames(pca_data$x) == sampleinfo_pancreas$sample
pca_plot <- ggplot(pca_data$x, aes(PC1, PC2)) +
  geom_point(aes(fill = factor(sampleinfo_pancreas$group)), 
             color = "black",  # 设置边框为黑色
             size = 4, 
             shape = 21) +     # 形状21是有边框的点
  geom_text(aes(label = sampleinfo_pancreas$sample), size = 3.5) +
  theme_cowplot(10) +
  theme(legend.position = "none") +
  stat_ellipse(aes(color = factor(sampleinfo_pancreas$group)), 
               type = "norm", level = 0.80, segments = 200) +
  scale_fill_manual(values = c('#ECECEC', '#244B84', '#DC716E', '#850E00')) + # 注意改为scale_fill_manual
  scale_color_manual(values = c('#ECECEC', '#244B84', '#DC716E', '#850E00')) + # 保留原来的颜色标度
  labs(x = "PC1", y = "PC2", title = "PCA analysis")
ggsave('pca_plot代谢组学.pdf',pca_plot,device='pdf',units='in',width = 6,height = 5)



### 差异分析
sampleinfo_pancreas$group<- factor(sampleinfo_pancreas$group,levels = c("ctrl","Na.VPA","CER.AP","Na.VPA.CER.AP"))
mm <- model.matrix(~0 + sampleinfo_pancreas$group)
colnames(mm) <- c('ctrl','Na.VPA','CER.AP','Na.VPA.CER.AP')

fit <- lmFit(df_pancreas, mm)
head(coef(fit))

VPA.CER.AP <- makeContrasts(Na.VPA.CER.AP-CER.AP, levels = colnames(mm))
VPA.CER.AP  <- contrasts.fit(fit, VPA.CER.AP) %>% eBayes()
VPA.CER.AP.table <- topTable(VPA.CER.AP, sort.by = "P", n = Inf)
head(VPA.CER.AP.table, 20)
length(which(VPA.CER.AP.table$P.Val < 0.05))##Na.VPA.CER.AP vs. CER.AP
VPA.CER.AP.table$metabolite <- rownames(VPA.CER.AP.table)
write.table(VPA.CER.AP.table, file = "VPA.CER.AP vs CER.AP.table.txt", row.names = F, sep = "\t", quote = F)

VPA.control <- makeContrasts(Na.VPA-ctrl, levels = colnames(mm))
VPA.control  <- contrasts.fit(fit, VPA.control) %>% eBayes()
VPA.control.table <- topTable(VPA.control, sort.by = "P", n = Inf)
head(VPA.control, 20)
length(which(VPA.control$p.value < 0.05)) ## Na.VPA vs. ctrl
VPA.control.table$metabolite <- rownames(VPA.control.table)
write.table(VPA.control.table, file = "VPA.control vs control.table.txt", row.names = F, sep = "\t", quote = F)

CER.AP.control <- makeContrasts(CER.AP-ctrl, levels = colnames(mm))
CER.AP.control  <- contrasts.fit(fit, CER.AP.control) %>% eBayes()
CER.AP.control.table <- topTable(CER.AP.control, sort.by = "P", n = Inf)
head(CER.AP.control, 20)
length(which(CER.AP.control$p.value < 0.05)) ## CER.AP vs ctrl
CER.AP.control.table$metabolite <- rownames(CER.AP.control.table)
write.table(CER.AP.control.table, file = "CER.AP vs control.table.txt", row.names = F, sep = "\t", quote = F)

### 火山图  ###
VPA.CER.AP.table <- VPA.CER.AP.table %>% mutate(group = 'VPA.CER.AP vs CER.AP')
VPA.control.table <- VPA.control.table  %>% mutate(group = 'VPA.control vs control')
CER.AP.control.table <- CER.AP.control.table  %>% mutate(group = 'CER.AP vs control')
df_volcaono <- rbind(VPA.CER.AP.table,VPA.control.table,CER.AP.control.table)

df_volcaono$group <- factor(df_volcaono$group, 
                            levels = c("VPA.control vs control",
                                       "CER.AP vs control", 
                                       "VPA.CER.AP vs CER.AP"))
##与之前绘制单组火山图一致，先根据设定阈值确定所有OTU的显著性
df_volcaono$group2<-as.factor(ifelse(df_volcaono$P.Value < 0.05, 
                            ifelse(df_volcaono$logFC>= 0 ,'up','down'),'ns'))
df_volcaono$group2 <- factor(df_volcaono$group2, levels = c("up", "down", "ns"))
##确定添加标签的数据，看后续是否要加
#df_volcaono$label<-ifelse(df$p_value<0.05&abs(df$log2FC)>=4,"Y","N")
#df_volcaono$label<-ifelse(df$label == 'Y', as.character(df$OTU), '')

##为了构建图形中所有数据点背景框，需要先确定每个组的最大值与最小值
df_bg <- df_volcaono %>%
  group_by(group) %>%
  summarize(max_logFC = max(logFC),min_logFC = min(logFC))

# 首先统计每组中up和down的基因数量
df_counts <- df_volcaono %>%
  filter(group2 %in% c("up", "down")) %>%  # 只统计up和down的基因
  count(group, group2, name = "count") %>%  # 按group和group2统计数量
  mutate(
    # 为up和down文本设置不同的y位置
    y_pos = ifelse(
      group2 == "up", 
      0.8 * max(df_volcaono$logFC),  # up文本放在接近顶部的位置
      0.8 * min(df_volcaono$logFC)   # down文本放在接近底部的位置
    )
  )

# 查看统计结果
df_counts

# 绘制火山图并添加数量标签
p1 <- ggplot() +
  geom_jitter(data = df_volcaono,
              mapping = aes(x = group, y = logFC, color = group2),
              size= 3, width = 0.4, alpha = 0.7) +
  scale_color_manual(values = c("up" = "#8E1B20",  # 莫兰迪红（带灰调）
                                "down" = "#0A3B59", # 莫兰迪蓝
                                "ns" = "#B8B8B8"), # 莫兰迪灰
                     name = "DE") +  # 设置legend标签为'DE'
  # 添加数量文本标签
  geom_text(
    data = df_counts,
    mapping = aes(x = group, y = y_pos, label = count, color = group2),
    fontface = "bold",
    size = 5,
    vjust = ifelse(df_counts$group2 == "up", -0.5, 1.5)  # 调整文本垂直位置
  ) +
  labs(x = "Contrast", y = "LogFoldChange", fill= NULL, color = NULL)  # 设置x轴标签为'contrast'

p1
ggsave('volcaoo_plot代谢组学.pdf',p1,device='pdf',units='in',width = 7,height = 5)



### 富集分析
table(VPA.control.table$metabolite %in% metadata_pancreas$Compound_ID)
VPA.control.table.enrich <- left_join(VPA.control.table, metadata_pancreas, 
                            by = c("metabolite" = "Compound_ID"))
VPA.control.table.enrich <- VPA.control.table.enrich[!is.na(VPA.control.table.enrich$`KEGG ID`),]
VPA.control.table.enrich.sig <- VPA.control.table.enrich[VPA.control.table.enrich$P.Value<0.05,] 
write.table(VPA.control.table.enrich.sig, file = "VPA.control.table.enrich.sig.txt", row.names = F, sep = "\t", quote = F)

table(VPA.CER.AP.table$metabolite %in% metadata_pancreas$Compound_ID)
VPA.CER.AP.table.enrich <- left_join(VPA.CER.AP.table, metadata_pancreas, 
                                      by = c("metabolite" = "Compound_ID"))
VPA.CER.AP.table.enrich <- VPA.CER.AP.table.enrich[!is.na(VPA.CER.AP.table.enrich$`KEGG ID`),]
VPA.CER.AP.table.enrich.sig <- VPA.CER.AP.table.enrich[VPA.CER.AP.table.enrich$P.Value<0.05,] 
write.table(VPA.CER.AP.table.enrich.sig, file = "VPA.CER.AP.table.enrich.sig.txt", row.names = F, sep = "\t", quote = F)

table(CER.AP.control.table$metabolite %in% metadata_pancreas$Compound_ID)
CER.AP.control.table.enrich <- left_join(CER.AP.control.table, metadata_pancreas, 
                                     by = c("metabolite" = "Compound_ID"))
CER.AP.control.table.enrich <- CER.AP.control.table.enrich[!is.na(CER.AP.control.table.enrich$`KEGG ID`),]
CER.AP.control.table.enrich.sig <- CER.AP.control.table.enrich[CER.AP.control.table.enrich$P.Value<0.05,] 
write.table(CER.AP.control.table.enrich.sig, file = "CER.AP.control.table.enrich.sig.txt", row.names = F, sep = "\t", quote = F)


bg_metabolism <- unique(metadata_pancreas$Name)
write.table(bg_metabolism, file = "bg_metabolism.txt", row.names = F, sep = "\t", quote = F)


co_metabolism <- intersect(VPA.control.table.enrich.sig$metabolite,VPA.CER.AP.table.enrich.sig$metabolite)

co.table.enrich.sig <- VPA.CER.AP.table.enrich.sig[VPA.CER.AP.table.enrich.sig$metabolite %in% co_metabolism,]
write.table(co.table.enrich.sig, file = "co.table.enrich.sig.txt", row.names = F, sep = "\t", quote = F)

###做venn图
all_contrast <- makeContrasts(Na.VPA.CER.AP-CER.AP,
                              Na.VPA-ctrl,
                              CER.AP-ctrl,levels = colnames(mm))
all_contrast.results <- contrasts.fit(fit, all_contrast) %>% eBayes()
results <- decideTests(all_contrast.results, adjust.method = "none")
vennDiagram(results, include=c("up","down"),
            counts.col=c("#8E1B20",  "#0A3B59"),
            circle.col = c("#E3D8A7", "#dbebfa", "#89B0A2"))

result_rows <- results[results[,1] == -1 & results[,2] == -1 & results[,3] == -1, ]
view(result_rows)
# 查看结果
print(result_rows)

kegg_meta <- readxl::read_xlsx('KEGG_all_metabolism_mmu.xlsx',sheet = 1) %>% as.data.frame() 
head(kegg_meta)

kegg_meta_no_comma <- kegg_meta %>%
  filter(!str_detect(CPDs, ",")) %>%
  mutate(CPDs = str_trim(CPDs))
kegg_meta_comma <- kegg_meta %>%
  filter(str_detect(CPDs, ",")) %>%
  separate_rows(CPDs, sep = ",") %>%
  mutate(CPDs = str_trim(CPDs))
kegg_meta_final <- bind_rows(kegg_meta_no_comma, kegg_meta_comma)
kegg_meta_final %>%
  arrange(PIDs, CPDs) %>%
  head(15)


kegg_meta_final$CPDs
bg_kegg <- unique(metadata_pancreas$`KEGG ID`)
bg_kegg <- bg_kegg[!is.na(bg_kegg)]
bg_kegg <- bg_kegg[bg_kegg %in% unique(kegg_meta_final$CPDs)]
filtered_kegg <- kegg_meta_final

filtered_kegg <- kegg_meta_final %>%
  filter(CPDs %in% bg_kegg)
saveRDS(kegg_meta_final,'filtered_kegg过滤后的自测背景库.rds')
filtered_kegg <- kegg_meta_final
##
library(clusterProfiler)

# 构造感兴趣的代谢物列表
metabolite_list1 <- VPA.control.table.enrich.sig$`KEGG ID`
metabolite_list2 <- VPA.CER.AP.table.enrich.sig$`KEGG ID`

# 构造代谢途径与代谢物的映射关系 (TERM2GENE)
term2gene = data.frame(
  Term = filtered_kegg$PIDs,
  Gene =filtered_kegg$CPDs
)

# 构造代谢途径描述信息 (TERM2NAME)
term2name = data.frame(
  Term = filtered_kegg$PIDs,
  Name = filtered_kegg$Description
)


ENRICH1 = enricher(gene = metabolite_list1,pvalueCutoff = 1,
                  pAdjustMethod = 'BH',qvalueCutoff = 1,
                  TERM2GENE = term2gene,TERM2NAME = term2name)
ENRICH2 = enricher(gene = metabolite_list2,pvalueCutoff = 1,
                   pAdjustMethod = 'BH',qvalueCutoff = 1,
                   TERM2GENE = term2gene,TERM2NAME = term2name)

sortdf_vpa <- ENRICH1@result[order(ENRICH1@result$RichFactor, decreasing = TRUE), ]
sortdf_vpa_ap <- ENRICH2@result[order(ENRICH2@result$RichFactor, decreasing = TRUE), ]
sortdf_vpa$Description <- factor(sortdf_vpa$Description, levels = sortdf_vpa$Description)
sortdf_vpa_ap$Description <- factor(sortdf_vpa_ap$Description, levels = sortdf_vpa_ap$Description)


pl <- ggplot(sortdf_vpa[1:10,], aes(RichFactor, fct_reorder(Description, RichFactor, .desc = F), 
                                    colour = pvalue)) + 
  geom_point(aes(size = Count)) + 
  scale_size_continuous(range = c(2,10)) + 
  scale_color_gradientn(colours = c("#1F77B4", "#9467BD")) +  # Nature蓝到紫渐变 
  theme_bw() + 
  ylab("") + 
  theme(legend.position = c(1,0), legend.justification = c(1, 0),
        legend.key.size = unit(0.5, "cm"),  # 减小legend图标大小
        legend.text = element_text(size = 8)) +  # 减小legend文字大小
  theme(legend.background = element_blank()) + 
  theme(legend.key = element_blank())

pr <- ggplot(sortdf_vpa_ap[1:10,], 
             aes(RichFactor,fct_reorder(Description, RichFactor, .desc = TRUE), 
                 colour = pvalue)) + 
  geom_point(aes(size=Count)) + 
  scale_color_gradientn(colours=c("#1F77B4", "#9467BD")) + 
  scale_size_continuous(range = c(2, 10)) + 
  ylab("") + 
  scale_y_discrete(position = "right") + #把term放到右侧 put term on the right 
  theme_bw() + 
  theme(legend.position=c(0, 0), legend.justification = c(0, 0),
        legend.key.size = unit(0.5, "cm"),  # 减小legend图标大小
        legend.text = element_text(size = 8)) +  # 减小legend文字大小
  theme(legend.background = element_blank()) + 
  theme(legend.key = element_blank())

library(cowplot)
plot_grid(pl, pr, labels = "")
ggsave('代谢组左vpa右vpa-ap富集分析图.pdf',width = 15,height = 5)

library(ComplexHeatmap) 
library(circlize) 
library(data.table) 

# Set heatmap color scheme
heatmap.BlWtRd <- c("#1F66AC", "#75AFD3", "grey90", "#FAB99B", "#B2192B")
brota_filter <- readRDS("E://Vincent Na-VPA/brota_filter.rds")
brota_filter


sinfo <- readxl::read_xlsx('Vicent代谢组学.xlsx',sheet = 2)
row_anno <- metadata_pancreas
metadata_pancreas$`Class I(English)`
metadata_pancreas$Name

# 准备热图数据：行是代谢物名称，列是样本名
# 移除Group列，获取表达矩阵
expr_matrix <- brota_filter %>% dplyr::select(-Group)

# 转置矩阵，使行变为代谢物，列变为样本
expr_matrix <- t(expr_matrix)

# 确保行名是代谢物名称
if (!all(rownames(expr_matrix) %in% metadata_pancreas$Name)) {
  # 如果不完全匹配，尝试更新行名
  common_metabolites <- intersect(rownames(expr_matrix), metadata_pancreas$Name)
  expr_matrix <- expr_matrix[common_metabolites, ]
  metadata_pancreas <- metadata_pancreas[metadata_pancreas$Name %in% common_metabolites, ]
}

# 根据expr_matrix的行顺序重新排序metadata_pancreas
metadata_pancreas <- metadata_pancreas[match(rownames(expr_matrix), metadata_pancreas$Name), ]

# 获取样本分组信息
brota_filter$Group <- factor(brota_filter$Group,levels =c(
  'ctrl', 'Na.VPA','CER.AP',  'Na.VPA.CER.AP'
) )
sample_groups <- brota_filter$Group
names(sample_groups) <- colnames(expr_matrix)

# 对样本按照分组因子水平排序
sorted_sample_order <- names(sample_groups)[order(sample_groups)]

# 对表达数据进行标准化处理（按行）
scaled_expr <- t(scale(t(expr_matrix)))

# 准备代谢物大类信息
metabolite_classes <- metadata_pancreas$`Class I(English)`
names(metabolite_classes) <- rownames(expr_matrix)

# 将代谢物大类转换为因子
grouped_metabolites <- data.frame(
  name = rownames(expr_matrix),
  class = metabolite_classes,
  stringsAsFactors = FALSE
)

# 按代谢物大类排序
grouped_metabolites <- grouped_metabolites %>%
  arrange(class)

# 获取排序后的行顺序
sorted_row_order <- grouped_metabolites$name

# 计算每个代谢物大类的数量
class_counts <- table(metabolite_classes)

# 创建颜色映射
# 为样本分组创建颜色
group_colors <- c(
  'ctrl' = '#ECECEC', 
  'Na.VPA' = '#244B84',
  'CER.AP' = '#DC716E', 
  'Na.VPA.CER.AP' = '#850E00'
)

# 为代谢物大类创建颜色
n_classes <- length(unique(metabolite_classes))
class_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_classes)
names(class_colors) <- unique(metabolite_classes)

# 基于排序后的行列顺序重新构建热图数据
sorted_scaled_expr <- scaled_expr[sorted_row_order, sorted_sample_order]
sorted_metabolite_classes <- metabolite_classes[sorted_row_order]
sorted_sample_groups <- sample_groups[sorted_sample_order]

rownames(brota_by_class) <- paste0(brota_by_class$Group, 
                                   ave(rep(1, nrow(brota_by_class)), 
                                       brota_by_class$Group, 
                                       FUN = seq_along))

# 准备基于brota_by_class的柱状图数据
# 1. 对brota_by_class进行排序，与热图的样本顺序一致
brota_by_class_sorted <- brota_by_class[match(sorted_sample_order, rownames(brota_by_class)), ]

# 2. 为每个代谢物创建对应的大类平均值数据
# 创建一个空的向量来存储每个代谢物对应的大类平均值
bar_values <- numeric(nrow(sorted_scaled_expr))

# 遍历每个代谢物大类
for (class_name in unique(sorted_metabolite_classes)) {
  # 获取该大类下的所有代谢物索引
  class_indices <- which(sorted_metabolite_classes == class_name)
  
  if (class_name %in% colnames(brota_by_class_sorted)) {
    # 计算该大类在所有样本中的平均值
    class_mean <- mean(brota_by_class_sorted[[class_name]])
    # 为该大类下的所有代谢物设置相同的平均值
    bar_values[class_indices] <- class_mean
  } else {
    # 如果该大类不在brota_by_class中，设置为0或NA
    bar_values[class_indices] <- 0
  }
}
bar_values <- scale(bar_values)


# 创建左侧的柱状图注释（基于brota_by_class数据）
row_anno <- rowAnnotation(
  # 左侧的柱状图，使用brota_by_class数据
  barplot = anno_barplot(
    x = bar_values,  # 使用brota_by_class中对应大类的平均值
    width = unit(1, "cm"),
    gp = gpar(fill = class_colors[sorted_metabolite_classes]),
    border = TRUE,
    axis = TRUE,
    axis_param = list(at = c(10, 12, 14, 16), labels = c(10, 12, 14, 16))
  ),
  # 代谢物大类的标签
  class = sorted_metabolite_classes,
  col = list(class = class_colors),
  annotation_width = c(barplot = unit(1, "cm"), class = unit(2, "cm")),
  gap = unit(0.1, "cm")
)

# 创建列注释：显示样本分组
col_anno <- HeatmapAnnotation(
  group = sorted_sample_groups,
  col = list(group = group_colors),
  height = unit(0.5, "cm")
)

# 创建热图

main_heatmap <- Heatmap(
  sorted_scaled_expr,  # 使用排序并标准化后的数据
  name = "Scaled Expression",
  show_row_names = FALSE,  # 不显示行名
  show_column_names = FALSE,  # 不显示列名
  column_names_gp = gpar(fontsize = 8),  # 设置列名字体大小
  col = colorRamp2(c(-2, 0, 2), c("#0A3B59", "white", "#8E1B20")),  # 设置热图颜色
  show_row_dend = FALSE,  # 不显示行聚类树
  show_column_dend = FALSE,  # 不显示列聚类树
  top_annotation = col_anno,  # 添加顶部样本分组注释
  left_annotation = row_anno,  # 添加左侧代谢物大类注释（基于brota_by_class数据）
  cluster_rows = T,  # 不对行进行聚类
  cluster_columns = FALSE  # 不对列进行聚类
)
main_heatmap 

kegg_results <- readxl::read_xlsx('代谢组学3组对比的KEGG富集分析.xlsx',sheet = 1) %>% as.data.frame() 
library(ggpubr)

#绘制富集分析棒棒糖图
colnames(kegg_results)
#pdf(file=outFile,width=7,height=6)
ggdotchart(kegg_results, x="pathwayName", y="ComponentRatio", color = "contrast",
           group = "contrast", 
           palette = "aaas",     #配色方案
           legend = "right",     #图例位置
           sorting = "descending",   #上升排序，区别于desc
           add = "segments",    #增加线段
           rotate = TRUE,       #横向显示
           dot.size = 5,        #圆圈大小
           label = round(kegg_results$pvalue),   #圆圈内数值
           font.label = list(color="white",size=9, vjust=0.5),   #圆圈内数值字体设置
           ggtheme = theme_pubr())

sorted_results <- kegg_results[order(kegg_results$contrast, kegg_results$pvalue), ]

# 提取每个 contrast 的前 5 条记录
top5_by_contrast <- do.call(
  rbind,
  by(sorted_results, sorted_results$contrast, head, n = 10)
)
top5_by_contrast$ComponentRatio
top5_by_contrast

split_ratios <- strsplit(top5_by_contrast$ComponentRatio, "/")
top5_by_contrast$metanumber <- map_chr(split_ratios, 1)  # 提取每个子列表的第1个元素
top5_by_contrast$metanumber <- as.numeric(top5_by_contrast$metanumber)



top5_by_contrast.CER.APvsctrl <- top5_by_contrast[top5_by_contrast$contrast == 'CER.AP-ctrl',]
top5_by_contrast.Na.VPAvsctrl <- top5_by_contrast[top5_by_contrast$contrast == 'Na.VPA-ctrl',]
top5_by_contrast.Na.VPA.CER.APvsCER.AP <- top5_by_contrast[top5_by_contrast$contrast == 'Na.VPA.CER.AP-CER.AP',]

p1 <- ggplot(top5_by_contrast.Na.VPA.CER.APvsCER.AP, 
             aes(x = reorder(pathwayName, -pvalue), y = pvalue, color = metanumber)) +
  geom_segment(aes(xend = pathwayName, yend = 0), color = "grey50") +
  geom_point(size = 5) +
  scale_color_viridis_c(name = "Number") +  #
  coord_flip() +
  labs(x = "Pathway", y = "P-value") +
  theme_cowplot() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12))+
  ggtitle('Na.VPA.CER.APvsCER.AP')
ggsave('Na.VPA.CER.APvsCER.AP.pdf',p1,device='pdf',units='in',width = 7,height = 5)

p2 <- ggplot(top5_by_contrast.Na.VPAvsctrl, 
             aes(x = reorder(pathwayName, -pvalue), y = pvalue, color = metanumber)) +
  geom_segment(aes(xend = pathwayName, yend = 0), color = "grey50") +
  geom_point(size = 5) +
  scale_color_viridis_c(name = "Number") +  #
  coord_flip() +
  labs(x = "Pathway", y = "P-value") +
  theme_cowplot() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12))+
  ggtitle('Na.VPAvsctrl')
ggsave('Na.VPAvsctrl.pdf',p2,device='pdf',units='in',width = 7,height = 5)


p3 <- ggplot(top5_by_contrast.CER.APvsctrl, 
             aes(x = reorder(pathwayName, -pvalue), y = pvalue, color = metanumber)) +
  geom_segment(aes(xend = pathwayName, yend = 0), color = "grey50") +
  geom_point(size = 5) +
  scale_color_viridis_c(name = "Number") +  #
  coord_flip() +
  labs(x = "Pathway", y = "P-value") +
  theme_cowplot() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12))+
  ggtitle('CER.APvsctrl')
ggsave('CER.APvsctrl.pdf',p3,device='pdf',units='in',width = 7,height = 5)

brota_filter <- readRDS("E://Vincent Na-VPA/brota_filter.rds")
brota_filter
colnames(brota_filter)
metadata_pancreas$`Class I(English)`

metabolites_in_data <- colnames(brota_filter)[-which(colnames(brota_filter) == "Group")]

# 3. 筛选metadata_pancreas中存在于brota_filter中的代谢物
metadata_filtered <- metadata_pancreas[metadata_pancreas$Name %in% metabolites_in_data, ]

# 4. 创建大类到代谢物的映射
class_to_metabolites <- split(metadata_pancreas$Name, metadata_filtered$`Class I(English)`)

brota_by_class <- data.frame(Group = brota_filter$Group)

for(class_name in names(class_to_metabolites)) {
  # 获取该大类下的所有代谢物
  mets_in_class <- class_to_metabolites[[class_name]]
  # 确保这些代谢物存在于brota_filter中
  mets_in_data <- mets_in_class[mets_in_class %in% colnames(brota_filter)]
  
  if(length(mets_in_data) > 0) {
    # 计算该大类的平均表达水平
    brota_by_class[[class_name]] <- rowMeans(brota_filter[, mets_in_data, drop = FALSE])
  }
}


correlation_matrix <- cor(brota_by_class[, -1])
brota_by_class_ctrl <- brota_by_class %>% filter(Group == 'ctrl')
brota_by_class_cerap <- brota_by_class %>% filter(Group == 'CER.AP')
brota_by_class_navpa <- brota_by_class %>% filter(Group == 'Na.VPA')
brota_by_class_navpacerap <- brota_by_class %>% filter(Group == 'Na.VPA.CER.AP')


#### 可视化
# 定义p值分组函数
categorize_pvalue <- function(p_values) {
  cut(p_values, 
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"),
      include.lowest = TRUE
  )
}

# 定义对应的颜色方案（可根据您的需求调整）
sig_colors <- c("<0.001" = "#E64B35FF",  # 红色
                "0.001-0.01" = "#DC0000FF",  # 橙色
                "0.01-0.05" = "#00A087FF",  # 金色
                ">=0.05" = "#8491B4FF")  # 灰色



correct_column <- which(colnames(brota_filter) == 'S-Adenosyl-L-methionine')
cor_data <- brota_filter[,correct_column]

mantel02 <- fortify_mantel(cor_data, brota_by_class[, -1], 
                           method = "pearson") %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf), 
                 labels = c("<0.25", "0.25-0.5", ">=0.5"), 
                 right = FALSE), 
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"), 
                       right = FALSE))
quickcor(brota_by_class[, -1], type = "upper") + geom_square() + 
  add_link(mantel02, mapping = aes(colour = p.value, size = r),
           diag.label = F) +
  scale_colour_manual(values = sig_colors) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  remove_axis("x")

correct_column <- which(colnames(brota_filter) == 'S-Adenosyl-L-methionine') #796
correct_column <- which(colnames(brota_filter) == 'Methionine') #623
correct_column <- which(colnames(brota_filter) %in% c('S-Adenosyl-L-methionine',
                                                    'Methionine')) #623

# 定义显著性水平的标签和对应的颜色
sig_labels <- c("p < 0.01", "0.01 ≤ p < 0.05", "p ≥ 0.05")
sig_colors <- c("#EE0000", "#FF8C00", "#999999")  # 红、橙、灰
names(sig_colors) <- sig_labels

# 创建一个函数用于统一处理p值分组
categorize_pvalue <- function(p) {
  cut(p, breaks = c(0, 0.01, 0.05, 1), 
      labels = sig_labels, 
      include.lowest = TRUE)
}


cor_data <- brota_filter[,correct_column]

mantel02 <- fortify_mantel(cor_data,brota_by_class[, -1], 
                            spec.select = list(spec01 = 1, spec02 = 2),
                           method = "pearson") %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf), 
                 labels = c("<0.25", "0.25-0.5", ">=0.5"), 
                 right = FALSE), 
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"), 
                       right = FALSE))
cor_metabo <- quickcor(brota_by_class[, -1], type = "upper") + geom_square() + 
  add_link(mantel02, mapping = aes(colour = p.value, size = r),
           diag.label = F) +  scale_colour_manual(values = sig_colors) +
  scale_fill_gradient2(low = "#0A3B59", mid = "white", high = "#8E1B20", 
                       midpoint = 0.5, limit = c(0, 1))+
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  remove_axis("x")
ggsave('cor_metabo.pdf',cor_metabo,width = 10,height = 5)


### 读取RNA-seq
rnaseq_withmetabo <- readRDS("E://Vincent Na-VPA/rnaseq_withmetabo.rds")
rnaseq_withmetabo 

cor_gene <-  c('Dhfr','Mthfr','Ggh','Shmt1','Shmt2','Tyms','Mthfd1','Mthfd2', 'Mthfd1l','Mthfd2l','Gart','Atic',
   'Mtr','Mtrr','Mat1a','Mat2a','Mat2b','Gnmt','Sardh','Chdh','Ahcy','Amd1','Odc1','Sms','Srm','Mtap',
   'Cbs','Cth','Gss','Gclc','Gclm','Gsr','Gpx1','Gpx3','Gpx4')


rnaseq_withmetabo_cor <- rnaseq_withmetabo[c('group',cor_gene)]
rownames(rnaseq_withmetabo_cor) <- gsub("\\.(\\d+)$", "\\1", rownames(rnaseq_withmetabo_cor))
rownames(rnaseq_withmetabo_cor)<- ifelse(grepl("\\d$", rownames(rnaseq_withmetabo_cor)), 
                       rownames(rnaseq_withmetabo_cor), 
                       paste0(rownames(rnaseq_withmetabo_cor), "3"))

cor_data2 <- cor_data[rownames(rnaseq_withmetabo_cor),]

mantel2 <- fortify_mantel(cor_data2,rnaseq_withmetabo_cor[, -1], 
                           spec.select = list(spec01 = 1, spec02 = 2),
                           method = "pearson") %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf), 
                 labels = c("<0.25", "0.25-0.5", ">=0.5"), 
                 right = FALSE), 
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"), 
                       right = FALSE))
cor_rnaseq <- quickcor(rnaseq_withmetabo_cor[, -1], type = "lower") + geom_square() + 
  add_link(mantel2, mapping = aes(colour = p.value, size = r),
           diag.label = F) +
  scale_colour_manual(values = sig_colors) +
  scale_fill_gradient2(low = "#0A3B59", mid = "white", high = "#8E1B20", 
                       midpoint = 0.5, limit = c(0, 1))+
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  remove_axis("x")
ggsave('cor_rnaseq.pdf',cor_rnaseq,width = 10,height = 5)
cor_metabo+cor_rnaseq

