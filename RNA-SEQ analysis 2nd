library(KEGG.db)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)
data(geneList,package="DOSE")
head(geneList)

#常规富集分析，接DEGs结果
#gene = co_DEG %>% 
#  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"), 
#        OrgDb = org.Mm.eg.db ) #%>%  # from 'symbol' to 'ENTREZID'

gene = genes_ctrl %>% 
  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"), 
        OrgDb = org.Mm.eg.db ) #%>%  # from 'symbol' to 'ENTREZID'

gene = genes_AP %>% 
  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"), 
        OrgDb = org.Mm.eg.db ) #%>%  # from 'symbol' to 'ENTREZID'
#GO enrichanalyse
#KEGG
kk <- enrichKEGG(gene = gene$ENTREZID,
                 #universe      = geneList$ENTREZID,
                 keyType = "kegg", pAdjustMethod = "none",
                 use_internal_data = T,
                 organism   = "mmu", # mmu is mouse，human is hsa
                 pvalueCutoff = 0.05 # filter via pvalue
)%>% DOSE::setReadable(OrgDb='org.Mm.eg.db',keyType='ENTREZID')

kk_enrich <- kk@result
metarelate_enrich <- kk_enrich[kk_enrich$category == "Metabolism",'Description' ]
metarelate_enrich <- sub(" - Mus musculus \\(house mouse\\)$",
                                      "", metarelate_enrich)
metarelate_enrich <- metarelate_enrich[!is.na(metarelate_enrich)]

filtered_kk_enrich <- kk_enrich[kk_enrich$category == "Metabolism" &kk_enrich$pvalue<0.05, ]
filtered_kk_enrich <- filtered_kk_enrich[!is.na(filtered_kk_enrich$Description),]
filtered_kk_enrich$Description <- sub(" - Mus musculus \\(house mouse\\)$",
                                      "", filtered_kk_enrich$Description)
#filtered_kk_enrich <- head(filtered_kk_enrich[order(filtered_kk_enrich$p.adjust), ], 15)
filtered_kk_enrich$Description <- str_wrap(filtered_kk_enrich$Description, width = 30)
meta_kegg_plot <-ggplot(filtered_kk_enrich[1:5,], 
                        aes(x = reorder(Description, -RichFactor), y = RichFactor, 
                            color = pvalue)) + 
  geom_point(aes(size=Count)) + 
  scale_color_gradientn(colours=c("#1F77B4", "#9467BD")) + 
  scale_size_continuous(range = c(2, 10)) + 
  coord_flip() + 
  labs(x = "Pathway Name", y = "RichFactor") + 
  # 移除scale_y_continuous(position = "right")这一行
  theme_bw() + 
  theme(legend.position = "right", 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10), 
        axis.title = element_text(size = 12))+ 
  theme(legend.key = element_blank())+ggtitle('RNA-SEQ enrichment analysis')
meta_kegg_plot 
write.csv(kk_enrich,'转录组富集分析结果.csv')
ggsave('RNA-SEQ enrichment analysis.pdf',,device='pdf',units='in',width = 12,height = 5)



## GSEA分析
alldiff <- res_AP[order(res_AP$stat,decreasing = T),] ###读取差异分析列表
alldiff <- res_ctrl[order(res_ctrl$stat,decreasing = T),] ###读取差异分析列表
genelist <- alldiff$stat
names(genelist) <- alldiff$ENTREZID
kk2 <- gseKEGG(geneList = genelist,
               eps = 0,
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               minGSSize = 10,
               maxGSSize = 500,seed = 42,
               organism     = 'mmu', use_internal_data =T,
               keyType       = "ENTREZID") %>% 
  DOSE::setReadable(OrgDb='org.Mm.eg.db',keyType='ENTREZID')
kk2@result$Description <- gsub(" - Mus musculus \\(house mouse\\)$", "", kk2@result$Description)


metarelate_enrich #之前做好的
kk2_metabolism <- kk2[kk2@result$Description %in% metarelate_enrich, ]
topPathwaysUp <- kk2_metabolism %>%
  filter(NES > 0) %>%           # 筛选NES > 0的行
  arrange(pvalue) %>%           # 按pvalue升序排序
  head(10) %>%                  # 取前10行
  pull(Description)             # 提取Description列的值
topPathwaysDown <- kk2_metabolism %>%
  filter(NES < 0) %>%           # 筛选NES > 0的行
  arrange(pvalue) %>%           # 按pvalue升序排序
  head(10) %>%                  # 取前10行
  pull(Description) 
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
length(topPathways)

kk2_metabolism[kk2_metabolism$Description %in%topPathways, ]

# 先创建过滤后的数据集
data_filtered <- kk2_metabolism[kk2_metabolism$Description %in% topPathways, ]

# 添加group变量：根据NES是否大于0来分组
data_filtered$group <- ifelse(data_filtered$NES > 0, "positive", "negative")

# 分别创建NES>0和NES<0的子集
data_positive <- subset(data_filtered, NES > 0)
data_positive <-data_positive %>%
  arrange(desc(abs(NES)))
data_positive$ID[1:5]

data_negative <- subset(data_filtered, NES < 0)
data_negative <-data_negative%>%
  arrange(desc(abs(NES)))
data_negative$ID[1:5]

gsea_kegg <- ggplot(data_filtered, aes(x = reorder(Description, NES), y = NES, fill = pvalue)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_gradientn(colours = c("#1F77B4", "#9467BD")) + 
  #scale_fill_manual(values = c("#0A3B59",  "#8E1B20"), guide = FALSE)  + 
  xlab("Top metabolism related KEGG pathways") + ylab("NES") + 
  # NES>0的文字显示在左侧（使用hjust=1右对齐） 
  geom_text(data = data_positive, 
            aes(x = reorder(Description, NES), y = 0, label = Description, color = group), 
            size = 3, 
            hjust = 1,  # 右对齐，文本显示在条形图左侧 
            nudge_x = -0.1) +  # 稍微向左调整位置 
  geom_text(data = data_negative, 
            aes(x = reorder(Description, NES), y = 0, label = Description, color = group), 
            size = 3, 
            hjust = 0,  # 左对齐，文本显示在条形图右侧 
            nudge_x = 0.1) +  # 稍微向右调整位置 
  theme_bw(8) + 
  theme(panel.grid = element_blank()) + 
  theme(panel.border = element_rect(size = 0.6)) + 
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  # 设置文本颜色与填充色一致 
  scale_color_manual(values = c("#0A3B59",  "#8E1B20"), guide = FALSE)

# 显示图形
print(gsea_kegg)

ggsave('gsea_kegg.pdf',gsea_kegg,height = 5,width = 7)

#terms <- c('mmu00790', 'mmu00670','mmu00270','mmu00260','mmu00230')
terms <- c('mmu00670','mmu00260')
#tmp_gesa <- kk2@result[kk2@result$ID %in% terms,c('NES','ID','Description','p.adjust')]
#gseaplot_rnaseq <- gseaplot2(kk2, geneSetID = terms, pvalue_table = FALSE,rel_heights = c(1.5, 0.5, 0.5),base_size = 15,subplots=1:2,
#          color = c("#4A475C", "#ACD4D6", "#A6BAAF",'#B98A82','#E4DBD2'), ES_geom = "line")
gseaplot_rnaseq <- gseaplot2(kk2, geneSetID = terms, pvalue_table = FALSE,rel_heights = c(1.5, 0.5, 0.5),base_size = 15,subplots=1:2,
                             color = c("#4A475C", "#ACD4D6"), ES_geom = "line")
ggsave('gseaplot_rnaseq.pdf',gseaplot_rnaseq,width = 5,height = 5,units = 'in')
# 创建原始gseaplot
data_positive$ID[1:5]
data_negative$ID[1:5]
gseaplot_pos_rnaseq <- gseaplot2(kk2, 
                             geneSetID = data_positive$ID[1:5], 
                             pvalue_table = FALSE, 
                             rel_heights = c(1.5, 0.5, 0.5), 
                             base_size = 15, 
                             subplots = 1, 
                             color = c("#4A475C", "#ACD4D6", "#A6BAAF",'#B98A82','#E4DBD2'), 
                             ES_geom = "dot") + 
  # 添加黑色边框
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8)) + 
  # 如果需要，也可以设置plot.border
  theme(plot.border = element_rect(color = "black", fill = NA, size = 1))
ggsave('gseaplot_pos_rnaseq.pdf',gseaplot_pos_rnaseq,width = 5,height = 5,units = 'in')

gseaplot_neg_rnaseq <- gseaplot2(kk2, 
                                 geneSetID = data_negative$ID[1:5], 
                                 pvalue_table = FALSE, 
                                 rel_heights = c(1.5, 0.5, 0.5), 
                                 base_size = 15, 
                                 subplots = 1, 
                                 color = c("#4A475C", "#ACD4D6", "#A6BAAF",'#B98A82','#E4DBD2'), 
                                 ES_geom = "dot") + 
  # 添加黑色边框
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8)) + 
  # 如果需要，也可以设置plot.border
  theme(plot.border = element_rect(color = "black", fill = NA, size = 1))
ggsave('gseaplot_neg_rnaseq.pdf',gseaplot_neg_rnaseq,width = 5,height = 5,units = 'in')

write.csv(kk2,'gsea_RNAseq.csv')
