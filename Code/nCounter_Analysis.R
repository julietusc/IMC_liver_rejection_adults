
#######################################################
# nCounter Transcriptomic and TCR Expression Analysis #
#######################################################

library(readr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

setwd('/users/sarahbangerth/Desktop/nCounter/Emamuallee_032023_final/NR vs TCMR Advanced Analysis 2023-03-30 22-59/results')


norm_data = read_csv('Normalization/mRNA_normalized_data_log2_counts.csv')
head(norm_data)

# Define markers for each cell type

b_cells = c('BLK-mRNA','TNFRSF17-mRNA','SPIB-mRNA','TCL1A-mRNA','FAM30A-mRNA','CD19-mRNA',
            'PNOC-mRNA','FCRL2-mRNA')

dc = c('CCL13-mRNA','HSD11B1-mRNA','CD209-mRNA')

nk = c('KIR3DL1/2-mRNA','IL21R-mRNA','XCL1/2-mRNA','NCR1-mRNA')

neutrophils = c('FPR1-mRNA','S100A12-mRNA','FCAR-mRNA','FCGR3A/B-mRNA','CEACAM3-mRNA','CSF3R-mRNA','SIGLEC5-mRNA')

mast_cells = c('TPSAB1/B2-mRNA','CPA3-mRNA','MS4A2-mRNA','HDC-mRNA')

macrophages = c('MS4A4A-mRNA','CD80-mRNA','IL12B-mRNA','IL23A-mRNA','TLR2-mRNA','TLR4-mRNA','IL27-mRNA','CXCL9-mRNA','CXCL10-mRNA',
                'CXCL11-mRNA','CXCL16-mRNA','TGFB1-mRNA','CCL17-mRNA','CCL22-mRNA','CCL24-mRNA','IL10-mRNA','CCL1-mRNA',
                'CCL18-mRNA')

t_cells = c('CD45RA-mRNA','CD45R0-mRNA','CD3E-mRNA','CD3D-mRNA',
            'TRAT1-mRNA','SH2D1A-mRNA','TBX21-mRNA','TNFRSF9-mRNA')

t_cell_activation = c('PDCD1-mRNA','CD274-mRNA','PDCD1LG2-mRNA','LAG3-mRNA','CD160-mRNA','CTLA4-mRNA','HAVCR2-mRNA',
                      'IL10-mRNA','IL6-mRNA','STAT3-mRNA')

cytotoxic = c('PRF1-mRNA','GZMA-mRNA','CTSW-mRNA',
              'GZMH-mRNA','GNLY-mRNA')

all_genes = c('BLK-mRNA','MS4A1-mRNA','TNFRSF17-mRNA','SPIB-mRNA','TCL1A-mRNA','FAM30A-mRNA','CD19-mRNA',
              'PNOC-mRNA','FCRL2-mRNA','CCL13-mRNA','HSD11B1-mRNA','CD209-mRNA','KIR3DL1/2-mRNA','IL21R-mRNA',
              'XCL1/2-mRNA','NCR1-mRNA','FPR1-mRNA','S100A12-mRNA','FCAR-mRNA','FCGR3A/B-mRNA','CEACAM3-mRNA',
              'CSF3R-mRNA','SIGLEC5-mRNA','TPSAB1/B2-mRNA','CPA3-mRNA','MS4A2-mRNA','HDC-mRNA',
              'CD68-mRNA','CD163-mRNA','MS4A4A-mRNA','CD84-mRNA',
              'CD80-mRNA','IL12A-mRNA','IL12B-mRNA','IL23A-mRNA','IL23R-mRNA','IL27-mRNA',
              'CXCL9-mRNA','CXCL16-mRNA','CXCL11-mRNA','CXCL10-mRNA',
              'CCL5-mRNA','TLR2-mRNA','TLR4-mRNA','TGFB1-mRNA','CCL17-mRNA','CCL22-mRNA','CCL24-mRNA','IL10-mRNA','CCL1-mRNA',
              'CCL18-mRNA','CXCL13-mRNA',
              'PTPRC-mRNA','CD45RA-mRNA','CD45R0-mRNA','CD3E-mRNA','CD3D-mRNA','CD3G-mRNA',
              'TRAT1-mRNA','SH2D1A-mRNA','CD6-mRNA','CD8A-mRNA','CD8B-mRNA','CD4-mRNA',
              'FOXP3-mRNA','CD274-mRNA','PDCD1LG2-mRNA',
              'LAG3-mRNA','CD244-mRNA','CD160-mRNA','CTLA4-mRNA','PDCD1-mRNA','HAVCR2-mRNA','TNFRSF9-mRNA',
              'BATF-mRNA','IL6-mRNA','STAT3-mRNA','KLRG1-mRNA','TBX21-mRNA','EOMES-mRNA',
              'GZMB-mRNA','PRF1-mRNA','KLRB1-mRNA','KLRK1-mRNA','GZMA-mRNA','CTSW-mRNA',
              'GZMH-mRNA','NKG7-mRNA','KLRD1-mRNA','GNLY-mRNA')

# scale and center expression
scaled_norm_data = norm_data
scaled_norm_data[,c(all_genes)] <- scale(norm_data[,c(all_genes)])

### B CELLS, DC, MAST CELLS,NEUTROPHILS

m <- as.matrix((scaled_norm_data[,c(b_cells,neutrophils,dc,mast_cells)]))
rownames(m) <- as.vector(scaled_norm_data[,'Sample'])$Sample

m = m[c('SP191689A1','SP197279A1','SP206085A1','SP218904A1',
        'SP203382A1','SP207572A1','SP213775A1','SP21414A1'),]

col_fun = colorRamp2(c(1,2,3,4,5,6,7,8),
                     c("#E0B0FF","#E0B0FF","#E0B0FF","#E0B0FF","#DA70D6","#DA70D6","#DA70D6","#DA70D6"))

ha_left = HeatmapAnnotation(Group = 1:8, col = list(Group = col_fun), annotation_name_rot=90, 
                            gp = gpar(col = "black"), which='row', show_legend = FALSE)

nCounter_samples1 = Heatmap(m, name = "Scaled expression",border='black',
                            cluster_columns = FALSE,
                            show_column_dend = FALSE,
                            cluster_rows = FALSE,
                            show_row_dend = FALSE,
                            column_names_rot = 90,
                            column_names_centered = FALSE,
                            show_column_names = TRUE,
                            row_title_gp = gpar(fontsize = 23),
                            left_annotation = ha_left,
                            rect_gp = gpar(col = "white", lwd = 2),
                            col = colorRamp2(c(-1, 0, 2), hcl_palette='viridis'),
                            heatmap_legend_param = list(at = c(-1:2), legend_width = unit(6,"cm"),
                                                        direction="horizontal", title_gp = gpar(fontsize=12),
                                                        title_position = "topcenter",
                                                        labels_gp = gpar(fontsize=12)),
                            column_names_side = "bottom",
                            row_names_side = 'left',
                            column_labels = c('BLK','CD269','SPIB','TCL1A','FAM30A','CD19',
                                              'PNOC','FCRL2','FPR1','S100A12','FCAR','CD16A/B',
                                              'CEACAM3','CSF3R','SIGLEC5','CCL13','HSD11B1','CD209',
                                              'TPSAB1/B2','CPA3','MS4A2','HDC'),
                            height = unit(8, "cm"),
                            width = unit(20,"cm"))

draw(nCounter_samples1, heatmap_legend_side = "top")

### MACROPHAGES, NK CELLS

m <- as.matrix((scaled_norm_data[,c(macrophages,nk)]))
rownames(m) <- as.vector(scaled_norm_data[,'Sample'])$Sample

m = m[c('SP191689A1','SP197279A1','SP206085A1','SP218904A1',
        'SP203382A1','SP207572A1','SP213775A1','SP21414A1'),]

col_fun = colorRamp2(c(1,2,3,4,5,6,7,8),
                     c("#E0B0FF","#E0B0FF","#E0B0FF","#E0B0FF","#DA70D6","#DA70D6","#DA70D6","#DA70D6"))

ha_left = HeatmapAnnotation(Group = 1:8, col = list(Group = col_fun), annotation_name_rot=90, 
                            gp = gpar(col = "black"), which='row', show_legend = FALSE)

nCounter_samples3 = Heatmap(m, name = "Scaled expression",border='black',
                            cluster_columns = FALSE,
                            show_column_dend = FALSE,
                            cluster_rows = FALSE,
                            show_row_dend = FALSE,
                            column_names_rot = 90,
                            column_names_centered = FALSE,
                            show_column_names = TRUE,
                            row_title_gp = gpar(fontsize = 23),
                            left_annotation = ha_left,
                            rect_gp = gpar(col = "white", lwd = 2),
                            col = colorRamp2(c(-1, 0, 2), hcl_palette='viridis'),
                            heatmap_legend_param = list(at = c(-1:2), legend_width = unit(6,"cm"),
                                                        direction="horizontal", title_gp = gpar(fontsize=12),
                                                        title_position = "topcenter",
                                                        labels_gp = gpar(fontsize=12)),
                            column_names_side = "bottom",
                            row_names_side = 'left',
                            column_labels = c('MS4A4A','CD80','IL12B','IL23A','TLR2','TLR4','IL27','CXCL9','CXCL10',
                                              'CXCL11','CXCL16','TGFB1','CCL17','CCL22','CCL24','IL10','CCL1','CCL18',
                                              'NKB1','IL21R','XCL1','NCR1'),
                            height = unit(8, "cm"),
                            width = unit(20,"cm"))

draw(nCounter_samples3, heatmap_legend_side = "top")

### T CELLS

m <- as.matrix((scaled_norm_data[,c(t_cells,cytotoxic,t_cell_activation)]))
rownames(m) <- as.vector(scaled_norm_data[,'Sample'])$Sample

m = m[c('SP191689A1','SP197279A1','SP206085A1','SP218904A1',
        'SP203382A1','SP207572A1','SP213775A1','SP21414A1'),]

col_fun = colorRamp2(c(1,2,3,4,5,6,7,8),
                     c("#E0B0FF","#E0B0FF","#E0B0FF","#E0B0FF","#DA70D6","#DA70D6","#DA70D6","#DA70D6"))

ha_left = HeatmapAnnotation(Group = 1:8, col = list(Group = col_fun), annotation_name_rot=90, 
                            gp = gpar(col = "black"), which='row', show_legend = FALSE)

nCounter_samples4 = Heatmap(m, name = "Scaled expression",border='black',
                            cluster_columns = FALSE,
                            show_column_dend = FALSE,
                            cluster_rows = FALSE,
                            show_row_dend = FALSE,
                            column_names_rot = 90,
                            column_names_centered = FALSE,
                            show_column_names = TRUE,
                            row_title_gp = gpar(fontsize = 23),
                            left_annotation = ha_left,
                            rect_gp = gpar(col = "white", lwd = 2),
                            col = colorRamp2(c(-1, 0, 2), hcl_palette='viridis'),
                            heatmap_legend_param = list(at = c(-1:2), legend_width = unit(6,"cm"),
                                                        direction="horizontal", title_gp = gpar(fontsize=12),
                                                        title_position = "topcenter",
                                                        labels_gp = gpar(fontsize=12)),
                            column_names_side = "bottom",
                            row_names_side = 'left',
                            column_labels = c('CD45RA','CD45R0','CD3E','CD3D','TRAT1','SH2D1A','TBX21','4-1BB',
                                              'PRF1','GZMA','CTSW','GZMH','GNLY','PD1','PD-L1','PD-L2',
                                              'LAG3','CD160','CTLA4','HAVCR2',
                                              'IL10','IL6','STAT3'),
                            height = unit(8, "cm"),
                            width = unit(20,"cm"))

draw(nCounter_samples4, heatmap_legend_side = "top")

################################################################
### Log2FoldChange for TCR in each TCMR sample vs the NR avg ###
################################################################

tcrA_genes = c('TRAV1-1-mRNA','TRAV1-2-mRNA','TRAV10-mRNA','TRAV11-mRNA','TRAV12-1-mRNA','TRAV12-2-mRNA','TRAV12-3-mRNA',
               'TRAV13-1-mRNA','TRAV13-2-mRNA','TRAV14-mRNA','TRAV16-mRNA','TRAV17-mRNA','TRAV18-mRNA','TRAV19-mRNA',
               'TRAV2-mRNA','TRAV20-mRNA','TRAV21-mRNA','TRAV22-mRNA','TRAV23-mRNA','TRAV24-mRNA','TRAV25-mRNA',
               'TRAV26-1-mRNA','TRAV26-2-mRNA','TRAV27-mRNA','TRAV29-mRNA','TRAV3-mRNA','TRAV30-mRNA','TRAV34-mRNA',
               'TRAV35-mRNA','TRAV36-mRNA','TRAV38-1-mRNA','TRAV38-2-mRNA','TRAV39-mRNA','TRAV4-mRNA','TRAV40-mRNA',
               'TRAV41-mRNA','TRAV5-mRNA','TRAV6-mRNA','TRAV7-mRNA','TRAV8-1-mRNA','TRAV8-2-mRNA','TRAV8-3-mRNA',
               'TRAV8-6-mRNA','TRAV9-1-mRNA','TRAV9-2-mRNA')

tcrB_genes = c('TRBV10-1-mRNA','TRBV10-2-mRNA','TRBV10-3-mRNA','TRBV11-1-mRNA',
               'TRBV11-2-mRNA','TRBV11-3-mRNA','TRBV12-3-mRNA','TRBV12-5-mRNA','TRBV13-mRNA','TRBV14-mRNA','TRBV15-mRNA',
               'TRBV16-mRNA','TRBV18-mRNA','TRBV19-mRNA','TRBV2-mRNA','TRBV20-1-mRNA','TRBV24-1-mRNA','TRBV25-1-mRNA',
               'TRBV27-mRNA','TRBV28-mRNA','TRBV29-1-mRNA','TRBV3-1-mRNA','TRBV30-mRNA','TRBV4-1-mRNA','TRBV4-2-mRNA',
               'TRBV4-3-mRNA','TRBV5-1-mRNA','TRBV5-4-mRNA','TRBV5-5-mRNA','TRBV5-6-mRNA','TRBV5-8-mRNA','TRBV6-1-mRNA',
               'TRBV6-2-mRNA','TRBV6-4-mRNA','TRBV6-5-mRNA','TRBV6-6-mRNA','TRBV6-8-mRNA','TRBV6-9-mRNA','TRBV7-2-mRNA',
               'TRBV7-3-mRNA','TRBV7-4-mRNA','TRBV7-6-mRNA','TRBV7-7-mRNA','TRBV7-8-mRNA','TRBV7-9-mRNA','TRBV9-mRNA')

tcrDG_genes = c('TRDV1-mRNA','TRDV2-mRNA','TRDV3-mRNA','TRGV2-mRNA','TRGV3/5-mRNA','TRGV4-mRNA','TRGV8-mRNA','TRGV9-mRNA')

# Calculate NR average for each gene (baseline)
head(norm_data)

avg_NR = apply(norm_data[norm_data$Sample %in% c('SP218904A1','SP206085A1','SP197279A1','SP191689A1'),-c(1)], 2, mean)

# Extract genes for each sample (function is not actually doing an average since we're only specifying one sample)
TCMR1 = apply(norm_data[norm_data$Sample %in% c('SP203382A1'),-c(1)], 2, mean)
TCMR2 = apply(norm_data[norm_data$Sample %in% c('SP207572A1'),-c(1)], 2, mean)
TCMR3 = apply(norm_data[norm_data$Sample %in% c('SP213775A1'),-c(1)], 2, mean)
TCMR4 = apply(norm_data[norm_data$Sample %in% c('SP21414A1'),-c(1)], 2, mean)

foldchange_df1 = as.data.frame(cbind(avg_NR, TCMR1))
foldchange_df2 = as.data.frame(cbind(foldchange_df1, TCMR2))
foldchange_df3 = as.data.frame(cbind(foldchange_df2, TCMR3))
foldchange_df = as.data.frame(cbind(foldchange_df3, TCMR4))

head(foldchange_df)

# Calculate fold change for each TCMR sample using the NR average as baseline
foldchange_df$log2FoldChange_TCMR1 = foldchange_df$TCMR1 - foldchange_df$avg_NR
foldchange_df$log2FoldChange_TCMR2 = foldchange_df$TCMR2 - foldchange_df$avg_NR
foldchange_df$log2FoldChange_TCMR3 = foldchange_df$TCMR3 - foldchange_df$avg_NR
foldchange_df$log2FoldChange_TCMR4 = foldchange_df$TCMR4 - foldchange_df$avg_NR

head(foldchange_df)

# Put together DF with all TCMR samples for visualization
tempdf1 = as.data.frame(foldchange_df$log2FoldChange_TCMR1)
colnames(tempdf1) = 'log2FoldChange'
tempdf1$Sample = 'SP203382A1'
tempdf1$Gene = rownames(foldchange_df)
head(tempdf1)

tempdf2 = as.data.frame(foldchange_df$log2FoldChange_TCMR2)
colnames(tempdf2) = 'log2FoldChange'
tempdf2$Sample = 'SP207572A1'
tempdf2$Gene = rownames(foldchange_df)

tempdf3 = as.data.frame(foldchange_df$log2FoldChange_TCMR3)
colnames(tempdf3) = 'log2FoldChange'
tempdf3$Sample = 'SP213775A1'
tempdf3$Gene = rownames(foldchange_df)

tempdf4 = as.data.frame(foldchange_df$log2FoldChange_TCMR4)
colnames(tempdf4) = 'log2FoldChange'
tempdf4$Sample = 'SP21414A1'
tempdf4$Gene = rownames(foldchange_df)

log2FoldChange_plot1 = rbind(tempdf1,tempdf2)
log2FoldChange_plot2 = rbind(log2FoldChange_plot1,tempdf3)
log2FoldChange_plot = rbind(log2FoldChange_plot2,tempdf4)
head(log2FoldChange_plot)

# Add UP/DOWN labels for plot color
log2FoldChange_plot = log2FoldChange_plot %>% mutate(Color = case_when(log2FoldChange > 0 ~ 'UP',
                                                                       log2FoldChange < 0 ~ 'DOWN',
                                                                       log2FoldChange == 0 ~ 'NONE'))

summary(log2FoldChange_plot$log2FoldChange)

# Visualize TCR alpha
samples_tcrA = log2FoldChange_plot[(log2FoldChange_plot$Gene %in% tcrA_genes)&(log2FoldChange_plot$Color=='UP'),] %>%
  group_by(Sample) %>% dplyr::summarize(Upregulated = n())
samples_tcrA$Total = 45
samples_tcrA$Prop_up = samples_tcrA$Upregulated / samples_tcrA$Total

ggplot(data = samples_tcrA,
       aes(x = Sample, y = Prop_up)) +
  geom_col(position = position_dodge(width = 0.8), width=0.8) +
  theme_bw() +
  ggtitle("") + xlab("") + ylab("% of upregulated\nTRAV genes") +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30),
        plot.title=element_text(size = 34, hjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=30),
        aspect.ratio=1,
        panel.grid=element_blank(),
        strip.text=element_text(size = 34, hjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) 



# Visualize TCR beta
samples_tcrB = log2FoldChange_plot[(log2FoldChange_plot$Gene %in% tcrB_genes)&(log2FoldChange_plot$Color=='UP'),] %>%
  group_by(Sample) %>% dplyr::summarize(Upregulated = n())
samples_tcrB$Total = 45
samples_tcrB$Prop_up = samples_tcrB$Upregulated / samples_tcrB$Total

ggplot(data = samples_tcrB,
       aes(x = Sample, y = Prop_up)) +
  geom_col(position = position_dodge(width = 0.8), width=0.8) +
  theme_bw() +
  ggtitle("") + xlab("") + ylab("% of upregulated\nTRBV genes") +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30),
        plot.title=element_text(size = 34, hjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=30),
        aspect.ratio=1,
        panel.grid=element_blank(),
        strip.text=element_text(size = 34, hjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) 


# Visualize TCR delta & gamma
samples_tcrDG = log2FoldChange_plot[(log2FoldChange_plot$Gene %in% tcrDG_genes)&(log2FoldChange_plot$Color=='UP'),] %>%
  group_by(Sample) %>% dplyr::summarize(Upregulated = n())
samples_tcrDG$Total = 45
samples_tcrDG$Prop_up = samples_tcrDG$Upregulated / samples_tcrDG$Total

ggplot(data = samples_tcrDG,
       aes(x = Sample, y = Prop_up)) +
  geom_col(position = position_dodge(width = 0.8), width=0.8) +
  theme_bw() +
  ggtitle("") + xlab("") + ylab("% of upregulated\nTRGV and TRDV genes") +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30),
        plot.title=element_text(size = 34, hjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=30),
        aspect.ratio=1,
        panel.grid=element_blank(),
        strip.text=element_text(size = 34, hjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,1)) 

#######################
#### External data ####
#######################

library(readxl)
fc_data = read_excel("/users/sarahbangerth/Downloads/mmc3.xls")
head(fc_data)
colnames(fc_data)[1] = 'Pt1 Baseline'

library(reshape2)
fc_data.melt = melt(fc_data, id.vars = c('gene','isoform'), 
                    measure.vars = c('Pt1 Baseline','Pt2 Baseline','Pt3 Baseline','Pt4 Baseline',
                                     'Pt5 Baseline','Pt6 Baseline','Pt1 W4','Pt2 W4',"Pt3 W4","Pt4 W4",
                                     'Pt5 W4','Pt6 W4'))

head(fc_data.melt)

# Add UP/DOWN labels for plot color
fc_data.melt = fc_data.melt %>% mutate(Color = case_when(value > 0 ~ 'UP',
                                                         value < 0 ~ 'DOWN',
                                                         value == 0 ~ 'NONE'))

all_genes = c('IL12B','IL23A','TLR2','TLR4','IL27','CD80','CXCL9','CXCL10',
              'CXCL11','CCL17','CCL22','CCL24','IL10','CCL1','CXCL16','CCL18',
              'MS4A4A','TGFB1','KIR3DL1','IL21R','XCL1','NCR1',
              'CD45RA','CD45R0','CD3E','CD3D','TRAT1','SH2D1A','TBX21','CD274','PDCD1LG2',
              'PRF1','GZMA','CTSW','GZMH','GNLY','LAG3','CD160','CTLA4','PDCD1','HAVCR2',
              'TNFRSF9','IL10','IL6','TGFB1','STAT3')

ggplot(fc_data.melt[fc_data.melt$gene %in% all_genes,], 
       aes(forcats::fct_rev(isoform),variable, fill = Color, size = value)) +
  geom_point(shape = 21, stroke = 0) +
  # geom_hline(yintercept = seq(.5, 4.5, 1), size = .2) +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(0, 13)) +
  # scale_fill_gradient(colours=viridis, breaks = c(-1, 0, 1, 2), labels = c("-1", "0", "1", "2"), limits = c(-1, 2)) +
  # scale_fill_viridis(breaks = c(-2,-1, 0, 1, 2), labels = c("-2","-1", "0", "1", "2")) +
  scale_fill_manual(values = c("#313695","#FFEA00")) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.text = element_text(size=15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25), 
                             label.position = "bottom",
                             title.position = "right", 
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  scale_y_discrete(limits=c('Pt6 W4','Pt5 W4',"Pt4 W4","Pt3 W4",'Pt2 W4','Pt1 W4',
                            'Pt6 Baseline','Pt5 Baseline','Pt4 Baseline','Pt3 Baseline','Pt2 Baseline','Pt1 Baseline')) +
  scale_x_discrete(limits=c('PRF1-001','GZMA-001',
                            'CD274-001','PDCD1LG2-001','HAVCR2-001','STAT3-001','STAT3-002','STAT3-012',
                            'MS4A4A-001','TLR2-001','TLR4-001','TLR4-002','IL27-001','CXCL9-001','CXCL10-001',
                            'CXCL16-001','CXCL16-002','CXCL16-003','CXCL16-005','TGFB1-001'),
                   labels=c('PRF1-001','GZMA-001',
                            'PD-L1-001','PD-L2-001','HAVCR2-001','STAT3-001','STAT3-002','STAT3-012',
                            'MS4A4A-001','TLR2-001','TLR4-001','TLR4-002','IL27-001','CXCL9-001','CXCL10-001',
                            'CXCL16-001','CXCL16-002','CXCL16-003','CXCL16-005','TGFB1-001'),
                   position='bottom') +
  labs(size = "log2 Fold Change", fill = "z-score:", x = NULL, y = NULL)


#######################################################
# CORRELATION BETWEEN nCOUNTER DATA AND EXTERNAL DATA #
#######################################################

# nCounter
head(norm_data)
avg_NR = apply(norm_data[norm_data$Sample %in% c('SP218904A1','SP206085A1','SP197279A1','SP191689A1'),-c(1)], 2, mean)
avg_TCMR = apply(norm_data[norm_data$Sample %in% c('SP21414A1','SP203382A1','SP207572A1','SP213775A1'),-c(1)], 2, mean)
foldchange_ncounter = as.data.frame(cbind(avg_NR, avg_TCMR))
foldchange_ncounter$FoldChange = foldchange_ncounter$avg_TCMR - foldchange_ncounter$avg_NR
foldchange_ncounter$gene = sapply(strsplit(rownames(foldchange_ncounter), "-"), "[", 1)
foldchange_ncounter$Data = 'nCounter'
head(foldchange_ncounter)

# External
head(fc_data.melt)
avg_baseline = fc_data.melt[fc_data.melt$variable %in% c('Pt1 Baseline','Pt2 Baseline','Pt3 Baseline','Pt4 Baseline'),] %>%
  group_by(gene) %>% dplyr::summarize(avg_baseline = mean(value))
avg_rj = fc_data.melt[fc_data.melt$variable %in% c('Pt1 W4','Pt2 W4','Pt3 W4','Pt4 W4'),] %>%
  group_by(gene) %>% dplyr::summarize(avg_rj = mean(value))
foldchange_ext = merge(avg_baseline,avg_rj, by='gene')
foldchange_ext$FoldChange = foldchange_ext$avg_rj - foldchange_ext$avg_baseline
foldchange_ext$Data = 'Ext'
head(foldchange_ext)

# Merge nCounter and External
foldchange_df = rbind(foldchange_ncounter[,c('gene','FoldChange','Data')],
                      foldchange_ext[,c('gene','FoldChange','Data')])

nk = c('KIR3DL1','IL21R','XCL1','NCR1')
macrophages = c('MS4A4A','CD80','IL12B','IL23A','TLR2','TLR4','IL27','CXCL9','CXCL10',
                'CXCL11','CXCL16','TGFB1','CCL17','CCL22','CCL24','IL10','CCL1','CCL18')
t_cells = c('CD45RA','CD45R0','CD3E','CD3D','TRAT1','SH2D1A','TBX21','TNFRSF9')
t_cell_activation = c('PDCD1','CD274','PDCD1LG2','LAG3','CD160','CTLA4','HAVCR2','IL10','IL6','STAT3')
cytotoxic = c('PRF1','GZMA','CTSW','GZMH','GNLY')

foldchange_df[foldchange_df$gene %in% c(t_cells,cytotoxic,t_cell_activation),]

ggplot(foldchange_df[foldchange_df$gene %in% c('STAT3','PRF1','PDCD1LG2','CD274','HAVCR2','GZMA'
                                               ,'MS4A4A','TLR2','TLR4','IL27','CXCL9','CXCL10','CXCL16','TGFB1'),], 
       aes(gene, FoldChange, fill=Data)) +
  geom_col(position = position_dodge(width = 0.8), width=0.8) +
  theme_bw() +
  ggtitle(NULL) + xlab("") + ylab("Fold Change") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=30),
        plot.title=element_text(size = 34, hjust = 0.5),
        # legend.position='none',
        legend.title=element_blank(),
        legend.text=element_text(size=30),
        aspect.ratio=1,
        panel.grid=element_blank(),
        strip.text=element_text(size = 34, hjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  coord_flip() +
  scale_fill_manual(values = c("gray","black")) +
  scale_x_discrete(limits=c('TGFB1','CXCL16','CXCL10','CXCL9','IL27','TLR4','TLR2','MS4A4A',
                            'STAT3','HAVCR2','PDCD1LG2','CD274','GZMA','PRF1'),
                   labels=c('TGFB1','CXCL16','CXCL10','CXCL9','IL27','TLR4','TLR2','MS4A4A',
                            'STAT3','HAVCR2','PD-L2','PD-L1','GZMA','PRF1'))


