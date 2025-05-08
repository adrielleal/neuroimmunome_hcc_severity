#Papper Main Figures
# Librarys ----------------------------------------------------------------
{
  library(tidyverse)
  library(dplyr)
  library(scales)
  library(ggnet)
  library(ggplot2)
  library(tidyverse)
  library(igraph)
  library(GGally)
  library(network)
  library(sna)
  library(RColorBrewer)
  library(intergraph)
  library(openxlsx)
  library(ggrepel)
  library(stats)
  library(rrvgo)
  library(stringr)
  library(shinydashboard)
  library(heatmaply)
  library(seriation)
  library(wordcloud)
  library(grid)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggVennDiagram)
  library(UpSetR)
  library(org.Rn.eg.db)
  library(DOSE)
  library(ComplexHeatmap)
  library(ggnewscale)
  library(enrichplot)
  library(clusterProfiler)
  library(readxl)
  library(rstatix)
  library(reshape2)
  library(ggpubr)
}

# MetaDEGs Characterization
# DivergiPlot (metaDEGs) ------------------------------------------------------------

data_metaDEGs <- read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/2 - Express Analysis/Metanalise - MetaVolcano/data.voltano_all_conditions.xlsx")

data_bar_metaDEGs <- data_metaDEGs %>%
  group_by(Group, diffexpressed) %>%
  summarise(count = n(), .groups = "drop")
data_bar_metaDEGs <- subset(data_bar_metaDEGs, diffexpressed != "Not significant")

result <- data_bar_metaDEGs %>%
  # Contar os valores "Up" e "Down" para cada grupo
  group_by(Group, diffexpressed) %>%
  summarise(total_count = sum(count)) %>%
  spread(key = diffexpressed, value = total_count, fill = 0) %>%  # Criar colunas separadas para 'Up' e 'Down'
  ungroup() %>%
  # Adicionar a coluna 'Names' com o formato desejado (por exemplo: "HBV | N = 2")
  mutate(Names = paste(Group, "| N =",
                       rowSums(select(., starts_with("Down") | starts_with("Up"))))) %>%
  # Adicionar a coluna 'Classification' conforme o tipo de grupo (exemplo: PBMC ou Liver)
  mutate(Classification = case_when(
    grepl("PBMC", Group) ~ "PBMC",
    grepl("Liver", Group) ~ "Liver",
    grepl("Invitro", Group) ~ "Invitro",
    grepl("G2", Group) ~ "TCGA - G2xG1",
    grepl("G3", Group) ~ "TCGA - G3xG1",
    grepl("G4", Group) ~ "TCGA - G4xG1",
    TRUE ~ "Other"
  )) %>%
  select(Names, all_UP = Up, all_DOWN = Down, Classification)  # Renomear as colunas para o formato desejado

result <- result %>%
  mutate(Names = sub(" .*", "", Names))  # Remove tudo após o primeiro espaço

tabela1 <- reshape::melt(as.data.frame(result), id.vars = c("Names", "Classification"), variable_name = "Type")
tabela1 <- tidyr::separate(tabela1, Type, sep = "_", into = c("Type", "UPDOWN"))

tabela1 <- tabela1 %>%
  mutate(Names = ifelse(Names == "HCC", "HCC (TCGA)", Names))
tabela1$Classification <- gsub("^TCGA - ", "", tabela1$Classification)

tabela2 <- tabela1
tabela2[c(5:7, 18:20), 1] <- tabela1[c(5:7, 18:20), 2]
tabela2[c(5:7, 18:20), 2] <- tabela1[c(5:7, 18:20), 1]

# Garantindo que neworder contém apenas os valores presentes em tabela1$Names
neworder <- rev(unique(tabela2$Names))  # Cria a ordem reversa
neworder <- neworder[neworder %in% tabela2$Names]  # Filtra apenas valores existentes em tabela1

# Transformando Names em fator e rearranjando
tabela_final <- dplyr::arrange(
  transform(
    tabela2,
    Names = factor(Names, levels = neworder)  # Define os níveis de acordo com neworder
  ),
  Names  # Reorganiza pela ordem dos fatores
)

# !!! plot for only total DEGs:
## skip to bottom to save
p1 <- ggplot(data = tabela_final) +
  # downregulated total DEGs, and labels
  geom_bar(data = subset(tabela_final, (UPDOWN == "DOWN")),
           position = position_dodge(), stat = "identity",
           aes(fill = paste(UPDOWN), y = -value, x = Names),
           alpha = 1, width = 0.3) +  # Ajuste o valor de width para juntar mais as barras
  geom_text(data = subset(tabela_final, (UPDOWN == "DOWN")),
            aes(y = -value, x = Names, label = value, hjust = 1.2),
            size = 4) +
  # upregulated total DEGs, and labels
  geom_bar(data = subset(tabela_final, (UPDOWN == "UP")),
           position = position_dodge(), stat = "identity",
           aes(fill = paste(UPDOWN), y = value, x = Names),
           alpha = 1, width = 0.3) +  # Ajuste o valor de width para juntar mais as barras
  geom_text(data = subset(tabela_final, (UPDOWN == "UP")),
            aes(y = value, x = Names, label = value, hjust = -0.2),
            size = 4) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),  # Ajuste o tamanho do texto aqui
        axis.text.y = element_text(size = 12, color = "black"),  # Ajuste o tamanho do texto aqui
        panel.background = element_rect(fill = "white", color = "grey"),
        panel.grid = element_blank(),
        strip.text.y = element_text(vjust = 1, hjust = 0.5, size = 12),  # Ajuste o tamanho do texto aqui
        strip.background = element_rect(colour = "grey", fill = "gray95"),
        text = element_text(size = 12)) +  # Ajuste o tamanho do texto aqui
  labs(x = element_blank(), y = element_blank()) +
  coord_flip() +
  facet_grid(rows = vars(Classification), scale = "free_y", space = "free_y") +
  scale_fill_manual(values = c("#93BFCF", "#D71313"),
                    labels = c("Down-regulated DEGs", "Up-regulated DEGs"))

#Salvando tabela
write.xlsx(tabela_final,
           "MetaDEGs_Counts_DivergPlot.xlsx")

ggsave(filename = "DivergiPlotALL.svg",
       plot = p1,
       device = "svg",  # Define o formato do arquivo
       width = 8,  # Largura do gráfico em polegadas
       height = 6,  # Altura do gráfico em polegadas
       units = "in", # Unidade de medida
       dpi = 300)  # qualidade

tabela_final <- read.xlsx("/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/2 - Express Analysis/Metanalise - MetaVolcano/MetaDEGs_Counts_DivergPlot.xlsx")

#BubbleHeatmap
tabela_final <- tabela_final %>%
  dplyr::filter(!Names %in% c("G4xG1", "G3xG1", "G2xG1"))

# Expandir o grid com todas as combinações de Names, Classification e UPDOWN
# tabela_final_expanded <- expand.grid(
#   Names = unique(tabela_final$Names),
#   Classification = unique(tabela_final$Classification),
#   UPDOWN = unique(tabela_final$UPDOWN)
# )
# 
# # Fazer o join com a tabela original para preencher os valores existentes
# tabela_final_expanded <- tabela_final_expanded %>%
#   left_join(tabela_final, by = c("Names", "Classification", "UPDOWN"))
# 
# # tabela_final_expanded <- tabela_final_expanded %>%
# #   mutate(value = ifelse(is.na(value), 0, value))
# 
# tabela_final_expanded$Type <- "all"
# 
# # tabela_final_expanded_up <- tabela_final_expanded %>%
# #   dplyr::filter(!UPDOWN %in% c("DOWN"))
# # 
# # tabela_final_expanded_down <- tabela_final_expanded %>%
# #   dplyr::filter(!UPDOWN %in% c("UP"))
# 
# # Reordenar a coluna Names como fator com níveis personalizados
# tabela_final_expanded <- tabela_final_expanded %>%
#   mutate(
#     Names = as.character(Names),
#     Classification = as.character(Classification),
#     Group = paste(Names, Classification, sep = " ")
#   )
# 
# tabela_final_expanded <- tabela_final_expanded %>%
#   mutate(
#     Classification = as.character(Classification),
#     UPDOWN = as.character(UPDOWN),
#     Group_diff = paste(Classification, UPDOWN, sep = " ")
#   )
# 
# tabela_final_expanded <- tabela_final_expanded %>%
#   mutate(Names = factor(Names, levels = c("HDV", "HCV", "HBV", "HAV", "HEV")))
# 
# tabela_final_expanded <- tabela_final_expanded %>%
#   mutate(Group_diff = factor(Group_diff, levels = c("PBMC UP", "PBMC DOWN", "Liver UP", "Liver DOWN", "Invitro UP", "Invitro DOWN")))
# 
# # Criar o gráfico com o eixo Y ordenado
# buble_heatmap_metadegs <- ggplot(tabela_final_expanded, aes(x = Group_diff, y = Names, size = value, fill = UPDOWN)) +
#   geom_point(shape = 21, color = "black") +  # Criar os pontos com borda preta
#   geom_text(aes(label = round(value, 1)), vjust = -1.4, size = 4) +  # Adicionar labels de "value"
#   # facet_wrap(~UPDOWN) +  # Criar facetas baseadas na variável UPDOWN
#   scale_size_continuous(
#     range = c(2, 10),  # Tamanhos dos pontos (ajuste conforme necessário)
#     name = "Number of meta-DEGs"     # Nome da escala de tamanho
#   ) +
#   scale_fill_manual(
#     values = c("UP" = "#D71313", "DOWN" = "#93BFCF"),  # Especificar cores para UP e DOWN
#     name = "Diff. expression"                         # Nome da escala de preenchimento
#   ) +
#   scale_x_discrete(position = "top") +  # Mover o eixo X para o topo
#   labs(
#     x = "Number of meta-DEGs",
#     y = "Hepatitis"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.y = element_text(size = 10, colour = "black"),
#     axis.title.y = element_text(size = 10, colour = "black"),  # Adicionar vírgula aqui
#     axis.text.x = element_blank(),
#     axis.title.x = element_blank(),
#     strip.text = element_blank(),   
#     plot.title = element_blank(),
#     panel.grid = element_blank(),
#     legend.position = "none",
#     legend.text = element_text(size = 10, colour = "black"),
#     legend.title = element_text(size = 10, colour = "black"),
#   ) + guides(fill = guide_legend(override.aes = list(size = 4)))
# 
# buble_heatmap_metadegs_labels <- ggplot(tabela_final_expanded, aes(x = Group_diff, y = Names, size = value, fill = UPDOWN)) +
#   geom_point(shape = 21, color = "black") +  # Criar os pontos com borda preta
#   geom_text(aes(label = round(value, 1)), vjust = -1.4, size = 4) +  # Adicionar labels de "value"
#   # facet_wrap(~UPDOWN) +  # Criar facetas baseadas na variável UPDOWN
#   scale_size_continuous(
#     range = c(2, 10),  # Tamanhos dos pontos (ajuste conforme necessário)
#     name = "Number of meta-DEGs"     # Nome da escala de tamanho
#   ) +
#   scale_fill_manual(
#     values = c("UP" = "#D71313", "DOWN" = "#93BFCF"),  # Especificar cores para UP e DOWN
#     name = "Diff. expression"                         # Nome da escala de preenchimento
#   ) +
#   scale_x_discrete(position = "top") +  # Mover o eixo X para o topo
#   labs(
#     x = "Number of meta-DEGs",
#     y = "Hepatitis"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.y = element_text(size = 10, colour = "black"),
#     axis.title.y = element_text(size = 10, colour = "black"),  # Adicionar vírgula aqui
#     axis.text.x = element_text(size = 10, colour = "black"),
#     axis.title.x = element_text(size = 10, colour = "black"),
#     strip.text = element_blank(),   
#     plot.title = element_blank(),
#     panel.grid = element_blank(),
#     legend.position = "right",
#     legend.text = element_text(size = 10, colour = "black"),
#     legend.title = element_text(size = 10, colour = "black"),
#   ) + guides(fill = guide_legend(override.aes = list(size = 4)))

tabela_final_barplot <- tabela_final %>%
  mutate(value = ifelse(UPDOWN == "DOWN", -abs(value), value))

# Definir a ordem desejada para o eixo Y
desired_order <- c("HEV", "HDV", "HCV", "HBV", "HAV")

# Reordenar a coluna 'ID' com base na ordem desejada
tabela_final_barplot$Names <- factor(tabela_final_barplot$Names,
                                     levels = desired_order)

metaDEGs_barplot <- ggplot(tabela_final_barplot, aes(x = value, y = Names, fill = UPDOWN)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +  # Barras com espessura consistente
  geom_text(
    aes(
      label = abs(value),
      hjust = ifelse(value > 0, -0.2, 1.2)  # Alinhar dentro ou fora com base no valor
    ),
    size = 2.5, 
    color = "black"
  ) +
  facet_wrap(~ Classification, scales = "free_y", ncol = 1, strip.position = "top") +  # Facetas por Tissue
  scale_fill_manual(
    values = c("UP" = "#D71313", "DOWN" = "#93BFCF"),  # Cores para "Up" e "Down"
    name = "Expression"
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(),  # Quebras automáticas para os valores
    name = "Number of Genes"
  ) +
  labs(
    y = "GSE"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    legend.title = element_text(size = 9, colour = "black"),
    legend.text = element_text(size = 9, colour = "black"),
    strip.text = element_blank(),
    panel.spacing = unit(0.6, "lines"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),  # Linhas no eixo Y
    panel.grid.minor.y = element_blank(),  # Desativa as linhas menores
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed")  # Linhas no eixo X
  ) +
  coord_cartesian(clip = "off")

# ggsave(filename = "MetaDEGs_Bubble.tiff",
#        plot = buble_heatmap_metadegs,
#        device = "tiff",  # Define o formato do arquivo
#        width = 4,  # Largura do gráfico em polegadas
#        height = 3,  # Altura do gráfico em polegadas
#        units = "in", # Unidade de medida
#        dpi = 300)  # qualidade
# 
# ggsave(filename = "MetaDEGs_Bubble_labels.tiff",
#        plot = buble_heatmap_metadegs_labels,
#        device = "tiff",  # Define o formato do arquivo
#        width = 4,  # Largura do gráfico em polegadas
#        height = 3,  # Altura do gráfico em polegadas
#        units = "in", # Unidade de medida
#        dpi = 300)  # qualidade

ggsave(filename = "MetaDEGs_BarPlot.svg",
       plot = metaDEGs_barplot,
       device = "svg",  # Define o formato do arquivo
       width = 3,  # Largura do gráfico em polegadas
       height = 4,  # Altura do gráfico em polegadas
       units = "in", # Unidade de medida
       dpi = 300)  # qualidade

write.xlsx(tabela_final,
           "Suppl.Table1b_barplot.xlsx")

# VolcanoPlot (metaDEGs) -------------------------------------------------------------

data.volcano <- data_metaDEGs

# Ajuste os limites dos eixos conforme necessário
min_value_x <- -4
max_value_x <- 4
min_value_y <- 0
max_value_y <- 15

# Filtrando os dados para garantir que eles estejam dentro dos limites desejados
data_filt <- data.volcano %>%
  dplyr::filter(metafc >= min_value_x & metafc <= max_value_x & -log10(metap) >= min_value_y & -log10(metap) <= max_value_y)

data_filt <- data_filt %>%
  dplyr::filter(!Group %in% c("HCC G2", "HCC G3", "HCC G4"))

#Ordenando
unique(data_filt$Group)
data_filt$Group <- factor(data_filt$Group,
                          levels = c("HAV Invitro", "HBV Invitro", "HCV Invitro", 
                                     "HDV Invitro", "HEV Invitro",
                                     "HBV Liver", "HBV PBMC", "HCV Liver", "HCV PBMC",
                                     "HDV Liver"))

data_filt$Group <- factor(data_filt$Group,
                          levels = c("HBV Liver", "HBV PBMC", "HCV Liver", "HCV PBMC",
                                     "HDV Liver", "HAV Invitro", "HBV Invitro", "HCV Invitro", 
                                     "HDV Invitro", "HEV Invitro"))
data_filt <- data_filt %>%
  mutate(
    Genes = as.character(Genes),
    group = as.character(group),
    cgene = paste(Genes, group, sep = " ")
  )

#Plotando
p2 <- ggplot(data = data_filt,
            aes(x = metafc,
                y = -log10(metap),
                col = diffexpressed,
                label = delabel)) +
  geom_point(size = 1) +
  # geom_text_repel(size = 3, color = "black") +
  facet_wrap(~Group, ncol = 5, scales = "free_x") +
  theme_bw(12) +
  scale_color_manual(values = c("#93BFCF", "#D8D8D8", "#D71313")) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 14),  # Ajuste o tamanho do texto no eixo X
        axis.text.y = element_text(size = 14),  # Ajuste o tamanho do texto no eixo Y
        strip.text = element_text(size = 14),   # Ajuste o tamanho do texto do facet_wrap
        # legend.text = element_text(size = 12),  # Aumenta o tamanho do texto da legenda
        # legend.title = element_text(size = 12), # Aumenta o tamanho do título da legenda
        legend.position = "none") +           # Muda a posição da legenda
  xlim(c(min_value_x, max_value_x)) +
  ylim(c(min_value_y, max_value_y)) +
  labs(x = "Log2 Fold Change", y = "-Log10 p-value")

#Salvando plot
ggsave(filename = "Volcano_All_2.tiff",
       plot = p2,
       device = "tiff",  # Define o formato do arquivo
       width = 10,  # Largura do gráfico em polegadas
       height = 6,  # Altura do gráfico em polegadas
       units = "in", # Unidade de medida
       dpi = 300)  # qualidade

# ShinyGO ------------------------------------------------------------------
setwd("/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/2 - Express Analysis/Metanalise - MetaVolcano/4 - ShinyGo - Genes caracterization/")

data_all <- read.xlsx("data_all_ShinyGO.xlsx")

data_all$Type <- factor(data_all$Type, levels = data_all %>%
                          group_by(Type) %>%
                          summarise(total_count = sum(count)) %>%
                          arrange(desc(total_count)) %>%
                          pull(Type))

data_all <- data_all %>%
  mutate(
    group = factor(group, levels = c(
      "HBV Liver", "HBV PBMC", "HCV Liver", "HCV PBMC", "HDV Liver",
      "HAV in vitro", "HBV in vitro", "HCV in vitro", "HDV in vitro", "HEV in vitro"))
  )

#Plotando
p3 <- ggplot(data_all, aes(x = Type, y = count, fill = differexpressed)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +  # Ajustando a espessura das barras para 0.7
  facet_wrap(~group, scales = "free", nrow = 2) +  # Facetas com a ordem definida
  scale_fill_manual(values = c("up" = "#D71313", "down" = "#93BFCF")) +  # Cores ajustadas manualmente
  theme_bw(base_line_size = 0) +
  theme(strip.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        # axis.line = element_line(linewidth = 0.5, colour = "black"),  # Linha dos eixos
        legend.position = "none"
  ) +
  labs(
    title = "Gene characterization by condition",
    x = "",
    y = ""
  )

# ggsave(filename = "Gene_characterization.svg",
#        plot = p1,
#        device = "svg", 
#        width = 16, 
#        height = 8, 
#        units = "in", 
#        dpi = 300)

combined_plot2 <- p2/p3

ggsave(filename = "Gene_characterization.tiff",
       plot = combined_plot2,
       device = "tiff", 
       width = 16, 
       height = 11, 
       units = "in", 
       dpi = 200)

ggsave(filename = "Gene_characterization.svg",
       plot = combined_plot2,
       device = "svg", 
       width = 16, 
       height = 14, 
       units = "in", 
       dpi = 300)

# Appyter Characterization (AppyterFinalCharacterization.RData)
load("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/2 - Appyter/AppyterFinalCharacterization.RData")
#--------------------- ScatterPlot ----
# import data
data_1 = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoHBV_liver.csv")
data_2 = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoHBV_pbmc.csv")
data_3 = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoHCV_liver.csv")
data_4 = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoHCV_pbmc.csv")
data_5 = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoHDV.csv")

#Referencia
data.cluster = read.xlsx("~/Documents/DataSets/MetaVolcano - ArrayxBulk/3 - Metanalise HBCDV - PBMC & Liver/Enrichment/HBV_Liver/ScatterPlot/HBV_liver_pathways.xlsx")

#Criando tabela com identificacao de GO e deletando de "Terms"
data.cluster <- data.cluster %>%
  mutate(
    GO = str_extract(term, "(?<=\\().*(?=\\))"))
data.cluster$term <- sub("\\s*\\(.*?\\)\\s*", "", data.cluster$term)
data_1 <- data_1 %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_1$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_1$Term)
data_2 <- data_2 %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_2$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_2$Term)
data_3 <- data_3 %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_3$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_3$Term)
data_4 <- data_4 %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_4$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_4$Term)
data_5 <- data_5 %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_5$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_5$Term)

#Unindo com referencial
data_1 <- left_join(data.cluster, data_1,
                    by = "GO", suffix = c(".cluster",
                                          "_1"),
                    keep = TRUE)
data_2 <- left_join(data.cluster, data_2,
                    by = "GO", suffix = c(".cluster",
                                          "_2"),
                    keep = TRUE)
data_3 <- left_join(data.cluster, data_3,
                    by = "GO", suffix = c(".cluster",
                                          "_3"),
                    keep = TRUE)
data_4 <- left_join(data.cluster, data_4,
                    by = "GO", suffix = c(".cluster",
                                          "_4"),
                    keep = TRUE)
data_5 <- left_join(data.cluster, data_5,
                    by = "GO", suffix = c(".cluster",
                                          "_5"),
                    keep = TRUE)

#Inserindo P-valores nao significantes
data_1$P.value[is.na(data_1$P.value)] <- 1
data_2$P.value[is.na(data_2$P.value)] <- 1
data_3$P.value[is.na(data_3$P.value)] <- 1
data_4$P.value[is.na(data_4$P.value)] <- 1
data_5$P.value[is.na(data_5$P.value)] <- 1

#Ajeitando tabela pro Plot
data_1 = data_1[,c(2:4,6,14,16)] 
data_2 = data_2[,c(2:4,6,14,16)] 
data_3 = data_3[,c(2:4,6,14,16)] 
data_4 = data_4[,c(2:4,6,14,16)] 
data_5 = data_5[,c(2:4,6,14,16)] 
View(data_1)
# grouping clusters by total number
data_1 = data_1 |>
  group_by(cluster) |>
  mutate(total_number = n())
data_2 = data_2 |>
  group_by(cluster) |>
  mutate(total_number = n())
data_3 = data_3 |>
  group_by(cluster) |>
  mutate(total_number = n())
data_4 = data_4 |>
  group_by(cluster) |>
  mutate(total_number = n())
data_5 = data_5 |>
  group_by(cluster) |>
  mutate(total_number = n()) 

# picking significant pvalues
data_1$sig = ifelse(data_1$P.value < 0.05, 1, 0)
data_2$sig = ifelse(data_2$P.value < 0.05, 1, 0)
data_3$sig = ifelse(data_3$P.value < 0.05, 1, 0)
data_4$sig = ifelse(data_4$P.value < 0.05, 1, 0)
data_5$sig = ifelse(data_5$P.value < 0.05, 1, 0)

#Contabilizando Significantes
data_1 = data_1 |>
  group_by(cluster) |>
  mutate(filtered = sum(sig))
data_2 = data_2 |>
  group_by(cluster) |>
  mutate(filtered = sum(sig))
data_3 = data_3 |>
  group_by(cluster) |>
  mutate(filtered = sum(sig))
data_4 = data_4 |>
  group_by(cluster) |>
  mutate(filtered = sum(sig))
data_5 = data_5 |>
  group_by(cluster) |>
  mutate(filtered = sum(sig))

# Defining groups names in data
data_1$group = 'HBV Liver'
data_2$group = 'HBV PBMC'
data_3$group = 'HCV Liver'
data_4$group = 'HCV PBMC'
data_5$group = 'HDV Liver'

# binding datas
data = rbind(data_1,
             data_2,
             data_3,
             data_4,
             data_5)
head(data)

clusters = paste("Cluster", 0:22, sep = " ")
data$cluster = factor(data$cluster, levels = clusters)
mid = mean(data$total_number)
data_all = subset(data, filtered >= 10)
data_all = data_all |>
  group_by(cluster, group) |>
  summarise(x_mean = mean(x),
            y_mean = mean(y))
data$alpha = ifelse(data$filtered < 5, 0.1, 1)

#Separando Neuro e Imune
neuro_keywords <- c('neuron', 'neuro', 'synapse', 'brain', 'nervous', "Anterograde",
                    "neural", "synaptic", "telencephalon", "neuronal", 'Neuroblast',
                    "axon", "neurotransmitter", "serotonin", "neurogenesis", 
                    "dopamine", "neuropeptide", "oligodendrocyte", 'Neurotransmitter',
                    "nerve impulse", "limbic", "Dendritic", 'Glial', 'Glutamate',
                    'Neuroinflammatory', 'Neuropeptide', 'Nerve', 'Sodium Ion Homeostasis',
                    'Behavior', 'Potassium Ion Homeostasis',
                    'Potassium Ion Import Across Plasma Membrane', 'Glutamine Family Amino Acid Catabolic Process',
                    'Glutamine Family Amino Acid Biosynthetic Process', 'Dendrite Morphogenesis',
                    'acetyl-CoA Metabolic Process', 'Sodium-Dependent Phosphate Transport',
                    'Nitric Oxide Biosynthetic Process', 'Cellular Glucuronidation', "Synaptic", "neuro", "Synapse",
                    "Perisynaptic", "Postsynaptic",
                    "Neurotransmitter", "Neuroinflammatory",
                    "Neuron", "Postsynapse", "Presynapse",
                    "Presynaptic", "Nervous", "Dopamine", "Axon",
                    "Neurogenesis", "Adrenergic", "Dendritic", "Neural",
                    "Neuromuscular", "Neuroblast", "Sensory", "Gliogenesis",
                    "Learning", "Brain")
immune_keywords <- c('immune', 'immunoglobulin', 'lymphocyte', 'cytokine', 
                     'B cell', 'Apoptotic', 'Leukocyte', 'Calcium', 'Inflammatory', 
                     'interleukin', 'T cell', 'monocyte', 'T-helper', 'interferon',
                     'chemokine', 'MHC', 'Neutrophil', 'Antigen', 'Regulation', 'Defense',
                     'ERK1', 'ERK2', 'MAPK', 'Natural Killer', 'Macrophage', 'Kinase',
                     'Humoral', 'Antibacterial', 'Modulation By Host Of Viral Process',
                     'Cellular Response To Tumor Necrosis Factor', 'Eosinophil Migration',
                     'Pathway-Restricted SMAD Protein Phosphorylation', 'Response To BMP',
                     'Endocytosis', 'Toll-Like Receptor 9 Signaling Pathway', 'Response To Tumor Necrosis Factor',
                     'Cellular Response To Decreased Oxygen Levels', 'TORC2 Signaling',
                     'Response To cGMP', 'Cellular Response To Hypoxia', 'Pathway-Restricted SMAD Protein Phosphorylation',
                     'Response To BMP', 'Mitotic Sister Chromatid Segregation', 'DNA Replication Checkpoint Signaling',
                     'G1/S Transition Of Mitotic Cell Cycle', 'G2/M Transition Of Mitotic Cell Cycle',
                     'DNA Repair', 'Acute-Phase Response', 'Response To Cadmium Ion',
                     'Response To Copper Ion', 'Cellular Response To Metal Ion', 'Cellular Response To Cadmium Ion',
                     'Cellular Response To Copper Ion', 'Cellular Response To Zinc Ion', 'Response To Fibroblast Growth Factor',"leukocyte", "T cell", "mononuclear",
                     "lymphocyte", "immune response", "B cell",
                     "interferon", "Inflammatory", "Immune", "Defense",
                     "Cytokine", "Kinase", "ERBB2", "G Protein", 
                     "Cysteine", "ERK1 And ERK2", "Interleukin",
                     "Macrophage", "Necrosis", "Tumor", "MAPK",
                     "Healing", "Humoral", "Neutrophil", "Leukocyte",
                     "Myeloid", "Response", "Negative Regulation Of Type II Interferon-Mediated Signaling Pathway",
                     "Positive Regulation Of Growth", "MHC", "Cell Growth",
                     "Toll-Like Receptor", "Monocyte", "Interferon", "STAT",
                     "Prostaglandin", "NF-kappaB","GTPase", "Apoptotic", "Viral",
                     "Modulation", "Immunoglobulin", "NIK/NF-kappaB", "B Cell",
                     "T Cell", "I-kappaB kinase/NF-kappaB", "kinase", "Death",
                     "Ubiquitin", "Ubiquitin-Dependent","Ubiquitination",
                     "Proliferation", "Growth", "Mast", "Regulation")
data <- data %>%
  mutate(pathway = case_when(
    str_detect(term, paste(neuro_keywords,
                           collapse = "|")) ~ "Nervous",
    str_detect(term, paste(immune_keywords,
                           collapse = "|")) ~ "Immune",
    TRUE ~ "Others"
  ))

data <- data %>%
  mutate(Cores = case_when(
    pathway == "Immune" ~ "#3182bd",
    pathway == "Nervous" ~ "tomato",
    pathway == "Others" ~ "#93BFCF",
    TRUE ~ NA_character_  # Se nenhum dos casos acima for atendido, definir como NA
  ))

data <- data %>%
  mutate(Shapes = case_when(
    pathway == "Immune" ~ 19,
    pathway == "Nervous" ~ 19,
    pathway == "Others" ~ 8,
    TRUE ~ NA_real_  # Definir como NA do tipo numeric
  ))

# Ordenando a variável group de acordo com Brain e PBMC
data <- data %>%
  mutate(group = factor(group, levels = c("HBV Liver",
                                          "HBV PBMC",
                                          "HCV Liver",
                                          "HCV PBMC",
                                          "HDV Liver")))
#Scatterplot
# ggplot(data, aes(x, y)) + 
#   geom_point(colour = data$Cores, 
#              shape = data$Shapes,
#              size = data$alpha,
#              aes(alpha = alpha)) + 
#   geom_label_repel(data = data_all, 
#                    mapping = aes(x = x_mean, y = y_mean, 
#                                  label = paste(cluster)), 
#                    colour = 'black', 
#                    size = 4, 
#                    nudge_x = -1,  # Ajuste manual de deslocamento horizontal
#                    nudge_y = 1.5) +  # Ajuste manual de deslocamento vertical
#   facet_wrap(~ group, nrow = 1, dir = "h") +  
#   # scale_x_continuous(breaks = scales::pretty_breaks(n = 4.5), labels = scales::number_format(accuracy = 1)) +  # Ajuste para números inteiros no eixo X
#   theme_bw() +  # Tema branco com tamanho de base ajustado
#   theme(axis.text.y = element_text(size = 12, colour = "black"),  # Ajusta o tamanho do texto dos eixos
#         axis.text.x = element_text(size = 12, colour = "black"),
#         strip.text = element_text(size = 12, colour = "black")) +  # Ajusta o tamanho do texto do facet_wrap
#   xlab("UMAP1") + 
#   ylab("UMAP2")

p1 <- ggplot(data, aes(x, y)) + 
  geom_point(colour = data$Cores, 
             shape = data$Shapes,
             aes( 
               alpha = ifelse(data$pathway == "Others", 0.01, data$alpha))) + 
  geom_label_repel(
    data = data_all, 
    mapping = aes(x = x_mean, y = y_mean, label = paste(cluster)), 
    colour = 'black', 
    size = 6, 
    nudge_x = -1,  # Ajuste manual de deslocamento horizontal
    nudge_y = 1.5  # Ajuste manual de deslocamento vertical
  ) +  
  facet_wrap(~ group, nrow = 1, dir = "h") +  
  labs(
    x = "UMAP1",  # Rótulo do eixo X
    y = "UMAP2",  # Rótulo do eixo Y
    colour = "Cluster Colour",  # Nome personalizado para a legenda de cor
    alpha = "Alpha"  # Nome personalizado para a legenda de alpha
  ) + 
  theme_bw() +  # Tema branco com tamanho de base ajustado
  theme(legend.position = "right",
        legend.text = element_text(size = 24, colour = "black"),
        legend.title = element_text(size = 24, colour = "black"),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 24, colour = "black"),
        strip.text = element_text(size = 24, colour = "black"),
        panel.spacing = unit(1.5, "lines")
  )

ggsave("sua_imagem.svg", plot = p1,
       width = 20, height = 5, units = "in")

ggsave("scatterplot.tiff", plot = p1,
       width = 20, height = 5, units = "in")

#Saving Tables
#Todos os Terms significantes
write.xlsx(data, "ScatterPlot_Final_Characterization.xlsx")

#--------------------- Grafico de Linhas ------------------
luvial = data[,c(3,5,6,10,12)]
# Selecionando valores menores que 0.05 na coluna P.value
luvial_filtered <- luvial[luvial$P.value < 0.05, ]
#Contabilizando
#Frequencia de system (Neuro, Immune e Others) por cluster
luvial_line <- luvial_filtered %>%
  group_by(pathway, cluster, group) %>%
  mutate(TermSys = ifelse(P.value < 0.05, 1, 0)) %>%
  mutate(TermSys = sum(TermSys))
View(luvial_line)
#Completando vales faltantes
all_combinations <- expand.grid(cluster = unique(luvial_line$cluster),
                                pathway = unique(luvial_line$pathway),
                                Term = unique(luvial_line$Term),
                                group = unique(luvial_line$group))
# Combinar com os dados existentes
merged_data <- merge(all_combinations, luvial_line,
                     by = c("cluster", "pathway", "Term", "group"),
                     all.x = TRUE)
# Preencher os valores ausentes com 0 em TermSys
merged_data[is.na(merged_data$TermSys), "TermSys"] <- 0
# Remover a coluna Size e Order, pois elas não são mais necessárias
# merged_data <- merged_data[, !(names(merged_data) %in%
#                                  c("Size", "Order"))]
merged_data_unique <- merged_data %>%
  distinct(Term, pathway, cluster, group, .keep_all = TRUE)

#Para o grafico
View(merged_data_unique)
line_luvial = merged_data_unique[,c(1,2,4,6)]
colnames(line_luvial)
View(line_luvial)
# Remover "Cluster " da coluna cluster
line_luvial$cluster <- gsub("Cluster ", "", line_luvial$cluster)
# Ordenando a variável group de acordo com grupos
unique(line_luvial$group)
line_luvial <- line_luvial %>%
  mutate(group = factor(group, levels = c("HBV Liver", "HBV PBMC",
                                          "HCV Liver", "HCV PBMC",
                                          "HDV Liver")))
# Convertendo para números
line_luvial$cluster <- as.numeric(line_luvial$cluster)
line_luvial <- line_luvial %>%
  mutate_all(~replace_na(., 0))
line_luvial <- line_luvial %>%
  mutate(cluster = factor(cluster, levels = 22:0))
View(line_luvial)
# Plotando
tiff("/Users/adrielnobile/LineSystem.tiff",
     width = 14,
     height = 6,
     res = 300, units = 'in') 
p2 <- ggplot(line_luvial, aes(y = cluster,
                              x = TermSys,
                              color = pathway)) +
  geom_ribbon(aes(xmin = 0, xmax = TermSys), alpha = 0.3) +  # Adicionando preenchimento
  geom_point(data = subset(line_luvial, TermSys != 0),
             aes(shape = pathway), size = 3) +
  scale_color_manual(values = c("Nervous" = "tomato",
                                "Immune" = "#3182bd",
                                "Others" = "#93BFCF")) +
  theme_bw() +  # Ajusta o tamanho da fonte global
  facet_grid(~group, scales = "free_y") +  # Definindo as facetas com escala livre no eixo y
  theme(axis.text.y = element_text(size = 24, colour = "black"),  # Ajusta o tamanho do texto no eixo y
        axis.text.x = element_text(size = 24, colour = "black"),  # Ajusta o tamanho do texto no eixo x
        axis.title.x = element_text(size = 24, colour = "black"),
        axis.title.y = element_text(size = 24, colour = "black"),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_blank(),  # Ajusta o tamanho do texto no facet_wrap
        legend.position = "right",
        legend.text = element_text(size = 24, colour = "black"),
        legend.title = element_text(size = 24, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Number of Biological Process (p-value < 0.05)",
       y = "Cluster",
       colour = "System", shape = "System") +
  scale_y_discrete(breaks = 0:22, labels = 0:22) +
  geom_vline(xintercept = 0, linetype = "twodash", color = "black") +
  geom_point(data = subset(line_luvial, TermSys == 0), shape = NA)

dev.off()
#Salvando Tabela
write.xlsx(line_luvial, "Appyter_Lineplot_Final.xlsx")

ggsave("line_plot.tiff", plot = p2,
       width = 20, height = 8, units = "in")

#--------------------  BarPlot.   ------------------------
bar_data <- luvial[luvial$P.value < 0.05, ]
View(bar_data)
#Contabilizando
bar_data = bar_data[,c(4:5)]
bar_data <- bar_data %>%
  group_by(pathway, group) %>%
  mutate(TermSys = n())
bar_data_unique <- bar_data %>%
  distinct(group, pathway, TermSys)
#Definindo Cores
cols = c("Nervous" = "tomato",
         "Immune" = "#3182bd",
         "Others" = "#93BFCF")
# Ordenando a variável group de acordo com Brain e PBMC
bar_data_unique <- bar_data_unique %>%
  mutate(group = factor(group, levels = c("HBV Liver", "HBV PBMC",
                                          "HCV Liver", "HCV PBMC",
                                          "HDV Liver")))
#Plotando
tiff("/Users/adrielnobile/BarPlotSystem.tiff",
     width = 30,
     height = 12,
     res = 100, units = 'in') 
p3 <- ggplot(bar_data_unique, aes(x = pathway,
                                  y = TermSys,
                                  fill = pathway,
                                  label = TermSys)) +
  geom_bar(stat = "identity", position = "dodge") +  # Cria o gráfico de barras
  geom_text(position = position_dodge(width = 0.9), vjust = 0.5, size = 8, color = "black") +  # Texto dentro das barras
  facet_wrap(~ group, nrow = 1, scales = "free_x") +  # Facetado por 'group' com escalas livres no eixo x
  scale_fill_manual(values = cols) +  # Cores personalizadas para o preenchimento
  theme_bw(base_size = 18) +  # Tema base com fonte maior
  labs(x = "Pathway", y = "Number of BP", fill = "Group") +  # Rótulos dos eixos e legenda
  theme(
    axis.text.x = element_blank(),  # Oculta o texto do eixo X
    axis.text.y = element_text(size = 22, colour = "black"),  # Ajusta tamanho e cor do texto do eixo Y
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, colour = "black"),
    strip.text = element_blank(),  # Remove os rótulos das facetas
    legend.position = "right",
    legend.text = element_text(size = 22, colour = "black"),
    legend.title = element_text(size = 22, colour = "black"),
    panel.grid.major = element_blank(),  # Remove as linhas de grade principais
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.5, "lines")
  )

# write.xlsx(bar_data_unique, "Appyter_BarPlot_Final.xlsx")

ggsave("barplot_plot.tiff", plot = p3,
       width = 20, height = 3, units = "in")

ggsave("barplot_plot.svg", plot = p3,
       width = 20, height = 3, units = "in")


# rrvGO -------------------------------------------------------------------
{
  library(rrvgo)
  library(stringr)
  library(shinydashboard)
  library(heatmaply)
  library(seriation)
  library(wordcloud)
  library(grid)
}

# import data
data_bl = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoHBV_liver.csv")
data_bp = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoHBV_pbmc.csv")
data_cl = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoHCV_liver.csv")
data_cp = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoHCV_pbmc.csv")
data_dl = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoHDV.csv")

data_liver <- rbind(data_bl, data_cl, data_dl)

data_liver <- data_liver %>%
  mutate(
    GOID = str_extract(Term, "(?<=\\()(GO:\\d+)(?=\\))")
  )

data_pbmc <- rbind(data_bp, data_cp)

data_pbmc <- data_pbmc %>%
  mutate(
    GOID = str_extract(Term, "(?<=\\()(GO:\\d+)(?=\\))")
  )

#LIVER
#Atribuir os IDs para simMatrix
simMatrix <- calculateSimMatrix(data_liver$GOID, #ele usa os IDs dos terms para fazer a análise
                                orgdb="org.Hs.eg.db",
                                ont="BP", #Selecione a Ontology
                                method="Rel")

#ele usa o -log10(qvalue) como base para calcular o score de cada term 
scores <- setNames(-log10(data_liver$Adjusted.P.value), data_liver$GOID) #Clusterizando
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

# Paleta de cores para o heatmap
cores_escala <- c("#F1F0E8", "#B3C8CF", "#89A8B2")
escala_cores <- colorRampPalette(cores_escala)
num_cores <- 100
heatmap_colors <- escala_cores(num_cores)

# Anotação com 'parent' e configuração das cores para cada categoria
annotationLabel <- "parentTerm"
ann <- data.frame(
  parent = factor(reducedTerms[match(rownames(simMatrix), reducedTerms$go), annotationLabel]), 
  row.names = rownames(simMatrix)
)

# Paleta de cores personalizada para a anotação
base_colors <- c("#493628", "#BC7C7C", "#867070", "#F6EFBD", "#A5B68D", "#9BB0C1")
annColors <- list(parent = colorRampPalette(base_colors)(length(levels(ann$parent)))) # Ajusta para 33 cores

# Nomeia as cores da paleta com os níveis da variável `parent`
names(annColors$parent) <- levels(ann$parent)

p7 <- pheatmap::pheatmap(
  simMatrix,                     # Matriz de similaridade
  annotation_row = ann,          # Anotação das linhas
  annotation_colors = annColors, # Cores personalizadas para a anotação
  annotation_legend = FALSE,
  col = heatmap_colors,          # Cores do heatmap
  show_rownames = FALSE, 
  show_colnames = FALSE,         # Não mostrar nomes das linhas e colunas
  clustering_distance_cols = "euclidean",
  clustering_method = "complete", # Método de clusterização
  fontsize = 6                   # Tamanho da fonte
)

p7

ggsave(filename = "rrvgo_LIVER_heatmap.tiff",
       plot = p7,
       device = "tiff", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300)  

p8 <- pheatmap::pheatmap(
  simMatrix,                     # Matriz de similaridade
  annotation_row = ann,          # Anotação das linhas
  annotation_colors = annColors, # Cores personalizadas para a anotação
  annotation_legend = TRUE,
  col = heatmap_colors,          # Cores do heatmap
  show_rownames = FALSE, 
  show_colnames = FALSE,         # Não mostrar nomes das linhas e colunas
  clustering_distance_cols = "euclidean",
  clustering_method = "complete", # Método de clusterização
  fontsize = 16                   # Tamanho da fonte
)

ggsave(filename = "rrvgo_LIVER_heatmap_labels.tiff",
       plot = p8,
       device = "tiff", 
       width = 12, 
       height = 23, 
       units = "in", 
       dpi = 300)  

#PBMC
simMatrix <- calculateSimMatrix(data_pbmc$GOID, #ele usa os IDs dos terms para fazer a análise
                                orgdb="org.Hs.eg.db",
                                ont="BP", #Selecione a Ontology
                                method="Rel")

#ele usa o -log10(qvalue) como base para calcular o score de cada term 
scores <- setNames(-log10(data_pbmc$Adjusted.P.value), data_pbmc$GOID) #Clusterizando
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

# Paleta de cores para o heatmap
cores_escala <- c("#F1F0E8", "#B3C8CF", "#89A8B2")
escala_cores <- colorRampPalette(cores_escala)
num_cores <- 100
heatmap_colors <- escala_cores(num_cores)

# Anotação com 'parent' e configuração das cores para cada categoria
annotationLabel <- "parentTerm"
ann <- data.frame(
  parent = factor(reducedTerms[match(rownames(simMatrix), reducedTerms$go), annotationLabel]), 
  row.names = rownames(simMatrix)
)

# Paleta de cores personalizada para a anotação
base_colors <- c("#493628", "#BC7C7C", "#867070", "#F6EFBD", "#A5B68D", "#9BB0C1")
annColors <- list(parent = colorRampPalette(base_colors)(length(levels(ann$parent)))) # Ajusta para 33 cores

# Nomeia as cores da paleta com os níveis da variável `parent`
names(annColors$parent) <- levels(ann$parent)

p9 <- pheatmap::pheatmap(
  simMatrix,                     # Matriz de similaridade
  annotation_row = ann,          # Anotação das linhas
  annotation_colors = annColors, # Cores personalizadas para a anotação
  annotation_legend = FALSE,
  col = heatmap_colors,          # Cores do heatmap
  show_rownames = FALSE, 
  show_colnames = FALSE,         # Não mostrar nomes das linhas e colunas
  clustering_distance_cols = "euclidean",
  clustering_method = "complete", # Método de clusterização
  fontsize = 6                   # Tamanho da fonte
)

ggsave(filename = "rrvgo_PBMC_heatmap.tiff",
       plot = p9,
       device = "tiff", 
       width = 8, 
       height = 8, 
       units = "in", 
       dpi = 300)  

p10 <- pheatmap::pheatmap(
  simMatrix,                     # Matriz de similaridade
  annotation_row = ann,          # Anotação das linhas
  annotation_colors = annColors, # Cores personalizadas para a anotação
  annotation_legend = TRUE,
  col = heatmap_colors,          # Cores do heatmap
  show_rownames = FALSE, 
  show_colnames = FALSE,         # Não mostrar nomes das linhas e colunas
  clustering_distance_cols = "euclidean",
  clustering_method = "complete", # Método de clusterização
  fontsize = 16                   # Tamanho da fonte
)

ggsave(filename = "rrvgo_PBMC_heatmap_labels.tiff",
       plot = p10,
       device = "tiff", 
       width = 12, 
       height = 12, 
       units = "in", 
       dpi = 300)  

# SynGO -------------------------------------------------------------------
setwd("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/5 - SynGO/")

data_syn_bl <- read.delim("HBV_Liver/HBV_Liver/HBV_Liver_Syn.txt")
data_syn_bp <- read.delim("HBV_PBMC/HBV_PBMC/HBV_PBMC_syn.txt")
data_syn_cl <- read.delim("HCV_Liver/HCV_Liver/HCV_Liver_syn.txt")
data_syn_cp <- read.delim("HCV_PBMC/HCV_PBMC/HCV_PBMC_syn.txt")
data_syn_dl <- read.delim("HDV_Liver/HDV_Liver/HDV_Liver_syn.txt")

data_syn_bl$group <- "HBV Liver"
data_syn_bp$group <- "HBV PBMC"
data_syn_cl$group <- "HCV Liver"
data_syn_cp$group <- "HCV PBMC"
data_syn_dl$group <- "HDV Liver"

data_syn <- rbind(data_syn_bl, data_syn_bp, data_syn_cl, data_syn_cp, data_syn_dl)

data_syn <- data_syn %>%
  mutate(
    GOID = str_extract(Term, "\\(GO:\\d+\\)"),
    Term_clean = str_remove(Term, "\\(GO:\\d+\\).*"),
    Ontology = str_extract(Term, "(?<=\\) ).*")
  ) %>%
  mutate(
    GOID = str_remove_all(GOID, "[()]"), # Remover os parênteses
    Term = str_trim(Term_clean)         # Substituir diretamente a coluna Term
  ) %>%
  dplyr::select(-Term_clean) # Usar dplyr::select para evitar conflitos

data_syn_filt <- data_syn %>%
  dplyr::filter(Ontology == "BP")

data_syn_filt <- data_syn_filt %>%
  group_by(group) %>%
  slice_min(order_by = Adjusted.P.value, n = 5, with_ties = FALSE) %>%  # Seleciona os 5 menores valores por grupo
  ungroup()

# p4 <- ggplot(data_syn_filt, aes(x = Term, y = Overlap,
#                                 size = as.numeric(str_extract(Overlap, "^\\d+")), color = Adjusted.P.value)) +
#   geom_point(alpha = 0.8) +
#   facet_wrap(~group, scales = "free", nrow = 1) +
#   scale_color_gradient(low ="#D71313", high = "#D8D8D8", trans = "log10") +
#   scale_size_continuous(name = "Overlap", range = c(2, 10)) + # Ajusta o tamanho dos pontos
#   labs(
#     title = "",
#     x = "Terms",
#     y = "Overlap",
#     color = "Adjusted p-value"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(size = 24, angle = 45, hjust = 1, color = "black"),
#     axis.title.x = element_text(size = 24, color = "black"),
#     axis.text.y = element_text(size = 24, colour = "black"),
#     axis.title.y = element_text(size = 24, color = "black"),
#     plot.title = element_text(hjust = 0.5),
#     legend.text = element_text(size = 24),
#     legend.title = element_text(size = 24),
#     strip.text = element_text(size = 24, colour = "black"),
#     panel.spacing = unit(1.5, "lines")
#   )

data_syn_filt <- data_syn_filt %>%
  mutate(Reduced_Term = str_replace_all(Term, 
                                        c("Regulation Of " = "Reg. ",
                                          "Synaptic " = "Syn. ",
                                          "Postsynapse " = "Postsyn. ",
                                          "Transmission" = "Transm.",
                                          "Receptor" = "Rec.",
                                          "Membrane" = "Mem.",
                                          "Activity" = "Act.",
                                          "Specialization" = "Spec.",
                                          "Process" = "Proc.",
                                          "Involved In" = "In.",
                                          "Modulating" = "Mod.",
                                          "Modulation" = "Mod.",
                                          "Between" = "",
                                          "Vesicle" = "Ves.",
                                          "Presynaptic" = "Presyn. ",
                                          "Active" = "Act",
                                          "Postsynaptic" = "Postsyn.",
                                          "Potential" = "Pot."))
  )

# Reordenar a variável 'Overlap' com base em 'Adjusted.P.value'
data_syn_filt <- data_syn_filt %>%
  mutate(Overlap_reordered = reorder(Overlap, Adjusted.P.value))

# Criar o gráfico
p4 <- ggplot(data_syn_filt, aes(x = reorder(Reduced_Term, Adjusted.P.value),
                                y = Overlap,
                                size = -log10(Adjusted.P.value),
                                color = Adjusted.P.value)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~group, scales = "free", nrow = 1) +
  scale_color_gradient(low ="#D71313", high = "#D8D8D8", trans = "log10") +
  scale_size_continuous(name = "-log10(Adj.P.value)", range = c(2, 10)) + # Ajusta o tamanho dos pontos
  labs(
    title = "",
    x = "Terms",
    y = "Overlap",
    color = "Adj. p-value"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 24, angle = 45, hjust = 1, color = "black"),
    axis.title.x = element_text(size = 24, color = "black"),
    axis.text.y = element_text(size = 26, colour = "black"),
    axis.title.y = element_text(size = 26, color = "black"),
    plot.title = element_blank(),
    legend.text = element_blank(),
    legend.title = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

ggsave("synapseGO.svg", plot = p4,
       width = 25, height = 10, units = "in")

#Network (SynGO)
data_net <- data_syn

data_net <- data_net %>%
  separate_rows(Genes, sep = ";")

data_net <- data_net %>%
  filter(Ontology == "BP")

data_net <- data_net %>%
  mutate(Reduced_Term = str_replace_all(Term, 
                                        c("Regulation Of " = "Reg. ",
                                          "Synaptic " = "Syn. ",
                                          "Postsynapse " = "Postsyn. ",
                                          "Transmission" = "Transm.",
                                          "Receptor" = "Rec.",
                                          "Membrane" = "Mem.",
                                          "Activity" = "Act.",
                                          "Specialization" = "Spec.",
                                          "Process" = "Proc.",
                                          "Involved In" = "In.",
                                          "Modulating" = "Mod.",
                                          "Modulation" = "Mod.",
                                          "Between" = "",
                                          "Vesicle" = "Ves.",
                                          "Presynaptic" = "Presyn. ",
                                          "Active" = "Act",
                                          "Postsynaptic" = "Postsyn.",
                                          "Potential" = "Pot."))
  )

data_net[1] <- NULL
names(data_net)[12] <- "Term"

#Concatenando
data_net <- data_net %>%
  mutate(
    Genes = as.character(Genes),
    group = as.character(group),
    cgene = paste(Genes, group, sep = " ")
  )

data_net <- data_net %>%
  mutate(
    Term = as.character(Term),
    group = as.character(group),
    cterm = paste(Term, group, sep = " ")
  )

names(data_net)[1] <- "Term_ID"
names(data_net)[8] <- "Genes_ID"
names(data_net)[13] <- "Genes"
names(data_net)[14] <- "Term"

# data_net_bl <- data_net %>%
#   filter(group == "HBV Liver")
# 
# data_net_bp <- data_net %>%
#   filter(group == "HBV PBMC")

data_net_b <- data_net %>%
  filter(group %in% c("HBV Liver", "HBV PBMC"))

# data_net_cl <- data_net %>%
#   filter(group == "HCV Liver")
# 
# data_net_cp <- data_net %>%
#   filter(group == "HCV PBMC")

data_net_c <- data_net %>%
  filter(group %in% c("HCV Liver", "HCV PBMC"))

data_net_dl <- data_net %>%
  filter(group == "HDV Liver")

# net <- network(data_net_b[, c("Term", "Genes")], directed = FALSE) #HBV
# net <- network(data_net_c[, c("Term", "Genes")], directed = FALSE) #HCV
net <- network(data_net_dl[, c("Term", "Genes")], directed = FALSE) #HDV

# Add diff attribute to network (COLOURS)
# vertex_color <- unique(data_net[, c("Genes", "group")])
vertex_color <- unique(data_net_dl[, c("Genes", "group")]) #TROCA O DATAFRAME
vertex_color <- setNames(vertex_color$group, vertex_color$Genes)

set.vertex.attribute(net, "group", vertex_color[network.vertex.names(net)])

##OBS: SEMPRE SUBSTITUA OS ARQUIVOS: data_net_b (HBV), data_net_c (HCV) data_net_dl (HDV)
# unique_systems <- unique(data_net$group)  # Analise quais são os nomes do system
df_color_dicio <- data.frame("group" = c("HBV Liver", "HBV PBMC", "HCV Liver",
                                         "HCV PBMC", "HDV Liver", "Synapse BP"), 
                             "color" = c("#DDF1F5", "#789DBC","#EDAEAE",
                                         "#BC7C7C", "#FFD09B", "#EEEDEB"))

color_palette <- setNames(df_color_dicio$color, df_color_dicio$group)  # O valor da cor

df_color <- data.frame("group" = get.vertex.attribute(net, "group"))
unique(df_color$group)

df_color %>% 
  mutate(group = ifelse(is.na(group), "Synapse BP", group)) -> df_color
set.vertex.attribute(net, "group", df_color$group)
print(df_color)

#--------  Adicionar atributo group ao objeto network (SHAPES). --------
# Criar uma lista de atributos
vertex_shape <- unique(data_net_dl[, c("Genes", "group")]) #TROCA O DATAFRAME
vertex_shape <- setNames(vertex_shape$group, vertex_shape$Genes)

# Adicionar atributos ao objeto network
set.vertex.attribute(net, "group", vertex_shape[network.vertex.names(net)]) #AQUI É SYSTEM!!!!!

# Definir os shapes manualmente
unique(data_net_dl$group)  # Analise quais são os nomes do system
# Crie um data.frame com esses nomes definindo as cores para eles
df_shape_dicio <- data.frame("group" = c("HBV Liver", "HBV PBMC", "HCV Liver",
                                         "HCV PBMC", "HDV Liver", "Synapse BP"), 
                             "shape" = c(19, 19, 19, 19, 19, 8))

# Defina os nomes em um data.frame
shape_palette <- setNames(df_shape_dicio$shape, df_shape_dicio$group)  # O valor da cor

# Obtenha os group correspondentes aos nós
df_shape <- data.frame("group" = get.vertex.attribute(net, "group"))
print(df_shape)
#Substitua os valores NA por "Other"/ou qualquer outra coisa (nesse caso Processos Biológicos/Genes)
df_shape %>% 
  mutate(group = ifelse(is.na(group), "Synapse BP", group)) -> df_shape
set.vertex.attribute(net, "group", df_shape$group)

#Extraindo nomes de nodes
node_names = as.data.frame(network.vertex.names(net))
names(node_names)[1] <- "names"
node_names <- node_names %>%
  mutate(
    names = str_remove_all(names, " HBV Liver| HBV PBMC| HCV Liver| HCV PBMC| HDV Liver")
  )
node_names = node_names$names

#Extraindo genes e BPs
#Genes
genes_synapse <- data_net %>%
  dplyr::select(Genes) %>%
  unique()
genes_synapse <- genes_synapse$Genes

HBV_liver_synapse_genes <- unique(data_net_bl$Genes)
HBV_pbmc_synapse_genes <- unique(data_net_bp$Genes)

HCV_liver_synapse_genes <- unique(data_net_cl$Genes)
HCV_pbmc_synapse_genes <- unique(data_net_cp$Genes)

HDV_synapse_genes <- unique(data_net_dl$Genes)

# Termos de interesse
terms_of_interest <- c(
  "Postsynapse",
  "Integral Component Of Postsynaptic Membrane",
  "Presynapse",
  "Integral Component Of Presynaptic Membrane",
  "Synapse",
  "Integral Component Of Postsynaptic Density Membrane"
)

genes_compartilhados <- data_net_dl %>%
  filter(Term %in% terms_of_interest) %>%
  select(Term, Genes) %>%
  distinct() %>%  # Remove repetições idênticas
  group_by(Genes) %>%
  summarise(n_terms = n_distinct(Term), .groups = "drop") %>%
  filter(n_terms > 1)  # Fica apenas com os genes compartilhados por 2 ou mais termos

# Filtrar genes correspondentes aos termos
genes_filtrados <- data_net_dl %>%
  filter(Term %in% terms_of_interest) %>%
  distinct(Genes) %>%
  pull(Genes)

# Criar um vetor único combinando termos e genes
vetor_final <- c(terms_of_interest, genes_filtrados)

#BPs
term_synapse <- data_net %>%
  dplyr::select(Term) %>%
  unique()
term_synapse <- term_synapse$Term

#HBV
# p11 <- ggnet2(
#   net,
#   size = "degree",
#   alpha = 0.9,
#   edge.color = "gray85",
#   shape.legend = TRUE,
#   legend.size = 18,
#   color = "group",
#   palette = color_palette,
#   shape = "group",
#   shape.palette = shape_palette
# ) +
#   geom_text(
#     aes(label = node_names),
#     size = ifelse(node_names %in% genes_synapse, 4, 3),  # Tamanho definido diretamente
#     fontface = "bold",
#     # max.overlaps = 60,
#     # box.padding = 0.5,
#     # point.padding = 0.2,
#     # min.segment.length = 0.5
#   ) +
#   theme(
#     legend.text = element_blank(),
#     legend.title = element_blank(),
#     legend.position = "none"
#   )
# 
# p11_1 <- ggnet2(
#   net,
#   size = "degree",
#   alpha = 0.9,
#   edge.color = "gray85",
#   shape.legend = TRUE,
#   legend.size = 18,
#   color = "group",
#   palette = color_palette,
#   shape = "group",
#   shape.palette = shape_palette
# ) +
#   geom_text(
#     aes(label = node_names),
#     size = ifelse(node_names %in% genes_synapse, 4, 3),  # Tamanho definido diretamente
#     fontface = "bold",
#     # max.overlaps = 60,
#     # box.padding = 0.5,
#     # point.padding = 0.2,
#     # min.segment.length = 0.5
#   ) +
#   theme(
#     legend.text = element_text(size = 12),
#     legend.title = element_text(size = 12),
#     legend.position = "right"
#   )

ggsave(filename = "network_hbv.svg",
       plot = p11,
       device = "svg", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300)  

ggsave(filename = "network_hbv_labels.svg",
       plot = p11_1,
       device = "svg", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300)  

#HCV
# p12 <- ggnet2(
#   net,
#   size = "degree",
#   alpha = 0.9,
#   edge.color = "gray85",
#   shape.legend = TRUE,
#   legend.size = 18,
#   color = "group",
#   palette = color_palette,
#   shape = "group",
#   shape.palette = shape_palette
# ) +
#   geom_text(
#     aes(label = node_names),
#     size = ifelse(node_names %in% genes_synapse, 4, 3),  # Tamanho definido diretamente
#     fontface = "bold",
#     # max.overlaps = 60,
#     # box.padding = 0.5,
#     # point.padding = 0.2,
#     # min.segment.length = 0.5
#   ) +
#   theme(
#     legend.text = element_blank(),
#     legend.title = element_blank(),
#     legend.position = "none"
#   )

# p12_1 <- ggnet2(
#   net,
#   size = "degree",
#   alpha = 0.9,
#   edge.color = "gray85",
#   shape.legend = TRUE,
#   legend.size = 18,
#   color = "group",
#   palette = color_palette,
#   shape = "group",
#   shape.palette = shape_palette
# ) +
#   geom_text(
#     aes(label = node_names),
#     size = ifelse(node_names %in% genes_synapse, 4, 3),  # Tamanho definido diretamente
#     fontface = "bold",
#     # max.overlaps = 60,
#     # box.padding = 0.5,
#     # point.padding = 0.2,
#     # min.segment.length = 0.5
#   ) +
#   theme(
#     legend.text = element_text(size = 12),
#     legend.title = element_text(size = 12),
#     legend.position = "right"
#   )

ggsave(filename = "network_hcv.svg",
       plot = p12,
       device = "svg", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300)  

ggsave(filename = "network_hcv_labels.svg",
       plot = p12_1,
       device = "svg", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300) 

p13 <- ggnet2(
  net,
  size = "degree",
  alpha = 0.9,
  edge.color = "gray85",
  shape.legend = TRUE,
  legend.size = 18,
  color = "group",
  palette = color_palette,
  shape = "group",
  shape.palette = shape_palette,
) +
  ggrepel::geom_text_repel(
    aes(label = node_names), bg.color = "white", bg.r = 0.1,
    size = ifelse(node_names %in% vetor_final, 4, 0),  # Tamanho definido diretamente
    fontface = "bold",
    # max.overlaps = 60,
    # box.padding = 0.5,
    # point.padding = 0.2,
    # min.segment.length = 0.5
  ) +
  theme(
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

p13_1 <- ggnet2(
  net,
  size = "degree",
  alpha = 0.9,
  edge.color = "gray85",
  shape.legend = TRUE,
  # legend.size = 18,
  color = "group",
  palette = color_palette,
  shape = "group",
  shape.palette = shape_palette,
) +
  geom_text(
    aes(label = node_names),
    size = ifelse(node_names %in% vetor_final, 4, 0),  # Tamanho definido diretamente
    fontface = "bold",
    # max.overlaps = 60,
    # box.padding = 0.5,
    # point.padding = 0.2,
    # min.segment.length = 0.5
  ) +
  theme(
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 12, colour = "black"),
    legend.position = "none"
  )

ggsave(filename = "network_hdv.svg",
       plot = p13,
       device = "svg", 
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300)  

ggsave(filename = "network_hdv_labels.svg",
       plot = p13_1,
       device = "svg",
       width = 6, 
       height = 6, 
       units = "in", 
       dpi = 300)  

# Dotplot (SynGO) ---------------------------------------------------------

data_heat_syn <- data_syn

data_heat_bl <- data_heat_syn %>%
  dplyr::filter(group == "HBV Liver")

data_heat_bp <- data_heat_syn %>%
  dplyr::filter(group == "HBV PBMC")

data_heat_b <- data_heat_syn %>%
  dplyr::filter(group %in% c("HBV Liver", "HBV PBMC"))

data_heat_cl <- data_heat_syn %>%
  dplyr::filter(group == "HCV Liver")

data_heat_cp <- data_heat_syn %>%
  dplyr::filter(group == "HCV PBMC")

data_heat_c <- data_heat_syn %>%
  dplyr::filter(group %in% c("HCV Liver", "HCV PBMC"))

data_heat_dl <- data_heat_syn %>%
  dplyr::filter(group == "HDV Liver")

data_heat_c$Hepatitis <- "HCV"
data_heat_b$Hepatitis <- "HBV"
data_heat_dl$Hepatitis <- "HDV"

data_heat_syn <- rbind(data_heat_c, data_heat_b, data_heat_dl)

# data_heat_syn <- data_heat_syn %>%
#   group_by(group) %>%
#   arrange(Adjusted.P.value) %>% # Ordenar os termos dentro de cada grupo pelo Adjusted.P.value (crescente)
#   slice_head(n = 5) %>%         # Selecionar os 5 primeiros termos mais significativos
#   ungroup()

# Definir os termos que devem ser mantidos em todos os grupos
common_terms <- data_heat_syn %>%
  group_by(Term) %>%
  dplyr::filter(n_distinct(group) == n_distinct(data_heat_syn$group)) %>%
  pull(Term)

specific_terms <- unique(common_terms)

# Selecionar os 5 termos mais significativos dentro de cada grupo, incluindo os termos específicos
data_heat_syn <- data_heat_syn %>%
  group_by(group) %>%
  # Selecionar os 5 termos mais significativos, mas garantir que os termos específicos sejam incluídos
  dplyr::filter(Term %in% specific_terms | Term %in% Term[order(Adjusted.P.value)[1:5]]) %>%
  ungroup()

# data_heat_syn <- data_heat_syn %>%
#   separate_rows(Genes, sep = ";")

data_heat_syn <- data_heat_syn %>%
  mutate(Overlap_numeric = as.numeric(sub("/.*", "", Overlap)) / 
           as.numeric(sub(".*/", "", Overlap)))

# data_heat_syn <- data_heat_syn %>%
#   mutate(Alpha = (Overlap_numeric - min(Overlap_numeric)) / 
#            (max(Overlap_numeric) - min(Overlap_numeric)))

data_heat_syn <- data_heat_syn %>%
  group_by(group) %>%
  mutate(alpha = Overlap_numeric / max(Overlap_numeric, na.rm = TRUE)) %>%
  ungroup()

group_colors <- c(
  "HBV Liver" = "#DDF1F5",
  "HBV PBMC" = "#789DBC",
  "HCV Liver" = "#EDAEAE",
  "HCV PBMC" = "#BC7C7C",
  "HDV Liver" = "#FFD09B"
)

group_shapes <- c(
  "HBV Liver" = 19,
  "HBV PBMC" = 18,
  "HCV Liver" = 19,
  "HCV PBMC" = 18,
  "HDV Liver" = 19 
)

data_heat_syn <- data_heat_syn %>%
  mutate(alpha = case_when(
    Adjusted.P.value < 0.0005 ~ 1,   # Para valores < 0.0005
    Adjusted.P.value < 0.005  ~ 0.75,   # Para valores < 0.005
    Adjusted.P.value < 0.05   ~ 0.5,    # Para valores < 0.01
    TRUE ~ 0.1                    # Para valores >= 0.05
  ))

p24 <- ggplot(data_heat_syn, aes(x = Overlap, y = reorder(Term, Adjusted.P.value))) +
  geom_point(aes(size = alpha, color = group,
                 shape = group,
                 alpha = alpha)) +  # Adicionando forma e ajustando a transparência com alpha
  scale_size_continuous(name = "Adj.p-value", 
                        breaks = c(0.1, 0.5, 0.75, 1), 
                        labels = c(">= 0.05", "< 0.05", "< 0.005", "< 0.0005"),
                        range = c(2, 10)) +
  scale_shape_manual(values = group_shapes, name = "Group") +  # Escala de formas manual
  scale_color_manual(values = group_colors, name = "Group") +  # Escala de cores manual
  scale_alpha_continuous(name = "Adj.p-value", 
                         breaks = c(0.1, 0.5, 0.75, 1), 
                         labels = c(">= 0.05", "< 0.05", "< 0.005", "< 0.0005"),
                         range = c(0.1, 1)) +  # Ajuste da escala de alpha com legendas customizadas
  theme_bw() +
  facet_wrap(~group, scales = "free_x", ncol = 5) +
  theme(axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 18, colour = "black"),
        legend.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.position = "right") +  # Mostrar legenda de grupo
  labs(x = "Genes Overlap", y = "", title = "") +
  guides(alpha = guide_legend(title = "Adj.p-value"),
         size = guide_legend(title = "Adj.p-value"),
         color = guide_legend(override.aes = list(size = 5)))

#ALTERNATIVO
# Adiciona uma coluna auxiliar para -log10(Adjusted.P.value)
data_heat_syn <- data_heat_syn %>%
  mutate(logAdjP = -log10(Adjusted.P.value))

# Organiza os termos dentro de cada grupo com base no valor -log10(Adj.P.value)
data_heat_syn <- data_heat_syn %>%
  mutate(logAdjP = -log10(Adjusted.P.value)) %>%
  group_by(group) %>%
  mutate(Term = fct_reorder(Term, logAdjP)) %>%
  ungroup()

# Gráfico
p24 <- ggplot(data_heat_syn, aes(x = logAdjP, y = Term)) +
  geom_point(aes(size = alpha, color = group,
                 shape = group,
                 alpha = alpha), stroke = 0.6) +
  scale_x_continuous(
    name = expression(-log[10](Adjusted~P~value)),
    limits = c(0, 22),  # ajuste se necessário
    breaks = seq(0, 22, by = 5)
  ) +
  scale_size_continuous(name = "Adj.p-value", 
                        breaks = c(0.1, 0.5, 0.75, 1), 
                        labels = c(">= 0.05", "< 0.05", "< 0.005", "< 0.0005"),
                        range = c(3, 10)) +
  scale_shape_manual(values = group_shapes, name = "Group") +
  scale_color_manual(values = group_colors, name = "Group") +
  scale_alpha_continuous(name = "Adj.p-value", 
                         breaks = c(0.1, 0.5, 0.75, 1), 
                         labels = c(">= 0.05", "< 0.05", "< 0.005", "< 0.0005"),
                         range = c(0.3, 1)) +
  theme_bw(base_size = 14) +
  facet_wrap(~group, scales = "fixed", ncol = 5) +  # escala fixa aqui
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  guides(
    alpha = guide_legend(title = "Adj.p-value"),
    size = guide_legend(title = "Adj.p-value"),
    color = guide_legend(override.aes = list(size = 5))
  )

ggsave(filename = "Dotplot_SynGO2024_ALTERNATIVO.tiff",
       plot = p24,
       device = "tiff", 
       width = 14, 
       height = 6, 
       units = "in", 
       dpi = 300)  

ggsave(filename = "Dotplot_SynGO2024_ALTERNATIVO.svg",
       plot = p24,
       device = "svg", 
       width = 14, 
       height = 6, 
       units = "in", 
       dpi = 300)  

write.xlsx(data_heat_syn, file.path("Suppl.Table2c_enrich_SynGO_commumBPs.xlsx"))
write.xlsx(data_heat_c, file.path("Suppl.Table2b_NetworkHCV.xlsx"))
write.xlsx(data_heat_b, file.path("Suppl.Table2b_NetworkHBV.xlsx"))
write.xlsx(data_heat_dl, file.path("Suppl.Table2b_NetworkHDV.xlsx"))

data_heat_syn <- read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/5 - SynGO/Suppl.Table2c_enrich_SynGO_commumBPs.xlsx")

# Enriquecimento em Comum -------------------------------------------------
enrich_liver <- read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoComumLiver.csv")
enrich_liver <- enrich_liver %>%
  group_by(Genes) %>%
  distinct(Term, .keep_all = TRUE) %>%  # Mantém a primeira ocorrência de Term por Genes
  ungroup()

enrich_pbmc <- read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/4 - Pathways Caracterization/Commum/EnriquecimentoComumPBMC.csv")
enrich_pbmc <- enrich_pbmc %>%
  group_by(Genes) %>%
  distinct(Term, .keep_all = TRUE) %>%  # Mantém a primeira ocorrência de Term por Genes
  ungroup()

# Filtrar o data_heat_syn original para manter apenas os termos em common_terms_df
enrich_synGO <- data_heat_syn %>%
  dplyr::filter(Term %in% common_terms)
enrich_synGO$sytem <- "Synapse"
names(enrich_synGO)[13] <- "hepatitis"

#BubbleHeatmap figado
colnames(enrich_liver)
enrich_liver <- enrich_liver %>%
  distinct(Term, Genes, hepatites, .keep_all = TRUE)

extrair_extremos <- function(df) {
  # Filtrar as categorias "Immune" e "Nervous"
  df_immune <- subset(df, sytem == "Immune")
  df_nervous <- subset(df, sytem == "Nervous")
  
  # Ordenar pelos menores valores de P.value
  df_immune <- df_immune[order(df_immune$P.value), ]
  df_nervous <- df_nervous[order(df_nervous$P.value), ]
  
  # Extrair os 10 elementos mais significantes de cada categoria
  top_immune <- head(df_immune, 10)
  top_nervous <- head(df_nervous, 10)
  
  # Combinar os resultados
  resultado <- rbind(top_immune, top_nervous)
  
  return(resultado)
}

top20_tems <- extrair_extremos(enrich_liver)
data_filt <- enrich_liver %>%
  dplyr::filter(Term %in% c("Nervous System Development (GO:0007399)",
                     "Mitotic Cytokinesis (GO:0000281)",
                     "Cytoskeleton-Dependent Cytokinesis (GO:0061640)",
                     "Acute Inflammatory Response (GO:0002526)",
                     "Cellular Response To Chemokine (GO:1990869)",
                     "Regulation Of Fibrinolysis (GO:0051917)",
                     "Cellular Response To Cytokine Stimulus (GO:0071345)", 
                     "Regulation Of Cell Cycle Process (GO:0010564)",
                     "Lymphocyte Chemotaxis (GO:0048247)",
                     "Regulation Of Microtubule Polymerization Or Depolymerization (GO:0031110)", 
                     "Acute-Phase Response (GO:0006953)",
                     "Cellular Glucuronidation (GO:0052695)",
                     "Positive Regulation Of Neural Precursor Cell Proliferation (GO:2000179)", 
                     "Negative Regulation Of Dendritic Cell Apoptotic Process (GO:2000669)",
                     "Regulation Of Dendritic Cell Apoptotic Process (GO:2000668)",
                     "Nitric Oxide Biosynthetic Process (GO:0006809)",
                     "Transport Across Blood-Brain Barrier (GO:0150104)",
                     "Nitric Oxide Biosynthetic Process (GO:0006809)",
                     "Cellular Glucuronidation (GO:0052695)",
                     "Anterograde Trans-Synaptic Signaling (GO:0098916)"))

data_filt <- data_filt %>%
  distinct(Term, hepatites, .keep_all = TRUE)

data_filt <- data_filt %>%
  distinct(Term, hepatites, .keep_all = TRUE)
# Criando a coluna "GeneRatio" a partir de "Overlap"
data_filt <- data_filt %>%
  mutate(GeneRatio = as.numeric(sub("/.*", "", Overlap)) / as.numeric(sub(".*/", "", Overlap)))
data_filt$GO <- str_extract(data_filt$Term, "\\((.*?)\\)")
data_filt$Term <- sub("\\s*\\([^\\)]+\\)", "",
                      data_filt$Term)
data_filt$Term <- sub("\\s*\\([^\\)]+\\)", "",
                      data_filt$Term)
# Defina a ordem desejada dos níveis na coluna sytem
data_filt$sytem <- factor(data_filt$sytem, levels = c("Nervous", "Immune"))

# Reordene a coluna Term com base na ordem da coluna sytem
data_filt <- data_filt[order(data_filt$sytem, data_filt$Term), ]

# Se necessário, você pode reordenar Term como um fator também
data_filt$Term <- factor(data_filt$Term, levels = unique(data_filt$Term))
names(data_filt)[13] <- "System"

p1 <- ggplot(data = data_filt,
             aes(x = hepatites,
                 y = Term)) +
  geom_point(aes(size = -log10(GeneRatio),
                 fill = P.value,
                 shape = System),
             color = "black") +
  scale_shape_manual(values = c("Immune" = 21,
                                "Nervous" = 23)) +  # Mapear shapes para valores da coluna sytem
  scale_fill_gradientn(colours = c("#93BFCF", "#D8D8D8", "#D71313"),
                       name = "P-value",
                       trans = "log10") +  # Configurando a escala de cores
  scale_size_continuous(name = "GeneRatio",
                        range = c(2, 10)) +  # Configurando o tamanho das bolhas
  theme_bw(base_size = 12) +  # Estilo do tema minimalista
  theme(strip.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 24),
        axis.text.x = element_text(size = 24),  # Rotacionando o texto do eixo x
        axis.text.y = element_text(size = 18),  # Tamanho do texto no eixo y
        panel.grid.major = element_blank(),  # Removendo grades principais
        panel.grid.minor = element_blank(), # Removendo grades secundárias
        legend.position = "right") + 
  labs(x = "", y = "", title = "Liver")

#SalvandoPlot
ggsave(filename = "Liver_PBMC.tiff",
       plot = p1,
       device = "tiff",  # Define o formato do arquivo
       width = 11,  # Largura do gráfico em polegadas
       height = 8,  # Altura do gráfico em polegadas
       units = "in", # Unidade de medida
       dpi = 300)  # qualidade

ggsave(filename = "Liver_PBMC.svg",
       plot = p1,
       device = "svg",  # Define o formato do arquivo
       width = 11,  # Largura do gráfico em polegadas
       height = 8,  # Altura do gráfico em polegadas
       units = "in", # Unidade de medida
       dpi = 300)  # qualidade

enrich_pbmc <- enrich_pbmc %>%
  distinct(Term, Genes, hepatites, .keep_all = TRUE)

top20_tems <- extrair_extremos(enrich_pbmc)

data_filt <- enrich_pbmc %>%
  dplyr::filter(Term %in% c(
    "Modulation Of Chemical Synaptic Transmission (GO:0050804)",
    "Nervous System Development (GO:0007399)",
    "Regulation Of Neuron Death (GO:1901214)",
    "Neuropeptide Signaling Pathway (GO:0007218)",
    "Neuron Differentiation (GO:0030182)",
    "Generation Of Neurons (GO:0048699)",
    "Neuron Projection Morphogenesis (GO:0048812)", 
    "Modulation Of Chemical Synaptic Transmission (GO:0050804)",
    "Regulation Of Apoptotic Process (GO:0042981)",
    "Positive Regulation Of Intracellular Signal Transduction (GO:1902533)", 
    "Positive Regulation Of Multicellular Organismal Process (GO:0051240)",
    "Inflammatory Response (GO:0006954)",
    "Regulation Of Cell Migration (GO:0030334)", 
    "Regulation Of GTPase Activity (GO:0043087)",
    "Cytokine-Mediated Signaling Pathway (GO:0019221)",
    "Regulation Of ERK1 And ERK2 Cascade (GO:0070372)",
    "Regulation Of Tumor Necrosis Factor Production (GO:0032680)",
    "Cellular Response To Interleukin-1 (GO:0071347)",
    "Response To Type II Interferon (GO:0034341)",
    "T Cell Receptor Signaling Pathway (GO:0050852)"
  ))

data_filt <- data_filt %>%
  distinct(Term, hepatites, .keep_all = TRUE)
# Criando a coluna "GeneRatio" a partir de "Overlap"
data_filt <- data_filt %>%
  mutate(GeneRatio = as.numeric(sub("/.*", "", Overlap)) / as.numeric(sub(".*/", "", Overlap)))
data_filt$GO <- str_extract(data_filt$Term, "\\((.*?)\\)")
data_filt$Term <- sub("\\s*\\([^\\)]+\\)", "",
                      data_filt$Term)
data_filt$Term <- sub("\\s*\\([^\\)]+\\)", "",
                      data_filt$Term)
# Defina a ordem desejada dos níveis na coluna sytem
data_filt$sytem <- factor(data_filt$sytem, levels = c("Nervous", "Immune"))

# Reordene a coluna Term com base na ordem da coluna sytem
data_filt <- data_filt[order(data_filt$sytem, data_filt$Term), ]

# Se necessário, você pode reordenar Term como um fator também
data_filt$Term <- factor(data_filt$Term, levels = unique(data_filt$Term))
names(data_filt)[13] <- "System"

p2 <- ggplot(data = data_filt,
             aes(x = hepatites,
                 y = Term)) +
  geom_point(aes(size = -log10(GeneRatio),
                 fill = P.value,
                 shape = System),
             color = "black") +
  scale_shape_manual(values = c("Immune" = 21,
                                "Nervous" = 23)) +  # Mapear shapes para valores da coluna sytem
  scale_fill_gradientn(colours = c("#93BFCF", "#D8D8D8", "#D71313"),
                       name = "P-value",
                       trans = "log10") +  # Configurando a escala de cores
  scale_size_continuous(name = "GeneRatio",
                        range = c(2, 10)) +  # Configurando o tamanho das bolhas
  theme_bw(base_size = 12) +  # Estilo do tema minimalista
  theme(strip.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 24),
        axis.text.x = element_text(size = 24),  # Rotacionando o texto do eixo x
        axis.text.y = element_text(size = 18),  # Tamanho do texto no eixo y
        panel.grid.major = element_blank(),  # Removendo grades principais
        panel.grid.minor = element_blank(), # Removendo grades secundárias
        legend.position = "right") + 
  labs(x = "", y = "", title = "PBMC")

# ggsave(filename = "Bubble_PBMC.tiff",
#        plot = p2,
#        device = "tiff",  # Define o formato do arquivo
#        width = 10,  # Largura do gráfico em polegadas
#        height = 8,  # Altura do gráfico em polegadas
#        units = "in", # Unidade de medida
#        dpi = 400)  # qualidade

#Salvando tabelas
write.xlsx(data_filt, file.path("SupplementaryFig4c_enrich_pbmc.xlsx"))
write.xlsx(data_filt, file.path("SupplementaryFig4c_enrich_liver.xlsx"))


enrich_liver <- enrich_liver %>%
  group_by(Genes) %>%
  distinct(Term, .keep_all = TRUE) %>%  # Mantém a primeira ocorrência de Term por Genes
  ungroup()

#Unindo tudo
enrich_unique_liver <- enrich_liver %>%
  group_by(hepatites) %>%
  distinct(Term, .keep_all = TRUE) %>%  # Mantém a primeira ocorrência de Term por Genes
  ungroup()

enrich_unique_pbmc <- enrich_pbmc %>%
  group_by(hepatites) %>%
  distinct(Term, .keep_all = TRUE) %>%  # Mantém a primeira ocorrência de Term por Genes
  ungroup()

enrich_unique_liver <- enrich_unique_liver %>%
  mutate(group = paste(hepatites, tissue, sep = " "))
enrich_unique_liver$diff <- NULL
enrich_unique_liver$tissue <- NULL
enrich_unique_liver$hepatites <- NULL

enrich_unique_pbmc <- enrich_unique_pbmc %>%
  mutate(group = paste(hepatites, tissue, sep = " "))
enrich_unique_pbmc$diff <- NULL
enrich_unique_pbmc$tissue <- NULL
enrich_unique_pbmc$hepatites <- NULL

enrich_synGO$Ontology <- NULL
enrich_synGO$alpha <- NULL
enrich_synGO$Overlap_numeric <- NULL
enrich_synGO$GOID <- NULL
enrich_synGO$hepatitis <- NULL

#unindo dfs
enrich_all <- rbind(enrich_unique_liver, enrich_unique_pbmc, enrich_synGO)

write.xlsx(enrich_synGO, file.path("enrich_SynGO_commumBPs.xlsx"))

# Intersection Analysis --------------------------------------------------------
library(UpSetR)
library(ComplexUpset)
library(dplyr)
library(tidyr)
# Carregando dados

data_HBVLup = read.table("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/1 - MetaDEGs/HBV_Liver/Up/GO_Biological_Process_2023_table.txt", sep = "\t", header = TRUE)
data_HBVLdown = read.table("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/1 - MetaDEGs/HBV_Liver/Down/GO_Biological_Process_2023_table (1).txt", sep = "\t", header = TRUE)
data_HBVPup = read.table("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/1 - MetaDEGs/HBV_PBMC/Up/GO_Biological_Process_2023_table.txt", sep = "\t", header = TRUE)
data_HBVPdown = read.table("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/1 - MetaDEGs/HBV_PBMC/Down/GO_Biological_Process_2023_table.txt", sep = "\t", header = TRUE)
data_HCVLup = read.table("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/1 - MetaDEGs/HCV_Liver/Up/GO_Biological_Process_2023_table.txt", sep = "\t", header = TRUE)
data_HCVLdown = read.table("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/1 - MetaDEGs/HCV_Liver/Down/GO_Biological_Process_2023_table.txt", sep = "\t", header = TRUE)
data_HCVPup = read.table("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/1 - MetaDEGs/HCV_PBMC/Up/GO_Biological_Process_2023_table.txt", sep = "\t", header = TRUE)
data_HCVPdown = read.table("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/1 - MetaDEGs/HCV_PBMC/Down/GO_Biological_Process_2023_table.txt", sep = "\t", header = TRUE)
data_HDVLup = read.table("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/1 - MetaDEGs/HDV_Liver/Up/GO_Biological_Process_2023_table.txt", sep = "\t", header = TRUE)
data_HDVLdown = read.table("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/1 - MetaDEGs/HDV_Liver/Down/GO_Biological_Process_2023_table.txt", sep = "\t", header = TRUE)
enrich_synGO <- read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/5 - SynGO/enrich_SynGO_commumBPs.xlsx")

#ou carregue os objetos já processados:
load("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/Pathwayscomuns.RData")

#Editando dados
data_HBVLup$group <- "HBV Liver"
data_HBVLdown$group <- "HBV Liver"
data_HBVPup$group <- "HBV PBMC"
data_HBVPdown$group <- "HBV PBMC"
data_HCVLup$group <- "HCV Liver"
data_HCVLdown$group <- "HCV Liver"
data_HCVPup$group <- "HCV PBMC"
data_HCVPdown$group <- "HCV PBMC"
data_HDVLup$group <- "HDV Liver"
data_HDVLdown$group <- "HDV Liver"

enrich_all <- rbind(data_HBVLup, data_HBVLdown, data_HBVPup, data_HBVPdown,
                    data_HCVLup, data_HCVLdown, data_HCVPup, data_HCVPdown,
                    data_HDVLup, data_HDVLdown)

#System
#Separando Neuro e Imune
neuro_keywords <- c('neuron', 'neuro', 'synapse', 'brain', 'nervous', "Anterograde",
                    "neural", "synaptic", "telencephalon", "neuronal", 'Neuroblast',
                    "axon", "neurotransmitter", "serotonin", "neurogenesis", 
                    "dopamine", "neuropeptide", "oligodendrocyte", 'Neurotransmitter',
                    "nerve impulse", "limbic", "Dendritic", 'Glial', 'Glutamate',
                    'Neuroinflammatory', 'Neuropeptide', 'Nerve', 'Sodium Ion Homeostasis',
                    'Behavior', 'Potassium Ion Homeostasis',
                    'Potassium Ion Import Across Plasma Membrane', 'Glutamine Family Amino Acid Catabolic Process',
                    'Glutamine Family Amino Acid Biosynthetic Process', 'Dendrite Morphogenesis',
                    'acetyl-CoA Metabolic Process', 'Sodium-Dependent Phosphate Transport',
                    'Nitric Oxide Biosynthetic Process', 'Cellular Glucuronidation', "Synaptic", "neuro", "Synapse",
                    "Perisynaptic", "Postsynaptic",
                    "Neurotransmitter", "Neuroinflammatory",
                    "Neuron", "Postsynapse", "Presynapse",
                    "Presynaptic", "Nervous", "Dopamine", "Axon",
                    "Neurogenesis", "Adrenergic", "Dendritic", "Neural",
                    "Neuromuscular", "Neuroblast", "Sensory", "Gliogenesis",
                    "Learning", "Brain")
immune_keywords <- c('immune', 'immunoglobulin', 'lymphocyte', 'cytokine', 
                     'B cell', 'Apoptotic', 'Leukocyte', 'Calcium', 'Inflammatory', 
                     'interleukin', 'T cell', 'monocyte', 'T-helper', 'interferon',
                     'chemokine', 'MHC', 'Neutrophil', 'Antigen', 'Regulation', 'Defense',
                     'ERK1', 'ERK2', 'MAPK', 'Natural Killer', 'Macrophage', 'Kinase',
                     'Humoral', 'Antibacterial', 'Modulation By Host Of Viral Process',
                     'Cellular Response To Tumor Necrosis Factor', 'Eosinophil Migration',
                     'Pathway-Restricted SMAD Protein Phosphorylation', 'Response To BMP',
                     'Endocytosis', 'Toll-Like Receptor 9 Signaling Pathway', 'Response To Tumor Necrosis Factor',
                     'Cellular Response To Decreased Oxygen Levels', 'TORC2 Signaling',
                     'Response To cGMP', 'Cellular Response To Hypoxia', 'Pathway-Restricted SMAD Protein Phosphorylation',
                     'Response To BMP', 'Mitotic Sister Chromatid Segregation', 'DNA Replication Checkpoint Signaling',
                     'G1/S Transition Of Mitotic Cell Cycle', 'G2/M Transition Of Mitotic Cell Cycle',
                     'DNA Repair', 'Acute-Phase Response', 'Response To Cadmium Ion',
                     'Response To Copper Ion', 'Cellular Response To Metal Ion', 'Cellular Response To Cadmium Ion',
                     'Cellular Response To Copper Ion', 'Cellular Response To Zinc Ion', 'Response To Fibroblast Growth Factor',"leukocyte", "T cell", "mononuclear",
                     "lymphocyte", "immune response", "B cell",
                     "interferon", "Inflammatory", "Immune", "Defense",
                     "Cytokine", "Kinase", "ERBB2", "G Protein", 
                     "Cysteine", "ERK1 And ERK2", "Interleukin",
                     "Macrophage", "Necrosis", "Tumor", "MAPK",
                     "Healing", "Humoral", "Neutrophil", "Leukocyte",
                     "Myeloid", "Response", "Negative Regulation Of Type II Interferon-Mediated Signaling Pathway",
                     "Positive Regulation Of Growth", "MHC", "Cell Growth",
                     "Toll-Like Receptor", "Monocyte", "Interferon", "STAT",
                     "Prostaglandin", "NF-kappaB","GTPase", "Apoptotic", "Viral",
                     "Modulation", "Immunoglobulin", "NIK/NF-kappaB", "B Cell",
                     "T Cell", "I-kappaB kinase/NF-kappaB", "kinase", "Death",
                     "Ubiquitin", "Ubiquitin-Dependent","Ubiquitination",
                     "Proliferation", "Growth", "Mast", "Regulation")
enrich_all <- enrich_all %>%
  mutate(sytem = case_when(
    str_detect(Term, paste(neuro_keywords,
                           collapse = "|")) ~ "Nervous",
    str_detect(Term, paste(immune_keywords,
                           collapse = "|")) ~ "Immune",
    TRUE ~ "Others"
  ))

enrich_all_together <- rbind(enrich_all, enrich_synGO)

# Filtrar para os grupos Figado
enrich_all_liver <- subset(enrich_all_together, group %in% c("HBV Liver", "HCV Liver", "HDV Liver"))

# Filtrar para os grupos PBMC
enrich_all_pbmc <- subset(enrich_all_together, group %in% c("HCV PBMC", "HBV PBMC"))

intersection_summary <- enrich_all_pbmc %>%
  group_by(Term) %>%
  summarize(
    Groups = paste(unique(group), collapse = ", "),
    Count = n()
  ) %>%
  arrange(desc(Count))

set_cols <- setdiff(colnames(upset_data), c("Term", "Count"))

# Preparar os dados para o Upset Plot
upset_data <- intersection_summary %>%
  separate_rows(Groups, sep = ", ") %>%   # Separar os grupos
  mutate(value = 1) %>%                   # Adicionar coluna para binarização
  pivot_wider(names_from = Groups, values_from = value, values_fill = 0) # Converter para formato largo

# Identificar as colunas dos conjuntos, excluindo 'Term' e 'Count'
set_cols <- setdiff(colnames(upset_data), c("Term", "Count"))

# Criar o Upset Plot sem o argumento set_sizes (usa o padrão)
upset(
  upset_data,
  intersect = set_cols,
  name = "Term Intersections"
)

intersection_summary_liver <- as.data.frame(intersection_summary)
intersection_summary_pbmc <- as.data.frame(intersection_summary)

write.xlsx(intersection_summary_liver, file.path("intersectionVenn_Liver.xlsx"))
write.xlsx(intersection_summary_pbmc, file.path("intersectionVenn_PBMC.xlsx"))

#Intersection NetWork
intersection_liver_filt <- subset(intersection_summary_liver, Groups == "HBV Liver, HCV Liver, HDV Liver")
intersection_pbmc_filt <- subset(intersection_summary_pbmc, Groups == "HBV PBMC, HCV PBMC")

intersection_liver_vector <- intersection_liver_filt$Term
intersection_pbmc_vector <- intersection_pbmc_filt$Term

filtered_enrich_all_liver <- subset(enrich_all_liver, Term %in% intersection_liver_vector)
filtered_enrich_all_pbmc <- subset(enrich_all_pbmc, Term %in% intersection_pbmc_vector)

#Criando Network
data_net_liver <- filtered_enrich_all_liver
data_net_pbmc <- filtered_enrich_all_pbmc

data_net_liver <- data_net_liver %>%
  separate_rows(Genes, sep = ";")
data_net_pbmc <- data_net_pbmc %>%
  separate_rows(Genes, sep = ";")

#Concatenando
data_net_liver <- data_net_liver %>%
  mutate(
    Genes = as.character(Genes),
    group = as.character(group),
    cgene = paste(Genes, group, sep = " ")
  )
data_net_liver <- data_net_liver %>%
  mutate(
    Term = as.character(Term),
    group = as.character(group),
    cterm = paste(Term, group, sep = " ")
  )

data_net_pbmc <- data_net_pbmc %>%
  mutate(
    Genes = as.character(Genes),
    group = as.character(group),
    cgene = paste(Genes, group, sep = " ")
  )
data_net_pbmc <- data_net_pbmc %>%
  mutate(
    Term = as.character(Term),
    group = as.character(group),
    cterm = paste(Term, group, sep = " ")
  )

net <- network(data_net_liver[, c("Term", "cgene")], directed = FALSE) #Liver
net <- network(data_net_pbmc[, c("Term", "cgene")], directed = FALSE) #PBMC

# Add diff attribute to network (COLOURS)
# vertex_color <- unique(data_net[, c("Genes", "group")])
vertex_color <- unique(data_net_pbmc[, c("cgene", "sytem")]) #TROCA O DATAFRAME: Liver ou PBMC
vertex_color <- setNames(vertex_color$sytem, vertex_color$cgene)

set.vertex.attribute(net, "sytem", vertex_color[network.vertex.names(net)])

unique(data_net_pbmc$sytem)
# unique_systems <- unique(data_net_liver$group)  # Analise quais são os nomes do system
df_color_dicio <- data.frame("sytem" = c("Immune", "Others",
                                         "Nervous", "Intersection BPs"), 
                             "color" = c("#3182bd", "#93BFCF",
                                         "tomato", "#EEEDEB"))

color_palette <- setNames(df_color_dicio$color, df_color_dicio$sytem)  # O valor da cor

df_color <- data.frame("sytem" = get.vertex.attribute(net, "sytem"))
unique(df_color$sytem)

df_color %>% 
  mutate(sytem = ifelse(is.na(sytem), "Intersection BPs", sytem)) -> df_color
set.vertex.attribute(net, "sytem", df_color$sytem)
print(df_color)

#Dotplot enrichment

#--------  Adicionar atributo group ao objeto network (SHAPES). --------
# Criar uma lista de atributos
vertex_shape <- unique(data_net_pbmc[, c("cgene", "group")]) #TROCA O DATAFRAME
vertex_shape <- setNames(vertex_shape$group, vertex_shape$cgene)

# Adicionar atributos ao objeto network
set.vertex.attribute(net, "group", vertex_shape[network.vertex.names(net)]) #AQUI É SYSTEM!!!!!

# Definir os shapes manualmente
unique(data_net_pbmc$group)  # Analise quais são os nomes do system
# Crie um data.frame com esses nomes definindo as cores para eles
# df_shape_dicio <- data.frame("group" = c("HBV Liver", "HCV Liver",
#                                          "HDV Liver", "Intersection BPs"), 
#                              "shape" = c(19, 19, 19, 8))

df_shape_dicio <- data.frame("group" = c("HBV PBMC", "HCV PBMC","Intersection BPs"), 
                             "shape" = c(19, 19, 8))
# Defina os nomes em um data.frame
shape_palette <- setNames(df_shape_dicio$shape, df_shape_dicio$group)  # O valor da cor

# Obtenha os group correspondentes aos nós
df_shape <- data.frame("group" = get.vertex.attribute(net, "group"))
print(df_shape)
#Substitua os valores NA por "Other"/ou qualquer outra coisa (nesse caso Processos Biológicos/Genes)
df_shape %>% 
  mutate(group = ifelse(is.na(group), "Intersection BPs", group)) -> df_shape
set.vertex.attribute(net, "group", df_shape$group)

#Liver
net_liver <- ggnet2(
  net,
  size = "degree",
  alpha = 0.9,
  edge.color = "gray90",
  shape.legend = TRUE,
  legend.size = 18,
  color = "sytem",
  palette = color_palette,
  shape = "group",
  shape.palette = shape_palette
) +
  theme(
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

#PBMC
net_pbmc <- ggnet2(
  net,
  size = "degree",
  alpha = 0.9,
  edge.color = "gray95",
  shape.legend = TRUE,
  legend.size = 18,
  color = "sytem",
  palette = color_palette,
  shape = "group",
  shape.palette = shape_palette
) +
  theme(
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

#Legends
net_legends <- ggnet2(
  net,
  size = "degree",
  alpha = 0.9,
  edge.color = "gray95",
  shape.legend = TRUE,
  legend.size = 18,
  color = "sytem",
  palette = color_palette,
  shape = "group",
  shape.palette = shape_palette
)

ggsave(filename = "Figure1f_NetworkLiver.tiff",
       plot = net_liver,
       device = "tiff", 
       width = 4, 
       height = 4, 
       units = "in", 
       dpi = 300)

ggsave(filename = "Figure1f_NetworkPBMC.tiff",
       plot = net_pbmc,
       device = "tiff", 
       width = 4, 
       height = 4, 
       units = "in", 
       dpi = 300)

ggsave(filename = "Figure1f_NetworkLegends.svg",
       plot = net_legends,
       device = "svg", 
       width = 8, 
       height = 8, 
       units = "in", 
       dpi = 300)

write.xlsx(data_net_liver, file.path("Suppl.Table1f_LiverNeuroimmuneNetwork.xlsx"))
write.xlsx(data_net_pbmc, file.path("Suppl.Table1f_PBMCNeuroimmuneNetwork.xlsx"))

#UpsetPlot Genes Intersection
intersection_genes_upset <- rbind(filtered_enrich_all_liver, filtered_enrich_all_pbmc)
intersection_genes_upset <- subset(intersection_genes_upset, sytem %in% c("Immune", "Nervous"))
intersection_genes_upset <- intersection_genes_upset %>%
  separate_rows(Genes, sep = ";")

intersection_summary <- intersection_genes_upset %>%
  group_by(Genes) %>%
  summarize(
    Groups = paste(unique(group), collapse = ", "),
    Count = n()
  ) %>%
  arrange(desc(Count))

set_cols <- setdiff(colnames(upset_data), c("Genes", "Count"))

# Preparar os dados para o Upset Plot
upset_data <- intersection_summary %>%
  separate_rows(Groups, sep = ", ") %>%   # Separar os grupos
  mutate(value = 1) %>%                   # Adicionar coluna para binarização
  pivot_wider(names_from = Groups, values_from = value, values_fill = 0) # Converter para formato largo

# Identificar as colunas dos conjuntos, excluindo 'Genes' e 'Count'
set_cols <- setdiff(colnames(upset_data), c("Genes", "Count"))

# Criar o Upset Plot sem o argumento set_sizes (usa o padrão)
upset_genes <- upset(
  upset_data,
  intersect = set_cols,
  name = "Genes Intersections"
)

intersection_summary_genes <- as.data.frame(intersection_summary)

write.xlsx(intersection_summary_genes, file.path("Suppl.Table1g_intersectionGenes_UpsetPlot.xlsx"))
write.xlsx(upset_data, file.path("Suppl.Table1g_intersectionGenes_position_UpsetPlot.xlsx"))

ggsave(filename = "Figure1g_upsetGenes.svg",
       plot = upset_genes,
       device = "svg", 
       width = 5, 
       height = 3, 
       units = "in", 
       dpi = 300)

# TCGA Analysis
# Enrichment Analysis (TCGA) ----------------------------------------------

#Carregando dados
enrich.data <- read.delim("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/Enrichment/GO_Biological_Process_2023_table.txt", header = TRUE, sep = "\t")
str(enrich.data)
syngo.data <- read.delim("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/Enrichment/SynGO_2024_table.txt", header = TRUE, sep = "\t")
str(syngo.data)
data.g14.lda <- read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/Enrichment/data_g14_lda_genes.xlsx")
str(data.g14.lda)

genes_filter = data.frame(genes = data.g14.lda$genes)

#Manipulando dados de Synapse
# Criar um novo data frame para armazenar as informações separadas
syngo.data_separated <- syngo.data %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"),            # Extrair o código GO
    Ontology = str_extract(Term, "\\b(CC|BP)\\b")             # Extrair a ontologia
  )
syngo.data_separated$Term <- sub("\\s*\\(.*?\\)\\s*", "", syngo.data_separated$Term)
syngo.data_separated$Term <- gsub("\\s*(BP|CC)\\s*", "", syngo.data_separated$Term)
#BubbleHeatmap
# Ordenar os dados por Combined.Score para garantir que os maiores círculos apareçam na frente
syngo.data <- syngo.data_separated[order(syngo.data_separated$Combined.Score, decreasing = TRUE), ]
# Criar o bubble heatmap

syngo_heatmap <- syngo.data

syngo_heatmap <- syngo_heatmap %>%
  dplyr::filter(P.value < 0.05)

syngo_heatmap <- syngo_heatmap %>%
  mutate(Reduced_Term = str_replace_all(Term, 
                                        c("Regulation Of " = "Reg. ",
                                          "Synaptic " = "Syn. ",
                                          "Postsynapse " = "Postsyn. ",
                                          "Transmission" = "Transm.",
                                          "Receptor" = "Rec.",
                                          "Membrane" = "Mem.",
                                          "Activity" = "Act.",
                                          "Specialization" = "Spec.",
                                          "Process" = "Proc.",
                                          "Involved In" = "In.",
                                          "Modulating" = "Mod.",
                                          "Modulation" = "Mod.",
                                          "Between" = "",
                                          "Vesicle" = "Ves.",
                                          "Presynaptic" = "Presyn. ",
                                          "Active" = "Act",
                                          "Postsynaptic" = "Postsyn.",
                                          "Potential" = "Pot."))
  )

# Ordenar o eixo y (Term) pelo menor P.value
syngo_heatmap <- syngo_heatmap %>%
  mutate(Reduced_Term = factor(Reduced_Term, levels = Reduced_Term[order(P.value)]))

all_combinations <- expand.grid(Genes = unique(syngo_heatmap$Genes),
                                Reduced_Term = unique(syngo_heatmap$Reduced_Term))
# Combinar com os dados existentes
syngo_heatmap <- merge(all_combinations,
                     syngo_heatmap, by = c("Genes",
                                          "Reduced_Term"),
                     all.x = TRUE)

p15 <- ggplot(syngo_heatmap, aes(x = Genes, y = Reduced_Term, fill = P.value)) +
  geom_tile(color = "black", linewidth = 0.1) +
  scale_fill_gradient2(
    low = "#D71313",
    mid = "#D71313",
    high = "#93BFCF",
    name = "p-value",
    na.value = "white" # Define cor para valores NA
  ) +
  labs(color = "p-value",
    x = NULL, y = NULL,
    title = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    legend.title = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, colour = "black")
  )

ggsave(filename = "syngo_TCGA_heatmap.tiff",
       plot = p15,
       device = "tiff", 
       width = 7.5, 
       height = 3, 
       units = "in", 
       dpi = 300)

#Salvando tabelas
write.xlsx(syngo_heatmap, "Suppl.Table4e_OutputSynGOHeatmap.xlsx")
write.xlsx(syngo.data, "Suppl.Table4e_InputSynGOHeatmap.xlsx")

#Unindo SynGo e EnrichR
#Manipulando dados do EnrichR
enrich.data_separated <- enrich.data %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
enrich.data_separated$Term <- sub("\\s*\\(.*?\\)\\s*", "", enrich.data_separated$Term)
enrich.data_separated$Ontology <- c("BP")
dim(enrich.data_separated)
dim(syngo.data)
#Identificando os dois df
syngo.data$df = c("SynGo")
enrich.data_separated$df = c("EnrichR")
#Unindo
enrich.total = rbind(syngo.data,
                     enrich.data_separated)

#Nomeando: neuro, imune, neuroimune
# Palavras-chave para cada categoria
neuro_keywords <- c("Synaptic", "neuro", "Synapse",
                    "Perisynaptic", "Postsynaptic",
                    "Neurotransmitter", "Neuroinflammatory",
                    "Neuron", "Postsynapse", "Presynapse",
                    "Presynaptic", "Nervous", "Dopamine", "Axon",
                    "Neurogenesis", "Adrenergic", "Dendritic", "Neural",
                    "Neuromuscular", "Neuroblast", "Sensory", "Gliogenesis",
                    "Learning", "Brain")
immune_keywords <- c("leukocyte", "T cell", "mononuclear",
                     "lymphocyte", "immune response", "B cell",
                     "interferon", "Inflammatory", "Immune", "Defense",
                     "Cytokine", "Kinase", "ERBB2", "G Protein", 
                     "Cysteine", "ERK1 And ERK2", "Interleukin",
                     "Macrophage", "Necrosis", "Tumor", "MAPK",
                     "Healing", "Humoral", "Neutrophil", "Leukocyte",
                     "Myeloid", "Response", "Negative Regulation Of Type II Interferon-Mediated Signaling Pathway",
                     "Positive Regulation Of Growth", "MHC", "Cell Growth",
                     "Toll-Like Receptor", "Monocyte", "Interferon", "STAT",
                     "Prostaglandin", "NF-kappaB","GTPase", "Apoptotic", "Viral",
                     "Modulation", "Immunoglobulin", "NIK/NF-kappaB", "B Cell",
                     "T Cell", "I-kappaB kinase/NF-kappaB", "kinase", "Death",
                     "Ubiquitin", "Ubiquitin-Dependent","Ubiquitination",
                     "Proliferation", "Growth")

# enrich.total.id <- enrich.total %>%
#   mutate(pathway = case_when(
#     str_detect(Term, paste(neuro_keywords, collapse = "|")) & 
#       !str_detect(Term, paste(immune_keywords, collapse = "|")) ~ "Nervous",
#     !str_detect(Term, paste(neuro_keywords, collapse = "|")) & 
#       str_detect(Term, paste(immune_keywords, collapse = "|")) ~ "Immune",
#     str_detect(Term, paste(neuro_keywords, collapse = "|")) & 
#       str_detect(Term, paste(immune_keywords, collapse = "|")) ~ "Neuroimmune",
#     TRUE ~ "Others"
#   ))

enrich.total.id <- enrich.total %>%
  mutate(pathway = case_when(
    str_detect(Term, paste(neuro_keywords,
                           collapse = "|")) ~ "Nervous",
    str_detect(Term, paste(immune_keywords,
                           collapse = "|")) ~ "Immune",
    TRUE ~ "Others"
  ))

enrich.total.id <- enrich.total.id %>%
  dplyr::filter(P.value < 0.05)
data.barplot.enrich <- enrich.total.id %>%
  separate_rows(Genes, sep = ";")
data.barplot.enrich <- data.barplot.enrich %>% 
  dplyr::filter(Genes %in% genes_filter$genes)

data.barplot.enrich <- data.barplot.enrich %>%
  mutate(pathway = ifelse(df == "SynGo", "Synapse", pathway))

data.barplot.enrich <- data.barplot.enrich %>%
  group_by(Genes, pathway) %>%
  summarize(SystemCount = n())

#Completando vales faltantes
all_combinations <- expand.grid(Genes = unique(data.barplot.enrich$Genes),
                                pathway = unique(data.barplot.enrich$pathway))
# Combinar com os dados existentes
merged_data <- merge(all_combinations, data.barplot.enrich, by = c("Genes", "pathway"), all.x = TRUE)
# Preencher os valores ausentes com 0 em SystemCount
merged_data[is.na(merged_data$SystemCount), "SystemCount"] <- 0
# Remover a coluna Size e Order, pois elas não são mais necessárias
merged_data <- merged_data[, !(names(merged_data) %in%
                                 c("Size", "Order"))]

#Plotando com Immune, Nervous e Synapse


# Filtrar apenas os elementos "Immune" e "Nervous" da coluna "pathway"
data.barplot.filtered <- data.barplot.enrich %>%
  dplyr::filter(pathway %in% c("Immune", "Nervous", "Synapse"))

gene_order <- c("HP", "CSF3R", "VSIG4",
                "C9", "LIFR", "TNFAIP6",
                "NLRP12", "EREG", #Immune
                "NRG1", "IL33", "GCH1", "PROK2", "OLFM1", "TRIM71",
                "SPOCK1", "DBH", "PLP1", "WDR62", "FER1L4")

# Converta a coluna 'Genes' para um fator com a ordem desejada
data.barplot.filtered$Genes <- factor(data.barplot.filtered$Genes, levels = gene_order)
# Defina os shapes desejados para cada pathway
shapes <- c("Immune" = 19, "Nervous" = 15)  # Por exemplo, círculos para Immune e triângulos para Nervous

# Plotando

p16 <- ggplot(data.barplot.filtered, aes(x = Genes, y = SystemCount, color = pathway, fill = pathway, group = pathway)) +
  geom_ribbon(aes(ymin = 0, ymax = SystemCount), alpha = 0.3) +  # Adicionando preenchimento
  geom_line() +
  geom_point(data = subset(data.barplot.filtered, SystemCount != 0),
             aes(shape = pathway),
             size = 2,
             alpha = 0.8) +  # Ajuste do alpha para um valor entre 0 e 1
  scale_color_manual(values = c("Immune" = "#3182bd", "Nervous" = "tomato", "Synapse" = "#BC7C7C")) +
  scale_fill_manual(values = c("Immune" = "#3182bd", "Nervous" = "tomato", "Synapse" = "#BC7C7C")) +
  scale_shape_manual(values = shapes) +  # Definindo os shapes manualmente
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Number of Biological Process") +
  xlab("Genes") +
  scale_x_discrete(drop = FALSE) +  # Garante que todos os níveis são exibidos, mesmo sem dados
  geom_hline(yintercept = 0, linetype = "twodash", color = "black") +
  geom_point(data = subset(data.barplot.filtered, SystemCount == 0), shape = NA)

ggsave(filename = "enrichment_TCGA_lineplot.tiff",
       plot = p16,
       device = "tiff", 
       width = 6.5, 
       height = 2.5, 
       units = "in", 
       dpi = 300)

#Barplot
# Criando todas as combinações possíveis de Genes e pathway
all_combinations <- expand.grid(Genes = unique(data.barplot.filtered$Genes),
                                pathway = unique(data.barplot.filtered$pathway))

# Mesclando com os dados existentes para garantir que todos os pares Genes e pathway tenham correspondentes
merged_data <- merge(all_combinations, data.barplot.filtered, by = c("Genes", "pathway"), all.x = TRUE)

# Se necessário, substituindo valores NA na coluna SystemCount por 0 ou outro valor adequado
merged_data$SystemCount[is.na(merged_data$SystemCount)] <- 0

gene_order <- c("C9", "CSF3R", "WDR62", "FER1L4", "GCH1", "SPOCK1",
                "PROK2", "HP", "PLP1", "OLFM1", "LIFR", "DBH", "VSIG4", "TNFAIP6", "IL33", "NRG1",
                "NLRP12", "EREG")

# Converta a coluna 'Genes' para um fator com a ordem desejada
merged_data$Genes <- factor(merged_data$Genes, levels = gene_order)

p17 <- ggplot(merged_data %>% dplyr::filter(!is.na(Genes)), aes(x = Genes, y = SystemCount, fill = pathway)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Criação do gráfico de barras
  scale_fill_manual(values = c("Immune" = "#3182bd", "Nervous" = "tomato", "Synapse" = "#BC7C7C")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("") +
  scale_x_discrete(drop = FALSE) +  # Garante que todos os níveis são exibidos, mesmo sem dados
  geom_hline(yintercept = 0, linetype = "twodash", color = "black")

ggsave(filename = "enrichment_TCGA_barplot.tiff",
       plot = p17,
       device = "tiff", 
       width = 6.5, 
       height = 2.5, 
       units = "in", 
       dpi = 300)

write.xlsx(merged_data, "Suppl.Table4d_OutBarPlot.xlsx")
write.xlsx(enrich.total.id, "Suppl.Table4d_InputBarPlot.xlsx")

#Enrich Heatmap
data_heatmap <- enrich.total.id
data_heatmap <- subset(data_heatmap, df != "SynGo")
data_heatmap <- data_heatmap %>%
  separate_rows(Genes, sep = ";")
data_heatmap <- data_heatmap %>%
  group_by(Genes) %>%
  slice_min(order_by = Adjusted.P.value, n = 1) %>%
  ungroup()

# Ordenar os termos pelo Adjusted.P.value para uma visualização mais intuitiva
data_heatmap <- data_heatmap %>%
  mutate(Term = reorder(Term, Adjusted.P.value))

# Extração da parte numérica da coluna Overlap
data_heatmap <- data_heatmap %>%
  mutate(Overlap_numeric = as.numeric(str_extract(Overlap, "^\\d+")))

#Dotplot
p18 <- ggplot(data_heatmap, aes(x = Genes, y = Term, 
                         size = Overlap_numeric, 
                         color = Adjusted.P.value)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(high = "#93BFCF", low = "#D71313", trans = "log10") + 
  scale_size_continuous(name = "Overlap", range = c(2, 10)) +  # Ajusta o tamanho dos pontos
  labs(
    title = "",
    x = "",
    y = "",
    color = "Adj.p-value"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text = element_text(size = 14, color = "black"),
    panel.spacing = unit(1.5, "lines")
  )

ggsave(filename = "enrichment_TCGA_WithoutSynapse_heatmap.tiff",
       plot = p18,
       device = "tiff", 
       width = 15, 
       height = 4, 
       units = "in", 
       dpi = 300)

write.xlsx(data_heatmap, "Suppl.Table4f_InputDotPlotEnrichment.xlsx")

# Boxplots (TCGA) ---------------------------------------------------------

library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(broom)
library(rstatix)
library(openxlsx)

setwd("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/")
# outDir <- "~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/Relative Effect/"
# dir.create(outDir)
my_data12 <- read.xlsx("Genes_GBM_G1&2.xlsx", rowNames = TRUE)
my_data13 <- read.xlsx("Genes_GBM_G1&3.xlsx", rowNames = TRUE)
my_data14 <- read.xlsx("Genes_GBM_G1&4.xlsx", rowNames = TRUE)

#Separando Graus
data_G1 <- my_data12 %>% dplyr::filter(neoplasmhistologicgrade == "G1")
data_G2 <- my_data12 %>% dplyr::filter(neoplasmhistologicgrade == "G2")
data_G3 <- my_data13 %>% dplyr::filter(neoplasmhistologicgrade == "G3")
data_G4 <- my_data14 %>% dplyr::filter(neoplasmhistologicgrade == "G4")

#Vetor de cores
group_color = c("G1" = '#FFAAA6',
                "G2" = '#DF6A6A',
                "G3" = '#BB5A5A',
                "G4" = 'darkred')
#Unindo
my_data2 = rbind(data_G1,
                 data_G2,
                 data_G3,
                 data_G4)
names(my_data2)[20] = "ID"

# Gather the data into long format
p_long <- my_data2 %>%
  pivot_longer(cols = -ID,
               names_to = "Genes",
               values_to = "Value")

names(p_long)[1] <- "Diagnostico"
names(p_long)[2] <- "Genes"
names(p_long)[3] <- "Expression"

# Gather the data into long format
p_long$Diagnostico <- factor(p_long$Diagnostico, levels = names(group_color))
p_long$Expression <- as.numeric(p_long$Expression)
p_long$Expression <- log2(p_long$Expression)
p_long$Expression <- ifelse(is.infinite(p_long$Expression) == TRUE, 0, p_long$Expression)
# Create boxplots using ggplot2 and facet_wrap
p_long <- na.omit(p_long)
#Ordenando Genes do Grafico
# p_long$Genes <- factor(p_long$Genes, levels = c("WDR62", "TRIM71", "HP", "PROK2", #SIG
#                                                 "EREG","FER1L4", "GCH1", "NLRP12", #SIG
#                                                 "OLFM1", "PLP1", "CSF3R", "NRG1",
#                                                 "C9", "IL33", "VSIG4", "LIFR",
#                                                 "DBH", "SPOCK1", "TNFAIP6"))

p_long$Genes <- factor(p_long$Genes, levels = c("TRIM71","WDR62", "HP", "FER1L4", "GCH1", "C9", "IL33", "VSIG4", "DBH", #Sig
                                                "OLFM1", "PLP1", "CSF3R", "NRG1", "PROK2", "NLRP12", "EREG", "LIFR", "SPOCK1", "TNFAIP6"))

#Teste t
# df_stat <- p_long %>%
#   group_by(Genes) %>%
#   t_test(Expression ~ Diagnostico) %>%
#   adjust_pvalue(method = 'fdr') %>%
#   add_significance() %>%
#   add_xy_position(dodge = .8, x = "Diagnostico", scales = "free")

#wilcox (analise nao parametrica entre os grupos)
df_stat <- p_long %>%
  group_by(Genes) %>%
  wilcox_test(Expression ~ Diagnostico) %>%  # Substitui t_test por wilcox_test
  adjust_pvalue(method = 'fdr') %>%         # Ajusta os valores de p usando FDR
  add_significance() %>%                    # Adiciona significância
  add_xy_position(dodge = .8, x = "Diagnostico", scales = "free")  # Define as posições no gráfico

#Kruskal
# df_stat_global <- p_long %>%
#   group_by(Genes) %>%
#   kruskal_test(Expression ~ Diagnostico) %>%
#   adjust_pvalue(method = "fdr") %>%
#   add_significance()

df_stat <- df_stat %>% 
  dplyr::filter(p.adj.signif != "ns")

df_stat <- df_stat %>%
  mutate(y.position = y.position + 0.8) 

# Criando boxplots usando ggplot2 e facet_wrap:
#Sem dots
# tiff("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/Plots/boxplot_tcga.tiff",
#      units="in", width = 18, height = 6, res = 300)
tiff("~/Boxplot_Legends.tiff", #Definir local
     units="in",
     width=10, height=8,
     res=300)
ggplot(p_long, aes(x = Diagnostico, y = Expression)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7,
               aes(fill = Diagnostico)) +
  scale_fill_manual(values = group_color) +
  stat_pvalue_manual(df_stat,
                     y.position = "y.position",
                     label = 'p.adj.signif',
                     tip.length = 0.01, # tamanho das "berinhas"
                     label.size = 4, # tamanho das estrelas
                     vjust = 0.65, # posição das estrelas
                     shape = 16) +
  facet_wrap(~ Genes, nrow = 3, scales = "free") +
  labs(x = "", y = "Expression") +
  theme_classic(18) +
  theme(axis.text.x = element_blank(), 
        strip.background = element_blank(), legend.position = "bottom")
dev.off()

write.xlsx(df_stat, file.path("Suppl.Table4g_OutputBoxPlotWilcoxTest.xlsx"))
write.xlsx(p_long, file.path("Suppl.Table4g_InputBoxPlotWilcoxTest.xlsx"))

df_stat <- read.xlsx("wilcox_test_TCGA.xlsx")

# df_stat$Genes <- factor(df_stat$Genes, levels = c("WDR62", "TRIM71", "HP", "PROK2", #SIG
#                                                   "EREG","FER1L4", "GCH1", "NLRP12", #SIG
#                                                   "OLFM1", "PLP1", "CSF3R", "NRG1",
#                                                   "C9", "IL33", "VSIG4", "LIFR",
#                                                   "DBH", "SPOCK1", "TNFAIP6"))
#Com dots
tiff("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/Plots/boxplot_tcga_dots.tiff",
     units="in", width = 10, height = 6, res = 300)

p14 <- ggplot(p_long, aes(x = Diagnostico, y = Expression)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7,
               aes(group = Diagnostico)) +
  geom_point(
    position = position_jitter(width = 0.2),
    size = 0.7,
    alpha = 0.4,
    aes(color = Diagnostico)) +  # Adding dotplots
  scale_color_manual(values = group_color)  + # Setting color manually
  stat_pvalue_manual(df_stat,
                     y.position = "y.position",
                     label = 'p.adj.signif',
                     tip.length = 0.01, #tamanho da berinha
                     label.size = 3.5, #tamanho da estrela
                     vjust = 0.65, #posicao da estrela
                     shape = 16) + 
  facet_wrap(~ Genes, ncol = 5, scales = "free") + #Aqui que se define a ordem (Levels de "Genes")
  labs(x = "", y = "Expression") +
  theme_classic(15) +
  theme(axis.text.x = element_blank(), 
        strip.background = element_blank())

p14 <- ggplot(p_long, aes(x = Diagnostico, y = Expression)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7,
               aes(group = Diagnostico),
               outlier.shape = NA) +  # Remove os dots de outliers do geom_boxplot
  geom_point(
    position = position_jitter(width = 0.2),
    size = 0.7,
    alpha = 0.4,
    aes(color = Diagnostico)) +  # Adiciona os pontos separados
  scale_color_manual(values = group_color) +  # Define manualmente as cores
  stat_pvalue_manual(df_stat,
                     y.position = "y.position",
                     label = 'p.adj.signif',
                     tip.length = 0.01,  # Tamanho das "berinhas"
                     label.size = 3.5,  # Tamanho das estrelas
                     vjust = 0.65,  # Posição das estrelas
                     shape = 16) + 
  facet_wrap(~ Genes, nrow = 2, ncol = 10, scales = "free") +  # Define a ordem (Levels de "Genes")
  labs(x = "", y = "Expression") +
  theme_classic(15) +
  theme(axis.text.x = element_blank(), 
        strip.background = element_blank())

p14 <- ggplot(p_long, aes(x = Diagnostico, y = Expression)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7,
               aes(group = Diagnostico)) +
  geom_point(
    position = position_jitter(width = 0.2),
    size = 0.5,
    alpha = 0.4,
    aes(color = Diagnostico)) +  # Adding dotplots
  scale_color_manual(values = group_color)  + # Setting color manually
  stat_pvalue_manual(df_stat,
                     y.position = "y.position",
                     label = 'p.adj.signif',
                     tip.length = 0.01, #tamanho da berinha
                     label.size = 2.5, #tamanho da estrela
                     vjust = 0.65, #posicao da estrela
                     shape = 16) + 
  facet_wrap(~ Genes, ncol = 10, scales = "free") + #Aqui que se define a ordem (Levels de "Genes")
  labs(x = "", y = "Expression") +
  theme_classic(15) +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        strip.text = element_text(size = 12, colour = "black"),
        plot.title = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

p14 <- ggplot(p_long, aes(x = Diagnostico, y = Expression)) +
  geom_point(
    position = position_jitter(width = 0.2),
    size = 0.7,
    alpha = 0.4,
    aes(color = Diagnostico)) +  # Primeiro os pontos
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7,
               aes(group = Diagnostico),
               outlier.shape = NA) +  # Depois as caixas, que vão sobrepor os pontos
  scale_color_manual(values = group_color) +
  stat_pvalue_manual(df_stat,
                     y.position = "y.position",
                     label = 'p.adj.signif',
                     tip.length = 0.01,
                     label.size = 3.5,
                     vjust = 0.65,
                     shape = 16) + 
  facet_wrap(~ Genes, nrow = 2, ncol = 10, scales = "free") +
  labs(x = "", y = "Expression") +
  theme_classic(15) +
  theme(axis.text.x = element_blank(), 
        strip.background = element_blank())

ggsave(filename = "All_Barplot_wilcox_TCGA.tiff",
       plot = p14,
       device = "tiff", 
       width = 12, 
       height = 3, 
       units = "in", 
       dpi = 300)  

ggsave(filename = "All_Barplot_wilcox_TCGA_V2.svg",
       plot = p14,
       device = "svg", 
       width = 12, 
       height = 3, 
       units = "in", 
       dpi = 300)

ggsave(filename = "All_Barplot_wilcox_TCGA_V3.svg",
       plot = p14,
       device = "svg", 
       width = 12, 
       height = 3, 
       units = "in", 
       dpi = 300)

#pvalue heatmap
data_heat <- read.xlsx("All_stat.xlsx")

# Supondo que seu data.frame 'data_heat' já esteja carregado
# Cria uma nova coluna com a transformação -log10(p.adj)
data_heat <- data_heat %>%
  mutate(p.adj.log = -log10(p.adj))

# Agrega os dados para que cada combinação única de gene e comparação tenha um único valor
# Aqui usamos a média, mas você pode escolher outra função (median, min, etc.)
data_agg <- data_heat %>%
  group_by(Genes, groups) %>%
  summarise(p.adj.log = mean(p.adj.log, na.rm = TRUE), .groups = "drop")

# Converte os dados de formato longo para wide, de forma que:
# - Cada linha seja um gene
# - Cada coluna seja uma comparação (definida na coluna 'groups')
data_wide <- pivot_wider(data_agg, names_from = groups, values_from = p.adj.log)

# Converte para data.frame e define os nomes das linhas a partir da coluna Genes
data_wide <- as.data.frame(data_wide)
rownames(data_wide) <- data_wide$Genes
data_wide$Genes <- NULL  # Remove a coluna Genes, pois ela já está nos rownames

# Define a ordem desejada para as colunas
data_wide$Control <- 0

# desired_order <- c("Control", "G3, G4", "G2, G4", "G1, G4", "G2, G3", "G1, G3", "G1, G2", "Control")
desired_order <- c("Control", "G1, G2", "G1, G3", "G2, G3", "G1, G4", "G2, G4", "G3, G4")

# Reordena as colunas do data_wide conforme essa ordem
data_wide <- data_wide[, desired_order]

colnames(data_wide) <- gsub(", ", "x", colnames(data_wide))

# Calcula a média de -log10(p.adj) para cada comparação (coluna) e para cada gene (linha)
col_means <- colMeans(data_wide, na.rm = TRUE)
row_means <- rowMeans(data_wide, na.rm = TRUE)

# Cria a anotação superior (top annotation) para as colunas
top_ann <- HeatmapAnnotation(
  Mean = anno_barplot(col_means, gp = gpar(fill = "gray80"), border = FALSE),
  annotation_name_side = "left"
)

# Cria a anotação lateral (left annotation) para as linhas
left_ann <- rowAnnotation(
  Mean = anno_barplot(row_means, gp = gpar(fill = "gray80"), border = FALSE),
  annotation_name_side = "top"
)

# Mapeando o valor mínimo para "blue" e o máximo para "red"
col_fun <- colorRamp2(c(min(as.matrix(data_wide), na.rm = FALSE), 
                        max(as.matrix(data_wide), na.rm = FALSE)), 
                      c("#93BFCF", "#D71313"))

# Plota o heatmap com as anotações superior e lateral de barras
tiff("pvalue_heatmap.tiff", #Definir local
     units="in",
     width=5, height=6,
     res=300)
Heatmap(as.matrix(data_wide),
        name = "-log10(p.adj)",
        col = col_fun,
        row_title = "Genes",
        cluster_columns = FALSE,
        top_annotation = top_ann,
        right_annotation = left_ann,
        heatmap_legend_param = list(title = "-log10(p.adj)"))
dev.off()

# Análises de receptores --------------------------------------------------
#BiocManager::install("BiocNeighbors")
devtools::install_github("jinworks/CellChat")

# library("CellChat")
library("ggalluvial")
library("writexl")
library(dplyr)
library(stringr)
library(openxlsx)
library(RColorBrewer)

# get CellChat list of ligand-receptors
# db = CellChat::CellChatDB.human
# db2 = db$interaction
db2 <- read.table("CellChatdb.tsv", header = TRUE)

db2 <- read.table("/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/CellChatdb.tsv", header = TRUE)

#write.table(db2, file = "./CellChatdb.tsv", sep = "\t", row.names = F)

# get SynGO list of synapse-related genes
# syngo = read_xlsx(path = "./../GSEbatch/syngo_genes.xlsx") #syngo
# syngo = syngo$hgnc_symbol #syngo 

G12 <- read.xlsx("G12_DESeq2Resultado.xlsx", rowNames = TRUE)
G13 <- read.xlsx("G13_DESeq2Resultado.xlsx", rowNames = TRUE)
G14 <- read.xlsx("G14_DESeq2Resultado.xlsx", rowNames = TRUE)

int_genes <- c("TRIM71",
               "WDR62",
               "HP",
               "FER1L4",
               "GCH1",
               "C9",
               "IL33",
               "VSIG4",
               "DBH",
               "NRG1",
               "OLFM1")

# data_all <- rbind(G12, G13, G14)

# data_all$Genes <- rownames(data_all)

# data_all <- subset(data_all, Genes != c("?|652919", "?|90288"))

# write.xlsx(data_all, "genes_tcga.xlsx")
# 
# write_xlsx(data_all,
#            "genes_tcga.xlsx")


# syngo = read_xlsx(path = "./../GSEbatch/syngo_genes.xlsx") #syngo
# syngo = syngo$hgnc_symbol #syngo 

# syngo <- data_all$Genes

# table for alluvial plot
tab_allu = db2[0,] 
# make columns to rbind
tab_allu = tab_allu %>% mutate(pres = NA) %>% mutate(disease = NA) %>% mutate(LR = NA)

# for (grupao in c("AD", "PD", "MS")) {
#   fname = paste0("./DEAs_",
#                  grupao,
#                  "_GSEbatch_megameta.tsv") #meta-analysis results
#   tab1 = read.table(file = fname, header = T, sep = "\t")
#   tab2 = tab1[,2:3]
#   colnames(tab2) = c("UP", "DOWN")
#   tab2 = tab2 %>%
#     pivot_longer(cols =  c("UP", "DOWN"),
#                  names_to = "direction", values_to = "Genes")
#   tab2 = tab2 %>% subset(Genes %in% syngo)

#join Genes into one string, separated by | (OR)
# genes = paste(data_all$Genes, collapse = "|")

genes <- c(HBV_synapse_genes, HCV_synapse_genes, HDV_synapse_genes)

pat = paste0("^(", genes, ")$")

# filter rows that have the genes from the meta DEgs SynGO table
tabL = db2 %>% #Pegando ligantes
  rowwise() %>% # by row
  dplyr::mutate(pres = ligand.symbol %>% # get ligand symbol (name) column
                  strsplit(split = ",") %>% # separate ligand names 
                  unlist() %>% # turn into vector to apply to separate names
                  gsub(pattern = "^ | $", replacement = "") %>%  # remove space in the beginning or end of string
                  grep(value = T, pattern = pat) %>% # get names which match the Genes from tab2 (meta DEGs and SynGO)
                  paste0(collapse = "; ")) %>% # join names, if there are multiple matches, to fit one column
  subset(pres != "") # filter table for matches only

tabR = db2 %>% #Pegando Receptores
  rowwise() %>% # by row
  dplyr::mutate(pres = receptor.symbol %>% # get receptor symbol (name) column
                  strsplit(split = ",") %>% # separate receptor names 
                  unlist() %>% # turn into vector to apply to separate names
                  gsub(pattern = "^ | $", replacement = "") %>%  # remove space in the beginning or end of string
                  grep(value = T, pattern = pat) %>% # get names which match the Genes from tab2 (meta DEGs and SynGO)
                  paste0(collapse = "; ")) %>% # join names, if there are multiple matches, to fit one column
  subset(pres != "") # filter table for matches only

# if there are results (has rows), add identification cols and append to tab_allu 
if (nrow(tabL) != 0) {
  # tabL$disease = grupao
  tabL$LR = "LIGAND"
  tab_allu = rbind(tab_allu, tabL)
}

if (nrow(tabR) != 0) {
  # tabR$disease = grupao
  tabR$LR = "RECEPTOR"
  tab_allu = rbind(tab_allu, tabR)
}

# reorder columns
tab_allu2 = tab_allu %>%
  relocate(pathway_name, ligand.symbol, receptor.symbol, 
           pres, 
           # disease,
           LR, is_neurotransmitter,
           .before = 1)
# add count for each line
tab_allu2$num_rows = 1

# # change names temporarily to avoid repeating in the different axes  
# {rm = which(tab_allu2$receptor.symbol == "NCAM2")
#   tab_allu2[rm, "receptor.symbol"] = "NCAM2!!!!!"
#   
#   rm = which(tab_allu2$pathway_name == "ADM")
#   tab_allu2[rm, "pathway_name"] = "ADM!!!!!"
#   
#   rm = which(tab_allu2$pathway_name == "MIF")
#   tab_allu2[rm, "pathway_name"] = "MIF!!!!!"
#   
#   rm = which(tab_allu2$pathway_name == "SELE")
#   tab_allu2[rm, "pathway_name"] = "SELE!!!!!"
#   
#   rm = which(tab_allu2$pathway_name == "SPP1")
#   tab_allu2[rm, "pathway_name"] = "SPP1!!!!!"
#   
#   rm = which(tab_allu2$pathway_name == "IAPP")
#   tab_allu2[rm, "pathway_name"] = "IAPP!!!!!"
#   
#   rm = which(tab_allu2$pathway_name == "FN1")
#   tab_allu2[rm, "pathway_name"] = "FN1!!!!!"
# }
# 
# # order receptors
# ar = levels(factor(tab_allu2$receptor.symbol))
# ar2 = c("HTR5A","CALCR, RAMP1","NR3C1","CD44","CD44, CD74",
#         "FPR1","FPR2","FPR2, FPR3","FPR3",
#         "EPB41L1","EPHA7", ar[30:44], #GABAA
#         "GABBR1","GABBR2", ar[46:70], #Glutamate
#         "GLRA2","NCAM2!!!!!","L1CAM","NGFR","NTRK2","SORT1",
#         ar[79:81], "LRP5, FZD3","LRP5, FZD7", ar[10:27] #WNT
# )
# 
# # order ligands
# al = levels(factor(tab_allu2$ligand.symbol))
# al2 = c("SLC18A1, TPH1", "SLC18A1, TPH2","SLC18A2, TPH1", "SLC18A2, TPH2","SLC6A4, TPH1","SLC6A4, TPH2", # 5-HT
#         "ADM","ANXA1", 
#         "CADM3", "CALCA", "CALCB","CYP11B1","CYP11B2", al[7:23],"DDC","IAPP", "FN1", al[27:31], al[33:42], al[44:55],"MIF", "SELE",
#         al[63:68], "SPP1", "SLC1A6, GLS","SLC1A6, GLS2",# Glutamate
#         "SHMT1, SLC6A9","SHMT2, SLC6A9", "SLC6A5, SHMT1", "SLC6A5, SHMT2", # Glycine
#         "NCAM1","NCAM2", "BDNF", "LRFN4","TAFA4","WNT3A" 
# )
# 
# # reorder
# tab_allu2 = tab_allu2 %>%
#   mutate(disease = factor(disease, levels = c("AD", "PD", "MS"))) %>%
#   mutate(LR = factor(LR, levels = c("LIGAND", "RECEPTOR"))) %>% 
#   mutate(receptor.symbol = factor(receptor.symbol, levels = ar2)) %>% 
#   mutate(ligand.symbol = factor(ligand.symbol, levels = al2)) 
# 
# # remove repetitive ligands
# rm = which(tab_allu2$ligand.symbol == "GAD2, SLC6A6")
# tab_allu2 = tab_allu2[-rm, ]
# rm = which(tab_allu2$disease == "PD" & tab_allu2$ligand.symbol == "SLC1A6, GLS2")
# tab_allu2 = tab_allu2[-rm, ]

# set color palette
myPalette <- colorRampPalette(rev(brewer.pal(7, "Accent")))

## plot
ggplot(data = tab_allu2, aes(axis1 = pathway_name,
                             axis2 = ligand.symbol,
                             axis3 = receptor.symbol,
                             y = num_rows,
                             fill = pathway_name)) +
  #geom_flow(aes(fill = LR), order = LR) +
  stat_alluvium(aes.bind = "flows",
                lode.guidance = "zagzig") +
  #stat_flow( aes.bind = "flows") + #stat_stratum() +
  #geom_stratum(aes.bind = T) + #geom_stratum(aes(fill = LR)) + #geom_flow() +
  geom_stratum(aes(color = LR),
               size = 0.3,
               alpha = 1,
               width = 1/2.5) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 2.8) +
  theme(panel.background = element_blank(), 
        axis.title = element_blank(),
        axis.text =  element_blank(),
        axis.ticks =  element_blank(), 
        strip.background =  element_blank() ) + 
  # facet_grid(vars(disease), scales = "free_y", space = "free") +
  scale_fill_manual(values = myPalette(3))

# save file (!!!)
ggsave(filename = "../../../figures/alluvial/cellchat5.svg", 
       width = 2500, height = 5000,  #width = 4400, height = 1500, 
       units = "px", dpi = 300, device = "svg", scale = 1)

#Com cores escolhidas
# Define manualmente as duas cores desejadas
minhasCores <- c("#8EACCD", "#C96868", "gray80")

## plot
p1 <- ggplot(data = tab_allu2, aes(axis1 = pathway_name,
                                   axis2 = ligand.symbol,
                                   axis3 = receptor.symbol,
                                   y = num_rows,
                                   fill = pathway_name)) +
  stat_alluvium(aes.bind = "flows",
                lode.guidance = "zagzig") +
  geom_stratum(color = "black",
               size = 0.3,
               alpha = 1,
               width = 1/2.5) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 2.8) +
  theme(panel.background = element_blank(), 
        axis.title = element_blank(),
        axis.text =  element_blank(),
        axis.ticks =  element_blank(), 
        strip.background =  element_blank()) + 
  scale_fill_manual(values = minhasCores)

ggsave(filename = "Alluvial_pathway.tiff",
       plot = p1,
       device = "tiff", 
       width = 7, 
       height = 4, 
       units = "in", 
       dpi = 300)

#Com genes da LDA
# Definir genes e seus grupos
genes_list <- list(
  "HBV Liver" = c("TRIM71", "C9", "TNFSF11", "GRIN2B", "AXL", "SPOCK1", "NLRP12", "GCH1", "LIFR", "VXN"),
  "HBV PBMC" = c("CCL18", "TNFAIP6", "CCR3", "INHBB", "LHX4", "OLFM1", "VSIG4", "PADI2", "B3GNT5", "FER1L4"),
  "HCV Liver" = c("ADAM22", "RIMBP3C", "CSF3R", "NRG1", "CCN2", "STMN3", "DBH", "PLP1", "BCAN", "EPHB2"),
  "HCV PBMC" = c("FAM83D", "CYBB", "ITGB1", "IRF1", "NOG", "TIAM2", "TPST1", "OLFM1", "PROK2", "CXCL16"),
  "HDV Liver" = c("WDR62", "EREG", "HP", "IRF8", "PTGER2", "PTGIS", "TASL", "XCL1", "MARCO", "IL33")
)

# Criar um dataframe unificado para int_genes
df_genes <- do.call(rbind, lapply(names(genes_list), function(grupo) {
  data.frame(Gene = genes_list[[grupo]], Group = grupo)
}))

# Criar padrão de busca para filtrar os genes desejados
genes <- paste(df_genes$Gene, collapse = "|")
pat <- paste0("^(", genes, ")$")

# Criar tabela alluvial vazia
tab_allu <- db2[0, ] %>% mutate(pres = NA, group = NA, LR = NA)

# Filtrar genes que aparecem como Ligantes
tabL <- db2 %>%
  rowwise() %>%
  mutate(pres = ligand.symbol %>%
           strsplit(split = ",") %>%
           unlist() %>%
           gsub(pattern = "^ | $", replacement = "") %>%
           grep(value = T, pattern = pat) %>%
           paste0(collapse = "; ")) %>%
  subset(pres != "")

# Filtrar genes que aparecem como Receptores
tabR <- db2 %>%
  rowwise() %>%
  mutate(pres = receptor.symbol %>%
           strsplit(split = ",") %>%
           unlist() %>%
           gsub(pattern = "^ | $", replacement = "") %>%
           grep(value = T, pattern = pat) %>%
           paste0(collapse = "; ")) %>%
  subset(pres != "")

# Adicionar coluna de grupo com base nos genes correspondentes
tabL$group <- df_genes$Group[match(tabL$pres, df_genes$Gene)]
tabR$group <- df_genes$Group[match(tabR$pres, df_genes$Gene)]

# Adicionar identificadores de Ligand e Receptor
if (nrow(tabL) != 0) {
  tabL$LR <- "LIGAND"
  tab_allu <- rbind(tab_allu, tabL)
}

if (nrow(tabR) != 0) {
  tabR$LR <- "RECEPTOR"
  tab_allu <- rbind(tab_allu, tabR)
}

# Ajustar colunas e reordenar
tab_allu2 <- tab_allu %>%
  relocate(group, pathway_name, ligand.symbol, receptor.symbol, pres, LR, is_neurotransmitter, .before = 1)

tab_allu2_expanded <- tab_allu2 %>%
  separate_rows(receptor.symbol, sep = ", ") %>%  # Separar receptores
  separate_rows(ligand.symbol, sep = ", ")   # Separar ligantes

# Adicionar contagem de linhas
tab_allu2_expanded$num_rows <- 1

# Definir cores para os grupos
# Remover redundâncias agrupando por grupo e pathway
tab_allu2_clean <- tab_allu2_expanded %>%
  distinct(group, pathway_name, ligand.symbol, receptor.symbol, .keep_all = TRUE) %>%
  mutate(num_rows = 1)  # Resetar a contagem

tab_allu2_clean <- tab_allu2_clean %>%
  filter(!is.na(group), !is.na(pathway_name), !is.na(ligand.symbol), !is.na(receptor.symbol))

# Definir cores mais contrastantes para os grupos
custom_palette <- c("HBV Liver" = "#DDF1F5", 
                    "HBV PBMC" = "#789DBC", 
                    "HCV Liver" = "#EDAEAE", 
                    "HCV PBMC" = "#BC7C7C", 
                    "HDV Liver" = "#FFD09B")

# Definir cores personalizadas para Ligand/Receptor
# LR_palette <- c("LIGAND" = "#D55E00",  # Vermelho
#                 "RECEPTOR" = "#0072B2",  # Azul
#                 "NA" = "black")  # Preto para valores ausentes

# Remover linhas onde `group`, `pathway_name`, `ligand.symbol` ou `receptor.symbol` sejam NA
tab_allu2_clean <- tab_allu2_clean %>%
  mutate(group = ifelse(is.na(group), "NA", group)) %>%  # Substituir NA por string "NA"
  mutate(pathway_name = ifelse(is.na(pathway_name), "NA", pathway_name)) %>%
  mutate(ligand.symbol = ifelse(is.na(ligand.symbol), "NA", ligand.symbol)) %>%
  mutate(receptor.symbol = ifelse(is.na(receptor.symbol), "NA", receptor.symbol))

custom_palette <- c("HBV Liver" = "#DDF1F5", 
                    "HBV PBMC" = "#789DBC", 
                    "HCV Liver" = "#EDAEAE", 
                    "HCV PBMC" = "#BC7C7C", 
                    "HDV Liver" = "#FFD09B",
                    "NA" = "#BC7C7C")

# Criar gráfico alluvial corrigido
p2 <- ggplot(data = tab_allu2_clean, aes(
  axis1 = pathway_name, 
  axis2 = ligand.symbol, 
  axis3 = receptor.symbol,
  y = num_rows,
  fill = group)) +
  stat_alluvium(aes(bind = "flows"), lode.guidance = "forward", na.rm = TRUE) +  # Melhor fluxo, removendo NA
  geom_stratum(aes(fill = group, color = group),  
               size = 0.3, alpha = 1, width = 1/2.5, na.rm = TRUE) +  # Barras mais visíveis, sem NA
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 2.8,
            check_overlap = TRUE, na.rm = TRUE) +  # Labels sem NA
  scale_fill_manual(values = custom_palette, na.value = "#BC7C7C", drop = TRUE) +  # Cor preta para NA
  scale_color_manual(values = custom_palette, na.value = "#BC7C7C", drop = TRUE) +  # Cor preta para NA
  theme_void() +  # Remover elementos desnecessários
  theme(legend.position = "right",
        legend.title = element_blank())  # Ajustar legenda

ggsave(filename = "Alluvial_pathway_LDA.tiff",
       plot = p2,
       device = "tiff", 
       width = 6, 
       height = 18, 
       units = "in", 
       dpi = 300)

#Alluvial filtrado
tab_allu2_clean_filt <- tab_allu2_clean[tab_allu2_clean$is_neurotransmitter == TRUE, ]

# ggplot(data = tab_allu2_clean_filt, aes(
#   axis1 = pathway_name, 
#   axis2 = ligand.symbol, 
#   axis3 = receptor.symbol,
#   y = num_rows, fill = group)) +
#   stat_alluvium(aes(bind = "flows"), lode.guidance = "forward", na.rm = TRUE) +  # Melhor fluxo, removendo NA
#   geom_stratum(aes(fill = group, color = group),  
#                size = 0.3, alpha = 1, width = 1/2.5, na.rm = TRUE) +  # Barras mais visíveis, sem NA
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum)),
#             size = 2.8,
#             check_overlap = TRUE, na.rm = TRUE) +  # Labels sem NA
#   scale_fill_manual(values = custom_palette, na.value = "#BC7C7C", drop = TRUE) +  # Cor preta para NA
#   scale_color_manual(values = custom_palette, na.value = "#BC7C7C", drop = TRUE) +  # Cor preta para NA
#   theme_void() +  # Remover elementos desnecessários
#   theme(legend.position = "right",
#         legend.title = element_blank())  # Ajustar legenda

custom_palette <- c("HBV Liver" = "#DDF1F5", 
                    "HBV PBMC" = "#789DBC", 
                    "HCV Liver" = "#EDAEAE", 
                    "HCV PBMC" = "#BC7C7C", 
                    "HDV Liver" = "#FFD09B")

p3 <- ggplot(data = tab_allu2_clean_filt, aes(axis1 = group,
  axis2 = pathway_name,
                             axis3 = ligand.symbol,
                             axis4 = receptor.symbol,
                             y = num_rows,
                             fill = group)) +
  stat_alluvium(aes.bind = "flows",
                lode.guidance = "zagzig") +
  geom_stratum(color = "black",
               size = 0.3,
               alpha = 1,
               width = 1/2.5) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 2.5) +
  theme(panel.background = element_blank(), 
        axis.title = element_blank(),
        axis.text =  element_blank(),
        axis.ticks =  element_blank(), 
        strip.background =  element_blank()) + 
  scale_fill_manual(values = custom_palette, na.value = "gray90", drop = TRUE)

ggsave(filename = "Alluvial_pathway_LDA_Neurotransmitter.tiff",
       plot = p3,
       device = "tiff", 
       width = 10, 
       height = 7, 
       units = "in", 
       dpi = 300)

#Barplot
# Contabilizar elementos únicos para cada categoria dentro de cada grupo
df_counts <- tab_allu2_clean %>%
  select(group, pathway_name, ligand.symbol, receptor.symbol) %>%
  pivot_longer(cols = c(pathway_name, ligand.symbol, receptor.symbol),
               names_to = "type", values_to = "element") %>%
  distinct(group, type, element) %>%
  count(group, type)

df_counts <- df_counts %>%
  mutate(type = recode(type, "receptor.symbol" = "Receptor",
                       "ligand.symbol" = "Ligand",
                       "pathway_name" = "Pathway"))

df_counts$type <- factor(df_counts$type, levels = c("Pathway", "Ligand", "Receptor"))

# Criar o barplot com a ordem corrigida no eixo X
p_barplot <- ggplot(df_counts, aes(x = type, y = n, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barras lado a lado
  geom_text(aes(label = n), vjust = -0.3, size = 4) +  # Rótulos nas barras
  facet_wrap(~ group, scales = "free_y") +  # Facetado por `group`
  scale_fill_manual(values = c("Pathway" = "gray80", 
                               "Ligand" = "#FFB4A2", 
                               "Receptor" = "#B5828C")) +  # Cores personalizadas
  labs(x = "Category", y = "Unique Count", fill = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "bottom")

# Remover duplicatas dentro de cada pathway
df_filtered_unique <- df_filtered %>%
  distinct(pathway_name, ligand, receptor, .keep_all = TRUE)

df_counts <- df_filtered_unique %>%
  pivot_longer(cols = c(ligand,
                        receptor),
               names_to = "type",
               values_to = "element") %>%
  count(pathway_name, type, element)

# Modificar os nomes na coluna "type"
df_counts <- df_counts %>%
  mutate(type = recode(type, "receptor" = "Receptor",
                       "ligand" = "Ligand"))

# Definir cores personalizadas para "ligand" e "receptor"
cores_personalizadas <- c("Ligand" = "#FFB4A2", "Receptor" = "#B5828C")  # Azul e Vermelho

# Definir a ordem desejada para o eixo X
df_counts <- df_counts %>%
  mutate(pathway_name = factor(pathway_name, levels = c("Noradrenaline", "NRG", "IL1")))

# Criar o gráfico com a nova ordem no eixo X
tiff("counting_ligand_receptor.tiff",
     width = 2500, height = 1700, res = 300)
ggplot(df_counts, aes(x = pathway_name, y = n, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barras lado a lado
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = 1.1, size = 8) +  # Rótulos nas barras
  facet_wrap(~ type, scales = "free_x") +  # Facetas separadas para ligand/receptor
  scale_fill_manual(values = cores_personalizadas) +  # Aplicar cores manuais
  scale_y_continuous(breaks = c(0, 5, 9)) +  # Definir quebras no eixo Y
  theme_bw() +
  labs(x = "", y = "Frequency", fill = "Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # Texto do eixo X
    axis.text.y = element_text(size = 20),  # Texto do eixo Y
    axis.title.y = element_text(size = 20),  # Título do eixo Y
    strip.text = element_text(size = 20)  # Tamanho do texto das facetas
  )
dev.off()

#CircosPlot

# SYNAPSE CIRCOSPLOT ------------------------------------------------------

library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(grid)
library(stringr)

HBV_liver_synapse_genes <- c("KCNC1", "FLRT2", "CACNA1A", "TRPV1", "ADGRL3", "SHH",
                             "SEPTIN5", "DCX", "SEPTIN3", "CPLX2", "GRM7", "ACAN",
                             "INSYN1", "ACP4", "TAFA4", "RAB3B", "PLPPR4")

HBV_pbmc_synapse_genes <- c("C1QB", "C1QC", "NRN1", "AMPH", "SYN2", "AKAP12", "NXPH4",
                            "LRRTM1", "CACNA1E", "MMP9", "IGSF21", "NSG1", "OLFM1",
                            "DPYSL3", "ZDHHC8")

HCV_liver_synapse_genes <- c("CLSTN2", "ADGRB3", "NRG1", "ADRA2A", "GPR158", "GRK3",
                             "NPTX2", "EPHA3", "VCAN", "CNTN1")

HCV_pbmc_synapse_genes <- c("MOB4", "FMNL2", "PLK2", "PDXP", "RAPGEF2", "PHACTR1",
                            "CPEB2", "ADD2", "RPS15", "FAM171B", "VLDLR", "YBX1",
                            "PABPC1", "RPL10A", "ITGB1", "SYNGAP1", "HOMER3", "RBMX",
                            "ACHE", "TMEM163")

HDV_synapse_genes <- c("CHRM2", "GRIA2", "NGFR", "GABBR2", "KCNJ6", "KCNC1", "CHRNA4",
                       "GABRA5", "NRXN1", "NPY1R", "GRIK3", "HTR2A", "SLC1A7", "GRIN2B",
                       "ANO1", "GRIN2A", "EFNB3", "P2RX3", "SCN11A", "PODXL", "FXYD6",
                       "NCAM1", "EPHB2", "GRIA3", "CADM2", "NTRK3", "ABHD6", "SEZ6L",
                       "GABRR3", "SLC6A8", "KCNQ3", "SUSD4", "CSPG5", "KCNN2", "GABRE",
                       "PTPRT", "ROBO2", "PLPPR4", "LRRC4", "ADAM22", "NRG1", "ADCY1",
                       "SHISA9", "ELFN2", "CELSR3", "LRRTM3", "IL1RAPL1", "ASIC2", "CSMD2",
                       "GPR158", "CYFIP2", "FGB", "GNAZ", "FGA", "FGG", "IL1RAP", "TDRD5",
                       "CLU", "ELAVL2", "HAPLN1", "TDRD1", "RAP2A", "VCAN", "DAB1", "NWD2",
                       "ADGRB2", "ABLIM3", "GNG4", "CDH11", "CDH23", "CDH13", "PABPC1",
                       "BCAS1", "RASGRF1", "ARHGAP39", "THY1", "COMT", "HTR4", "CPLX2",
                       "ENAH", "GNAL", "SNX27", "MAP2", "PRKAR2B", "MAP1B", "TSPAN7",
                       "AGO2", "ABL2", "TRIM47", "KIF2C", "ROR2", "ARHGEF2", "CPEB3",
                       "GPM6A", "LINGO1", "SEMA7A", "NPY5R", "ITGA2", "GAD1", "PRUNE2",
                       "SEPTIN3", "SIPA1L2", "SEPTIN8", "GAP43", "SYT10", "DCX", "SPTBN2",
                       "TUBB2B", "PLCB4", "EPHB1", "DLGAP2", "CHRNA3", "NRXN3", "APBA1",
                       "DGKI", "NR3C2", "CAMK2B", "MAP2K1", "HOMER1", "HOMER2", "CTNND2",
                       "AURKA", "DLG5", "TACC3", "SPOCK1", "KPNA2", "SOS1", "CD200",
                       "EIF4G2", "DBN1", "TENM3", "SPARCL1", "IL1RAPL2", "YWHAZ",
                       "CTTNBP2", "SLITRK4", "CDH8", "ACTN2", "CTNNA2", "MMP9", "C1QB",
                       "SNAP25", "C1QC", "FRRS1L", "KIF21B", "LGI1", "LAMA4", "LAMC1",
                       "CBLN2", "BCAN", "ITPKA", "RIMS2", "ERC2", "LHFPL4", "SCAMP5",
                       "SLC30A3", "SLC17A8", "SYT9", "TMEM163", "GNAO1", "IGF1", "PLCB1",
                       "HTR1D", "STXBP5L", "KCND3", "SCRIB", "RAB3B", "RAB3C", "RAB26",
                       "SORT1", "PLXNA4", "ACAN", "NTM", "CACNA2D2", "UNC13A", "LINGO2",
                       "PCDH9", "TAFA4", "VPS45", "DNAJC6", "SIPA1L3", "PLXNC1", "ATP6V1C1")

genes_list <- list(
  "HBV_liver_synapse_genes" = c(HBV_liver_synapse_genes),
  "HBV_pbmc_synapse_genes" = c(HBV_pbmc_synapse_genes),
  "HCV_liver_synapse_genes" = c(HCV_liver_synapse_genes),
  "HCV_pbmc_synapse_genes" = c(HCV_pbmc_synapse_genes),
  "HDV_synapse_genes" = c(HDV_synapse_genes)
)

# Neurotransmissores ------------------------------------------------------
# ======================
# 0. Preparação
# ======================

# Lista branca de pathways neuroativos
neuro_pathways <- c("GABA-B", "Glutamate",
                    "Dopamine","NRG", "NRXN")

# Paleta de cores associada
cores_pathways_base <- c(
  "GABA-B" = "#F08080",
  "Glutamate" = "#93BFCF",
  "Dopamine" = "#e5c185",
  "NRG" = "gray70",
  "NRXN" = "#008585"
)

# Filtrar base para pathways relevantes
db3 <- db2 %>% filter(pathway_name %in% neuro_pathways)

# ======================
# 1. Loop por grupo de genes
# ======================
for (grupo in names(genes_list)) {
  
  cat("⏳ Processando:", grupo, "\n")
  
  my_genes <- unique(genes_list[[grupo]])
  
  tab_allu <- db3[0, ] %>%
    mutate(
      LR = character(),
      ligand.symbol = character(),
      receptor.symbol = character(),
      interaction_name = character(),
      pathway_name = character(),
      is_neurotransmitter = logical()
    )
  pres <- data.frame(pres = character())
  
  # Ligantes
  dbL <- db3 %>%
    separate_rows(ligand.symbol, sep = ",\\s*") %>%
    filter(ligand.symbol %in% my_genes)
  presL <- dbL %>% distinct(ligand.symbol) %>% rename(pres = ligand.symbol)
  tabL <- db3 %>% filter(interaction_name %in% dbL$interaction_name)
  
  # Receptores
  dbR <- db3 %>%
    separate_rows(receptor.symbol, sep = ",\\s*") %>%
    filter(receptor.symbol %in% my_genes)
  presR <- dbR %>% distinct(receptor.symbol) %>% rename(pres = receptor.symbol)
  tabR <- db3 %>% filter(interaction_name %in% dbR$interaction_name)
  
  if (nrow(tabL) > 0) {
    tabL$LR <- "LIGAND"
    tab_allu <- bind_rows(tab_allu, tabL)
    pres <- bind_rows(pres, presL)
  }
  if (nrow(tabR) > 0) {
    tabR$LR <- "RECEPTOR"
    tab_allu <- bind_rows(tab_allu, tabR)
    pres <- bind_rows(pres, presR)
  }
  
  if (nrow(tab_allu) == 0) {
    cat("⚠️  Nenhuma interação encontrada para", grupo, "- pulando...\n")
    next
  }
  
  # ======================
  # 2. Preparar dados
  # ======================
  tab_allu2 <- tab_allu %>%
    relocate(pathway_name, ligand.symbol, receptor.symbol, LR, is_neurotransmitter) %>%
    mutate(num_rows = 1) %>%
    separate_rows(ligand.symbol, sep = ",\\s*") %>%
    separate_rows(receptor.symbol, sep = ",\\s*") %>%
    distinct(pathway_name, ligand.symbol, receptor.symbol, LR, num_rows)
  
  # Agrupar genes em famílias (colapsar nomes)
  collapse_gene <- function(gene) {
    case_when(
      str_detect(gene, "^GAB[AR]") ~ "GABA",
      str_detect(gene, "^HTR") ~ "5-HT",
      str_detect(gene, "^GRIA|^GRIN") ~ "Glutamate",
      str_detect(gene, "^GRM[1-8]$") ~ "GRM",
      str_detect(gene, "^GRID[1-2]$") ~ "GRID",
      
      str_detect(gene, "^SLC6A|^SLC18A|^CHAT") ~ "Ach",
      
      str_detect(gene, "^GLS2?$") ~ "GLS",
      str_detect(gene, "^CBLN") ~ "CBLN",
      
      str_detect(gene, "^NRXN") ~ "NRXN",
      
      str_detect(gene, "^NLGN[123]$") ~ "NLGN",
      str_detect(gene, "^CLSTN[123]$") ~ "CLSTN",
      str_detect(gene, "^LRRTM[1234]$") ~ "LRRTM",
      
      str_detect(gene, "^ERBB[234]$") ~ "ERBB",
      str_detect(gene, "^ITGA") ~ "ITGA",
      str_detect(gene, "^ITGB") ~ "ITGB",
      
      str_detect(gene, "^DRD|^DOPA") ~ "Dopamine",
      TRUE ~ gene
    )
  }
  
  tab_allu3 <- tab_allu2 %>%
    mutate(
      ligand = collapse_gene(ligand.symbol),
      receptor = collapse_gene(receptor.symbol)
    ) %>%
    group_by(ligand, receptor, pathway_name) %>%
    summarise(n = n(), .groups = "drop")
  
  # ======================
  # 3. Filtrar para pathways reais do grupo e definir cores
  # ======================
  pathways_presentes <- unique(tab_allu3$pathway_name)
  cores_pathways <- cores_pathways_base[pathways_presentes]
  
  tab_allu3 <- tab_allu3 %>%
    mutate(cores = cores_pathways[pathway_name])
  
  # ======================
  # 4. Cores por nó e ordem
  # ======================
  grid.col <- tab_allu3 %>%
    pivot_longer(cols = c(ligand, receptor)) %>%
    distinct(value, cores) %>%
    deframe()
  
  orderL <- tab_allu3 %>% arrange(pathway_name) %>% distinct(ligand)
  orderR <- tab_allu3 %>% arrange(pathway_name) %>% distinct(receptor)
  order1 <- c(rev(orderL$ligand), orderR$receptor)
  
  # ======================
  # 5. Salvar gráfico
  # ======================
  output_file <- paste0(
    "/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/5 - SynGO/",
    "circos_neuroimmune_", gsub(" ", "_", grupo), ".pdf"
  )
  
  pdf(output_file, width = 6.5, height = 6.5)
  
  circos.clear()
  chordDiagram(tab_allu3,
               grid.col = grid.col,
               preAllocateTracks = 1,
               order = order1,
               annotationTrack = "grid",
               annotationTrackHeight = mm_h(2))
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    name <- get.cell.meta.data("sector.index")
    
    circos.text(mean(xlim), ylim[1] + 0.1, name,
                cex = 0.8, facing = "clockwise", niceFacing = TRUE,
                font = ifelse(name %in% pres$pres, 2, 1),
                adj = c(0, 0.5))
  }, bg.border = NA)
  
  # ======================
  # 6. Legenda
  # ======================
  lgd <- Legend(
    at = names(cores_pathways),
    legend_gp = gpar(fill = cores_pathways),
    title = "Pathway",
    title_position = "topleft"
  )
  draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
  
  dev.off()
  
  cat("✅ Circos salvo:", output_file, "\n")
}

# Neurotransmissores e outros ---------------------------------------------

# ================================
# 1. Definir listas de genes
# ================================
HBV_liver_synapse_genes <- c("KCNC1", "FLRT2", "CACNA1A", "TRPV1", "ADGRL3", "SHH",
                             "SEPTIN5", "DCX", "SEPTIN3", "CPLX2", "GRM7", "ACAN",
                             "INSYN1", "ACP4", "TAFA4", "RAB3B", "PLPPR4")

HBV_pbmc_synapse_genes <- c("C1QB", "C1QC", "NRN1", "AMPH", "SYN2", "AKAP12", "NXPH4",
                            "LRRTM1", "CACNA1E", "MMP9", "IGSF21", "NSG1", "OLFM1",
                            "DPYSL3", "ZDHHC8")

HCV_liver_synapse_genes <- c("CLSTN2", "ADGRB3", "NRG1", "ADRA2A", "GPR158", "GRK3",
                             "NPTX2", "EPHA3", "VCAN", "CNTN1")

HCV_pbmc_synapse_genes <- c("MOB4", "FMNL2", "PLK2", "PDXP", "RAPGEF2", "PHACTR1",
                            "CPEB2", "ADD2", "RPS15", "FAM171B", "VLDLR", "YBX1",
                            "PABPC1", "RPL10A", "ITGB1", "SYNGAP1", "HOMER3", "RBMX",
                            "ACHE", "TMEM163")

HDV_synapse_genes <- c("CHRM2", "GRIA2", "NGFR", "GABBR2", "KCNJ6", "KCNC1", "CHRNA4",
                       "GABRA5", "NRXN1", "NPY1R", "GRIK3", "HTR2A", "SLC1A7", "GRIN2B",
                       "ANO1", "GRIN2A", "EFNB3", "P2RX3", "SCN11A", "PODXL", "FXYD6",
                       "NCAM1", "EPHB2", "GRIA3", "CADM2", "NTRK3", "ABHD6", "SEZ6L",
                       "GABRR3", "SLC6A8", "KCNQ3", "SUSD4", "CSPG5", "KCNN2", "GABRE",
                       "PTPRT", "ROBO2", "PLPPR4", "LRRC4", "ADAM22", "NRG1", "ADCY1",
                       "SHISA9", "ELFN2", "CELSR3", "LRRTM3", "IL1RAPL1", "ASIC2", "CSMD2",
                       "GPR158", "CYFIP2", "FGB", "GNAZ", "FGA", "FGG", "IL1RAP", "TDRD5",
                       "CLU", "ELAVL2", "HAPLN1", "TDRD1", "RAP2A", "VCAN", "DAB1", "NWD2",
                       "ADGRB2", "ABLIM3", "GNG4", "CDH11", "CDH23", "CDH13", "PABPC1",
                       "BCAS1", "RASGRF1", "ARHGAP39", "THY1", "COMT", "HTR4", "CPLX2",
                       "ENAH", "GNAL", "SNX27", "MAP2", "PRKAR2B", "MAP1B", "TSPAN7",
                       "AGO2", "ABL2", "TRIM47", "KIF2C", "ROR2", "ARHGEF2", "CPEB3",
                       "GPM6A", "LINGO1", "SEMA7A", "NPY5R", "ITGA2", "GAD1", "PRUNE2",
                       "SEPTIN3", "SIPA1L2", "SEPTIN8", "GAP43", "SYT10", "DCX", "SPTBN2",
                       "TUBB2B", "PLCB4", "EPHB1", "DLGAP2", "CHRNA3", "NRXN3", "APBA1",
                       "DGKI", "NR3C2", "CAMK2B", "MAP2K1", "HOMER1", "HOMER2", "CTNND2",
                       "AURKA", "DLG5", "TACC3", "SPOCK1", "KPNA2", "SOS1", "CD200",
                       "EIF4G2", "DBN1", "TENM3", "SPARCL1", "IL1RAPL2", "YWHAZ",
                       "CTTNBP2", "SLITRK4", "CDH8", "ACTN2", "CTNNA2", "MMP9", "C1QB",
                       "SNAP25", "C1QC", "FRRS1L", "KIF21B", "LGI1", "LAMA4", "LAMC1",
                       "CBLN2", "BCAN", "ITPKA", "RIMS2", "ERC2", "LHFPL4", "SCAMP5",
                       "SLC30A3", "SLC17A8", "SYT9", "TMEM163", "GNAO1", "IGF1", "PLCB1",
                       "HTR1D", "STXBP5L", "KCND3", "SCRIB", "RAB3B", "RAB3C", "RAB26",
                       "SORT1", "PLXNA4", "ACAN", "NTM", "CACNA2D2", "UNC13A", "LINGO2",
                       "PCDH9", "TAFA4", "VPS45", "DNAJC6", "SIPA1L3", "PLXNC1", "ATP6V1C1")

# ================================
# 2. Filtrar db2 por genes presentes em ligand ou receptor
# ================================

# Função para filtrar base por uma lista de genes
filtrar_por_genes <- function(db, genes) {
  db %>%
    separate_rows(ligand.symbol, sep = ",\\s*") %>%
    separate_rows(receptor.symbol, sep = ",\\s*") %>%
    filter(ligand.symbol %in% genes | receptor.symbol %in% genes) %>%
    distinct()
}

# Aplicar para cada condição
db2_HBV_liver <- filtrar_por_genes(db2, HBV_liver_synapse_genes)
db2_HBV_pbmc <- filtrar_por_genes(db2, HBV_pbmc_synapse_genes)
db2_HCV_liver <- filtrar_por_genes(db2, HCV_liver_synapse_genes)
db2_HCV_pbmc <- filtrar_por_genes(db2, HCV_pbmc_synapse_genes)
db2_HDV <- filtrar_por_genes(db2, HDV_synapse_genes)


# Circos 2 ----------------------------------------------------------------
# Adaptado para priorizar NRXN, NRG, Glutamate e GABA-B nos circos plots

# Paleta de cores base
cores_pathways_base <- c(
  "GABA-B" = "#008585",
  "Glutamate" = "#e5c185",
  "Dopamine" = "#057dcd",
  "NRG" = "#9370DB",
  "NRXN" = "#F08080"
)

prioritarios <- c("NRXN", "NRG", "Glutamate", "GABA-B", "Dopamine")

# Lista dos genes sinápticos por grupo
genes_sinapticos_por_grupo <- list(
  HBV_liver = c("KCNC1", "FLRT2", "CACNA1A", "TRPV1", "ADGRL3", "SHH", "SEPTIN5",
                "DCX", "SEPTIN3", "CPLX2", "GRM7", "ACAN", "INSYN1", "ACP4", 
                "TAFA4", "RAB3B", "PLPPR4"),
  HBV_pbmc = c("C1QB", "C1QC", "NRN1", "AMPH", "SYN2", "AKAP12", "NXPH4",
               "LRRTM1", "CACNA1E", "MMP9", "IGSF21", "NSG1", "OLFM1",
               "DPYSL3", "ZDHHC8"),
  HCV_liver = c("CLSTN2", "ADGRB3", "NRG1", "ADRA2A", "GPR158", "GRK3",
                "NPTX2", "EPHA3", "VCAN", "CNTN1"),
  HCV_pbmc = c("MOB4", "FMNL2", "PLK2", "PDXP", "RAPGEF2", "PHACTR1", "CPEB2",
               "ADD2", "RPS15", "FAM171B", "VLDLR", "YBX1", "PABPC1", "RPL10A",
               "ITGB1", "SYNGAP1", "HOMER3", "RBMX", "ACHE", "TMEM163"),
  HDV = c("CHRM2", "GRIA2", "NGFR", "GABBR2", "KCNJ6", "KCNC1", "CHRNA4",
          "GABRA5", "NRXN1", "NPY1R", "GRIK3", "HTR2A", "SLC1A7", "GRIN2B",
          "ANO1", "GRIN2A", "EFNB3", "P2RX3", "SCN11A", "PODXL", "FXYD6",
          "NCAM1", "EPHB2", "GRIA3", "CADM2", "NTRK3", "ABHD6", "SEZ6L",
          "GABRR3", "SLC6A8", "KCNQ3", "SUSD4", "CSPG5", "KCNN2", "GABRE",
          "PTPRT", "ROBO2", "PLPPR4", "LRRC4", "ADAM22", "NRG1", "ADCY1",
          "SHISA9", "ELFN2", "CELSR3", "LRRTM3", "IL1RAPL1", "ASIC2", "CSMD2",
          "GPR158", "CYFIP2", "FGB", "GNAZ", "FGA", "FGG", "IL1RAP", "TDRD5",
          "CLU", "ELAVL2", "HAPLN1", "TDRD1", "RAP2A", "VCAN", "DAB1", "NWD2",
          "ADGRB2", "ABLIM3", "GNG4", "CDH11", "CDH23", "CDH13", "PABPC1",
          "BCAS1", "RASGRF1", "ARHGAP39", "THY1", "COMT", "HTR4", "CPLX2",
          "ENAH", "GNAL", "SNX27", "MAP2", "PRKAR2B", "MAP1B", "TSPAN7",
          "AGO2", "ABL2", "TRIM47", "KIF2C", "ROR2", "ARHGEF2", "CPEB3",
          "GPM6A", "LINGO1", "SEMA7A", "NPY5R", "ITGA2", "GAD1", "PRUNE2",
          "SEPTIN3", "SIPA1L2", "SEPTIN8", "GAP43", "SYT10", "DCX", "SPTBN2",
          "TUBB2B", "PLCB4", "EPHB1", "DLGAP2", "CHRNA3", "NRXN3", "APBA1",
          "DGKI", "NR3C2", "CAMK2B", "MAP2K1", "HOMER1", "HOMER2", "CTNND2",
          "AURKA", "DLG5", "TACC3", "SPOCK1", "KPNA2", "SOS1", "CD200",
          "EIF4G2", "DBN1", "TENM3", "SPARCL1", "IL1RAPL2", "YWHAZ",
          "CTTNBP2", "SLITRK4", "CDH8", "ACTN2", "CTNNA2", "MMP9", "C1QB",
          "SNAP25", "C1QC", "FRRS1L", "KIF21B", "LGI1", "LAMA4", "LAMC1",
          "CBLN2", "BCAN", "ITPKA", "RIMS2", "ERC2", "LHFPL4", "SCAMP5",
          "SLC30A3", "SLC17A8", "SYT9", "TMEM163", "GNAO1", "IGF1", "PLCB1",
          "HTR1D", "STXBP5L", "KCND3", "SCRIB", "RAB3B", "RAB3C", "RAB26",
          "SORT1", "PLXNA4", "ACAN", "NTM", "CACNA2D2", "UNC13A", "LINGO2",
          "PCDH9", "TAFA4", "VPS45", "DNAJC6", "SIPA1L3", "PLXNC1", "ATP6V1C1")
)

# Lista dos data.frames filtrados
bases_filtradas <- list(
  HBV_liver = db2_HBV_liver,
  HBV_pbmc = db2_HBV_pbmc,
  HCV_liver = db2_HCV_liver,
  HCV_pbmc = db2_HCV_pbmc,
  HDV = db2_HDV
)

for (grupo in names(bases_filtradas)) {
  
  cat("⏳ Processando:", grupo, "\n")
  
  db3 <- bases_filtradas[[grupo]]
  
  if (nrow(db3) == 0) {
    cat("⚠️  Nenhuma interação encontrada para", grupo, "- pulando...\n")
    next
  }
  
  tab_allu2 <- db3 %>%
    relocate(pathway_name, ligand.symbol, receptor.symbol) %>%
    mutate(num_rows = 1) %>%
    separate_rows(ligand.symbol, sep = ",\\s*") %>%
    separate_rows(receptor.symbol, sep = ",\\s*") %>%
    distinct(pathway_name, ligand.symbol, receptor.symbol, num_rows)
  
  collapse_gene <- function(gene) {
    case_when(
      str_detect(gene, "^DRD|^DOPA") ~ "Dopamine", #Comentados sao pra HDV
      # str_detect(gene, "^GRM[1-8]$") ~ "GRM",
      # str_detect(gene, "^GRID[1-2]$") ~ "GRID",
      # str_detect(gene, "^GRIK[1-8]$") ~ "GRIK",
      # str_detect(gene, "^GRIA[1-8]$") ~ "GRIA",
      # str_detect(gene, "^GRIN[1-2]$") ~ "GRIN",
      # str_detect(gene, "^ERBB[234]$") ~ "ERBB",
      # str_detect(gene, "^ITGA") ~ "ITGA",
      # str_detect(gene, "^ITGB") ~ "ITGB",
      # str_detect(gene, "^NLGN[123]$") ~ "NLGN",
      # str_detect(gene, "^CLSTN[123]$") ~ "CLSTN",
      TRUE ~ gene
    )
  }
  
  tab_allu3 <- tab_allu2 %>%
    mutate(
      ligand = collapse_gene(ligand.symbol),
      receptor = collapse_gene(receptor.symbol)
    ) %>%
    group_by(ligand, receptor, pathway_name) %>%
    summarise(n = n(), .groups = "drop")
  
  pathways_disponiveis <- unique(tab_allu3$pathway_name)
  presentes_prioritarios <- intersect(prioritarios, pathways_disponiveis)
  
  if (length(presentes_prioritarios) >= 4) {
    pathways_selecionados <- presentes_prioritarios[1:4]
  } else {
    extras <- setdiff(pathways_disponiveis, presentes_prioritarios)
    pathways_selecionados <- c(presentes_prioritarios, head(extras, 4 - length(presentes_prioritarios)))
  }
  
  tab_allu3 <- tab_allu3 %>% filter(pathway_name %in% pathways_selecionados)
  
  cores_pathways <- setNames(
    sapply(pathways_selecionados, function(pw) {
      if (pw %in% names(cores_pathways_base)) cores_pathways_base[pw] else "gray80"
    }),
    pathways_selecionados
  )
  
  tab_allu3 <- tab_allu3 %>% mutate(cores = cores_pathways[pathway_name])
  
  grid.col <- tab_allu3 %>%
    pivot_longer(cols = c(ligand, receptor)) %>%
    distinct(value, cores) %>%
    deframe()
  
  orderL <- tab_allu3 %>% arrange(pathway_name) %>% distinct(ligand)
  orderR <- tab_allu3 %>% arrange(pathway_name) %>% distinct(receptor)
  order1 <- c(rev(orderL$ligand), orderR$receptor)
  
  output_file <- paste0("/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/5 - SynGO/",
                        "circos_prioritized_", grupo, ".pdf")
  
  pdf(output_file, width = 6.5, height = 6.5)
  
  circos.clear()
  chordDiagram(tab_allu3,
               grid.col = grid.col,
               preAllocateTracks = 1,
               order = order1,
               annotationTrack = "grid",
               annotationTrackHeight = mm_h(2))
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    name <- get.cell.meta.data("sector.index")
    
    genes_sinapticos <- genes_sinapticos_por_grupo[[grupo]]
    fonte <- if (name %in% genes_sinapticos) 2 else 1
    
    circos.text(mean(xlim), ylim[1] + 0.1, name,
                cex = 0.8,
                facing = "clockwise", niceFacing = TRUE,
                font = fonte,
                adj = c(0, 0.5))
  }, bg.border = NA)
  
  # lgd <- Legend(
  #   at = names(cores_pathways),
  #   legend_gp = gpar(fill = cores_pathways),
  #   title = "Pathway",
  #   title_position = "topleft"
  # )
  # draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
  # 
  dev.off()
  
  cat("✅ Circos salvo:", output_file, "\n")
}

#salvando tabelas
write.xlsx(db2, "Suppl.Table_CellChart.xlsx", rowNames = FALSE)
write.xlsx(db2_HBV_liver, "Suppl.Table_SynapseCircosHBVL.xlsx", rowNames = FALSE)
write.xlsx(db2_HBV_pbmc, "Suppl.Table_SynapseCircosHBVP.xlsx", rowNames = FALSE)
write.xlsx(db2_HCV_liver, "Suppl.Table_SynapseCircosHCVL.xlsx", rowNames = FALSE)
write.xlsx(db2_HCV_pbmc, "Suppl.Table_SynapseCircosHCVP.xlsx", rowNames = FALSE)
write.xlsx(db2_HDV, "Suppl.Table_SynapseCircosHDVL.xlsx", rowNames = FALSE)

#Synapse EnrichR
data_synapse_lr_BP <- read.delim("/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/5 - SynGO/Ligand_Synapse_GENES_BP.txt", header = TRUE, sep = "\t")
data_synapse_lr_CC <- read.delim("/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/5 - SynGO/Ligand_Synapse_GENES_CC.txt", header = TRUE, sep = "\t")
data_synapse_lr_MF <- read.delim("/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/5 - SynGO/Ligand_Synapse_GENES_MF.txt", header = TRUE, sep = "\t")

data_synapse_lr_BP$Ontology <- "BP"
data_synapse_lr_CC$Ontology <- "CC"
data_synapse_lr_MF$Ontology <- "MF"

genes_synapse_lr <- c("GRM7", "FLRT2", "LRRTM1", "GABBR2", "GAD1", "SLC6A8", "NRXN1",
                      "SLC17A8", "SLC1A7", "NRXN3", "LRRTM3", "GRIN2B", "GRIN2A",
                      "ADGRB", "CNTN1", "NRG1", "CLSTN2", "ITGB1")
genes_HBV_liver <- c("GRM7", "FLRT2") # HBV Liver
genes_HBV_pbmc <- c("LRRTM1") # HBV PBMC
genes_HDV_liver <- c("GABBR2", "GAD1", "SLC6A8", "NRXN1", "SLC17A8", "SLC1A7",
                     "NRXN3", "LRRTM3", "GRIN2B", "GRIN2A") # HDV Liver
genes_HCV_liver <- c("ADGRB", "CNTN1", "NRG1", "CLSTN2") # HCV Liver
genes_HCV_pbmc <- c("ITGB1") # HCV PBMC


{
data_synapse_lr_BP$Term <- sub("\\s*\\([^\\)]+\\)", "",
                               data_synapse_lr_BP$Term)
data_synapse_lr_BP <- data_synapse_lr_BP %>%
  separate_rows(Genes, sep = ";")

data_synapse_lr_CC$Term <- sub("\\s*\\([^\\)]+\\)", "",
                               data_synapse_lr_CC$Term)
data_synapse_lr_CC <- data_synapse_lr_CC %>%
  separate_rows(Genes, sep = ";")

data_synapse_lr_MF$Term <- sub("\\s*\\([^\\)]+\\)", "",
                               data_synapse_lr_MF$Term)
data_synapse_lr_MF <- data_synapse_lr_MF %>%
  separate_rows(Genes, sep = ";")
}


# Função para identificar grupo
identificar_grupo <- function(genes) {
  genes_split <- str_split(genes, pattern = ",|;|\\s+", simplify = FALSE)[[1]]
  grupos <- character(0)
  
  if (any(genes_split %in% genes_HBV_liver)) grupos <- c(grupos, "HBV Liver")
  if (any(genes_split %in% genes_HBV_pbmc)) grupos <- c(grupos, "HBV PBMC")
  if (any(genes_split %in% genes_HDV_liver)) grupos <- c(grupos, "HDV Liver")
  if (any(genes_split %in% genes_HCV_liver)) grupos <- c(grupos, "HCV Liver")
  if (any(genes_split %in% genes_HCV_pbmc)) grupos <- c(grupos, "HCV PBMC")
  
  if (length(grupos) == 0) {
    return(NA_character_)
  } else {
    return(paste(grupos, collapse = ";"))
  }
}

# Aplicar a função à coluna "Genes"
{
data_synapse_lr_BP <- data_synapse_lr_BP %>%
  mutate(Group = sapply(Genes, identificar_grupo))
data_synapse_lr_BP <- data_synapse_lr_BP %>%
  filter(Adjusted.P.value < 0.05)

data_synapse_lr_CC <- data_synapse_lr_CC %>%
  mutate(Group = sapply(Genes, identificar_grupo))
data_synapse_lr_CC <- data_synapse_lr_CC %>%
  filter(Adjusted.P.value < 0.05)

data_synapse_lr_MF <- data_synapse_lr_MF %>%
  mutate(Group = sapply(Genes, identificar_grupo))
data_synapse_lr_MF <- data_synapse_lr_MF %>%
  filter(Adjusted.P.value < 0.05)
}

data_synapse_lr <- rbind(data_synapse_lr_BP, 
                         data_synapse_lr_CC,
                         data_synapse_lr_MF)

# Etapa 1: Expandir os genes em múltiplas linhas
data_expanded <- data_synapse_lr %>%
  separate_rows(Genes, sep = ",\\s*") %>%
  mutate(Genes = str_trim(Genes))

# Etapa 2: Calcular número de genes únicos por Term e Ontology
genes_per_term <- data_expanded %>%
  group_by(Ontology, Term) %>%
  summarise(n_genes = n_distinct(Genes), .groups = "drop")

# Etapa 3: Selecionar os 10 termos com mais genes únicos por Ontology
top_terms_by_ontology <- genes_per_term %>%
  group_by(Ontology) %>%
  slice_max(order_by = n_genes, n = 10, with_ties = FALSE) %>%
  ungroup()

# Etapa 4: Filtrar o data frame original para manter apenas os termos selecionados
data_top_terms <- data_expanded %>%
  semi_join(top_terms_by_ontology, by = c("Ontology", "Term"))

# ✅ Resultado final
print(data_top_terms %>% distinct(Ontology, Term) %>% arrange(Ontology))

#DOTPLOT
{
  # Expandir genes
  data_plot <- data_top_terms %>%
    separate_rows(Genes, sep = ",\\s*") %>%
    mutate(
      Genes = str_trim(Genes),
      padj_log = -log10(Adjusted.P.value),
      alpha = case_when(
        Adjusted.P.value < 0.0005 ~ 1,
        Adjusted.P.value < 0.005  ~ 0.75,
        Adjusted.P.value < 0.05   ~ 0.5,
        TRUE                      ~ 0.1
      )
    )
  
  # Ordenar fatores
  data_plot$Genes <- factor(data_plot$Genes, levels = sort(unique(data_plot$Genes)))
  data_plot$Term <- factor(data_plot$Term, levels = rev(unique(data_plot$Term)))
  
  # Cores e formas
  group_colors <- c(
    "HBV Liver" = "#DDF1F5",
    "HBV PBMC" = "#789DBC",
    "HCV Liver" = "#EDAEAE",
    "HCV PBMC" = "#BC7C7C",
    "HDV Liver" = "#FFD09B"
  )
  
  group_shapes <- c(
    "HBV Liver" = 19,
    "HBV PBMC" = 18,
    "HCV Liver" = 19,
    "HCV PBMC" = 18,
    "HDV Liver" = 19 
  )
  
  # Loop para salvar um gráfico por ontologia
  for (ontology in unique(data_plot$Ontology)) {
    
    data_subset <- data_plot %>% filter(Ontology == ontology)
    
    output_path <- paste0("/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/5 - SynGO/Dotplot_Synapse_", ontology, ".tiff")
    
    tiff(output_path, width = 8.5, height = 5, res = 300, units = 'in')
    
    p <- ggplot(data_subset, aes(x = Genes, y = Term)) +
      geom_point(aes(
        size = alpha,
        color = Group,
        shape = Group,
        alpha = alpha
      )) +
      scale_size_continuous(
        name = "Adj.p-value", 
        breaks = c(0.1, 0.5, 0.75, 1), 
        labels = c(">= 0.05", "< 0.05", "< 0.005", "< 0.0005"),
        range = c(2, 10)
      ) +
      scale_shape_manual(values = group_shapes, name = "Group") +
      scale_color_manual(values = group_colors, name = "Group") +
      scale_alpha_continuous(
        name = "Adj.p-value",
        breaks = c(0.1, 0.5, 0.75, 1),
        labels = c(">= 0.05", "< 0.05", "< 0.005", "< 0.0005"),
        range = c(0.1, 1)
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.position = "none"
      ) +
      labs(
        x = "Genes",
        y = "",
        title = paste("", ontology)
      ) +
      guides(
        alpha = guide_legend(title = "Adj.p-value"),
        size = guide_legend(title = "Adj.p-value"),
        color = guide_legend(override.aes = list(size = 5))
      )
    
    print(p)
    dev.off()
  }
  
  
  # Expandir genes
  data_plot <- data_top_terms %>%
    separate_rows(Genes, sep = ",\\s*") %>%
    mutate(
      Genes = str_trim(Genes),
      padj_log = -log10(Adjusted.P.value),
      
      # Alpha suavemente variável
      alpha = case_when(
        Adjusted.P.value < 0.0005 ~ 0.95,
        Adjusted.P.value < 0.005  ~ 0.85,
        Adjusted.P.value < 0.05   ~ 0.75,
        TRUE                      ~ 0.6
      )
    )
  
  # Ordenar fatores
  data_plot$Genes <- factor(data_plot$Genes, levels = sort(unique(data_plot$Genes)))
  data_plot$Term <- factor(data_plot$Term, levels = rev(unique(data_plot$Term)))
  
  # Cores e formas
  group_colors <- c(
    "HBV Liver" = "#DDF1F5",
    "HBV PBMC" = "#789DBC",
    "HCV Liver" = "#EDAEAE",
    "HCV PBMC" = "#BC7C7C",
    "HDV Liver" = "#FFD09B"
  )
  
  group_shapes <- c(
    "HBV Liver" = 19,
    "HBV PBMC" = 18,
    "HCV Liver" = 19,
    "HCV PBMC" = 18,
    "HDV Liver" = 19 
  )
  
  # Filtrar apenas MF
  data_subset <- data_plot %>% filter(Ontology == "MF")
  
  # Salvar o gráfico para MF
  output_path <- "/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/3 - Enrichment/5 - SynGO/Dotplot_Synapse_MF.tiff"
  
  tiff(output_path, width = 8.5, height = 5, res = 300, units = 'in')
  
  p <- ggplot(data_subset, aes(x = Genes, y = Term)) +
    geom_point(aes(
      size = alpha,
      color = Group,
      shape = Group,
      alpha = alpha
    )) +
    scale_size_continuous(
      name = "Adj.p-value", 
      breaks = c(0.6, 0.75, 0.85, 0.95), 
      labels = c(">= 0.05", "< 0.05", "< 0.005", "< 0.0005"),
      range = c(2, 10)
    ) +
    scale_shape_manual(values = group_shapes, name = "Group") +
    scale_color_manual(values = group_colors, name = "Group") +
    scale_alpha_continuous(
      name = "Adj.p-value",
      breaks = c(0.6, 0.75, 0.85, 0.95),
      labels = c(">= 0.05", "< 0.05", "< 0.005", "< 0.0005"),
      range = c(0.6, 0.95)
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 15),
      strip.text = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.position = "none"
    ) +
    labs(
      x = "Genes",
      y = "",
      title = "MF"
    ) +
    guides(
      alpha = guide_legend(title = "Adj.p-value"),
      size = guide_legend(title = "Adj.p-value"),
      color = guide_legend(override.aes = list(size = 5))
    )
  
  print(p)
  dev.off()
  }

#Salvar tabelas
write.xlsx(data_top_terms, "Suppl.Table_30SigSynapseDotPlot.xlsx", rowNames = FALSE)
write.xlsx(data_synapse_lr, "Suppl.Table_AllSigSynapseDotPlot.xlsx", rowNames = FALSE)

# TCGA CIRCOSPLOT ---------------------------------------------------------
# **Carregar os dados**
# data_path <- "ligand-receptor.txt"
# df <- read.delim(data_path, header = TRUE, sep = "\t")

df <- tab_allu2

#**Separar múltiplos receptores e reorganizar os dados**

df_long <- df %>%
  separate_rows(receptor, sep = "_") %>%
  select(pathway_name, ligand, receptor)

df_long <- df_long %>%
  mutate(ligand = gsub("_", " ", ligand))

length(unique(df_long$ligand)) #Ver quantidade de ligantes

# **Selecionar apenas 10 ligantes**
set.seed(123)  
ligands_to_show <- sample(unique(df_long$ligand), size = 5)  

df_filtered <- df_long %>% filter(ligand %in% ligands_to_show)

df_filtered <- df_filtered %>%
  mutate(ligand = gsub("_", " ", ligand))

# **Filtrar pathways e receptores relacionados aos ligantes selecionados**
pathways_filtered <- unique(df_filtered$pathway_name)
receptors_filtered <- unique(df_filtered$receptor)

# **Criar lista de nós filtrados garantindo unicidade**
all_nodes <- unique(c(df_filtered$pathway_name,
                      df_filtered$ligand,
                      df_filtered$receptor))

# **Criar matrizes de adjacência garantindo mesmas colunas**
pathway_ligand_matrix <- table(factor(df_filtered$pathway_name,
                                      levels = all_nodes),
                               factor(df_filtered$ligand,
                                      levels = all_nodes))

ligand_receptor_matrix <- table(factor(df_filtered$ligand,
                                       levels = all_nodes),
                                factor(df_filtered$receptor,
                                       levels = all_nodes))

#  **Concatenar as matrizes ajustadas**
adj_matrix <- pathway_ligand_matrix + ligand_receptor_matrix

# **Verificar se a matriz contém conexões**
print(sum(adj_matrix))  # Deve ser maior que 0

# **Definir grupos para os elementos do diagrama**
group <- structure(
  c(rep("A", length(pathways_filtered)),  
    rep("B", length(ligands_to_show)),  
    rep("C", length(receptors_filtered))  
  ),
  names = c(pathways_filtered, ligands_to_show, receptors_filtered)
)

# **Definir cores personalizadas para Pathway, Ligand e Receptor**
grid.col <- setNames(
  c(rep("#FFCDB2", length(pathways_filtered)),  
    rep("#FFB4A2", length(ligands_to_show)),  
    rep("#B5828C", length(receptors_filtered))  
  ),
  c(pathways_filtered, ligands_to_show, receptors_filtered)
)

# minhasCores <- c("#8EACCD", "#C96868", "gray80")

# **Criar o diagrama com mais espaço entre os setores**
# Abrir dispositivo gráfico para salvar como PNG
tiff("diagram_ligand_receptor.tiff", width = 1400, height = 1400, res = 300)

# Criar o diagrama com mais espaço entre os setores
circos.clear()

# Criar o diagrama de corda com trilha alocada
chordDiagram(adj_matrix, group = group, grid.col = grid.col,
             annotationTrack = c("grid", "axis"),
             preAllocateTracks = list(
               track.height = mm_h(4),
               track.margin = c(mm_h(4), 0)
             ))

# Ajustar rótulos nos setores
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.25, niceFacing = TRUE)
}, bg.border = NA)

# Destacar setores específicos
highlight.sector(unique(df_filtered$pathway_name), track.index = 1, col = "#FFCDB2",
                 text = "Pathway", cex = 0.4, text.col = "black", niceFacing = TRUE)
highlight.sector(unique(df_filtered$ligand), track.index = 1, col = "#FFB4A2",
                 text = "Ligand", cex = 0.4, text.col = "black", niceFacing = TRUE)
highlight.sector(unique(df_filtered$receptor), track.index = 1, col = "#B5828C",
                 text = "Receptor", cex = 0.4, text.col = "black", niceFacing = TRUE)

# Fechar o dispositivo gráfico e salvar o arquivo
dev.off()

#Contagem
# Carregar pacotes necessários
library(ggplot2)
library(dplyr)
library(tidyr)

# Remover duplicatas dentro de cada pathway
df_filtered_unique <- df_filtered %>%
  distinct(pathway_name, ligand, receptor, .keep_all = TRUE)

df_counts <- df_filtered_unique %>%
  pivot_longer(cols = c(ligand,
                        receptor),
               names_to = "type",
               values_to = "element") %>%
  count(pathway_name, type, element)

# Modificar os nomes na coluna "type"
df_counts <- df_counts %>%
  mutate(type = recode(type, "receptor" = "Receptor",
                       "ligand" = "Ligand"))

# Definir cores personalizadas para "ligand" e "receptor"
cores_personalizadas <- c("Ligand" = "#B5828C", "Receptor" = "#FFB4A2")  # Azul e Vermelho

# Definir a ordem desejada para o eixo X
df_counts <- df_counts %>%
  mutate(pathway_name = factor(pathway_name, levels = c("Noradrenaline", "NRG", "IL1")))

# Definir os rótulos personalizados
custom_labels <- c("Receptor" = "Ligand", "Ligand" = "Receptor")

df_counts$type <- factor(df_counts$type, levels = c("Receptor", "Ligand"))

# Criar o gráfico com a nova ordem no eixo X
tiff("counting_ligand_receptor.tiff",
     width = 2500, height = 1700, res = 300)
ggplot(df_counts, aes(x = pathway_name, y = n, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barras lado a lado
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = 1.1, size = 8) +  # Rótulos nas barras
  facet_wrap(~ type, scales = "free_x", labeller = labeller(type = custom_labels)) +  # Facetas personalizadas
  scale_fill_manual(values = cores_personalizadas) +  # Aplicar cores manuais
  scale_y_continuous(breaks = c(0, 5, 9)) +  # Definir quebras no eixo Y
  theme_bw() +
  labs(x = "", y = "Frequency", fill = "Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # Texto do eixo X
    axis.text.y = element_text(size = 20),  # Texto do eixo Y
    axis.title.y = element_text(size = 20),  # Título do eixo Y
    strip.text = element_text(size = 20)  # Tamanho do texto das facetas
  )
dev.off()

#Contagem
# Carregar pacotes necessários
library(ggplot2)
library(dplyr)
library(tidyr)

# Remover duplicatas dentro de cada pathway
df_filtered_unique <- df_filtered %>%
  distinct(pathway_name, ligand, receptor, .keep_all = TRUE)

df_counts <- df_filtered_unique %>%
  pivot_longer(cols = c(ligand,
                        receptor),
               names_to = "type",
               values_to = "element") %>%
  count(pathway_name, type, element)

# Modificar os nomes na coluna "type"
df_counts <- df_counts %>%
  mutate(type = recode(type, "receptor" = "Receptor",
                       "ligand" = "Ligand"))

# Definir cores personalizadas para "ligand" e "receptor"
cores_personalizadas <- c("Ligand" = "#FFB4A2", "Receptor" = "#B5828C")  # Azul e Vermelho

# Definir a ordem desejada para o eixo X
df_counts <- df_counts %>%
  mutate(pathway_name = factor(pathway_name, levels = c("Noradrenaline", "NRG", "IL1")))

# Criar o gráfico com a nova ordem no eixo X
tiff("counting_ligand_receptor.tiff",
     width = 2500, height = 1700, res = 300)
ggplot(df_counts, aes(x = pathway_name, y = n, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barras lado a lado
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = 1.1, size = 8) +  # Rótulos nas barras
  facet_wrap(~ type, scales = "free_x") +  # Facetas separadas para ligand/receptor
  scale_fill_manual(values = cores_personalizadas) +  # Aplicar cores manuais
  scale_y_continuous(breaks = c(0, 5, 9)) +  # Definir quebras no eixo Y
  theme_bw() +
  labs(x = "", y = "Frequency", fill = "Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # Texto do eixo X
    axis.text.y = element_text(size = 20),  # Texto do eixo Y
    axis.title.y = element_text(size = 20),  # Título do eixo Y
    strip.text = element_text(size = 20)  # Tamanho do texto das facetas
  )
dev.off()

# Correloplots (TCGA) -----------------------------------------------------
library(corrplot)
library(GGally)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(scales)
library(cluster)
library(dplyr)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
#Correloplots
data_cor <- my_data2

# Lista de genes desejados
genes_of_interest <- c("WDR62", "TRIM71", "HP", "FER1L4", "GCH1", "C9", "IL33", "VSIG4", "DBH", "OLFM1", "NRG1")

# Selecionando as colunas de genes no data frame my_data2
data_cor <- data_cor[, c(genes_of_interest, "ID")]

write.xlsx(data_cor, "Suppl.Table5b_inputcorreloplot.xlsx", rowNames = TRUE)

# Separar os dados pelo grupo ID
data_cor_grouped <- data_cor %>% 
  split(.$ID)  # Cria uma lista de data frames, um para cada grupo

# Função para calcular a correlação e gerar o heatmap
plot_correlation_heatmap <- function(data_group, order_clusters = NULL) {
  # Selecionar apenas as colunas numéricas para calcular a correlação
  data_group <- data_group[sapply(data_group, is.numeric)]
  
  # Calcular a correlação de Spearman
  cor_matrix <- cor(data_group, method = "spearman")
  
  # Reorganizar a matriz de correlação, se necessário
  if (!is.null(order_clusters)) {
    cor_matrix <- cor_matrix[order_clusters, order_clusters]
  }
  
  # Definir paleta de cores
  col_func <- colorRampPalette(c("#93BFCF", "white", "#D71313"))
  
  # Gerar o heatmap
  corrplot(cor_matrix,
           method = "circle",
           addrect = 2,
           rect.col = "black",
           rect.lwd = 3,
           type = "full",
           tl.pos = "lt",  # Ajusta a posição dos rótulos
           col = col_func(200),
           tl.col = "black")  # Cor dos rótulos
  
  # Retornar a matriz de correlação para salvar depois
  return(cor_matrix)
}

# Inicializa a variável que irá armazenar a ordem hierárquica
order_g1 <- NULL

# Loop para cada grupo em data_cor_grouped
for (group in names(data_cor_grouped)) {
  # Subconjunto dos dados para o grupo específico
  data_group <- data_cor_grouped[[group]]
  
  # Para G1, calcular a ordem hierárquica
  if (group == "G1") {
    cor_matrix_g1 <- cor(data_group[sapply(data_group, is.numeric)], method = "spearman")
    order_g1 <- hclust(dist(cor_matrix_g1))$order
  }
  
  # Nome do arquivo de saída para o heatmap
  output_filename <- paste0("heatmap_", group, ".tiff")
  
  # Nome do arquivo de saída para a tabela de correlação
  cor_table_filename <- paste0("correlation_table_", group, ".xlsx")
  
  # Salvar a tabela de correlação como CSV
  cor_matrix <- plot_correlation_heatmap(data_group, order_clusters = order_g1)
  write.xlsx(cor_matrix, file = cor_table_filename)  # Salva a tabela de correlação
  
  # Salvar o heatmap como TIFF
  tiff(output_filename, width = 4, height = 4, units = "in", res = 300)
  
  # Plotar o heatmap para o grupo atual
  cat("Plotando heatmap para o grupo:", group, "\n")
  plot_correlation_heatmap(data_group, order_clusters = order_g1)
  
  # Fechar o dispositivo gráfico
  dev.off()
}

# Manova - Relative Effect (TCGA) -----------------------------------------
{
  library(npmv)
  library(reshape2)
  library(ggplot2)
  library(openxlsx)
  library(dplyr)
}

# setwd("~/Análise IC/Colaboração/Adriel/Relative Effect/")
setwd("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/")
outDir <- "~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/Relative Effect/"
dir.create(outDir)

my_data12 <- read.xlsx("Genes_GBM_G1&2.xlsx", rowNames = TRUE)
my_data13 <- read.xlsx("Genes_GBM_G1&3.xlsx", rowNames = TRUE)
my_data14 <- read.xlsx("Genes_GBM_G1&4.xlsx", rowNames = TRUE)

#Manipulando as tabelas
data_G1 <- my_data12 %>% filter(neoplasmhistologicgrade == "G1")
data_G2 <- my_data12 %>% filter(neoplasmhistologicgrade == "G2")
data_G3 <- my_data13 %>% filter(neoplasmhistologicgrade == "G3")
data_G4 <- my_data14 %>% filter(neoplasmhistologicgrade == "G4")

#Unindo novamento
data_G12 = rbind(data_G1, data_G2) #G2 com G3
data_G13 = rbind(data_G1, data_G3) #G2 com G3
data_G14 = rbind(data_G1, data_G4) #G3 com G4
# data_G14 = rbind(data_G1, data_G4) #G1 com G4

#Salvando como "my_data"
my_data = data_G14

# my_data <- my_data[,-1] #only acute

# my_data <- my_data %>% 
#   rename(groups = neoplasmhistologicgrade) #Renomeando coluna com grupos

colnames(my_data)[colnames(my_data) == "neoplasmhistologicgrade"] <- "groups"

#Verificando quantos elementos tem do grupo
sum(grepl("G4",my_data$groups))
View(my_data)
#Trocando os identificadores por 1 e 0
my_data$groups <- c(rep(0, 60), #G1
                    rep(1, 205))  #G2

my_data$groups <- c(rep(0, 60), #G1
                    rep(1, 137))  #G3

my_data$groups <- c(rep(0, 60), #G1
                    rep(1, 11))  #G4

# my_data$groups <- c(rep(0, 60), #G1
#                     rep(1, 11))  #G4

dplyr::glimpse(my_data)

n_aab = nrow(my_data)
list_boot_aab = list()

n_boot = 1000

# Manova with bootstrap 
for (i in 0:n_boot) {
  if (i == 0) {
    df_sample = 1:n_aab
  } else {
    df_sample = sample(1:n_aab, n_aab, replace = TRUE)
  }
  
  manova_aab = nonpartest(VSIG4 | IL33 | HP | DBH | GCH1 |
                            FER1L4 | C9 | NRG1 | OLFM1 | WDR62 | TRIM71 ~ groups,
                          data = my_data[df_sample, ], permtest = FALSE, plots = FALSE,
                          permreps = 1000, tests = c(1, 0, 0, 0)
  )
  
  if (i == 0) obs_aab = manova_aab
  
  list_boot_aab[[i + 1]] = manova_aab
  
  cat("ANOVA/ boot/ aab/ i = ", i, "\n")
}

obs_manova = list_boot_aab[[1]]$twogroupreleffects
row_manova = rownames(obs_manova)

head(obs_manova)

releffects = lapply(list_boot_aab[-1], "[[", 2)

#intervalo de confiança
df_perc_boot = lapply(seq_len(nrow(obs_manova)), function(i) {
  q <- sapply(seq_len(ncol(obs_manova)), function(j) {
    v <- sapply(releffects, function(r) {
      row <- rownames(r)
      if(row[1] != row_manova[1]) r <- r[2:1,]
      r[i, j]
    })
    quantile(v, probs = c(0.025, 0.975))
  })
  colnames(q) <- colnames(obs_manova)
  df <- data.frame(Var2 = colnames(q), t(q),
                   group = ifelse(i == 1, row_manova[1], row_manova[2]))
  df$Var2 <- factor(df$Var2, levels = sort(unique(df$Var2)), ordered = TRUE)
  return(df)
})

#
df_perc_boot = lapply(seq_len(nrow(obs_manova)), function(i) {
  q <- sapply(seq_len(ncol(obs_manova)), function(j) {
    v <- sapply(releffects, function(r) {
      row <- rownames(r)
      if(row[1] != row_manova[1]) r <- r[2:1,]
      r[i, j]
    })
    quantile(v, probs = c(0.025, 0.975))
  })
  colnames(q) <- colnames(obs_manova)
  df <- data.frame(Var2 = colnames(q), t(q),
                   group = ifelse(i == 1, row_manova[1], row_manova[2]))
  df$Var2 <- factor(df$Var2, levels = sort(unique(df$Var2)), ordered = TRUE)
  return(df)
})

names(df_perc_boot) = row_manova

head(df_perc_boot)

df = melt(
  cbind(obs_manova, group = row_manova),
  id.vars = "group"
)
df$Var1 = "point_est"
df = df[, c(4, 2, 3, 1)]
colnames(df)[2] = "Var2"

df = dplyr::arrange(df, group, desc(value), Var2)
df$Var2 = as.character(df$Var2)
df$Var2 = factor(df$Var2, levels = unique(df$Var2), ordered = TRUE)

for (k in seq_len(2)) {
  df_perc_boot[[k]]$Var2 <- factor(df_perc_boot[[k]]$Var2, levels = levels(df$Var2), ordered = TRUE)
}

# Agora fazer o plot com ggplot2
p19 <- ggplot(data = df, aes(x = Var2, y = value, group = group)) +
  geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[1]],
              aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group), fill = "#FFAAA6",
              alpha = 0.3) +
  geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[2]],
              aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group), fill = "#DF6A6A",
              alpha = 0.3) +
  geom_line(aes(color = ifelse(group == 0, "#FFAAA6", "#DF6A6A"))) +
  geom_point(aes(size = value, color = ifelse(group == 0, "#FFAAA6", "#DF6A6A"))) +  # Cor direta no geom_point
  scale_colour_identity() +  # Mantém as cores sem mapear legendas automaticamente
  theme_classic() +
  theme(
    legend.text = element_text(size = 12),
    axis.text.x = element_text(size = 15, angle = 90),
    axis.text.y = element_text(size = 20, angle = 90, hjust = 0.5),
    axis.line.y = element_blank(),
    axis.title = element_text(size = 20),
    panel.grid.major = element_line(),
    legend.position = "none",
    legend.box.background = element_rect()
  ) +
  labs(y = "", x = '')

ggsave(filename = "Relative_effect_TCGA_G12.tiff",
       plot = p19,
       device = "tiff", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300)  

p20 <- ggplot(data = df, aes(x = Var2, y = value, group = group)) +
  geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[1]],
              aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group), fill = "#FFAAA6",
              alpha = 0.3) +
  geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[2]],
              aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group), fill = '#BB5A5A',
              alpha = 0.3) +
  geom_line(aes(color = ifelse(group == 0, "#FFAAA6", '#BB5A5A'))) +
  geom_point(aes(size = value, color = ifelse(group == 0, "#FFAAA6", '#BB5A5A'))) +  # Cor direta no geom_point
  scale_colour_identity() +  # Mantém as cores sem mapear legendas automaticamente
  theme_classic() +
  theme(
    legend.text = element_text(size = 12),
    axis.text.x = element_text(size = 15, angle = 90),
    axis.text.y = element_text(size = 20, angle = 90, hjust = 0.5),
    axis.line.y = element_blank(),
    axis.title = element_text(size = 20),
    panel.grid.major = element_line(),
    legend.position = "none",
    legend.box.background = element_rect()
  ) +
  labs(y = "", x = '')

ggsave(filename = "Relative_effect_TCGA_G13.tiff",
       plot = p20,
       device = "tiff", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300)  

p21 <- ggplot(data = df, aes(x = Var2, y = value, group = group)) +
  geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[1]],
              aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group), fill = "#FFAAA6",
              alpha = 0.3) +
  geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[2]],
              aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group), fill = "#DF6A6A",
              alpha = 0.3) +
  geom_line(aes(color = ifelse(group == 0, "#FFAAA6", "darkred"))) +
  geom_point(aes(size = value, color = ifelse(group == 0, "#FFAAA6", "darkred"))) +  # Cor direta no geom_point
  scale_colour_identity() +  # Mantém as cores sem mapear legendas automaticamente
  theme_classic() +
  theme(
    legend.text = element_text(size = 12),
    axis.text.x = element_text(size = 15, angle = 90),
    axis.text.y = element_text(size = 20, angle = 90, hjust = 0.5),
    axis.line.y = element_blank(),
    axis.title = element_text(size = 20),
    panel.grid.major = element_line(),
    legend.position = "none",
    legend.box.background = element_rect()
  ) +
  labs(y = "", x = '')

ggsave(filename = "Relative_effect_TCGA_G14.tiff",
       plot = p21,
       device = "tiff", 
       width = 6, 
       height = 4, 
       units = "in", 
       dpi = 300)  

write.xlsx(my_data, "Suppl.Table5a_G2XG1_Inputcorreloplot.xlsx", rowNames = TRUE)
write.xlsx(df, "Suppl.Table5a_G2XG1_Outputcorreloplot.xlsx", rowNames = TRUE)

write.xlsx(my_data, "Suppl.Table5a_G3XG1_Inputcorreloplot.xlsx", rowNames = TRUE)
write.xlsx(df, "Suppl.Table5a_G3XG1_Outputcorreloplot.xlsx", rowNames = TRUE)

write.xlsx(my_data, "Suppl.Table5a_G4XG1_Inputcorreloplot.xlsx", rowNames = TRUE)
write.xlsx(df, "Suppl.Table5a_G4XG1_Outputcorreloplot.xlsx", rowNames = TRUE)

# Bubble Heatmap - Genes Expression --------------------------------------

#VolcanoPlot tudo em todo lugar, ao mesmo tempo
setwd("/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/2 - Express Analysis/Metanalise - MetaVolcano/")
#Carregando Dados (In Vitro)
HAV = read.xlsx("1 - MetaVolcano OutPut/InVitro/HAV/meta_DEGs.xlsx")
HBV = read.xlsx("1 - MetaVolcano OutPut/InVitro/HBV/meta_DEGs.xlsx")
HCV = read.xlsx("1 - MetaVolcano OutPut/InVitro/HCV/meta_DEGs.xlsx")
HDV = read.xlsx("1 - MetaVolcano OutPut/InVitro/HDV/BDGSE112118/BDGSE112118.xlsx")
#Editando HDV
HDV = HDV[,c(6, 5, 2, 4)]
colnames(HDV)[colnames(HDV) == "pvalue"] <- "metap"
colnames(HDV)[colnames(HDV) == "log2FoldChange"] <- "metafc"
colnames(HDV)[colnames(HDV) == "stat"] <- "idx"
HEV = read.xlsx("1 - MetaVolcano OutPut/InVitro/HEV/meta_DEGs.xlsx")

#Carregando Dados (In Vivo)
HBV_liver <- read.xlsx("1 - MetaVolcano OutPut/HBV_Liver/meta_DEGs.xlsx")
HBV_PBMC <- read.xlsx("1 - MetaVolcano OutPut/HBV_PBMC/meta_DEGs.xlsx")
HCV_liver <- read.xlsx("1 - MetaVolcano OutPut/HCV_Liver/meta_DEGs.xlsx")
HCV_PBMC <- read.xlsx("1 - MetaVolcano OutPut/HCV_PBMC/meta_DEGs.xlsx")
HDV_liver <- read.xlsx("1 - MetaVolcano OutPut/HDV_Liver/meta_DEGs.xlsx")

#TCGA
data.g12 = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/G12deseqTotal.xlsx", rowNames = TRUE)
data.g13 = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/G13deseqTotal.xlsx", rowNames = TRUE)
data.g14 = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/G14deseqTotal.xlsx", rowNames = TRUE)

#Single-Cell
Bexpb <- read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/9 - Single-cell/Results/HBV/blood_treat.xlsx") #Expressao Total
Bexpbt <- read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/9 - Single-cell/Results/HBV/blood.xlsx") #Expressao Total
Bexpl <- read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/9 - Single-cell/Results/HBV/liver_treat.xlsx") #Expressao Total
Bexplt <- read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/9 - Single-cell/Results/HBV/liver.xlsx") #Expressao Total
Cexpb <- read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/9 - Single-cell/Results/HCV/hcvblood_diff.xlsx") #Expressao Total
Cexpl <- read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/9 - Single-cell/Results/HCV/hcvliver_diff.xlsx") #Expressao Total

colnames(HBV_liver)
#Editando
data.g12$gene_symbol <- rownames(data.g12)
data.g13$gene_symbol <- rownames(data.g13)
data.g14$gene_symbol <- rownames(data.g14)

data.g12 = data.g12[,c(1:2,6,7)]
data.g13 = data.g13[,c(1:2,6,7)]
data.g14 = data.g14[,c(1:2,6,7)]

#names
names(data.g12)[1] <- "idx"
names(data.g12)[2] <- "metafc"
names(data.g12)[3] <- "metap"
names(data.g13)[1] <- "idx"
names(data.g13)[2] <- "metafc"
names(data.g13)[3] <- "metap"
names(data.g14)[1] <- "idx"
names(data.g14)[2] <- "metafc"
names(data.g14)[3] <- "metap"

#Nomeando Hepatites
HAV$Hepatitis <- c("HAV")
HBV$Hepatitis <- c("HBV")
HCV$Hepatitis <- c("HCV")
HDV$Hepatitis <- c("HDV")
HEV$Hepatitis <- c("HEV")
HBV_liver$Hepatitis <- c("HBV")
HBV_PBMC$Hepatitis <- c("HBV")
HCV_liver$Hepatitis <- c("HCV")
HCV_PBMC$Hepatitis <- c("HCV")
HDV_liver$Hepatitis <- c("HDV")
data.g12$Hepatitis <- c("G2")
data.g13$Hepatitis <- c("G3")
data.g14$Hepatitis <- c("G4")

#Nomeando Tecido
HAV$Tissue <- c("In vitro")
HBV$Tissue <- c("In vitro")
HCV$Tissue <- c("In vitro")
HDV$Tissue <- c("In vitro")
HEV$Tissue <- c("In vitro")
HBV_liver$Tissue <- c("Liver")
HBV_PBMC$Tissue <- c("PBMC")
HCV_liver$Tissue <- c("Liver")
HCV_PBMC$Tissue <- c("PBMC")
HDV_liver$Tissue <- c("Liver")
data.g12$Tissue <- c("Liver")
data.g13$Tissue <- c("Liver")
data.g14$Tissue <- c("Liver")

#Nomeando Grupo
HAV$Group <- c("HAV Invitro")
HBV$Group <- c("HBV Invitro")
HCV$Group <- c("HCV Invitro")
HDV$Group <- c("HDV Invitro")
HEV$Group <- c("HEV Invitro")
HBV_liver$Group <- c("HBV Liver")
HBV_PBMC$Group <- c("HBV PBMC")
HCV_liver$Group <- c("HCV Liver")
HCV_PBMC$Group <- c("HCV PBMC")
HDV_liver$Group <- c("HDV Liver")
data.g12$Group <- c("HCC G2")
data.g13$Group <- c("HCC G3")
data.g14$Group <- c("HCC G4")

#Garantindo que os nomes sejam os mesmos
names(HBV_liver)[1] <- "gene_symbol"
names(HBV_PBMC)[1] <- "gene_symbol"
names(HCV_liver)[1] <- "gene_symbol"
names(HCV_PBMC)[1] <- "gene_symbol"
names(HDV_liver)[1] <- "gene_symbol"

#Salvando como .xlsx
# Lista dos data frames
{data_list <- list(
  HAV = HAV,
  HBV = HBV,
  HCV = HCV,
  HDV = HDV,
  HEV = HEV,
  HBV_liver = HBV_liver,
  HBV_PBMC = HBV_PBMC,
  HCV_liver = HCV_liver,
  HCV_PBMC = HCV_PBMC,
  HDV_liver = HDV_liver,
  data_g12 = data.g12,
  data_g13 = data.g13,
  data_g14 = data.g14
)

# Nome do arquivo Excel de saída
output_file <- "data_frames_volcano.xlsx"

# Criar um arquivo Excel e adicionar cada data frame como uma aba
wb <- createWorkbook()

for (name in names(data_list)) {
  addWorksheet(wb, name)                      # Criar uma aba com o nome do data frame
  writeData(wb, name, data_list[[name]])      # Escrever o data frame na aba correspondente
}

# Salvar o arquivo Excel
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("Data frames salvos no arquivo:", output_file, "\n")}

#Unindo Dados
data.volcano = rbind(HAV, HBV, HCV, HDV, HEV,
                     HBV_liver, HBV_PBMC, HCV_liver,HCV_PBMC,
                     HDV_liver, data.g12, data.g13, data.g14)

#Criando Volcano
#Criando coluna de expressao (Up, Down e Not Significant)
data.volcano$diffexpressed <- "NO"
data.volcano$diffexpressed <- ifelse(data.volcano$metafc > 1 & data.volcano$metap < 0.05,
                                     "Up", ifelse(data.volcano$metafc < -1 & data.volcano$metap < 0.05,
                                                  "Down", "Not significant"))
# Criando os identificadores significantes
data.volcano$delabel <- NA
mask <- !is.na(data.volcano$diffexpressed) & data.volcano$diffexpressed != "Not significant" &
  !is.na(data.volcano$gene_symbol) & data.volcano$gene_symbol != "Not significant"
data.volcano$delabel[mask] <- data.volcano$gene_symbol[mask]

#Plotando
colnames(data.volcano)

#Salvando Tabelas
data.volcano <- data.volcano[!is.na(data.volcano$diffexpressed), ]

# write.xlsx(data.volcano,
#            "data.voltano_all_conditions.xlsx")

# Ajuste os limites dos eixos conforme necessário
min_value_x <- -4
max_value_x <- 4
min_value_y <- 0
max_value_y <- 15

# Filtrando os dados para garantir que eles estejam dentro dos limites desejados
data_filt <- data.volcano %>%
  filter(metafc >= min_value_x & metafc <= max_value_x & -log10(metap) >= min_value_y & -log10(metap) <= max_value_y)

#Ordenando
unique(data_filt$Group)
data_filt$Group <- factor(data_filt$Group,
                          levels = c("HAV Invitro", "HBV Invitro", "HCV Invitro", 
                                     "HDV Invitro", "HEV Invitro",
                                     "HBV Liver", "HBV PBMC", "HCV Liver", "HCV PBMC",
                                     "HDV Liver", "HCC G2", "HCC G3", "HCC G4"))
data_filt$Group <- factor(data_filt$Group,
                          levels = c("HAV Invitro", "HBV Invitro", "HCV Invitro", 
                                     "HDV Invitro", "HEV Invitro",
                                     "HBV PBMC", "HCV PBMC", "HBV Liver", "HCV Liver",
                                     "HDV Liver", "HCC G2", "HCC G3", "HCC G4"))

# Plotando - Volcano com os dados filtrados

# p <- ggplot(data = data_filt,
#             aes(x = metafc,
#                 y = -log10(metap),
#                 col = diffexpressed,
#                 label = delabel)) +
#   geom_point(size = 1) +
#   # geom_text_repel(size = 3, color = "black") +
#   facet_wrap(~Group, nrow = 2, scales = "free_x") +
#   theme_bw(12) +
#   scale_color_manual(values = c("#93BFCF", "#D8D8D8", "#D71313")) +
#   theme(panel.grid = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.text.x = element_text(size = 14),  # Ajuste o tamanho do texto no eixo X
#         axis.text.y = element_text(size = 14),  # Ajuste o tamanho do texto no eixo Y
#         strip.text = element_text(size = 14),   # Ajuste o tamanho do texto do facet_wrap
#         legend.text = element_text(size = 14),  # Aumenta o tamanho do texto da legenda
#         legend.title = element_text(size = 14), # Aumenta o tamanho do título da legenda
#         legend.position = "right") +           # Muda a posição da legenda
#   xlim(c(min_value_x, max_value_x)) +
#   ylim(c(min_value_y, max_value_y)) +
#   labs(x = "Log2 Fold Change", y = "-Log10 p-value")

#Bubble Heatmap - Expressao Diferencial
genes_of_interest <- c("WDR62", "TRIM71", "HP", "FER1L4",
                       "GCH1", "C9", "IL33", "VSIG4", "DBH", 
                       "OLFM1", "NRG1")

data_bubble <- data.volcano %>%
  filter(gene_symbol %in% genes_of_interest)

order_vector <- c("OLFM1", "IL33", "VSIG4", "C9", "HP", 
                  "GCH1", "NRG1", "DBH", "TRIM71",
                  "WDR62", "FER1L4")

data_bubble <- data.volcano %>%
  filter(gene_symbol %in% genes_of_interest) %>%
  mutate(log_metap = -log10(metap),  # Transformação para facilitar visualização
         gene_symbol = factor(gene_symbol, levels = order_vector))  # Ordenar gene_symbol pelo vetor order_vector

data_bubble$metafc <- ifelse(data_bubble$metafc > 2, 2, data_bubble$metafc)
data_bubble$metafc <- ifelse(data_bubble$metafc < -2, -2, data_bubble$metafc)

data_bubble$diffexpressed <- ifelse(data_bubble$metafc > 0 & data_bubble$metap < 0.05,
                                    "Up", ifelse(data_bubble$metafc < -0 & data_bubble$metap < 0.05,
                                                 "Down", "Not significant"))

# p23 <- ggplot(data_bubble, aes(x = Hepatitis, y = gene_symbol)) +
#   geom_point(aes(size = log_metap, color = diffexpressed)) +  # Tamanho -> p-valor, Cor -> log fold-change
#   scale_size_continuous(name = "-log10(p-value)", range = c(2, 10)) +  # Escala de tamanho ajustada
#   scale_color_manual(values = c("Down" = "#93BFCF", "Up" = "#D71313", "Not significant" = "#D8D8D8")) +
#   theme_minimal(base_size = 12) +
#   facet_wrap(~Tissue, scales = "free") +
#   theme(axis.text.x = element_text(size = 10),  # Rotacionar os labels do eixo X
#         axis.text.y = element_text(size = 10),
#         legend.position = "right") +
#   labs(x = "Group", y = "Gene Symbol", title = "Bubble Heatmap of Meta Analysis Results")

# Criar o bubble heatmap
p23 <- ggplot(data_bubble, aes(x = Hepatitis, y = gene_symbol)) +
  geom_point(aes(size = log_metap, color = metafc)) +  # Tamanho -> p-valor, Cor -> log fold-change
  scale_size_continuous(name = "-log10(p-value)", range = c(2, 10)) +  # Escala de tamanho ajustada
  scale_color_gradient2(low = "#93BFCF", mid = "#D8D8D8", high = "#D71313", midpoint = 0, 
                        name = "Log2FC") +  # Escala de cor ajustada
  theme_bw() +
  facet_wrap(~Tissue, scales = "free_x") +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y =  element_blank(),
        strip.text = element_text(size = 12,colour = "black"),
        legend.position = "none") +
  labs(x = "", y = "", title = "")

p23_1 <- ggplot(data_bubble, aes(x = Hepatitis, y = gene_symbol)) +
  geom_point(aes(size = log_metap, color = metafc)) +  # Tamanho -> p-valor, Cor -> log fold-change
  scale_size_continuous(name = "-log10(p-value)", range = c(2, 10)) +  # Escala de tamanho ajustada
  scale_color_gradient2(low = "#93BFCF", mid = "#D8D8D8", high = "#D71313", midpoint = 0, 
                        name = "Log2FC") +  # Escala de cor ajustada
  theme_void() +  # Remove todos os elementos visuais do gráfico
  theme(
    legend.position = "right",  # Posição da legenda
    legend.title = element_text(size = 12, colour = "black"),  # Título da legenda
    legend.text = element_text(size = 12, colour = "black")  # Texto da legenda
  ) +
  labs(title = "")

ggsave(filename = "All_diff_Bubble_Genes.tiff",
       plot = p23,
       device = "tiff", 
       width = 10, 
       height = 4, 
       units = "in", 
       dpi = 300) 

ggsave(filename = "All_diff_Bubble_Genes_labels.tiff",
       plot = p23_1,
       device = "tiff", 
       width = 10, 
       height = 4, 
       units = "in", 
       dpi = 300) 

write.xlsx(data_bubble, "Suppl.Table5c_InputBubbleHeatmap.xlsx", rowNames = TRUE)

# Diseasome (TCGA) --------------------------------------------------------
#-------------------      Networking.
{library(ggplot2)
  library(tidyverse)
  library(ggnet)
  library(igraph)
  library(GGally)
  library(ggnet)
  library(network)
  library(sna)
  library(RColorBrewer)
  library(intergraph)
  library(openxlsx)
  library(dplyr)
  library(ggrepel)}
#Carregando e editando Tabelas
setwd("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/")

data_disease = read.delim("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/Diseasome/DisGeNET_table.txt", header = TRUE, sep = "\t")
data_deseq = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/results/G14_DESeq2Resultado.xlsx")
names(data_deseq)[1] = "Genes"

#Filtrando Significantes
data <- data_disease[data_disease$Adjusted.P.value < 0.05, ]
data <- data %>%
  separate_rows(Genes, sep = ";")

#Unindo: Genes/Doencas com up/down
data_enrich_genes <- left_join(data,
                               data_deseq,
                               by = "Genes",
                               suffix = c("", "_deseq"))
data_enrich_genes = na.omit(data_enrich_genes)

#Filtrando Genes de Interesse
genes_of_interest <- c("WDR62", "TRIM71", "HP", "FER1L4", "GCH1", "C9", "IL33", "VSIG4", "DBH", "OLFM1", "NRG1")

# Selecionando as colunas de genes no data frame my_data2
data_enrich_genes <- data_enrich_genes %>%
  dplyr::filter(Genes %in% genes_of_interest)

#Criando objeto Networkig
net_data <- data_enrich_genes[,c(1,9,16)]
table(unique(net_data$Genes))
table(unique(net_data$Term))
unique(net_data$Genes)
unique(net_data$Term)

dfnet_total <- net_data %>% distinct(Term,
                                     Genes, 
                                     .keep_all = TRUE)
net <- network(dfnet_total[, c("Term", "Genes")])

#Separando tipos celulares como grupos
u <- data.frame(Term = network.vertex.names(net))

bp_df <- net_data[c(1,3)] #Term e diff
gene_df <- net_data[c(2,3)] #Genes e diff
colnames(gene_df)[1] <- colnames(bp_df)[1]

u <- unique(left_join(u, rbind(bp_df,gene_df)))
# # Crie um vetor de palavras-chave para genes e doenças
# genes_keywords <- c("GCH1", "LIFR", "IL33", "TNFAIP6", "HP", "VSIG4", "EREG", "PLP1", "PROK2", "NRG1", "DBH", "WDR62", "OLFM1", "FER1L4", "SPOCK1", "C9")
# diseases_keywords <- c("Intermittent fever", "Arthritis", "Skin toxicity", "Major Depressive Disorder", "Visceral Pain", "Inflammation", "Cerebral Palsy", "Unipolar Depression", "Familial Dystonia", "Myoclonic dystonia", "Dystonia Musculorum Deformans", "Congenital Abnormality", "Colitis", "Schizophrenia, Childhood", "Osteopenia", "Mood Disorders", "Dysautonomia", "Depressive disorder", "Cardiovascular Diseases", "Deglutition Disorders", "Autonomic nervous system disorders", "Choreoathetoid movements", "Amyloidosis", "Stomach Carcinoma", "Malignant neoplasm of stomach", "Psychotic symptom", "Chronic schizophrenia", "Pain", "Pain, Crushing", "Pain, Migratory", "Pain, Splitting", "Tardive Dyskinesia", "Suffering, Physical", "Ache", "Drug-induced tardive dyskinesia", "Radiating pain", "Choreoathetosis", "Pain, Burning", "NEUROTICISM", "Congestive heart failure", "Exanthema", "Extrapyramidal Disorders", "Alzheimer's Disease", "Dystonia", "Systemic Inflammatory Response Syndrome", "Mental Depression", "NLRP12")
# 
# # Criar uma nova coluna chamada "Type" na estrutura de dados "u" com base nas palavras-chave
# u <- mutate(u, Type = ifelse(Term %in% genes_keywords,
#                              "Genes", 
#                              ifelse(Term %in% diseases_keywords,
#                                     "Diseases",
#                                     NA)))

u %>% mutate(Down = ifelse(diff == "Down", 1, 0),
             Up = ifelse(diff == "Up", 1, 0)) %>%
  group_by(Term) %>%
  summarise(Down = sum(Down),
            Up = sum(Up)) %>%
  mutate(condition = case_when(Down == 1 & Up == 0 ~"Down",
                               Up == 1 & Down == 0 ~"Up",
                               .default = "Both")) -> u

net %v% "Groups" <- as.character(u$condition)

unique(u$condition)
tiff("/Users/adrielnobile/networking_diseasome_original.tiff",
     width = 12,
     height = 8,
     res = 300, units = 'in') 
ggnet2(net, 
       size.cut = 50,
       # size.min = 1.1,
       # size.max = 10,
       size = "degree",
       shape = "Groups",
       label = TRUE,
       label.size = 3,
       alpha = 1,
       color = "Groups",
       size.legend.title = "Degree",  
       size.legend.position = "bottom",
       shape.mapping = aes(shape = degree),  # Use degree para shape diretamente
       shape.palette = c("Down" = 19,
                         "Up" = 19,
                         "Both" = 19),
       shape.legend = TRUE,
       palette = c("Down" = "#93BFCF",
                   "Up" = "#D71313",
                   "Both" = "gray50"),
       legend.size = 18)
dev.off()

#salvando como svg
ggsave("networking_diseasome.svg",  # Nome do arquivo e formato
       ggnet2(net, 
              size.cut = 50,
              # size.min = 1.1,
              # size.max = 10,
              size = "degree",
              shape = "Groups",
              label = TRUE,
              label.size = 4,
              alpha = 1,
              color = "Groups",
              size.legend.title = "Degree",  
              size.legend.position = "bottom",
              shape.mapping = aes(shape = degree),  # Use degree para shape diretamente
              shape.palette = c("Down" = 19,
                                "Up" = 19,
                                "Both" = 19),
              shape.legend = TRUE,
              palette = c("Down" = "#93BFCF",
                          "Up" = "#D71313",
                          "Both" = "grey50"),
              legend.size = 18),
       device = "svg",  # Define o formato do arquivo
       width = 10,  # Largura do gráfico em polegadas
       height = 6,  # Altura do gráfico em polegadas
       units = "in")  # Unidade de medida
View(u)
#Salvando Tabela de interacoes
write.xlsx(u,
           "Interacoes da Networking.xlsx")

#Criando o Heatmap Diseasome
#Completando vales faltantes
genes_of_interest <- c("WDR62", "TRIM71", "HP", "FER1L4", "GCH1", "C9", "IL33", "VSIG4", "DBH", "OLFM1", "NRG1")

# Selecionando as colunas de genes no data frame my_data2
data <- data %>%
  dplyr::filter(Genes %in% genes_of_interest)

#Transformando maiusculos
data$Term <- as.factor(sapply(as.character(data$Term), function(x) {
  # Deixar a string inteira em minúsculas e depois capitalizar a primeira letra
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}))

all_combinations <- expand.grid(Genes = unique(data$Genes),
                                Term = unique(data$Term))
# Combinar com os dados existentes
merged_data <- merge(all_combinations, data, by = c("Genes", "Term"), all.x = TRUE)
# Preencher os valores ausentes com 1 (nao é necessario)
# merged_data[is.na(merged_data$Adjusted.P.value), "Adjusted.P.value"] <- 1

# Ordenar os dados de acordo com o "Adjusted P-value"
merged_data <- merged_data %>%
  arrange(Adjusted.P.value) %>%
  mutate(Term = factor(Term, levels = unique(Term)))

order_vector <- c("OLFM1", "IL33", "VSIG4", "C9", "HP", "GCH1", "NRG1", "DBH", "TRIM71", "WDR62", "FER1L4")

# Ordenar a coluna Genes com base no vetor order_vector
merged_data$Genes <- factor(merged_data$Genes, levels = order_vector)

# Plot do heatmap
p22 <- ggplot(merged_data, aes(x = Genes,  # A ordem será respeitada conforme definida no fator
                               y = Term,
                               fill = Adjusted.P.value)) +
  geom_tile(color = "black", na.rm = FALSE) +  # Adiciona contornos pretos, incluindo valores NA
  coord_flip() +
  scale_fill_gradient(low = "#D71313",
                      high = "#93BFCF",
                      na.value = "white") +  # Define a cor para valores NA como branco
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, colour = "black"),  # Define o tamanho da fonte para o eixo X
        axis.text.y = element_text(size = 12, colour = "black")) +  # Define o tamanho da fonte para o eixo Y
  labs(x = "Genes", y = "Term", fill = "Adj.p-value")

ggsave(filename = "Diseasome_TCGA_Genes.tiff",
       plot = p22,
       device = "tiff", 
       width = 12, 
       height = 5, 
       units = "in", 
       dpi = 300) 

write.xlsx(merged_data, "Suppl.Table5d_InputDiseasomeHeatmap.xlsx", rowNames = FALSE)

# LDA ---------------------------------------------------------------------
#Barplot
library(ggplot2)

# Criar os dados para a tabela
data <- data.frame(
  Tissue = c("Liver", "PBMC", "Liver", "PBMC", "Liver"),
  Hepatitis = c("HBV", "HBV", "HCV", "HCV", "HDV"),
  Group = c("HBV Liver", "HBV PBMC", "HCV Liver", "HCV PBMC", "HDV Liver"),
  Accuracy = c(88, 60, 80, 73, 60)
)

# Criar o gráfico
# Alterar a ordem dos níveis da variável Hepatitis
data$Group <- factor(data$Group, levels = c("HDV Liver", "HCV PBMC", "HCV Liver",
                                            "HBV PBMC", "HBV Liver"))

# Criar o gráfico com a nova ordem
LDAbarplot <- ggplot(data, aes(x = Accuracy, y = Group, fill = Hepatitis, label = paste0(Accuracy, "%"))) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
  geom_text(hjust = -0.2, size = 5) + # Adiciona os valores de precisão nas barras
  scale_fill_manual(values = c("HBV" = "#DDF1F5", "HCV" = "#EDAEAE", "HDV" = "#FFD09B")) +
  theme_minimal() +
  labs(
    title = "Barplot with model Accuracy value",
    x = "Model Accuracy (%)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12, colour = "black"),
    axis.title.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12, colour = "black"),
    legend.text = element_text(size = 12, colour = "black"),
    strip.text = element_blank(),
    panel.spacing = unit(0.6, "lines"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),  # Linhas no eixo Y
    panel.grid.minor.y = element_blank(),  # Desativa as linhas menores
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed")  # Linhas no eixo X
  ) +
  coord_cartesian(xlim = c(0, 100)) # Define o limite do eixo x para 0-100

ggsave(filename = "LDAbarplot.tiff",
       plot = LDAbarplot,
       device = "tiff", 
       width = 5, 
       height = 5, 
       units = "in", 
       dpi = 300)


# Supplementary Fig. 4c ---------------------------------------------------
# Combine data
df <- bind_rows(EnriquecimentoLiver, EnriquecimentoPBMC)

# Clean Term names for display
df$Term <- str_to_title(gsub(".*?\\) ", "", df$Term))

# Select top 10 terms per system and tissue based on lowest Adjusted.P.value
top_terms <- df %>%
  group_by(tissue, sytem) %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 10) %>%
  ungroup()

# Filter original dataframe to only include top terms
df_filtered <- df %>%
  filter(Term %in% top_terms$Term)

# Calculate GeneRatio
top_terms$GeneRatio <- as.numeric(sub("/.*", "", top_terms$Overlap)) /
  as.numeric(sub(".*/", "", top_terms$Overlap))

# Prepare plot
ggplot(top_terms, aes(x = hepatites, y = Term)) +
  geom_point(aes(size = GeneRatio, color = Adjusted.P.value, shape = sytem)) +
  scale_color_gradient(low = "red", high = "lightblue", name = "P-value") +
  scale_shape_manual(values = c("Immune" = 16, "Nervous" = 18), name = "System") +
  facet_wrap(~ tissue, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 10, face = "bold")) +
  labs(size = "GeneRatio", x = NULL, y = NULL)

