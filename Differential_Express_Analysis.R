# Script for GEO/NCBI Microarray datasets Differential Expression --------
# Referencia: https://sbc.shef.ac.uk/geo_tutorial/tutorial.nb.html
# Instalar pacotes
install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery", force = T)
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")
install.packages("xml2")
#Buscar dataset direto do GEO----
library(GEOquery)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(limma)
library(dplyr)

# com essa função voce consegue importar todos os dados disponibilizados pelo estudo. Dentre eles as anotaçoes e metadados de IDs
my_id = "GSE107170"
gse = getGEO(my_id)
##quantas plataformas de sequenciamento foram usadas----
# Alguns conjuntos de dados no GEO podem ser derivados de diferentes plataformas de microarray
# aqui voce pode ver quantas plataformas de sequenciamento foram usadas e selecionar uma
length(gse) 
gse = gse[[1]]
seriesmatrix = pData(gse) # informação da amostra
data = exprs(gse) ## salvando dados de expressão

# data = as.data.frame(data)
# data[104:107] = NULL
# data = as.matrix(data)

#Verificar normalização
#Usando as informaçoes de metadados
seriesmatrix = pData(gse)
View(seriesmatrix)
colnames(seriesmatrix)
print(seriesmatrix$`subject condition:ch1`)
print(unique(seriesmatrix$`infection outcome:ch1`))

print((seriesmatrix$`subject condition:ch1`)) #HBV liner #46
print((seriesmatrix$source_name_ch1))
print(unique(seriesmatrix$source_name_ch1))

sampleInfo <- data.frame(id = 1:307)  # Substitua isto com os detalhes reais dos seus dados
sampleInfo$condition <- c(rep("HBV_Infected", 39), rep("HBV_Control", 81),
                          rep("HBV_Infected", 10), rep("HBV_Control", 10),
                          rep("HDV_Infected", 12), rep("HDV_Control", 23),
                          rep("HDV_Infected", 29), rep("HDV_Infected", 5),
                          rep("HDV_Control", 5), rep("HCV_Control", 3),
                          rep("HCV_Infected", 5), "HCV_Control",
                          rep("HCV_Infected", 6), "HCV_Control",
                          "HCV_Infected", rep("HCV_Control", 5),
                          rep("HCV_Infected", 5), rep("HCV_Control", 4),
                          rep("HCV_Infected", 5), rep("HCV_Control", 8),
                          rep("HCV_Infected", 10), rep("HCV_Control", 2),
                          rep("HCV_Infected", 5), rep("HCV_Control", 3),
                          rep("HCV_Infected", 5), rep("HCV_Control", 2),
                          rep("HCV_Control", 2), rep("HCV_Control", 2),
                          "HCV_Infected", "HCV_Control", "HCV_Infected",
                          "HCV_Control", "HCV_Infected", "HCV_Control",
                          "HCV_Infected", "HCV_Control", "HCV_Infected",
                          "HCV_Control", "HCV_Infected", "HCV_Control",
                          "HCV_Infected", "HCV_Control", "HCV_Infected",
                          "HCV_Control", "HCV_Infected", "HCV_Control")
sampleInfo$id = NULL
View(sampleInfo)
print(unique(sampleInfo$condition))
condition = as.data.frame(factor(sampleInfo$condition),
                          levels = c("HBV_Infected",
                                     "HBV_Control",
                                     "HDV_Infected",
                                     "HDV_Control",
                                     "HCV_Control",
                                     "HCV_Infected"))
#Montando sampleTable
condition <- as.data.frame(factor(seriesmatrix[[38]],
                                  levels = c("Infected", #Inserir Grupo Infectado
                                             "Control"))) #Inserir Grupo Controle
View(condition)
condition <- as.data.frame(factor(c(rep("Infected", 20),
                                    rep("Control", 8),
                                    rep("Infected", 8),
                                    rep("Control", 4),
                                    rep("Infected", 8),
                                    rep("Control", 8)),
                                  levels = c("Infected",
                                             "Control")))

# condition <- as.data.frame(condition[-c(104:107), ])
colnames(condition) = "condition"
rownames(condition) <- colnames(data)
View(condition)
#Batch
batch = as.data.frame((factor((seriesmatrix[[8]]))))
# batch <- as.data.frame(batch[-c(104:107), ])
colnames(batch) = "batch"
rownames(batch) <- colnames(data)
View(batch)

#other
age = as.data.frame((factor((seriesmatrix[[33]]))))
colnames(age) = "age"
rownames(age) <- colnames(data)

View(batch)
#Sample Table padrao
sampleTable = data.frame(condition = as.factor(condition[[1]]))
View(sampleTable)
# sampleTable$condition[7] = "tumoral"

#Batch
sampleTable = data.frame(condition = as.factor(condition[[1]]),
                         batch = as.factor(batch[[1]]))
View(sampleTable)

#####componentes principais----
pca = prcomp(t(exprs(gse)))
raw_thing = as.matrix(data)
seobj_adj1 <- SummarizedExperiment(assays=raw_thing,
                                   colData=sampleTable)
pca_obj_adj1 <- plotPCA(DESeqTransform(seobj_adj1),
                        intgroup=c("condition"))
#Com Batch
pca_obj_adj1 <- plotPCA(DESeqTransform(seobj_adj1),
                        intgroup=c("condition",
                                   "batch"))

PCA_rawdata <- ggplot(pca_obj_adj1$data,
                      aes(x = PC1,
                          y = PC2,
                          color = condition,
                          shape = condition)) +
  geom_point(size = 8, aes(alpha = 0.7)) +  # Adjust alpha for transparency
  stat_ellipse(aes(x = PC1, y = PC2,
                   color = condition),
               type = "t", level = 0.95) +
  labs(
    x = sprintf("PC1: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[1])),
    y = sprintf("PC2: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[2])),
    title = "RAW DATA"
  ) +
  theme_bw(18)

ggsave("PCA_bruta_Counts.svg", 
       plot=PCA_rawdata, 
       device = "svg", 
       width = 10, 
       height = 8)

# ---------------   DensityPlot    ------------------
conditions <- colData(seobj_adj1)$condition
# Criando um density plot com ggplot2
print(unique(conditions))
Raw_Density = ggplot(data.frame(condition = conditions),
                     aes(x = condition,
                         fill = condition)) +
  geom_density(alpha = 0.7) +
  labs(
    title = "Density Plot",
    x = "Condition",
    y = "Density"
  ) +
  scale_fill_manual(values = c("HBV_Infected" = "tomato", #Infectado
                               "HBV_Control" = "skyblue",
                               "HDV_Infected" = "red",
                               "HDV_Control" = "blue",
                               "HCV_Control" = "#4F709C",
                               "HCV_Infected" = "#9A3B3B")) +  # Controle
  theme_bw(18)

#Observaca: caso haja mais de duas condicoes, usar "#86A789" como terceira cor

ggsave("Density_Counts.svg", 
       plot=Raw_Density, 
       device = "svg", 
       width = 10, 
       height = 8)

# ---------------   BoxPlot.    ------------------
gs = conditions
print(unique(gs))
#Cor com DUAS condicoes
colors_by_group = c("dead" = "tomato", #Infectado
                    "survivial" = "skyblue")
#Cor com TRES condicoes
colors_by_group <- c("HBV_Infected" = "tomato", #Infectado
                     "HBV_Control" = "skyblue",
                     "HDV_Infected" = "red",
                     "HDV_Control" = "blue",
                     "HCV_Control" = "#4F709C",
                     "HCV_Infected" = "#9A3B3B")
colors_for_plot = colors_by_group[as.character(conditions)]

ord <- order(gs) # Ordenando as amostras por grupo
par(mar=c(7,4,2,1)) # Configurando as margens do gráfico
boxplot(data[, ord],
        boxwex = 0.6,
        notch = TRUE,
        outline = FALSE,
        las = 2, col = colors_for_plot[ord])

#Salvar os dados de expressao----
counts = as.data.frame(data)
write.xlsx(counts,
           "DataCounts_Original_Batch.xlsx",
           rowNames = TRUE)

# ANALISE DE EXPRESSAO DIFERENCIAL ----
#manipule os dados de acordo com a comparação que voce quer fazer
#i.e., selecione linhas e colunas de seu interesse, mas tome cuidado :) 
# agora vamos criar uma matriz de desing a partir das colunas do arquivo
# aqui eu sugiro que independente da comparação que voce faça voce nao mexa nas colunas, no entanto choices!
# ps: no exemplo eu nao vou mexer
#   para entendimento: Se voce usa uma comparação controle com outra condiçao X que nao é a unica do estudo, 
#     seus resultados podem ser um viés
## Renomear os nomes de sampleInfo com sampleTable
print(unique(sampleTable$condition))
View(sampleTable)
sampleInfo <- as.data.frame(c(rep("Liver biopsy in inactive carrier of chronic hepatitis B patients", 11),
                              rep("Liver biopsy in chronic hepatitis B patients in immune clearance phase", 48),
                              rep("Liver biopsy in chronic hepatitis B patients in immune tolerant phase", 24)))
row.names(sampleInfo) <- row.names(sampleTable)
colnames(sampleInfo) <- colnames(sampleTable)
colnames(sampleInfo) <- "condition"
sampleInfo = sampleTable
## Criando o Design
design = model.matrix(~0+sampleInfo$condition) #Padrao (Sem Batch)
design <- model.matrix(~0 + condition:batch, data = sampleTable) #Com Batch

print(colnames(design))
gp <- as.list(c("Infected_Myeloid", #Definindo novos nomes pra colunas de Design
                "Control_Myeloid",
                "Infected_Plasmocytoid",
                "Control_Plasmocytoid"))
gp = as.list(c(colnames(design)))
print(gp)
colnames(design) <- gp #Renomeando Colunas de Design
print(colnames(design))

fit <- lmFit(counts, design) #Criando o modelo Linear
exprMatrix_noBatch <- removeBatchEffect(counts, batch$batch) #remocao de batch

print(colnames(design))
cont.matrix <- makeContrasts(Infected_Plasmocytoid - Control_Plasmocytoid, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit, 0.01) #semcontraste

fit2 <- eBayes(fit2, 0.01)
limao <- topTable(fit2,
                  adjust="fdr", #p-values corrigidos para controle de FDR
                  sort.by="B", # Tamanho do ajuste pelo modelo bayesiano.
                  number=Inf) # extraindo todos os resultados 
View(limao)
write.xlsx(limao,
           "ExpressionResults.xlsx",
           row.names = T)

##### PCA sem Batch----
raw_thing2 = as.matrix(exprMatrix_noBatch)
seobj_adj2 <- SummarizedExperiment(assays=raw_thing2,
                                   colData=sampleTable)
pca_obj_adj2 <- plotPCA(DESeqTransform(seobj_adj2),
                        intgroup=c("condition",
                                   "batch"))

PCA_WithoutBatch <- ggplot(pca_obj_adj2$data,
                           aes(x = PC1,
                               y = PC2,
                               color = condition,
                               shape = condition)) +
  geom_point(size = 8, aes(alpha = 0.7)) +  # Adjust alpha for transparency
  stat_ellipse(aes(x = PC1, y = PC2,
                   color = condition),
               type = "t", level = 0.95) +
  labs(
    x = sprintf("PC1: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[1])),
    y = sprintf("PC2: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[2])),
    title = "RAW DATA"
  ) +
  theme_bw(18)

ggsave("PCA_WithoutBatch.svg", 
       plot=PCA_WithoutBatch, 
       device = "svg", 
       width = 10, 
       height = 8)

#Salvando tabela sem Batch
NobatchCounts = as.data.frame(exprMatrix_noBatch)
write.xlsx(NobatchCounts,
           "exprMatrix_noBatch.xlsx",
           row.names = T)


# DESeq2 used on GEO/NCBI Datasets (Bulk-RNAseq) ----------------------------------------
#Análise de expressão diferencial (DESeq2) - By Adriel Leal
# ============ Bibliotecas ===========

#Carregar bibliotecas
{
  library(openxlsx)
  library(org.Hs.eg.db)
  library(DESeq2)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GEOquery)
  library(ggplot2)
  library(ggrepel)
  library(rrvgo)
  library(tidyr)
  library(scales)
  library(sva)
}
{
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("sva")
}
#Criando Diretorio
{
  outDir <- "atl career"
  dir.create(outDir)
}

# ============ Importando Dados ===========
##Importar dado do PC
data = openxlsx::read.xlsx("Samples atl career.xlsx", colNames = T, rowNames = T)
data = GSE164266

##Importar do GEO
#Counts Table
urld = "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path = paste(urld, "acc=GSE173897", #Substitua pelo seu GSE
             "file=GSE173897_raw_counts_GRCh38.p13_NCBI.tsv.gz", #Substitua pelo seu GSE
             sep="&");
#Salvando GSE em um dataframe
GSE173897 = as.data.frame(data.table::fread(path,
                                            header=T,
                                            colClasses="integer"),
                          rownames="GeneID")
data = GSE173897
rownames(data) = data$GeneID
data[1] = NULL
View(data)
#Baixando Series Matrix
gse_code = "GSE173897" #Substituir pelo seu GSE
geo_url = paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
                gse_code,
                sep="")
seriesmatrix = as.data.frame(getGEO(gse_code,
                                    destdir = ".",
                                    getGPL = FALSE))
View(seriesmatrix)
#Definir quais são os grupos com base no seriesMatrix 
#Levels importa para definir qual será UP (positivo na expressão diferencial)
#Definindo Condicoes
print(colnames(seriesmatrix))
print(unique(seriesmatrix$GSE173897_series_matrix.txt.gz.hbv.status.ch1))
condition = as.data.frame((factor((seriesmatrix[[48]]),
                                  levels = c("Active", #Inserir Grupo "Infectado"
                                             "Inactive")))) #Inserir Grupo "Controle"
colnames(condition) = "condition"
rownames(condition) <- colnames(data)
dplyr::glimpse(condition) #Verificar propriedade das variaveis

# ============ Se necessario - Definindo Batch ===========
#Clones | Nao precida de Levels
clones = as.data.frame((factor((seriesmatrix[[39]]))))
colnames(clones) = "clones"
rownames(clones) <- colnames(data)

#Batch
batch = as.data.frame((factor((seriesmatrix[[47]]))))
colnames(batch) = "batch"
rownames(batch) <- colnames(data)

#Sex
sex = as.data.frame((factor((seriesmatrix[[47]]))))
colnames(sex) = "sex"
rownames(sex) <- colnames(data)

dplyr::glimpse(clones) #Verificar propriedade das variaveis

#Carregando dados
data.counts = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/TCGA counts.xlsx", rowNames = TRUE)
data.series = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/TCGA final filtered.xlsx")
data.counts.trans = as.data.frame(t(data.counts))
View(data.counts.trans)
data.counts.trans <- rownames_to_column(data.counts.trans, var = "samples")
data.merge = merge(data.counts.trans, data.series,
                   by.x = "samples", by.y = "sampleID")
data.series = data.merge[,c(1,20532:20546)]
data.counts.trans = data.merge[,c(1:20531)]
row.names(data.counts.trans) = data.counts.trans$samples
samples = data.counts.trans$samples
data.counts.trans = data.counts.trans[,c(2:20531)]
row.names(data.counts.trans) = data.merge$samples
data.normal = as.data.frame(t(data.counts.trans))
data.normal.save = data.normal
data.normal.save <- rownames_to_column(data.normal.save,
                                       var = "genes")

#Salvando
write.xlsx(data.normal.save,
           "countsmatriz Editada.xlsx")
data.counts = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/countsmatriz Editada.xlsx", rowNames = TRUE)
data.series = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/TCGA final filtered.xlsx")

colnames(data.series)
# ============ Criando SampleTable ===========
#Sample Table padrao
sampleTable = data.frame(condition = as.factor(condition[[1]]))
#Se houver batch (clones)
sampleTable = data.frame(condition = as.factor(condition[[1]]),
                         batch = as.factor(batch[[1]]))
#Mais de um batch
sampleTable = data.frame(condition = as.factor(condition[[1]]),
                         clones = as.factor(clones[[1]]),
                         batch = as.factor(batch[[1]]))

rownames(sampleTable) = colnames(data)
glimpse(data)
View(sampleTable)
##============== Combat-seq ===========
count_matrix = as.matrix(data)
batch = batch
condition = condition
adjusted = ComBat_seq(counts = count_matrix,
                      batch = batch,
                      group = NULL,
                      shrink = F)

##========  Vizualizacao de dados  ===========
# ---------------   PCA.    ------------------
#raw_thing <- scale(count_matrix) #Escalonando dados
raw_thing <- count_matrix #Sem Escalonar dados
#ComBatch:
# col_data <- data.frame(Batch = factor(batch), Group = factor(condition))
# rownames(col_data) <- colnames(count_matrix)
#SemBatch
# col_data <- data.frame(Group = factor(condition))
# col_data = data.frame(Group = as.factor(condition))
seobj_adj1 <- SummarizedExperiment(assays=raw_thing,
                                   colData=sampleTable)
#Com Batch
# pca_obj_adj1 <- plotPCA(DESeqTransform(seobj_adj1),
#                         intgroup=c("Batch", "Group"))
#Sem Batch
pca_obj_adj1 <- plotPCA(DESeqTransform(seobj_adj1),
                        intgroup=c("condition"))

PCA_rawdata <- ggplot(pca_obj_adj1$data,
                      aes(x = PC1,
                          y = PC2,
                          color = condition,
                          shape = condition)) +
  geom_point(size = 8, aes(alpha = 0.7)) +  # Adjust alpha for transparency
  stat_ellipse(aes(x = PC1, y = PC2,
                   color = condition),
               type = "t", level = 0.95) +
  labs(
    x = sprintf("PC1: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[1])),
    y = sprintf("PC2: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[2])),
    title = "RAW DATA"
  ) +
  theme_bw(18)

ggsave("PCA_bruta_Counts.svg", 
       plot=PCA_rawdata, 
       device = "svg", 
       width = 10, 
       height = 8)

# ---------------   DensityPlot    ------------------
conditions <- colData(seobj_adj1)$condition
# Criando um density plot com ggplot2
print(unique(conditions))
Raw_Density = ggplot(data.frame(condition = conditions),
                     aes(x = condition,
                         fill = condition)) +
  geom_density(alpha = 0.7) +
  labs(
    title = "Density Plot",
    x = "Condition",
    y = "Density"
  ) +
  scale_fill_manual(values = c("Tumor tissue" = "tomato", #Infectado
                               "Non-neoplastic liver tissue" = "skyblue")) +  # Controle
  theme_bw(18)

ggsave("Density_Counts.svg", 
       plot=Raw_Density, 
       device = "svg", 
       width = 10, 
       height = 8)

# ---------------   BoxPlot.    ------------------
gs = conditions
print(gs)
colors_by_group <- ifelse(gs == "Tumor tissue",
                          "tomato",
                          "skyblue")
#Com mais de uma condicao
colors_for_plot <- colors_by_group[as.character(conditions)]
colors_by_group <- c("Tumor tissue" = "tomato", #Infectado
                     "Non-neoplastic liver tissue" = "skyblue")
colors_for_plot <- colors_by_group[as.character(conditions)]

ord <- order(gs) # Ordenando as amostras por grupo
par(mar=c(7,4,2,1)) # Configurando as margens do gráfico
boxplot(data[, ord],
        boxwex = 0.6,
        notch = TRUE,
        outline = FALSE,
        las = 2, col = colors_for_plot[ord])
#Salve Manualmente como pdf

# ---------------   DESeq2    ------------------
### Criar o DESeqDataSet
levels(sampleTable$condition) <- make.names(levels(sampleTable$condition))
deseq = DESeqDataSetFromMatrix(countData = data,
                               colData = sampleTable,
                               design = ~condition)
d.deseq = DESeq(deseq)
res_deseq = results(d.deseq) #Salvando Resultado do DESeq
resOrdered = as.data.frame(res_deseq[order(res_deseq$pvalue),]) #Ordenando por pvalor
View(resOrdered)
#Salvando tabela
resOrdered$geneID = rownames(resOrdered)
write.xlsx(resOrdered,
           "DESeq2Resultado.xlsx",
           rowNames = TRUE)
write.xlsx(data,
           "CountsWITHBatch.xlsx",
           rowNames = TRUE)

##PCA (Com Batch Effect)
vst_before = vst(deseq)
pca_before = plotPCA(vst_before, intgroup = "condition")
PCA_Deseq2 <- ggplot(pca_before$data,
                     aes(x = PC1,
                         y = PC2,
                         color = condition,
                         shape = condition)) +
  geom_point(size = 8, aes(alpha = 0.7)) +  # Adjust alpha for transparency
  stat_ellipse(aes(x = PC1, y = PC2,
                   color = condition),
               type = "t", level = 0.95) +
  labs(
    x = sprintf("PC1: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[1])),
    y = sprintf("PC2: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[2])),
    title = "RAW DATA"
  ) +
  theme_bw(18)

ggsave("PCA_DESeq2.svg", 
       plot=PCA_Deseq2, 
       device = "svg", 
       width = 10, 
       height = 8)

# ---------------   DESeq sem Batch Effect    ------------------
## Removendo o efeito de Batch (Combat)
colData(deseq)$batch <- batch$batch #Adicionando informação do batch ao DESeqDataSet
count_matrix <- as.matrix(counts(deseq)) #Convertendo o DESeqDataSet para uma matriz de expressão
combat_result <- ComBat(count_matrix, batch = colData(deseq)$batch) # Remoção do efeito de batch usando combat
# Ajusta a matriz resultante para garantir que os valores sejam inteiros não negativos
count_matrix_no_batch <- round(combat_result)
count_matrix_no_batch[count_matrix_no_batch < 0] <- 0
View(count_matrix_no_batch)
count_matrix_no_batch = as.data.frame(count_matrix_no_batch)
count_matrix_no_batch$geneID = rownames(count_matrix_no_batch)
write.xlsx(count_matrix_no_batch,
           "DataCounts_Without_Batch.xlsx",
           rowNames = TRUE)

## Criando um novo DESeqDataSet com a matriz sem o efeito de batch
count_matrix_no_batch$geneID = NULL
deseq_no_batch <- DESeqDataSetFromMatrix(
  countData = count_matrix_no_batch,
  colData = colData(deseq),
  design = ~condition + batch
)

#PCA (Sem Batch Effect)
vst_after <- vst(deseq_no_batch)
pca_after <- plotPCA(vst_after, intgroup = c("condition",
                                             "batch"))

PCA_AfterBatch_Deseq2 = ggplot(pca_after$data,
                               aes(x = PC1,
                                   y = PC2,
                                   color = batch,
                                   shape = condition)) +
  geom_point(size = 8, aes(alpha = 0.7)) +  # Adjust alpha for transparency
  stat_ellipse(aes(x = PC1, y = PC2,
                   color = condition),
               type = "t", level = 0.95) +
  labs(
    x = sprintf("PC1: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[1])),
    y = sprintf("PC2: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[2])),
    title = "RAW DATA"
  ) +
  theme_bw(18)

ggsave("PCA_DESeq2_NoBatch.svg", 
       plot=PCA_AfterBatch_Deseq2, 
       device = "svg", 
       width = 10, 
       height = 8)

#Salvando resultados sem Batch Effect
d.deseq_batch = DESeq(deseq_no_batch)
res_deseq_batched = results(d.deseq_batch) #Salvando Resultado do DESeq
resOrdered_nobatch = as.data.frame(res_deseq_batched[order(res_deseq_batched$pvalue),]) #Ordenando por pvalor
View(resOrdered_nobatch)

write.xlsx(resOrdered_nobatch,
           "DESeq2NoBatchResult.xlsx",
           rowNames = TRUE)
# ---------------   Inserindo GeneSymbol    ------------------
#salvando a coluna de genes para usar de referencia
resOrdered$genes = row.names(resOrdered)
gene_ids <- as.character(c(resOrdered$genes))
# Carregar a base de dados de anotação para o organismo de interesse
org <- "org.Hs.eg.db"
db <- org.Hs.eg.db
keytypes(org.Hs.eg.db) #consultar os identificadores
gene_names <- mapIds(db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID")
resOrdered$genes_symbol = gene_names
resOrdered = na.omit(resOrdered)

#Filtrando Dados (Sig, Up e Down)
resSig <- as.data.frame(subset(resOrdered, pvalue < 0.05))
resl2fc <- as.data.frame(resSig[order(resSig$log2FoldChange),])
resl2fc_up <- as.data.frame(subset(resl2fc, log2FoldChange > 1))
resl2fc_down <- as.data.frame(subset(resl2fc, log2FoldChange < 1))

#Salvando Tabelas
write.xlsx(resOrdered, "DEGs.xlsx")
write.xlsx(resSig, "DEGs_Sig.xlsx")
write.xlsx(resl2fc, "DEGs_fcSig.xlsx")
write.xlsx(resl2fc_up, "DEGs_Up.xlsx")
write.xlsx(resl2fc_down, "DEGs_Down.xlsx")
write.xlsx(data, "Counts_table.xlsx", rowNames = TRUE)

#Alternativa ao Annotation: biomart mas que retorna a mesma proporção
#columns(org.Hs.ed.db)
table$entrez <- mapIds(org.Hs.eg.db, 
                       keys = rownames(table), 
                       keytype = "SYMBOL",
                       column = "ENTREZID")
length(na.exclude(table$entrez))
openxlsx::write.xlsx(as.data.frame(table$entrez),file.path(outDir, "entrez_up.xlsx"))

#down
table <- openxlsx::read.xlsx(file.path(outDir,"transcripts_down.xlsx"), rowNames = T)
table = as.data.frame(resl2fc_down) #de 704 for pra 581...

transcript_ids <- rownames(table)

table$entrez <- mapIds(org.Hs.eg.db, 
                       keys = rownames(table), 
                       keytype = "SYMBOL",
                       column = "ENTREZID")
length(na.exclude(table$entrez))
openxlsx::write.xlsx(as.data.frame(table$entrez),file.path(outDir,"entrez_down.xlsx"))



# DESeq2 TCGA -------------------------------------------------------------
#Análise de expressão diferencial (DESeq2) - By Adriel Leal
# ============ Bibliotecas ===========

#Carregar bibliotecas
{
  library(openxlsx)
  library(org.Hs.eg.db)
  library(DESeq2)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GEOquery)
  library(ggplot2)
  library(ggrepel)
  library(rrvgo)
  library(tidyr)
  library(scales)
  library(sva)
  library(ComplexHeatmap)
  library(circlize)
}
{
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("sva")
}

setwd("/Users/adrielnobile/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/")
#Carregando dados
data.counts = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/countsmatriz Editada.xlsx", rowNames = TRUE)
data.series = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/TCGA final filtered.xlsx")
series.total = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/TCGA clinical.xlsx")
XenaData = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/00_tblXenaHubInfo.csv")
#Filtrando tabela de series
names(series.total)
varClinKeep = c("sampleID","cohort", "primarydisease", "gender",
                "additionalpharmaceuticaltherapy", "additionalradiationtherapy",
                "ageatinitialpathologicdiagnosis", "daystodeath",
                "histhepatocarcfact", "neoplasmhistologicgrade",
                "pathologicstage", "personneoplasmcancerstatus",
                "radiationtherapy", "sampletype", "viralhepatitisserology",
                "vitalstatus", "days_to_last_followup", "days_to_death");
clinDF01 = as.data.frame(do.call(cbind, filterTCGA03));
clinFinal = clinDF01[varClinKeep];

NA -> clinFinal[clinFinal == ""];
clinFinal$neoplasmhistologicgrade <- ifelse(is.na(clinFinal$neoplasmhistologicgrade),
                                            "GX",
                                            clinFinal$neoplasmhistologicgrade)

colnames(data.series)
data = data.counts

#Separando grupos
data.trans = as.data.frame(t(data.counts))
sampletable = data.series[,c(1,10)]
data.trans <- rownames_to_column(data.trans, var = "sampleID")
data.info <- left_join(data.trans, sampletable, by = "sampleID")
unique(data.info$neoplasmhistologicgrade)

# Filtrando os dados para G1 e G2
countsG12 <- data.info[data.info$neoplasmhistologicgrade %in% c("G1", "G2"), ]
# Filtrando os dados para G1 e G3
countsG13 <- data.info[data.info$neoplasmhistologicgrade %in% c("G1", "G3"), ]
# Filtrando os dados para G1 e G4
countsG14 <- data.info[data.info$neoplasmhistologicgrade %in% c("G1", "G4"), ]
#Testando
unique(countsG12$neoplasmhistologicgrade)
unique(countsG13$neoplasmhistologicgrade)
unique(countsG14$neoplasmhistologicgrade)

#Criando sampletables
sample12 = countsG12[20532]
sample13 = countsG13[20532]
sample14 = countsG14[20532]

countsG12 <- countsG12[, !(names(countsG12) == "neoplasmhistologicgrade")]
countsG13 <- countsG13[, !(names(countsG13) == "neoplasmhistologicgrade")]
countsG14 <- countsG14[, !(names(countsG14) == "neoplasmhistologicgrade")]

sample12$samples = countsG12$sampleID
sample13$samples = countsG13$sampleID
sample14$samples = countsG14$sampleID

#Definindo como row.names
rownames(countsG12) = sample12$samples
rownames(countsG13) = sample13$samples
rownames(countsG14) = sample14$samples

#Removendo colunas
countsG12 <- countsG12[, !(names(countsG12) == "sampleID")]
countsG13 <- countsG13[, !(names(countsG13) == "sampleID")]
countsG14 <- countsG14[, !(names(countsG14) == "sampleID")]

#Transpondo colunas
countsG12 = as.data.frame(t(countsG12))
countsG13 = as.data.frame(t(countsG13))
countsG14 = as.data.frame(t(countsG14))

#Salvando tabelas DESeq2
write.xlsx(countsG12,
           "countsG12.xlsx", rowNames = TRUE)
write.xlsx(countsG13,
           "countsG13.xlsx", rowNames = TRUE)
write.xlsx(countsG14,
           "countsG14.xlsx", rowNames = TRUE)
#Salvando sampleTable
write.xlsx(sample12,
           "samplesG12.xlsx")
write.xlsx(sample13,
           "samplesG13.xlsx")
write.xlsx(sample14,
           "samplesG14.xlsx")

#Carregando tabelas e samplestables
#counts
countsG12 = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/countsG12.xlsx", rowNames = TRUE)
countsG13 = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/countsG13.xlsx", rowNames = TRUE)
countsG14 = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/countsG14.xlsx", rowNames = TRUE)

#sampletables
sample12 = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/samplesG12.xlsx")
sample13 = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/samplesG13.xlsx")
sample14 = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/10 - TCGA/samplesG14.xlsx")

#Retirando log2
original_counts <- 2^countsG13 - 1
original_counts = round(original_counts)
data = original_counts
# ============ Criando SampleTable ===========
#Sample Table padrao
sampleTable = data.frame(condition = as.factor(sample13[[1]]))
condition = data.frame(condition = (sample13[[1]]))
#Se houver batch (clones)
sampleTable = data.frame(condition = as.factor(condition[[1]]),
                         batch = as.factor(batch[[1]]))
#Mais de um batch
sampleTable = data.frame(condition = as.factor(condition[[1]]),
                         clones = as.factor(clones[[1]]),
                         batch = as.factor(batch[[1]]))
glimpse(data)
View(sampleTable)
##============== Combat-seq ===========
count_matrix = as.matrix(data)
batch = batch
condition = condition
# adjusted = ComBat_seq(counts = count_matrix,
#                       batch = batch,
#                       group = NULL,
#                       shrink = F)

##========  Vizualizacao de dados  ===========
# ---------------   PCA.    ------------------
#raw_thing <- scale(count_matrix) #Escalonando dados
raw_thing <- data #Sem Escalonar dados
#ComBatch:
# col_data <- data.frame(Batch = factor(batch), Group = factor(condition))
# rownames(col_data) <- colnames(count_matrix)
#SemBatch
# col_data <- data.frame(Group = factor(condition))
# col_data = data.frame(Group = as.factor(condition))
seobj_adj1 <- SummarizedExperiment(assays=as.matrix(raw_thing),
                                   colData=sampleTable)
#Com Batch
# pca_obj_adj1 <- plotPCA(DESeqTransform(seobj_adj1),
#                         intgroup=c("Batch", "Group"))
#Sem Batch
pca_obj_adj1 <- plotPCA(DESeqTransform(seobj_adj1),
                        intgroup=c("condition"))
ggplot(pca_obj_adj1$data,
       aes(x = PC1,
           y = PC2,
           color = condition,
           shape = condition)) +
  geom_point(size = 8, aes(alpha = 0.7)) +  # Adjust alpha for transparency
  stat_ellipse(aes(x = PC1, y = PC2,
                   color = condition),
               type = "t", level = 0.95) +
  labs(
    x = sprintf("PC1: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[1])),
    y = sprintf("PC2: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[2])),
    title = "RAW DATA"
  ) +
  theme_bw(18)

ggsave("PCA_bruta_Counts.svg", 
       plot=PCA_rawdata, 
       device = "svg", 
       width = 10, 
       height = 8)

# ---------------   DensityPlot    ------------------
conditions <- colData(seobj_adj1)$condition
# Criando um density plot com ggplot2
print(unique(conditions))
ggplot(data.frame(condition = conditions),
       aes(x = condition,
           fill = condition)) +
  geom_density(alpha = 0.7) +
  labs(
    title = "Density Plot",
    x = "Condition",
    y = "Density"
  ) +
  scale_fill_manual(values = c("GX" = "tomato", #Infectado
                               "G4" = "skyblue",
                               "G3" = "steelblue",
                               "G2" = "darkred",
                               "G1" = "grey")) +  # Controle
  theme_bw(18)

ggsave("Density_Counts.svg", 
       plot=Raw_Density, 
       device = "svg", 
       width = 10, 
       height = 8)

# ---------------   BoxPlot.    ------------------
gs = conditions
print(gs)
# colors_by_group <- ifelse(gs == "Tumor tissue",
#                           "tomato",
#                           "skyblue")
#Com mais de uma condicao
colors_for_plot <- colors_by_group[as.character(conditions)]
colors_by_group <- c("Tumor tissue" = "tomato", #Infectado
                     "Non-neoplastic liver tissue" = "skyblue")
colors_for_plot <- colors_by_group[as.character(conditions)]

ord <- order(gs) # Ordenando as amostras por grupo
par(mar=c(7,4,2,1)) # Configurando as margens do gráfico
boxplot(data[, ord],
        boxwex = 0.6,
        notch = TRUE,
        outline = FALSE,
        las = 2, col = colors_for_plot[ord])
#Salve Manualmente como pdf

# ---------------   DESeq2    ------------------
### Criar o DESeqDataSet
levels(sampleTable$condition) <- make.names(levels(sampleTable$condition))
rownames(sampleTable) = colnames(data)
deseq = DESeqDataSetFromMatrix(countData = data,
                               colData = sampleTable,
                               design = ~condition)
d.deseq = DESeq(deseq)
res_deseq = results(d.deseq) #Salvando Resultado do DESeq
resOrdered = as.data.frame(res_deseq[order(res_deseq$pvalue),]) #Ordenando por pvalor
View(resOrdered)
colnames(resOrdered)
res_filtrados <- subset(resOrdered, log2FoldChange < -1 |
                          log2FoldChange > 1 &
                          padj < 0.05)
res_filtrados$diff =ifelse(res_filtrados$log2FoldChange > 0, 'Up',
                           'Down')
table(res_filtrados$diff)
#Salvando tabela
write.xlsx(res_filtrados,
           "G13_DESeq2Resultado.xlsx",
           rowNames = TRUE)
write.xlsx(resOrdered,
           "G13deseqTotal.xlsx",
           rowNames = TRUE)
##---------    PCA.   -------------
#Caso necessário: reunir resultados do objeto deseq


##PCA
vst_before = vst(deseq)
pca_before = plotPCA(vst_before, intgroup = "condition")
raw_thing <- data #Sem Escalonar dados
seobj_adj1 <- SummarizedExperiment(assays=as.matrix(raw_thing),
                                   colData=sampleTable)
pca_obj_adj1 <- plotPCA(DESeqTransform(seobj_adj1),
                        intgroup=c("condition"))
ggplot(pca_before$data,
       aes(x = PC1,
           y = PC2,
           color = condition,
           shape = condition)) +
  geom_point(size = 8, aes(alpha = 0.7)) +  # Adjust alpha for transparency
  stat_ellipse(aes(x = PC1, y = PC2,
                   color = condition),
               type = "t", level = 0.95) +
  labs(
    x = sprintf("PC1: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[1])),
    y = sprintf("PC2: %s Variance",
                percent(pca_obj_adj1$plot_env$percentVar[2])),
    title = "RAW DATA"
  ) +
  theme_bw(18)
ggsave("PCA_DESeq2.svg", 
       plot=PCA_Deseq2, 
       device = "svg", 
       width = 10, 
       height = 8)

