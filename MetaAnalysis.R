# Librarys ----------------------------------------------------------------

BiocManager::install("MetaVolcanoR", eval = FALSE)

library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(MetaVolcanoR)
library(UpSetR)

# In Vitro Metanalysis ----------------------------------------------------------------
setwd("~adrielnobile/Documents/DataSets/Metanalise array_bulk/InVitro/MetaVolcano/")

#A
AGSE13 <- read_excel("AAGSE130460.xlsx")
AGSE11 <- read_excel("BAGSE114916.xlsx")
AGSE23 <- read_excel("BAGSE234478.xlsx")
#B
BGSE11 <- read_excel("ABGSE118295.xlsx")
BGSE72 <- read_excel("ABGSE72068.xlsx")
BGSE12 <- read_excel("BBGSE126831.xlsx")
BGSE13 <- read_excel("BBGSE135860.xlsx")
#C
CGSE14 <- read_excel("ACGSE140114.xlsx")
CGSE29 <- read_excel("ACGSE29889.xlsx")
CGSE21 <- read_excel("BCGSE211161.xlsx")
#D
DGSE10 <- read_excel("ADGSE109824.xlsx")
DGSE11 <- read_excel("BDGSE112118.xlsx")
#E
EGSE53 <- read_excel("AEGSE53731.xlsx")
EGSE22 <- read_excel("BEGSE224795.xlsx")

#C
genes.22 <- as.data.frame(EGSE22)
genes.53 <- as.data.frame(EGSE53)
genes.53[2:5] = NULL
genes.10 = as.data.frame(na.omit(genes.10))

genes.10 <- as.data.frame(na.omit(DGSE10))
genes.11 <- as.data.frame(na.omit(DGSE11))

#=======================AnnotationDbi ------

# Vetor de IDs de gene para os quais você deseja obter os nomes externos
gene_ids <- c(genes.10$ENTREZ_GENE_ID)

# Carregar a base de dados de anotação para o organismo de interesse
org <- org.Hs.eg.db
db <- org.Hs.eg.db

keytypes(org.Hs.eg.db)# Realizar a transformação dos IDs de gene para nomes externos
gene_names.13 <- mapIds(db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBLTRANS")

gene_names.13 <- mapIds(db, keys = gene_ids, column = "ENSEMBL", keytype = "ENTREZID")

# Exibir o resultado
print(gene_names.11)
gene.ensembl <- as.data.frame(gene_names.11)
gene.ensembl$...1 <- rownames(gene.ensembl)
GSE_11 <- merge(gene.ensembl, GSE11, by = "...1")
GSE_11 <- na.omit(GSE_11)

head(GSE_11)
GSE_11 <- GSE_11[, -1]
colnames(genes.22)[colnames(genes.22) == "P.Value"] <- "pvalue"
colnames(genes.53)[colnames(genes.53) == "P.Value"] <- "pvalue"
print(colnames(genes.11))
colnames(genes.22)[colnames(genes.22) == "log2FoldChange"] <- "logFC"
#=======================Confidence Interval===============================
# Set the desired confidence level (e.g., 95%)
confidence_level <- 0.95

# Calculate the degrees of freedom (typically n - 1 for t-distribution)
df.22 <- nrow(genes.22) - 1

# Calculate the margin of error
margin_of_error.22 <- qt((1 + confidence_level) / 2, df.22) * genes.22$logFC / sqrt(nrow(genes.22))

# Calculate the lower and upper confidence limits
genes.22$lower_CI <- genes.22$logFC - margin_of_error.22
genes.22$upper_CI <- genes.22$logFC + margin_of_error.22

# Calculate the degrees of freedom (typically n - 1 for t-distribution)
df.53 <- nrow(genes.53) - 1

# Calculate the margin of error
margin_of_error.53 <- qt((1 + confidence_level) / 2, df.53) * genes.53$logFC / sqrt(nrow(genes.53))

# Calculate the lower and upper confidence limits
genes.53$lower_CI <- genes.53$logFC - margin_of_error.53
genes.53$upper_CI <- genes.53$logFC + margin_of_error.53

#Renomar para identificadores serem iguais: pvalue e log2FoldChange
print(colnames(genes.29))
colnames(genes.29)[colnames(genes.29) == "P.Value"] <- "pvalue"

print(colnames(genes.11))
colnames(genes.11)[colnames(genes.11) == "logFC"] <- "log2FoldChange"
#========================MetaVolc===================================

row.names(genes.11) = genes.11$Gene.Symbol
genes.11 <- subset(genes.11, !grepl(" ", genes.11) & !grepl(" ", as.character(Coluna2)))
backup = genes.11
genes.11$Gene.Symbol <- subset(genes.11$Gene.Symbol, !grepl(" ", genes.11$Gene.Symbol))
rownames(genes.11) = make.names(genes.11$Gene.Symbol, unique = T)
rownames(genes.12) = make.names(genes.12$Gene.Symbol, unique = T)
rownames(genes.13) = make.names(genes.13$Gene.Symbol, unique = T)
rownames(genes.72) = make.names(genes.72$Gene.Symbol, unique = T)

data_list = list(
  EGSE22 = genes.22,
  EGSE53 = genes.53)

duplicated(genes.11$gene_symbol)
data <- genes.11[!duplicated(genes.11$Gene.Symbol),]

data_list <- list(
  BGSE11 = na.omit(genes.11),
  BGSE12 = na.omit(genes.12),
  BGSE13 = na.omit(genes.13),
  BGSE72 = na.omit(genes.72))

head(data_list[[1]])
class(data_list)

diffexplist <- data_list

?rem_mv
meta_degs_rem <- rem_mv(diffexp = data_list,
                        pcriteria = "pvalue",
                        foldchangecol = 'logFC', 
                        genenamecol = "gene_symbol",
                        collaps = TRUE,
                        llcol = "lower_CI",
                        rlcol = "upper_CI",
                        vcol = NULL,
                        cvar = TRUE,
                        metathr = 0.1,
                        jobname = "MetaVolcano",
                        outputfolder = ".", 
                        draw = 'HTML',
                        ncores = 1)

head(meta_degs_rem@metaresult, 3)
meta_degs_rem@MetaVolcano

duplicated(genes.11$ID)

z <- list(
  "EGSE22" = sample(genes.22$gene_symbol),
  "EGSE53" = sample(genes.53$gene_symbol))

dvenn = ggvenn::ggvenn(z, fill_color = c("#4682B4", "#FFD700", "#CD534CFF", "#A9A9A9"),
                       show_elements = F,
                       stroke_size = 0.5, set_name_size = 4)
print(dvenn)

meta_degs_vote <- votecount_mv(diffexp = data_list,
                               pcriteria='pvalue',
                               foldchangecol='logFC',
                               genenamecol='gene_symbol',
                               geneidcol=NULL,
                               pvalue=0.05,
                               foldchange=1, 
                               metathr=0.1,
                               collaps=T,
                               jobname="MetaVolcano", 
                               outputfolder=".",
                               draw='HTML')

head(meta_degs_vote@metaresult, 3)
meta_degs_vote@degfreq
meta_degs_vote@MetaVolcano

meta_degs_comb <- combining_mv(diffexp=data_list,
                               pcriteria='pvalue', 
                               foldchangecol='logFC',
                               genenamecol='gene_symbol',
                               geneidcol=NULL,
                               metafc='Mean',
                               metathr=0.1, 
                               collaps=TRUE,
                               jobname="MetaVolcano",
                               outputfolder=".",
                               draw='HTML')

meta_degs_comb@MetaVolcano
metavolcanolist <- as.data.frame(meta_degs_comb@metaresult)
head(meta_degs_comb@metaresult, 3)
meta_degs_comb@MetaVolcano

#Lista de MetaGenes Gerais
meta_DEGs <- as.data.frame(meta_degs_comb@metaresult)
#Verificando Duplicatas
duplicated(meta_DEGs$gene_symbol)
meta_DEGs <- meta_DEGs[!duplicated(meta_DEGs$gene_symbol),]

#Verificando significantes, up e down
meta_DEGsUp <- meta_DEGs[meta_DEGs$metafc > 1 & meta_DEGs$metap < 0.05, ]
meta_DEGsDown <- meta_DEGs[meta_DEGs$metafc < -1 & meta_DEGs$metap < 0.05, ]
meta_DEGsSig <- meta_DEGs[meta_DEGs$metap < 0.05, ]

#Salvando como Excel e csv
write.xlsx(meta_DEGs, "meta_DEGs.xlsx")
write.xlsx(meta_DEGsUp, "meta_DEGsUp.xlsx")
write.xlsx(meta_DEGsDown, "meta_DEGsDown.xlsx")
write.xlsx(meta_DEGsSig, "meta_DEGsSig.xlsx")

# Salvando gráficos com o tamanho das fontes aumentado
#Venn
ggsave(filename = "dvenn.png",  # Nome do arquivo de saída
       plot = dvenn,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)

#SummaryFoldChange
ggsave(filename = "SummaryFoldChange.png",  # Nome do arquivo de saída
       plot = meta_degs_rem@MetaVolcano,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)

meta_degs_vote@metaresult
#Volcano
ggsave(filename = "Volcano.png",  # Nome do arquivo de saída
       plot = meta_degs_comb@MetaVolcano,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)
#MetaResults
ggsave(filename = "MetaResults.png",  # Nome do arquivo de saída
       plot = meta_degs_vote@MetaVolcano,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)

# PBMC MetaAnalysis -------------------------------------------------------

getwd()
setwd("~adrielnobile/Documents/DataSets/MetaVolcano - ArrayxBulk/PBMC_InVivo/")
setwd("~adrielnobile/Documents/DataSets/MetaVolcano - ArrayxBulk/PBMC_InVivo/MetaVolcano/")

#B
BGSE168 <- read_excel("ABGSE168048.xlsx")
BGSE164 <- read_excel("BBGSE164266.xlsx")

#C
CGSE65 <- read_excel("ACGSE65123.xlsx")
CGSE93 <- read_excel("ACGSE93711.xlsx")
CGSE21 <- read_excel("BCGSE212871.xlsx")

#E
EGSE36 <- read_excel("AEGSE36539.xlsx")

#C
genes.21 = as.data.frame(na.omit(CGSE21))
genes.65 = as.data.frame(na.omit(CGSE65))
genes.93 = as.data.frame(na.omit(CGSE93))

#=======================AnnotationDbi ------
# Vetor de IDs de gene para os quais você deseja obter os nomes externos
gene_ids <- c(genes.10$ENTREZ_GENE_ID)

# Carregar a base de dados de anotação para o organismo de interesse
org <- org.Hs.eg.db
db <- org.Hs.eg.db

keytypes(org.Hs.eg.db)# Realizar a transformação dos IDs de gene para nomes externos
gene_names.13 <- mapIds(db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBLTRANS")

gene_names.13 <- mapIds(db, keys = gene_ids, column = "ENSEMBL", keytype = "ENTREZID")

#Renomar para identificadores serem iguais: pvalue, log2FoldChange e genesymbol
print(colnames(genes.93))
colnames(genes.93)[colnames(genes.93) == "P.Value"] <- "pvalue"
colnames(genes.11)[colnames(genes.11) == "logFC"] <- "log2FoldChange"
colnames(genes.65)[colnames(genes.65) == "GENE_SYMBOL"] <- "Gene.Symbol"
colnames(genes.15)[colnames(genes.15) == "ID"] <- "Gene.Symbol"
colnames(genes.21)[colnames(genes.21) == "log2FoldChange"] <- "logFC"

#Duplicados
duplicated(genes.168$Gene.Symbol)
genes.168 <- genes.168[!duplicated(genes.168$Gene.Symbol),]

#=======================Confidence Interval===============================
# Set the desired confidence level (e.g., 95%)
confidence_level <- 0.95

#genes.21
# Calculate the degrees of freedom (typically n - 1 for t-distribution)
df.21 <- nrow(genes.21) - 1
# Calculate the margin of error
margin_of_error.21 <- qt((1 + confidence_level) / 2, df.21) * genes.21$logFC / sqrt(nrow(genes.21))
# Calculate the lower and upper confidence limits
genes.21$lower_CI <- genes.21$logFC - margin_of_error.21
genes.21$upper_CI <- genes.21$logFC + margin_of_error.21

#genes.65
# Calculate the degrees of freedom (typically n - 1 for t-distribution)
df.65 <- nrow(genes.65) - 1
# Calculate the margin of error
margin_of_error.65 <- qt((1 + confidence_level) / 2, df.65) * genes.65$logFC / sqrt(nrow(genes.65))
# Calculate the lower and upper confidence limits
genes.65$lower_CI <- genes.65$logFC - margin_of_error.65
genes.65$upper_CI <- genes.65$logFC + margin_of_error.65

#genes.93
# Calculate the degrees of freedom (typically n - 1 for t-distribution)
df.93 <- nrow(genes.93) - 1
# Calculate the margin of error
margin_of_error.93 <- qt((1 + confidence_level) / 2, df.93) * genes.93$logFC / sqrt(nrow(genes.93))
# Calculate the lower and upper confidence limits
genes.93$lower_CI <- genes.93$logFC - margin_of_error.93
genes.93$upper_CI <- genes.93$logFC + margin_of_error.93

#========================MetaVolc===================================
#HBV
data_list = list(
  BGSE164 = genes.164,
  BGSE168 = genes.168)

#HCV
data_list = list(
  CGSE21 = genes.21,
  CGSE65 = genes.65,
  CGSE93 = genes.93)

#Microarray
data_list = list(
  CGSE93 = genes.93,
  BGSE168 = genes.168,
  CGSE65 = genes.65)

#Bulk
data_list = list(
  BGSE164 = genes.164,
  CGSE21 = genes.21)

#MicroarrayxBulk
data_list = list(
  CGSE93 = genes.93,
  BGSE168 = genes.168,
  CGSE65 = genes.65,
  BGSE164 = genes.164,
  CGSE21 = genes.21)

meta_degs_rem <- rem_mv(diffexp = data_list,
                        pcriteria = "pvalue",
                        foldchangecol = 'logFC', 
                        genenamecol = "Gene.Symbol",
                        collaps = TRUE,
                        llcol = "lower_CI",
                        rlcol = "upper_CI",
                        vcol = NULL,
                        cvar = TRUE,
                        metathr = 0.1,
                        jobname = "MetaVolcano",
                        outputfolder = ".", 
                        draw = 'HTML',
                        ncores = 1)
meta_degs_rem@MetaVolcano

#HBV
z <- list(
  "BGSE164" = sample(genes.164$Gene.Symbol),
  "BGSE164" = sample(genes.168$Gene.Symbol))

#HCV
z <- list(
  "CGSE21" = sample(genes.21$Gene.Symbol),
  "CGSE65" = sample(genes.65$Gene.Symbol),
  "CGSE93" = sample(genes.93$Gene.Symbol))

#Microarray
z <- list(
  "BGSE168" = sample(genes.168$Gene.Symbol),
  "CGSE65" = sample(genes.65$Gene.Symbol),
  "CGSE93" = sample(genes.93$Gene.Symbol))

#Bulk
z <- list(
  "BGSE164" = sample(genes.164$Gene.Symbol),
  "CGSE21" = sample(genes.21$Gene.Symbol))

#MicroarrayxBulk
z <- list(
  "BGSE168" = sample(genes.168$Gene.Symbol),
  "CGSE65" = sample(genes.65$Gene.Symbol),
  "CGSE93" = sample(genes.93$Gene.Symbol),
  "BGSE164" = sample(genes.164$Gene.Symbol),
  "CGSE21" = sample(genes.21$Gene.Symbol))

#UpsetPlot
upset_plot <- UpSetR::upset(fromList(z), 
                            order.by = "freq", 
                            sets.bar.color = "#4682B4",
                            main.bar.color = "#A9A9A9",
                            matrix.color = "#FFD700",
                            text.scale = 1.0,
                            nsets = 10)
print(upset_plot)

#Venn
dvenn = ggvenn::ggvenn(z, fill_color = c("#4682B4", "#FFD700", "#CD534CFF", "#A9A9A9"),
                       show_elements = F,
                       stroke_size = 0.5, set_name_size = 4)
print(dvenn)

#
meta_degs_vote <- votecount_mv(diffexp = data_list,
                               pcriteria='pvalue',
                               foldchangecol='logFC',
                               genenamecol='Gene.Symbol',
                               geneidcol=NULL,
                               pvalue=0.05,
                               foldchange=0.5, 
                               metathr=0.1,
                               collaps=T,
                               jobname="MetaVolcano", 
                               outputfolder=".",
                               draw='HTML')
meta_degs_vote@degfreq
meta_degs_vote@MetaVolcano

#
meta_degs_comb <- combining_mv(diffexp=data_list,
                               pcriteria='pvalue', 
                               foldchangecol='logFC',
                               genenamecol='Gene.Symbol',
                               geneidcol=NULL,
                               metafc='Mean',
                               metathr=0.1, 
                               collaps=TRUE,
                               jobname="MetaVolcano",
                               outputfolder=".",
                               draw='HTML')

meta_degs_comb@MetaVolcano

#Lista de MetaGenes Gerais
meta_DEGs <- as.data.frame(meta_degs_comb@metaresult)
#Verificando Duplicatas
duplicated(meta_DEGs$Gene.Symbol)
meta_DEGs <- meta_DEGs[!duplicated(meta_DEGs$Gene.Symbol),]
print(colnames(meta_DEGs))
#Verificando significantes, up e down
meta_DEGsUp = subset(meta_DEGs, metafc > 1 & metap < 0.05)
meta_DEGsDown = subset(meta_DEGs, metafc < -1 & metap < 0.05)
meta_DEGsSig = subset(meta_DEGs, metap < 0.05)

meta_DEGsUp = subset(meta_DEGs, metafc > 0.5 & metap < 0.05)
meta_DEGsDown = subset(meta_DEGs, metafc < -0.5 & metap < 0.05)
meta_DEGsSig = subset(meta_DEGs, metap < 0.05)

#Salvando como Excel e csv
write.xlsx(meta_DEGs, "meta_DEGs.xlsx")
write.xlsx(meta_DEGsUp, "meta_DEGsUp.xlsx")
write.xlsx(meta_DEGsDown, "meta_DEGsDown.xlsx")
write.xlsx(meta_DEGsSig, "meta_DEGsSig.xlsx")

# Salvando gráficos com o tamanho das fontes aumentado
#Venn
ggsave(filename = "dvenn.png",  # Nome do arquivo de saída
       plot = dvenn,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)

#SummaryFoldChange
ggsave(filename = "SummaryFoldChange.png",  # Nome do arquivo de saída
       plot = meta_degs_rem@MetaVolcano,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)

meta_degs_vote@metaresult
#Volcano
ggsave(filename = "Volcano.png",  # Nome do arquivo de saída
       plot = meta_degs_comb@MetaVolcano,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)
#MetaResults
ggsave(filename = "MetaResults.png",  # Nome do arquivo de saída
       plot = meta_degs_vote@MetaVolcano,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)

#Frequencia de Degs
ggsave(filename = "DegFreq.png",  # Nome do arquivo de saída
       plot = meta_degs_vote@degfreq,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)


# Liver Metanalysis -------------------------------------------------------

setwd("~adrielnobile/Documents/DataSets/MetaVolcano - ArrayxBulk/Liver_InVivo/")

#B
BGSE10B <- read_excel("ABGSE107170.xlsx")
BGSE14 <- read_excel("ABGSE14668.xlsx")
BGSE38 <- read_excel("ABGSE38941.xlsx")
BGSE47 <- read_excel("ABGSE47197.xlsx")
BGSE55 <- read_excel("ABGSE55092.xlsx")
BGSE65 <- read_excel("ABGSE65359.xlsx")
BGSE94 <- read_excel("BBGSE94660.xlsx")
BGSE16 = as.data.frame(read_excel("DEGs Log2Modificado.xlsx"))

#C
CGSE10C <- read_excel("ACGSE107170.xlsx")
CGSE78 <- read_excel("ACGSE78737.xlsx")
CGSE15 <- read_excel("BCGSE154211.xlsx")

#D
DGSE10D <- read_excel("MetaVolcano/ADGSE107170.xlsx")
DGSE98 = as.data.frame(read_tsv("GSE98383.tsv"))

genes.10D = as.data.frame(na.omit(DGSE10D))
genes.98 = as.data.frame(na.omit(DGSE98))
#B
genes.10B = as.data.frame(na.omit(BGSE10B))
genes.14 = as.data.frame(na.omit(BGSE14))
genes.38 = as.data.frame(na.omit(BGSE38))
genes.47 = as.data.frame(na.omit(BGSE47))
genes.55 = as.data.frame(na.omit(BGSE55))
genes.65 = as.data.frame(na.omit(BGSE65))
genes.94 = as.data.frame(na.omit(BGSE94))
genes.16 = as.data.frame(na.omit(BGSE16))

#=======================AnnotationDbi ------

# Vetor de IDs de gene para os quais você deseja obter os nomes externos
gene_ids <- c(genes.10$ENTREZ_GENE_ID)

# Carregar a base de dados de anotação para o organismo de interesse
org <- org.Hs.eg.db
db <- org.Hs.eg.db

keytypes(org.Hs.eg.db)# Realizar a transformação dos IDs de gene para nomes externos
gene_names.13 <- mapIds(db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBLTRANS")

gene_names.13 <- mapIds(db, keys = gene_ids, column = "ENSEMBL", keytype = "ENTREZID")

#Renomar para identificadores serem iguais: pvalue, log2FoldChange e genesymbol
print(colnames(BGSE16))
View(BGSE16)
colnames(DGSE98)[colnames(DGSE98) == "P.Value"] <- "pvalue"
colnames(genes.11)[colnames(genes.11) == "logFC"] <- "log2FoldChange"
colnames(BGSE16)[colnames(BGSE16) == "genes_symbol"] <- "Gene.Symbol"
colnames(DGSE98)[colnames(DGSE98) == "ID"] <- "Gene.Symbol"
print(colnames(genes.11))
colnames(BGSE16)[colnames(BGSE16) == "log2FoldChange"] <- "logFC"
colnames(genes.15)[colnames(genes.15) == "padj"] <- "adj.P.Val"

#Duplicados
duplicated(DGSE98$ID)
DGSE98 <- DGSE98[!duplicated(DGSE98$ID),]

#=======================Confidence Interval===============================
# Set the desired confidence level (e.g., 95%)
confidence_level <- 0.95

#Microarray
#genes.10B
# Calculate the degrees of freedom (typically n - 1 for t-distribution)
df.10B <- nrow(genes.10B) - 1
# Calculate the margin of error
margin_of_error.10B <- qt((1 + confidence_level) / 2, df.10B) * genes.10B$logFC / sqrt(nrow(genes.10B))
# Calculate the lower and upper confidence limits
genes.10B$lower_CI <- genes.10B$logFC - margin_of_error.10B
genes.10B$upper_CI <- genes.10B$logFC + margin_of_error.10B

#genes.14
# Calculate the degrees of freedom (typically n - 1 for t-distribution)
df.14 <- nrow(genes.14) - 1
# Calculate the margin of error
margin_of_error.14 <- qt((1 + confidence_level) / 2, df.14) * genes.14$logFC / sqrt(nrow(genes.14))
# Calculate the lower and upper confidence limits
genes.14$lower_CI <- genes.14$logFC - margin_of_error.14
genes.14$upper_CI <- genes.14$logFC + margin_of_error.14

#genes.38
df.38 <- nrow(genes.38) - 1
# Calculate the margin of error
margin_of_error.38 <- qt((1 + confidence_level) / 2, df.38) * genes.38$logFC / sqrt(nrow(genes.38))
# Calculate the lower and upper confidence limits
genes.38$lower_CI <- genes.38$logFC - margin_of_error.38
genes.38$upper_CI <- genes.38$logFC + margin_of_error.38

#genes.47
df.47 <- nrow(genes.47) - 1
# Calculate the margin of error
margin_of_error.47 <- qt((1 + confidence_level) / 2, df.47) * genes.47$logFC / sqrt(nrow(genes.47))
# Calculate the lower and upper confidence limits
genes.47$lower_CI <- genes.47$logFC - margin_of_error.47
genes.47$upper_CI <- genes.47$logFC + margin_of_error.47

#genes.55
df.55 <- nrow(genes.55) - 1
# Calculate the margin of error
margin_of_error.55 <- qt((1 + confidence_level) / 2, df.55) * genes.55$logFC / sqrt(nrow(genes.55))
# Calculate the lower and upper confidence limits
genes.55$lower_CI <- genes.55$logFC - margin_of_error.55
genes.55$upper_CI <- genes.55$logFC + margin_of_error.55

#genes.65
df.65 <- nrow(genes.65) - 1
# Calculate the margin of error
margin_of_error.65 <- qt((1 + confidence_level) / 2, df.65) * genes.65$logFC / sqrt(nrow(genes.65))
# Calculate the lower and upper confidence limits
genes.65$lower_CI <- genes.65$logFC - margin_of_error.65
genes.65$upper_CI <- genes.65$logFC + margin_of_error.65

#genes.78
df.78 <- nrow(genes.78) - 1
# Calculate the margin of error
margin_of_error.78 <- qt((1 + confidence_level) / 2, df.78) * genes.78$logFC / sqrt(nrow(genes.78))
# Calculate the lower and upper confidence limits
genes.78$lower_CI <- genes.78$logFC - margin_of_error.78
genes.78$upper_CI <- genes.78$logFC + margin_of_error.78

#genes.65
df.10C <- nrow(genes.10C) - 1
# Calculate the margin of error
margin_of_error.10C <- qt((1 + confidence_level) / 2, df.10C) * genes.10C$logFC / sqrt(nrow(genes.10C))
# Calculate the lower and upper confidence limits
genes.10C$lower_CI <- genes.10C$logFC - margin_of_error.10C
genes.10C$upper_CI <- genes.10C$logFC + margin_of_error.10C

#genes.10D
df.10D <- nrow(genes.10D) - 1
# Calculate the margin of error
margin_of_error.10D <- qt((1 + confidence_level) / 2, df.10D) * genes.10D$logFC / sqrt(nrow(genes.10D))
# Calculate the lower and upper confidence limits
genes.10D$lower_CI <- genes.10D$logFC - margin_of_error.10D
genes.10D$upper_CI <- genes.10D$logFC + margin_of_error.10D

#genes.98
df.98 <- nrow(genes.98) - 1
# Calculate the margin of error
margin_of_error.98 <- qt((1 + confidence_level) / 2, df.98) * genes.98$logFC / sqrt(nrow(genes.98))
# Calculate the lower and upper confidence limits
genes.98$lower_CI <- genes.98$logFC - margin_of_error.98
genes.98$upper_CI <- genes.98$logFC + margin_of_error.98

#Bulk
#genes.94
df.94 <- nrow(genes.94) - 1
# Calculate the margin of error
margin_of_error.94 <- qt((1 + confidence_level) / 2, df.94) * genes.94$logFC / sqrt(nrow(genes.94))
# Calculate the lower and upper confidence limits
genes.94$lower_CI <- genes.94$logFC - margin_of_error.94
genes.94$upper_CI <- genes.94$logFC + margin_of_error.94

#genes.15
df.15 <- nrow(genes.15) - 1
# Calculate the margin of error
margin_of_error.15 <- qt((1 + confidence_level) / 2, df.15) * genes.15$logFC / sqrt(nrow(genes.15))
# Calculate the lower and upper confidence limits
genes.15$lower_CI <- genes.15$logFC - margin_of_error.15
genes.15$upper_CI <- genes.15$logFC + margin_of_error.15

#genes.16
# Calculate the degrees of freedom (typically n - 1 for t-distribution)
df.16 <- nrow(genes.16) - 1
# Calculate the margin of error
margin_of_error.16 <- qt((1 + confidence_level) / 2, df.16) * genes.16$logFC / sqrt(nrow(genes.16))
# Calculate the lower and upper confidence limits
genes.16$lower_CI <- genes.16$logFC - margin_of_error.16
genes.16$upper_CI <- genes.16$logFC + margin_of_error.16

#========================MetaVolc===================================
data <- subset(data, !grepl(" ", data) & !grepl(" ", as.character()))
data = merged_sig
data = na.omit(data)
backup = data

#HBV
data_list = list(
  BGSE10B = genes.10B,
  BGSE14 = genes.14,
  BGSE38 = genes.38,
  BGSE47 = genes.47,
  BGSE55 = genes.55,
  BGSE65 = genes.65,
  BGSE94 = genes.94)

#HCV
data_list = list(
  CGSE10C = genes.10C,
  CGSE15 = genes.15,
  CGSE78 = genes.78)

#HDV
data_list = list(
  DGSE10D = genes.10D,
  DGSE98 = genes.98)

#Microarray
data_list = list(
  BGSE10B = genes.10B,
  BGSE14 = genes.14,
  BGSE38 = genes.38,
  BGSE47 = genes.47,
  BGSE55 = genes.55,
  BGSE65 = genes.65,
  CGSE10C = genes.10C,
  CGSE78 = genes.78,
  DGSE10D = genes.10D)

#Bulk
data_list = list(
  BGSE94 = genes.94,
  CGSE15 = genes.15)

#MicroarrayXBulk
data_list = list(
  BGSE10B = genes.10B,
  BGSE14 = genes.14,
  BGSE38 = genes.38,
  BGSE47 = genes.47,
  BGSE55 = genes.55,
  BGSE65 = genes.65,
  CGSE10C = genes.10C,
  CGSE78 = genes.78,
  DGSE10D = genes.10D,
  BGSE94 = genes.94,
  CGSE15 = genes.15)

#Oncogenic
data_list = list(
  BGSE10B = genes.10B,
  BGSE14 = genes.14,
  BGSE38 = genes.38,
  BGSE47 = genes.47,
  BGSE55 = genes.55,
  BGSE65 = genes.65,
  CGSE10C = genes.10C,
  CGSE78 = genes.78,
  DGSE10D = genes.10D,
  BGSE94 = genes.94,
  BGSE16 = genes.16,
  CGSE15 = genes.15)

print(colnames(genes.10B))
meta_degs_rem <- rem_mv(diffexp = data_list,
                        pcriteria = "pvalue",
                        foldchangecol = 'logFC', 
                        genenamecol = "Gene.Symbol",
                        collaps = TRUE,
                        llcol = "lower_CI",
                        rlcol = "upper_CI",
                        vcol = NULL,
                        cvar = TRUE,
                        metathr = 0.1,
                        jobname = "MetaVolcano",
                        outputfolder = ".", 
                        draw = 'HTML',
                        ncores = 1)
meta_degs_rem@MetaVolcano

#MicroarrayXBulk
#HBV
z <- list(
  "BGSE10B" = sample(genes.10B$Gene.Symbol),
  "BGSE14" = sample(genes.14$Gene.Symbol),
  "BGSE38" = sample(genes.38$Gene.Symbol),
  "BGSE47" = sample(genes.47$Gene.Symbol),
  "BGSE55" = sample(genes.55$Gene.Symbol),
  "BGSE65" = sample(genes.65$Gene.Symbol),
  "BGSE94" = sample(genes.94$Gene.Symbol))

z <- list(
  "BGSE10B" = sample(genes.10B$Gene.Symbol),
  "BGSE14" = sample(genes.14$Gene.Symbol),
  "BGSE38" = sample(genes.38$Gene.Symbol),
  "BGSE47" = sample(genes.47$Gene.Symbol),
  "BGSE55" = sample(genes.55$Gene.Symbol),
  "BGSE65" = sample(genes.65$Gene.Symbol),
  "CGSE10C" = sample(genes.10C$Gene.Symbol),
  "CGSE78" = sample(genes.78$Gene.Symbol),
  "DGSE10D" = sample(genes.10D$Gene.Symbol),
  "BGSE94" = sample(genes.94$Gene.Symbol),
  "CGSE15" = sample(genes.15$Gene.Symbol))

#HCV
z <- list(
  "CGSE10C" = sample(genes.10C$Gene.Symbol),
  "CGSE78" = sample(genes.78$Gene.Symbol),
  "CGSE15" = sample(genes.15$Gene.Symbol))

z = list(
  "DGSE10D" = sample(genes.10D$Gene.Symbol),
  "DGSE98" = sample(genes.98$Gene.Symbol))

#Microarray
z <- list(
  "BGSE10B" = sample(genes.10B$Gene.Symbol),
  "BGSE14" = sample(genes.14$Gene.Symbol),
  "BGSE38" = sample(genes.38$Gene.Symbol),
  "BGSE47" = sample(genes.47$Gene.Symbol),
  "BGSE55" = sample(genes.55$Gene.Symbol),
  "BGSE65" = sample(genes.65$Gene.Symbol),
  "CGSE10C" = sample(genes.10C$Gene.Symbol),
  "CGSE78" = sample(genes.78$Gene.Symbol),
  "DGSE10D" = sample(genes.10D$Gene.Symbol))

#Oncogenic
z <- list(
  "BGSE10B" = sample(genes.10B$Gene.Symbol),
  "BGSE14" = sample(genes.14$Gene.Symbol),
  "BGSE38" = sample(genes.38$Gene.Symbol),
  "BGSE47" = sample(genes.47$Gene.Symbol),
  "BGSE55" = sample(genes.55$Gene.Symbol),
  "BGSE65" = sample(genes.65$Gene.Symbol),
  "CGSE10C" = sample(genes.10C$Gene.Symbol),
  "CGSE78" = sample(genes.78$Gene.Symbol),
  "DGSE10D" = sample(genes.10D$Gene.Symbol),
  "DGSE98" = sample(genes.98$Gene.Symbol),
  "BGSE94" = sample(genes.94$Gene.Symbol),
  "BGSE16" = sample(genes.16$Gene.Symbol),
  "CGSE15" = sample(genes.15$Gene.Symbol))

#Bulk
z <- list(
  "BGSE94" = sample(genes.94$Gene.Symbol),
  "CGSE15" = sample(genes.15$Gene.Symbol))

#UpsetPlot
upset_plot <- UpSetR::upset(fromList(z), 
                            order.by = "freq", 
                            sets.bar.color = "#4682B4",
                            main.bar.color = "#A9A9A9",
                            matrix.color = "#FFD700",
                            text.scale = 1.0,
                            nsets = 14)
print(upset_plot)

#Venn
dvenn = ggvenn::ggvenn(z, fill_color = c("#4682B4", "#FFD700", "#CD534CFF", "#A9A9A9"),
                       show_elements = F,
                       stroke_size = 0.5, set_name_size = 4)
print(dvenn)

#
meta_degs_vote <- votecount_mv(diffexp = data_list,
                               pcriteria='pvalue',
                               foldchangecol='logFC',
                               genenamecol='Gene.Symbol',
                               geneidcol=NULL,
                               pvalue=0.05,
                               foldchange=1, 
                               metathr=0.1,
                               collaps=T,
                               jobname="MetaVolcano", 
                               outputfolder=".",
                               draw='HTML')
meta_degs_vote@degfreq
meta_degs_vote@MetaVolcano

#
meta_degs_comb <- combining_mv(diffexp=data_list,
                               pcriteria='pvalue', 
                               foldchangecol='logFC',
                               genenamecol='Gene.Symbol',
                               geneidcol=NULL,
                               metafc='Mean',
                               metathr=0.1, 
                               collaps=TRUE,
                               jobname="MetaVolcano",
                               outputfolder=".",
                               draw='HTML')

meta_degs_comb@MetaVolcano

#Lista de MetaGenes Gerais
meta_DEGs <- as.data.frame(meta_degs_comb@metaresult)
#Verificando Duplicatas
duplicated(meta_DEGs$Gene.Symbol)
meta_DEGs <- meta_DEGs[!duplicated(meta_DEGs$Gene.Symbol),]

#Verificando significantes, up e down
meta_DEGsUp = subset(meta_DEGs, metafc > 1 & metap < 0.05)
meta_DEGsDown = subset(meta_DEGs, metafc < -1 & metap < 0.05)
meta_DEGsSig = subset(meta_DEGs, metap < 0.05)
#Salvando como Excel e csv
write.xlsx(meta_DEGs, "meta_DEGs.xlsx")
write.xlsx(meta_DEGsUp, "meta_DEGsUp.xlsx")
write.xlsx(meta_DEGsDown, "meta_DEGsDown.xlsx")
write.xlsx(meta_DEGsSig, "meta_DEGsSig.xlsx")

# Salvando gráficos com o tamanho das fontes aumentado
#Venn
ggsave(filename = "dvenn.png",  # Nome do arquivo de saída
       plot = dvenn,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)

#SummaryFoldChange
ggsave(filename = "SummaryFoldChange.png",  # Nome do arquivo de saída
       plot = meta_degs_rem@MetaVolcano,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)

meta_degs_vote@metaresult
#Volcano
ggsave(filename = "Volcano.png",  # Nome do arquivo de saída
       plot = meta_degs_comb@MetaVolcano,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)
#MetaResults
ggsave(filename = "MetaResults.png",  # Nome do arquivo de saída
       plot = meta_degs_vote@MetaVolcano,        # Objeto do gráfico
       width = 8,                            # Largura da imagem
       height = 6,                           # Altura da imagem
       units = "in")                         # Unidades da imagem (polegadas)
