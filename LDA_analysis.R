#Linear Discriminant Analysis (LDA) - By Adriel Leal Nóbile
install.packages('devtools')
devtools::install_github('fawda123/ggord')
{library(MASS)
library(ggplot2)
library(ggord)
library(dplyr)
library(GEOquery)
library(tidyr)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(klaR)
library(psych)
library(MASS)
library(ggord)
library(devtools)
library(GEOquery)
library(curl)}

#--------------------------  Executando LDA ---------------------------------
#====================.  Passo 1 - Carregando os dados ======================
merge_HBV_Liver = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/merge_HBV_Liver.csv")
merge_HCV_Liver = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/merge_HCV_Liver.csv")
merge_HDV_Liver = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/merge_HDV_Liver.csv")
merge_HBV_PBMC = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/merge_HBV_PBMC.csv")
merge_HCV_PBMC = read.csv("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/merge_HCV_PBMC.csv")

## Manipulando data.frame
#HBV Liver
df <- merge_HBV_Liver
rownames(df) = df$ENSEMBL #Definindo variaveis|genes como rownames
#GeneInfo - Separando informacoes de grupos dos genes
geneInfo <- subset(df, select = sytem)
geneInfo$genes <- rownames(geneInfo)
df = df[4:55] #Separando counts
df <- as.data.frame(lapply(df,
                           as.numeric)) #Certificar que os counts sao numericos
rownames(df) = rownames(geneInfo) #Genes como rownames
df <- as.data.frame(t(df))
df
df$Groups = c(rep(c('infected', 'control'),
                  length.out = 42),
              rep("control", 5), rep("infected", 5))

write.xlsx(df, "Suppl.Table3a_HBVL_LDACounts.xlsx", rowNames = TRUE)

#HBV_PBMC
df <- merge_HBV_PBMC
rownames(df) = df$geneID #Definindo variaveis|genes como rownames
#GeneInfo - Separando informacoes de grupos dos genes
geneInfo <- subset(df, select = sytem)
geneInfo$genes <- rownames(geneInfo)
df = df[3:97] #Separando counts
df <- as.data.frame(lapply(df,
                           as.numeric)) #Certificar que os counts sao numericos
rownames(df) = rownames(geneInfo) #Genes como rownames
df <- as.data.frame(t(df))
df
gse_code = "GSE173897" #Substituir pelo seu GSE
geo_url = paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
                gse_code,
                sep="")
seriesmatrix = as.data.frame(getGEO(gse_code,
                                    destdir = ".",
                                    getGPL = FALSE))

condition = as.data.frame((factor((seriesmatrix[[48]]),
                                  levels = c("Active", #Inserir Grupo "Infectado"
                                             "Inactive")))) #Inserir Grupo "Controle"
groups = condition
groups = ifelse(groups == "Inactive", "control",
                "infected")
colnames(groups) = c("Groups")
df$Groups = groups[,c(1)]

df <- as.data.frame(t(df))

write.xlsx(df, "Suppl.Table3a_HBVP_LDACounts.xlsx", rowNames = TRUE)
#HCV_Liver
df <- merge_HCV_Liver
rownames(df) = df$geneID #Definindo variaveis|genes como rownames
#GeneInfo - Separando informacoes de grupos dos genes
geneInfo <- subset(df, select = sytem)
geneInfo$genes <- rownames(geneInfo)
df = df[3:12] #Separando counts
df <- as.data.frame(lapply(df,
                           as.numeric)) #Certificar que os counts sao numericos
rownames(df) = rownames(geneInfo) #Genes como rownames
df <- as.data.frame(t(df))
df
df$Groups = c(rep(c('control', 'infected'),
                  length.out = 10))
# df <- as.data.frame(t(df))
write.xlsx(df, "Suppl.Table3a_HCVL_LDACounts.xlsx", rowNames = TRUE)

#HCV_PBMC
df = merge_HCV_PBMC
rownames(df) = df$geneID #Definindo variaveis|genes como rownames
#GeneInfo - Separando informacoes de grupos dos genes
geneInfo <- subset(df, select = sytem)
geneInfo$genes <- rownames(geneInfo)
df = df[3:130] #Separando counts
df <- as.data.frame(lapply(df,
                           as.numeric)) #Certificar que os counts sao numericos
rownames(df) = rownames(geneInfo) #Genes como rownames
df <- as.data.frame(t(df))
df$Groups = c(rep("infected", 36),
              rep("control", 70),
              rep("infected", 1),
              rep("control", 1),
              rep("control", 1),
              rep("infected", 1),
              rep("infected", 1),
              rep("control", 1),
              rep("control", 1),
              rep("infected", 1),
              rep("infected", 1),
              rep("control", 1),
              rep("infected", 5),
              rep("control", 5),
              rep("infected", 2))

write.xlsx(df, "Suppl.Table3a_HCVP_LDACounts.xlsx", rowNames = TRUE)

#HDV_Liver
df = merge_HDV_Liver
df = df[!duplicated(df$external_gene_name),] #Remover duplicata
rownames(df) = df$external_gene_name #Definindo variaveis|genes como rownames
#GeneInfo - Separando informacoes de grupos dos genes
geneInfo <- subset(df, select = sytem)
geneInfo$genes <- rownames(geneInfo)
df = df[4:77] #Separando counts
df <- as.data.frame(lapply(df,
                           as.numeric)) #Certificar que os counts sao numericos
rownames(df) = rownames(geneInfo) #Genes como rownames
df <- as.data.frame(t(df))
df$Groups = c(rep("infected", 12),
              rep("control", 23),
              rep("infected", 29),
              rep("infected", 5),
              rep("control", 5))
write.xlsx(df, "Suppl.Table3a_HDVL_LDACounts.xlsx", rowNames = TRUE)

#====================.  Passo 2 - Escalonar os dados ======================
print(c(nrow(df),  #Verificando numero de linhas
        ncol(df))) #Verificando numero de colunas
num_colunas <- ncol(df) #Salvando numero de colunas
df[, 1:(num_colunas - 1)] <- scale(df[, 1:(num_colunas - 1)]) # Escala todas as colunas, exceto a última
apply(df[, 1:(num_colunas - 1)], 2, mean) #encontre a média de cada variável preditora
apply(df[, 1:(num_colunas - 1)], 2, sd) #encontre o desvio padrão de cada variável preditora

#Se necessario, filtrar os dados
#Alternativo
# backup_df = df
# df = backup_df
# col_zero <- colSums(df < 1)
# row_zero <- rowSums(df < 1)
# excl_row <- which(row_zero > 500)
# excl_col <- which(col_zero >= 50)
# 
# table(col_zero)
# table(row_zero)
# table(excl_row)
# table(excl_col)
# table(df$Groups)
# 
# df <- df[-excl_row, -excl_col]
# View(df)
#====================.  Passo 3 - criar modelos de treinamento e teste =====
set.seed(2)

sample <- sample(c(TRUE,
                   FALSE),
                 nrow(df),
                 replace=TRUE,
                 prob=c(0.7, #70% do conjunto de dados como conjunto treino
                        0.3)) #30% restantes como conjunto de teste
train <- df[sample, ] #Separando conjunto de treino
test <- df[!sample, ] #Separando conjunto de teste

#====================.  Passo 4 - ajuste o modelo LDA  ====================
model <- lda(Groups~.,# fit modelo LDA
             data=train)
cor(df[1:61])
model

#====================.  Passo 5 - use o modelo para fazer previsões =======
predicted <- predict(model,
                     test) #use o modelo LDA para fazer previsões sobre dados de teste
names(predicted)
head(predicted$class) #classe prevista para as primeiras seis observações no conjunto de teste
head(predicted$x) #visualizar discriminantes lineares para as primeiras seis observações no conjunto de teste
mean(predicted$class==test$Groups) #encontrar precisão do modelo

#====================.  Salvando dados de Interesse ====================
#Plots do pacote
ldahist(data = predicted$x[,1], g = train$Groups)
dev.off()
partimat(Groups~., data = train, method = "lda")

# Obtenha as cargas para as variáveis (genes)
cargas_lda <- model$scaling
# Crie um dataframe para visualizar as cargas
cargas_df <- data.frame(Variable = colnames(train[, -1]),
                        LD1 = cargas_lda[, 1])
cargas_df <- cargas_df[order(abs(cargas_df$LD1),
                             decreasing = TRUE), ]
# Exiba as principais variáveis (genes) que contribuem para LD1
print(cargas_df)
# Adicione rótulos aos pontos para as 5 principais variáveis com maiores cargas absolutas
top_genes <- head(cargas_df[order(abs(cargas_df$LD1),
                                  decreasing = TRUE), ], 50)

#Salvando em Excel
accuracy <- mean(predicted$class == test$Groups)
df_accuracy <- data.frame( # Criar um data frame com a precisão
  LD1 = accuracy
)
#Separar o resultado da analise em data frames
#Class
class_info <- as.data.frame(predicted[["class"]])
colnames(class_info) <- c("Class")
#InfectedxControl
df_posterior <- as.data.frame(predicted$posterior)
#LD1
df_x <- as.data.frame(predicted$x) #LD1
#Combinando
df_combinado = cbind(class_info,
                     df_posterior,
                     df_x)
df_combinado$accuracy = c(accuracy)
#Merge geneInfo X cargas_df
cargas_df = merge(geneInfo, cargas_df,
                  by.x = "genes", by.y = "Variable")
cargas_df$hepatites = c("HBV PBMC")
#Merge geneInfo X Df_Combinado
geneInfo_50 = merge(geneInfo, top_genes,
                    by.x = "genes", by.y = "Variable")
geneInfo_50$hepatites = c("HBV PBMC")
#Excel
openxlsx::write.xlsx(df_combinado, "LDA_results.xlsx") #Resultados LDA
openxlsx::write.xlsx(geneInfo_50, "Top_50_Genes.xlsx") #Genes de Interesse
openxlsx::write.xlsx(cargas_df, "Informacoes totais de Genes.xlsx") #Genes e Systems

#====================.  Passo 6 - Vizualizando os dados ====================
#Histograma
ggplot(data.frame(LD1 = predicted$x[, 1],
                  Groups = test$Groups),
       aes(x = LD1,
           fill = Groups)) +
  geom_histogram(binwidth = 0.5,
                 position = "identity",
                 alpha = 0.7) +
  labs(title = "Histograma LDA", x = "LD1")
#Density Plot
#DataFrame Base
data_df <- data.frame(LD1 = predicted$x[, 1],
                      Groups = as.factor(predicted$class))
#Cor dos Plots
cores_manual <- c("control" = "#93BFCF",
                  "infected" = "#D71313")

#DensityPlot
densidade = ggplot(data_df, aes(x = LD1, fill = Groups)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = cores_manual) +
  labs(title = "Density Plot LDA", x = "LD1", fill = "Groups") +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 16),     # Tamanho da letra da legenda
    legend.title = element_text(size = 14),    # Tamanho do título da legenda
    plot.title = element_text(size = 16),      # Tamanho do título do plot
    axis.text.x = element_text(size = 14),     # Tamanho da letra do rótulo do eixo x
    axis.text.y = element_text(size = 14),     # Tamanho da letra do rótulo do eixo y
    axis.title.x = element_text(size = 16),    # Tamanho do título do eixo x
    axis.title.y = element_text(size = 16)     # Tamanho do título do eixo y
  )
ggsave("DensityPlot.svg", 
       plot=densidade, 
       device = "svg", 
       width = 12, 
       height = 8)
#ScatterPlot
scatterplot <- ggplot(data_df,
                      aes(x = LD1, y = 1,
                          color = Groups)) +
  geom_point(position = position_jitter(width = 0.1),
             size = 5) +
  scale_color_manual(values = cores_manual) +
  labs(title = "Scatterplot LDA",
       x = "LD1",
       y = "") +
  theme_bw() +
  theme(
    legend.text = element_text(size = 16),     # Tamanho da letra da legenda
    legend.title = element_text(size = 14),    # Tamanho do título da legenda
    plot.title = element_text(size = 16),      # Tamanho do título do plot
    axis.text.x = element_text(size = 14),     # Tamanho da letra do rótulo do eixo x
    axis.text.y = element_text(size = 14),     # Tamanho da letra do rótulo do eixo y
    axis.title.x = element_text(size = 16),    # Tamanho do título do eixo x
    axis.title.y = element_text(size = 16)     # Tamanho do título do eixo y
  )
ggsave("ScatterPlot.svg", 
       plot=scatterplot, 
       device = "svg", 
       width = 8, 
       height = 6)

#=============.    Unindo resultados de LDA ====================
#Todos os predicteds
LDA_HBV_Liver = predicted
LDA_HBV_PBMC = predicted
LDA_HCV_Liver = predicted
LDA_HCV_PBMC = predicted
LDA_HDV_Liver = predicted

#Prototipos: Unindo varios resultados de LDA
combine_lda_results <- function(lda_list) {
  # Verifica se a lista não está vazia
  if (length(lda_list) == 0) {
    stop("A lista de resultados do LDA está vazia.")
  }
  
  # Verifica se todos os elementos da lista têm a mesma estrutura
  if (!all(sapply(lda_list, function(lda) identical(names(lda), names(lda_list[[1]]))))) {
    stop("Os data.frames na lista não têm a mesma estrutura.")
  }
  
  # Combinação das informações usando rbind
  combined_lda <- list(
    class = factor(unlist(lapply(lda_list, function(lda) lda$class))),
    posterior = do.call(rbind, lapply(lda_list, function(lda) lda$posterior)),
    x = do.call(rbind, lapply(lda_list, function(lda) lda$x))
  )
  
  return(combined_lda)
}

# Exemplo de uso da função com predicted e predicted2
combined_predicted = combine_lda_results(list(LDA_HBV_Liver,
                                              LDA_HBV_PBMC,
                                              LDA_HCV_Liver,
                                              LDA_HCV_PBMC,
                                              LDA_HDV_Liver))

# Adicione uma coluna 'source' para identificar a origem dos dados
table(LDA_HBV_Liver$x) #17
table(LDA_HBV_PBMC$x) #63
table(LDA_HCV_Liver$x) #5
table(LDA_HCV_PBMC$x) #11
table(LDA_HDV_Liver$x) #48

combined_predicted$source = c(rep("HBV Liver", 17),
                              rep("HBV PBMC", 63),
                              rep("HCV Liver", 5),
                              rep("HCV PBMC", 11),
                              rep("HDV Liver", 48))
View(combined_predicted)
# Crie um dataframe com as informações combinadas
combined_df <- data.frame(
  LD1 = combined_predicted$x[, 1],
  source = combined_predicted$source
)
#Nomeando Infectados e Controle
combined_df$Groups <- ifelse(combined_df$LD1 < 0,
                             "control",
                             "infected")
View(combined_df)
#Carregando os dados de LDA combinados
combined_df = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/LDA Results Combined/LDA Results Combined.xlsx")
#Criando vetor de cores
cores_manual <- c("control" = "#93BFCF",
                  "infected" = "#D71313")
#Criando Vetor de Formas
formas <- c("HBV Liver" = 16,
            "HBV PBMC" = 17,
            "HCV Liver" = 15,
            "HCV PBMC" = 3,
            "HDV Liver" = 7)

# Crie um scatterplot com todos os LDA
tiff("/Users/adrielnobile/ScatterLDAJuntos.tiff",
     width = 8,
     height = 5,
     res = 400, units = 'in') 
ggplot(combined_df,
       aes(x = LD1, y = 1,
           color = Groups,
           shape = source)) +
  geom_point(position = position_jitter(width = 0.1),
             size = 5) +
  scale_color_manual(values = cores_manual) +
  scale_shape_manual(values = formas) +
  labs(title = "Scatterplot LDA",
       x = "LD1",
       y = "") +
  theme_bw(18) +
  theme(
    legend.text = element_text(size = 16),     # Tamanho da letra da legenda
    legend.title = element_text(size = 18),    # Tamanho do título da legenda
    plot.title = element_text(size = 16),      # Tamanho do título do plot
    axis.text.x = element_text(size = 14),     # Tamanho da letra do rótulo do eixo x
    axis.text.y = element_text(size = 14),     # Tamanho da letra do rótulo do eixo y
    axis.title.x = element_text(size = 16),    # Tamanho do título do eixo x
    axis.title.y = element_text(size = 16)     # Tamanho do título do eixo y
  )
dev.off()

#Alternativo:
# ggsave("ScatterPlot Todos LDA.svg", 
#        plot=scatterall, 
#        device = "svg", 
#        width = 8, 
#        height = 6)
# Crie um DensitiyPlot com todos os LDA
tiff("/Users/adrielnobile/DensityLDAJuntos.tiff",
     width = 15,
     height = 5,
     res = 400, units = 'in') 
ggplot(combined_df,
       aes(x = LD1, fill = Groups)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~source, nrow = 1) +
  scale_fill_manual(values = cores_manual) +
  labs(title = "Density Plot LDA", x = "LD1", fill = "Groups") +
  theme_bw(18) +
  theme(
    legend.text = element_text(size = 16),     # Tamanho da letra da legenda
    legend.title = element_text(size = 18),    # Tamanho do título da legenda
    plot.title = element_text(size = 16),      # Tamanho do título do plot
    axis.text.x = element_text(size = 14),     # Tamanho da letra do rótulo do eixo x
    axis.text.y = element_text(size = 14),     # Tamanho da letra do rótulo do eixo y
    axis.title.x = element_text(size = 16),    # Tamanho do título do eixo x
    axis.title.y = element_text(size = 16)     # Tamanho do título do eixo y
  )
dev.off()
#Alternativo:
# ggsave("DensityPlot Todos LDA.svg", 
#        plot=densityall, 
#        device = "svg", 
#        width = 6, 
#        height = 20)
#Salvando tabela
openxlsx::write.xlsx(HDV_Liver_genes,
                     "Top_50_Genes.xlsx") #Genes e Systems
#BarPlot com Genes
HBV_Liver_genes = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/HBV_Liver/Top_50_Genes.xlsx")
HBV_PBMC_genes = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/HBV_PBMC/Top_50_Genes.xlsx")
HCV_Liver_genes = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/HCV_Liver/Top_50_Genes.xlsx")
HCV_PBMC_genes = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/HCV_PBMC/Top_50_Genes.xlsx")
HDV_Liver_genes = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/HDV_Liver/Top_50_Genes.xlsx")
str(HBV_Liver_genes)
head(HBV_Liver_genes)
#Venn
#Listando Itens
x = list(HBV_Liver = HBV_Liver_genes$genes,
         HBV_PBMC = HBV_PBMC_genes$genes,
         HCV_Liver = HCV_Liver_genes$genes,
         HCV_PBMC = HCV_PBMC_genes$genes,
         HDV_Liver = HDV_Liver_genes$genes)
#Diagrama de Venn
liverson = ggVennDiagram(x,
                         color = ,
                         label = "count",
                         label_size = 12,
                         set_size = 12,
                         edge_size = 0.5) +
  scale_fill_gradient(low = "#F4FAFE",
                      high = "#4981BF") +
  theme(legend.position = "none")

#Funcao para extrair 5 maiores e 5 menores e salvar em outro df
extrair_extremos <- function(df) {
  # Selecionar os 5 maiores valores de LD1
  maiores <- df %>%
    slice_max(order_by = LD1, n = 10)
  
  # Selecionar os 5 menores valores de LD1
  menores <- df %>%
    slice_min(order_by = LD1, n = 10)
  
  # Combinar os resultados em um novo dataframe
  resultados <- bind_rows(maiores, menores)
  
  return(resultados)
}

#Aplicar a funcao
HBV_Liver_extremos <- extrair_extremos(HBV_Liver_genes)
HBV_PBMC_extremos <- extrair_extremos(HBV_PBMC_genes)
HCV_Liver_extremos <- extrair_extremos(HCV_Liver_genes)
HCV_PBMC_extremos <- extrair_extremos(HCV_PBMC_genes)
HDV_Liver_extremos <- extrair_extremos(HDV_Liver_genes)

All_Top_Genes = rbind(HBV_Liver_extremos,
                      HBV_PBMC_extremos,
                      HCV_Liver_extremos,
                      HCV_PBMC_extremos,
                      HDV_Liver_extremos)
All_Top_Genes$Condition = ifelse(All_Top_Genes$LD1 > 0, "Infected", "Control")
All_Top_Genes$LD1 <- abs(All_Top_Genes$LD1) #Transformando valores neg em pos
# Criar o gráfico de barras
cores <- c("Infected" = "#D71313", "Control" = "#93BFCF")

#Plotando. e salvando
tiff("/Users/adrielnobile/GenesBarPlot.tiff",
     width = 15,
     height = 5,
     res = 400, units = 'in') 
ggplot(All_Top_Genes, aes(x = reorder(genes, LD1, decreasing = TRUE) ,
                          y = LD1,
                          fill = Condition)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cores) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 20),  # Ajuste o tamanho do texto no eixo X
        axis.text.y = element_text(size = 18),  # Ajuste o tamanho do texto no eixo Y
        axis.title = element_text(size = 18),  # Ajuste o tamanho do título do eixo
        strip.text = element_text(size = 18)) +  # Ajuste o tamanho do texto na legenda
  labs(title = "Most representative genes for Control and Infected condition in LDA",
       x = "Genes",
       y = "LD1") +
  coord_flip() +
  facet_wrap(~hepatites,
             nrow = 1,
             scales = "free") +
  theme_bw(18)
dev.off()
#BarPlots Separados
#Carregar dados
HBV_Liver_genes = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/HBV_Liver/Top_50_Genes.xlsx")
HBV_PBMC_genes = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/HBV_PBMC/Top_50_Genes.xlsx")
HCV_Liver_genes = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/HCV_Liver/Top_50_Genes.xlsx")
HCV_PBMC_genes = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/HCV_PBMC/Top_50_Genes.xlsx")
HDV_Liver_genes = read.xlsx("~/Documents/DataSets/Artigo 1 - HCC NeuroImmune Dyregulation/5 - LDA/HDV_Liver/Top_50_Genes.xlsx")

#Um metodo
# obter_top_genes <- function(df, n = 20) {
#   return(head(df[order(abs(df$LD1), decreasing = TRUE), ], n))
# }
# HBV_Liver_extremos <- obter_top_genes(HBV_Liver_genes)
# HBV_PBMC_extremos <- obter_top_genes(HBV_PBMC_genes)
# HCV_Liver_extremos <- obter_top_genes(HCV_Liver_genes)
# HCV_PBMC_extremos <- obter_top_genes(HCV_PBMC_genes)
# HDV_Liver_extremos <- obter_top_genes(HDV_Liver_genes)

#Alternativo
extrair_extremos <- function(df) {
  # Selecionar os 5 maiores valores de LD1
  maiores <- df %>%
    slice_max(order_by = LD1, n = 10)
  
  # Selecionar os 5 menores valores de LD1
  menores <- df %>%
    slice_min(order_by = LD1, n = 10)
  
  # Combinar os resultados em um novo dataframe
  resultados <- bind_rows(maiores, menores)
  
  return(resultados)
}
HBV_Liver_extremos <- extrair_extremos(HBV_Liver_genes)
HBV_PBMC_extremos <- extrair_extremos(HBV_PBMC_genes)
HCV_Liver_extremos <- extrair_extremos(HCV_Liver_genes)
HCV_PBMC_extremos <- extrair_extremos(HCV_PBMC_genes)
HDV_Liver_extremos <- extrair_extremos(HDV_Liver_genes)
#Vetor de cores
cores <- c("Infected" = "#D71313", "Control" = "#93BFCF")
HBV_Liver_extremos$Condition = ifelse(HBV_Liver_extremos$LD1 > 0,
                                      "Infected", "Control")
HBV_PBMC_extremos$Condition = ifelse(HBV_PBMC_extremos$LD1 > 0,
                                     "Infected", "Control")
HCV_Liver_extremos$Condition = ifelse(HCV_Liver_extremos$LD1 > 0,
                                      "Infected", "Control")
HCV_PBMC_extremos$Condition = ifelse(HCV_PBMC_extremos$LD1 > 0,
                                     "Infected", "Control")
HDV_Liver_extremos$Condition = ifelse(HDV_Liver_extremos$LD1 > 0,
                                      "Infected", "Control")
All_Top_Genes = rbind(HBV_Liver_extremos,
                      HBV_PBMC_extremos,
                      HCV_Liver_extremos,
                      HCV_PBMC_extremos,
                      HDV_Liver_extremos)
#BarPlot
tiff("/Users/adrielnobile/GenesBarPlot.tiff",
     width = 5,
     height = 6,
     res = 300, units = 'in') 
All_Top_Genes %>%
  mutate(genes = reorder(genes, LD1)) %>%
  ggplot(aes(x = genes,
             y = LD1,
             fill = Condition)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cores) +
  theme(axis.text.x = element_blank(),  # Remove rótulos do eixo X
        axis.text.y = element_blank(),  # Remove rótulos do eixo Y
        axis.title = element_blank(),   # Remove títulos dos eixos
        strip.text = element_blank(),   # Remove texto na legenda
        legend.position = "none") +     # Remove a legenda
  coord_flip() +
  facet_wrap(~hepatites,
             nrow = 1,
             scales = "free") +
  theme_bw(18)
dev.off()

#Salvando
openxlsx::write.xlsx(All_Top_Genes,
                     "All_Top_Genes.xlsx") #Genes e Systems
#Alternativo
# ggsave("barplotGenes.svg", 
#        plot=barplotGenes, 
#        device = "svg", 
#        width = 10, 
#        height = 8)