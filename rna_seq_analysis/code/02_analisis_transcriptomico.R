#---------------------------------------------------------------------------------------------------------------------------------------------------------
# Script name:  02_analisis_transcriptomico.R
#
# Description:  Este script realiza el procesamiento secundario de los datos de RNA-seq correspondientes al dataset OSD-635 
#               del repositorio OSDR de la NASA. Toma como entrada los archivos .sf con las cuantificaciones a nivel de transcrito 
#               generadas por Salmon y los metadatos del experimento. Ejecuta la importación, el filtrado, 
#               la normalización y el análisis de expresión génica diferencial (DGE) mediante edgeR. 
#               Posteriormente, realiza un análisis de enriquecimiento funcional (ORA y GSEA) sobre las bases 
#               de datos de Gene Ontology (GO) y KEGG utilizando clusterProfiler, y exporta los resultados y gráficos.
#
# Author:       Laurentiu Jalba Staver
# Institution:  Universidad Internacional de Valencia
# Date:         07/2026
# Version:      1.0
# Dependencies: tximport, edgeR, EnsDb.Hsapiens.v86, org.Hs.eg.db, ggplot2, pheatmap, clusterProfiler, ggrepel, enrichplot, pathview
#---------------------------------------------------------------------------------------------------------------------------------------------------------

# 0. CARGA DE PAQUETES -----------------------------------------------------------------------------------------------------------------------------------

library(tximport)
library(edgeR)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(ggrepel)
library(enrichplot)
library(pathview)

# Creación de directorio de salida para los resultados
dir_out <- "results/R_results"
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

# 1.IMPORTACION DE METADATOS -------------------------------------------------------------------------------------------------------------------------------

# 1.1 Se importa el archivo de metadatos
metadata_raw <- read.delim("data/metadata/OSD-635_metadata_OSD-635-ISA/s_OSD-635.txt", sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# 1.2 Se seleccionan las columnas con los datos de interés
metadata <- metadata_raw[, c("Sample Name", "Factor Value[Spaceflight]")]

# 1.3 Se cambian los nombres de las columnas para mayor claridad
colnames(metadata) <- c("Muestra", "Condicion")

# 1.4 Se convierte la columna muestra en índice
rownames(metadata) <- metadata$Muestra
metadata$Muestra <- NULL

# 1.5 Se modifica el nombre de "Condicion" para sustituir el espacio en blanco por "_"
metadata$Condicion <- gsub(" ", "_", metadata$Condicion)

# 1.6 La columna "Condicion" se convierte en factor
metadata$Condicion <- as.factor(metadata$Condicion)
grupo <- metadata$Condicion

# 2. IMPORTACION DE ARCHIVOS DE CONTEOS -------------------------------------------------------------------------------------------------------------------

# 2.1 Se obtiene la ruta de los archivos de conteos
archivos_conteos_raw <- list.files("results", pattern = "_quant\\.sf$", full.names = TRUE)

# 2.2 Se asigna el nombre de muestra a su correspondiente archivo
archivos_conteos <- sapply(rownames(metadata), function(Muestra) {
  grep(pattern = Muestra, x = archivos_conteos_raw, value = TRUE)
})

# 2.3 Asignación de los transcritos a genes
tran_df <- transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
tran_to_gene <- as.data.frame(tran_df[,c("tx_id", "gene_id")])

# Se eliminan los sufijos de versión de los identificadores de ENSEMBL para garantizar la coincidencia exacta con las
# cuantificaciones de Salmon
tran_to_gene$tx_id <- sub("\\..*", "", tran_to_gene$tx_id)     
tran_to_gene$gene_id <- sub("\\..*", "", tran_to_gene$gene_id)

# 2.4 Importación con tximport
txi <- tximport(archivos_conteos, 
                type = "salmon", 
                tx2gene = tran_to_gene, 
                countsFromAbundance = "lengthScaledTPM",
                ignoreTxVersion = TRUE)

cat("Genes cuantificados:", nrow(txi$counts))

# 3. ANOTACIÓN DE LOS GENES --------------------------------------------------------------------------------------------------------------------------------

genes_ids <- rownames(txi$counts)
anotacion <- AnnotationDbi::select(org.Hs.eg.db, 
                    keys = genes_ids, 
                    columns = c("ENTREZID", "SYMBOL", "GENENAME"), 
                    keytype = "ENSEMBL")

# 3.1 Manejo de nombres duplicados
anotacion <- anotacion[order(anotacion$ENSEMBL, is.na(anotacion$ENTREZID)), ]
anotacion <- anotacion[!duplicated(anotacion$ENSEMBL), ]
rownames(anotacion) <- anotacion$ENSEMBL

# 4. CREACIÓN DEL OBJETO DGEList Y FILTRADO DE GENES DE BAJA EXPRESIÓN O NULA -------------------------------------------------------------------------------

y <- DGEList(counts = txi$counts, group = grupo, genes = anotacion)
genes_filtrados <- filterByExpr(y)

# Se indica que se recalcule el tamaño efectivo de las bibliotecas con los genes de baja expresión o nula excluidos para
# mejorar la precisión estadística.
y_filtrado <- y[genes_filtrados, , keep.lib.sizes = FALSE]
genes_mantenidos <- (sum(genes_filtrados) / length(genes_filtrados)) * 100

cat("Genes que se mantienen tras el filtrado:", sum(genes_filtrados),"-",round(genes_mantenidos, 2), "%\n")
cat("Genes anotados con ENTREZID:", sum(!is.na(y_filtrado$genes$ENTREZID)))

# 5. NORMALIZACIÓN TMM Y GRÁFICO MDS ------------------------------------------------------------------------------------------------------------------------

# 5.1. Normalización
y_filtrado <- calcNormFactors(y_filtrado)

# 5.2. Gráfico MDS
colores <- c("darkblue", "orange")
png(file.path(dir_out, "Grafico_MDS.png"), width = 800, height = 600, res = 120)
plotMDS(y_filtrado, col = colores[as.numeric(grupo)], pch = 16, cex = 2)
legend("topleft", legend = levels(grupo), pch = 16, col = colores)
dev.off()

# 6. ESTIMACIÓN DE LA DISPERSIÓN Y MODELO ESTADÍSTICO ------------------------------------------------------------------------------------------------------

# 6.1 Matriz de diseño
design <- model.matrix(~ 0 + grupo)
colnames(design) <- levels(grupo)

# 6.2 Estimación de la dispersión
y_filtrado <- estimateDisp(y_filtrado, design, robust=TRUE)
cat("Dispersión común estimada:", y_filtrado$common.dispersion, "\n")
png(file.path(dir_out, "Grafico_BCV.png"), width = 800, height = 600, res = 120)
plotBCV(y_filtrado)
dev.off()

# 6.3 Ajuste del modelo lineal de regresión generalizada (GLM)
# Se hace uso del método de cuasi-verosimilitud con la opción robust = TRUE para reducir la influencia de valores atípicos
# y mejorar la precisión estadística, especialmente en experimentos con pocas réplicas biológicas, como es el caso.
fit <- glmQLFit(y_filtrado, design, robust=TRUE)
png(file.path(dir_out, "Grafico_QLDisp.png"), width = 800, height = 600, res = 120)
plotQLDisp(fit)
dev.off()

# 7. EXPRESIÓN DIFERENCIAL ---------------------------------------------------------------------------------------------------------------------------------

# 7.1 Se define el contraste
espacio_vs_tierra <- makeContrasts(Space_Flight - Ground_Control, levels = design)

# 7.2 Se aplica el test estadístico (QLF)
test_hn <- glmQLFTest(fit, contrast = espacio_vs_tierra)

# 7.3 Corrección de falsos positivos (FDR)
test_hn_corrected <- topTags(test_hn, n = Inf)

# 7.4 Filtrado por FDR y logFC
genes_DE <- decideTests(test_hn, adjust.method = "BH", p.value = 0.05, lfc = 1)
cat("Resultados de los genes diferencialmente expresados:\n")
summary(genes_DE)

# 8. VISUALIZACIÓN DE RESULTADOS ----------------------------------------------------------------------------------------------------------------------------

# 8.1 Gráfico de volcán
# Preparación de datos
data_volcano <- test_hn_corrected$table
data_volcano$DE <- "Not significant"
data_volcano$DE[data_volcano$logFC > 1 & data_volcano$FDR < 0.05] <- "Up"
data_volcano$DE[data_volcano$logFC < -1 & data_volcano$FDR < 0.05] <- "Down"

# Número de genes en leyenda
n_up <- sum(data_volcano$DE == "Up")
n_down <- sum(data_volcano$DE == "Down")
n_ns <- sum(data_volcano$DE == "Not significant")

# Creación de las etiquetas para la leyenda
etiquetas_leyenda <- c(
  "Down" = paste0("Down (", n_down, ")"),
  "Not significant" = paste0("Not significant (", n_ns, ")"),
  "Up" = paste0("Up (", n_up, ")")
)

# Preparación de datos para etiquetar genes en el gráfico
up_sig <- data_volcano[data_volcano$DE == "Up", ]
down_sig <- data_volcano[data_volcano$DE == "Down", ]
up_sig_clean <- up_sig[!is.na(up_sig$SYMBOL), ]
down_sig_clean <- down_sig[!is.na(down_sig$SYMBOL), ]
top_up <- head(up_sig_clean[order(-up_sig_clean$logFC, up_sig_clean$FDR), ], 5)
top_down <- head(down_sig_clean[order(down_sig_clean$logFC, down_sig_clean$FDR), ], 5)
genes_etiqueta <- rbind(top_up, top_down)

grafico_volcano <- ggplot(data_volcano, aes(x = logFC, y = -log10(FDR), col = DE)) + 
  geom_point(size = 0.5, alpha = 0.7) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("Down" = "blue", 
                                "Not significant" = "grey80", 
                                "Up" = "red"),
                     labels = etiquetas_leyenda) +
  geom_label_repel(data = genes_etiqueta, 
                   aes(label = SYMBOL), 
                   size = 3.5, 
                   color = "black", 
                   fontface = "bold", 
                   box.padding = 0.5,
                   force = 2,
                   max.overlaps = Inf, 
                   show.legend = FALSE) + 
  theme_minimal() + 
  labs(x = "Log2 Fold Change", y = "-Log10(FDR)")
print(grafico_volcano)
ggsave(file.path(dir_out, "Grafico_Volcano.png"), plot = grafico_volcano, width = 8, height = 6, dpi = 300)

# 8.2 Matriz de correlación y Heatmap de top 50 genes diferencialmente expresados
# Transformación a logCPM
logcpm <- cpm(y_filtrado, log = TRUE)
rownames(logcpm) <- make.unique(y_filtrado$genes$SYMBOL)
colnames(logcpm) <- paste(y_filtrado$samples$group, 1:3, sep = "_")

# Matriz de correlación
cor_matrix <- cor(logcpm)
pheatmap(cor_matrix, 
         display_numbers = TRUE, 
         clustering_method = "ward.D2",
         filename = file.path(dir_out, "Grafico_Correlacion.png"), width = 7, height = 6)

# Selección de genes diferencialmente expresados para el mapa de calor
DEG <- test_hn_corrected$table[test_hn_corrected$table$FDR <= 0.05 & abs(test_hn_corrected$table$logFC) >= 1 ,]
DEG_ordenado <- DEG[order(DEG$FDR), ]
DEG_ordenado <- DEG_ordenado[!is.na(DEG_ordenado$SYMBOL), ]
DEG_top50 <- head(DEG_ordenado, 50)
DEG_selection <- logcpm[DEG_top50$SYMBOL, ] 
cat("Genes diferencialmente expresados:", nrow(DEG))
cat("Porcentaje de genes diferencialmente expresados:", nrow(DEG)/nrow(y_filtrado)*100)

# Mapa de calor
pheatmap(DEG_selection, scale = "row", 
         cluster_rows = T, 
         cluster_cols = T, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2", 
         cutree_cols = 2,
         cutree_rows = 2,
         display_numbers = F, 
         fontsize_row = 6, 
         border_color = NA,
         filename = file.path(dir_out, "Grafico_Heatmap_Top50.png"), width = 8, height = 10)

# 9. ANÁLISIS FUNCIONAL -------------------------------------------------------------------------------------------------------------------------------------

# 9.1 Preparación de los datos
DEG_limpio <- DEG[!is.na(DEG$ENTREZID), ]
genes_up <- unique(as.character(DEG_limpio$ENTREZID[DEG_limpio$logFC >= 1]))
genes_down <- unique(as.character(DEG_limpio$ENTREZID[DEG_limpio$logFC <= -1]))
universo <- unique(na.omit(y_filtrado$genes$ENTREZID))
cat("DEGs anotados:", nrow(DEG_limpio))
cat("Cobertura de DEGs anotados", round(100*nrow(DEG_limpio)/nrow(DEG),2))

# Función para ontología de genes
go <- function(lista_genes, uni, ontologia, titulo) {
  res <- enrichGO(gene = lista_genes, 
                  universe = uni, 
                  OrgDb = org.Hs.eg.db, 
                  keyType = "ENTREZID", 
                  ont = ontologia, 
                  pAdjustMethod = "BH", 
                  pvalueCutoff = 0.05)
  grafico <- dotplot(res, showCategory = 10)
  print(grafico)
  nombre_archivo <- file.path(dir_out, paste0("Grafico_", gsub(" ", "_", titulo), ".png"))
  ggsave(filename = nombre_archivo, plot = grafico, width = 8, height = 6, dpi = 300)
  return(res) 
}

# Función para rutas metabólicas kegg
kegg <- function(lista_genes, uni, titulo) {
  res <- enrichKEGG(gene = lista_genes, 
                    universe = uni, 
                    organism = "hsa", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05)
  grafico <- barplot(res, showCategory = 10)
  print(grafico)
  nombre_archivo <- file.path(dir_out, paste0("Grafico_KEGG_", gsub(" ", "_", titulo), ".png"))
  ggsave(filename = nombre_archivo, plot = grafico, width = 8, height = 6, dpi = 300)
  return(res)
}

# 9.2 GENE ONTOLOGY
# 9.2.1 GO genes sobreexpresados
go_up_bp <- go(genes_up, universo, "BP", "GO-BP  genes sobreexpresados")
tabla_go_up <- as.data.frame(go_up_bp)
head(tabla_go_up[, c("ID", "Description", "p.adjust", "Count")], 10)

go_up_cc <- go(genes_up, universo, "CC", "GO-CC  genes sobreexpresados")
tabla_go_up_cc <- as.data.frame(go_up_cc)
head(tabla_go_up_cc[, c("ID", "Description", "p.adjust", "Count")], 10)

# 9.2.2 GO genes infraexpresados
go_down_bp <- go(genes_down, universo, "BP", "GO-BP  genes infraexpresados")
tabla_go_down <- as.data.frame(go_down_bp)
head(tabla_go_down[, c("ID", "Description", "p.adjust", "Count")], 10)

go_down_cc <- go(genes_down, universo, "CC", "GO-CC  genes infraexpresados")
tabla_go_down_cc <- as.data.frame(go_down_cc)
head(tabla_go_down_cc[, c("ID", "Description", "p.adjust", "Count")], 10)


# 9.3 RUTAS METABÓLICAS KEGG
# 9.3.1 Kegg genes sobreexpresados
kegg_up <- kegg(genes_up, universo, "Rutas KEGG enriquecidas en genes sobreexpresados")
tabla_kegg_up <- as.data.frame(kegg_up)
print(tabla_kegg_up)

# 9.3.2 Kegg genes infraexpresados
kegg_down <- kegg(genes_down, universo, "Rutas KEGG enriquecidas en genes infraexpresados")
tabla_kegg_down <- as.data.frame(kegg_down)
print(tabla_kegg_down)

# 9.4a Mapa ruta ciclo celular
vector_fc <- DEG_limpio$logFC
names(vector_fc) <- DEG_limpio$ENTREZID
ruta_met <- pathview(gene.data = vector_fc, 
                     pathway.id = "hsa04110", 
                     species = "hsa", 
                     limit = list(gene = max(abs(vector_fc)), cpd = 1))

archivo_pathview <- "hsa04110.pathview.png"
if (file.exists(archivo_pathview)) {
  file.rename(from = archivo_pathview, to = file.path(dir_out, "Grafico_Pathview_ciclo_celular.png"))
} else {
  cat("Pathview no generó la imagen")
}

# 9.4b Mapa ruta cadherinas
vector_fc <- DEG_limpio$logFC
names(vector_fc) <- DEG_limpio$ENTREZID
ruta_met <- pathview(gene.data = vector_fc, 
                     pathway.id = "hsa04519", 
                     species = "hsa", 
                     limit = list(gene = max(abs(vector_fc)), cpd = 1))

archivo_pathview <- "hsa04519.pathview.png"
if (file.exists(archivo_pathview)) {
  file.rename(from = archivo_pathview, to = file.path(dir_out, "Grafico_Pathview_cadherinas.png"))
} else {
  cat("Pathview no generó la imagen")
}

# 9.5 GENE SET ENRICHMENT ANALYSIS (GSEA)
set.seed(123)

# Extracción de la tabla con los genes
tabla_genes <- test_hn_corrected$table[!is.na(test_hn_corrected$table$ENTREZID), ]
tabla_genes <- tabla_genes[!duplicated(tabla_genes$ENTREZID), ]

# Se extraen de la tabla los valores de logFC
tabla_genes_logfc <- tabla_genes$logFC

# Se asigna el ENTREZID a cada valor de logFC
names(tabla_genes_logfc) <- tabla_genes$ENTREZID

# Se ordena el vector que almacena el logFC de mayor a menor
tabla_genes_logfc <- sort(tabla_genes_logfc, decreasing = T)

# Aplicación de GSEA para procesos biológicos
gsea_go <- gseGO(geneList = tabla_genes_logfc,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 eps = 0)

# Se grafican los resultados para procesos biológicos
grafico_gsea_go <- dotplot(gsea_go, showCategory = 10, split = ".sign", label_format = 45) + 
  facet_grid(.~.sign) +
  labs(x = "Gene Ratio", 
       y = "Biological Process (BP)") +
  theme(strip.text = element_text(face = "bold", size = 12),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.title.x = element_text(margin = margin(t = 20))) 
print(grafico_gsea_go)
ggsave(file.path(dir_out, "Grafico_GSEA_GO.png"), plot = grafico_gsea_go, width = 10, height = 9, dpi = 300)

print(grafico_gsea_go$data)

# Aplicación de GSEA para rutas KEGG
gsea_kegg <- gseKEGG(geneList = tabla_genes_logfc,
                     organism = "hsa",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     eps = 0)

tabla_gsea_kegg <- as.data.frame(gsea_kegg)
head(tabla_gsea_kegg, 10)

# Se grafican los resultados para rutas KEGG
grafico_gsea_kegg <- dotplot(gsea_kegg, showCategory = 10, split = ".sign", label_format = 45) + 
  facet_grid(.~.sign) +
  labs(x = "Gene Ratio", 
       y = "Rutas metabólicas y de señalización (KEGG)") +
  theme(strip.text = element_text(face = "bold", size = 12),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.title.x = element_text(margin = margin(t = 20)))
print(grafico_gsea_kegg)
ggsave(file.path(dir_out, "Grafico_GSEA_KEGG.png"), plot = grafico_gsea_kegg, width = 10, height = 9, dpi = 300)

print(grafico_gsea_kegg$data)

# 10. EXPORTACIÓN DE RESULTADOS -----------------------------------------------------------------------------------------------------------------------------

# 10.1 Tablas de expresión de genes
write.csv(test_hn_corrected$table, file = file.path(dir_out, "Resultados_01_todos_genes.csv"), row.names = F)
write.csv(DEG, file = file.path(dir_out, "Resultados_02_DEGs.csv"), row.names = F)
write.csv(DEG_top50, file = file.path(dir_out, "Resultados_03_top50_FDR.csv"), row.names = F)
write.csv(as.data.frame(logcpm), file = file.path(dir_out, "Resultados_04_logCPM.csv"), row.names = T)

# 10.2 Ontología de genes
write.csv(as.data.frame(go_up_bp), file = file.path(dir_out, "Resultados_05_GO_BP_UP.csv"), row.names = F)
write.csv(as.data.frame(go_up_cc), file = file.path(dir_out, "Resultados_06_GO_CC_UP.csv"), row.names = F)
write.csv(as.data.frame(go_down_bp), file = file.path(dir_out, "Resultados_07_GO_BP_DOWN.csv"), row.names = F)
write.csv(as.data.frame(go_down_cc), file = file.path(dir_out, "Resultados_08_GO_CC_DOWN.csv"), row.names = F)

# 10.3 Rutas metabólicas KEGG
write.csv(as.data.frame(kegg_up), file = file.path(dir_out, "Resultados_09_KEGG_UP.csv"), row.names = F)
write.csv(as.data.frame(kegg_down), file = file.path(dir_out, "Resultados_10_KEGG_DOWN.csv"), row.names = F)

# 10.4 Análisis GSEA
write.csv(as.data.frame(gsea_go), file = file.path(dir_out, "Resultados_11_GSEA_GO.csv"), row.names = F)
write.csv(as.data.frame(gsea_kegg), file = file.path(dir_out, "Resultados_12_GSEA_KEGG.csv"), row.names = F)

# 10.5 Información de la sesión
writeLines(
  capture.output(sessionInfo()),
  file.path(dir_out, "Resultados_sessioninfo.txt")
)







