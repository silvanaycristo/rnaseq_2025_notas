## ----model.matrix---------------------------------------------
## ?model.matrix
## Crea un modelo estadístico
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat
colnames(mat)


## ----lm_example-----------------------------------------------
summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))

## En la tabla de coeficientes nos dice el estimado de beta 0, 1 y 2

## ----EMM_example1---------------------------------------------
## Datos de ejemplo
(sampleData <- data.frame(
  genotype = rep(c("A", "B"), each = 4),
  treatment = rep(c("ctrl", "trt"), 4)
))

## Creemos las imágenes usando ExploreModelMatrix
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment,
  textSizeFitted = 4
)


## Veamos las imágenes
cowplot::plot_grid(plotlist = vd$plotlist)

mat2 <- with(
  sampleData,
  model.matrix(~ genotype + treatment)
)
mat2
sampleData

cbind(mat2, sampleData)
colnames(mat2)

## ----Ejemplo 2 adicional-------------------

(sampleData <- data.frame(
  Response = rep(c("Resistant", "Sensitive"), c(12, 18)),
  Patient = factor(rep(c(1:6, 8, 11:18), each = 2)),
  Treatment = factor(rep(c("pre","post"), 15)),
  ind.n = factor(rep(c(1:6, 2, 5:12), each = 2))))
#>     Response Patient Treatment ind.n
#> 1  Resistant       1       pre     1
#> 2  Resistant       1      post     1
#> 3  Resistant       2       pre     2
#> 4  Resistant       2      post     2
#> 5  Resistant       3       pre     3
#> 6  Resistant       3      post     3
#> 7  Resistant       4       pre     4
#> 8  Resistant       4      post     4
#> 9  Resistant       5       pre     5
#> 10 Resistant       5      post     5
#> 11 Resistant       6       pre     6
#> 12 Resistant       6      post     6
#> 13 Sensitive       8       pre     2
#> 14 Sensitive       8      post     2
#> 15 Sensitive      11       pre     5
#> 16 Sensitive      11      post     5
#> 17 Sensitive      12       pre     6
#> 18 Sensitive      12      post     6
#> 19 Sensitive      13       pre     7
#> 20 Sensitive      13      post     7
#> 21 Sensitive      14       pre     8
#> 22 Sensitive      14      post     8
#> 23 Sensitive      15       pre     9
#> 24 Sensitive      15      post     9
#> 25 Sensitive      16       pre    10
#> 26 Sensitive      16      post    10
#> 27 Sensitive      17       pre    11
#> 28 Sensitive      17      post    11
#> 29 Sensitive      18       pre    12
#> 30 Sensitive      18      post    12
vd <- VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ Response + Response:ind.n + Response:Treatment,
  textSizeFitted = 3
)
cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

## Interpreta ResponseResistant.Treatmentpre del ejercicio 2.

## 1) Pre v.s post da el orden donde un valor positivo indica mayor expresión en el pre.
## 2) Individuo constante.
## 3) La respuesta es resistente.


## ----Ejemplo 3 (adicional)-------------------

(sampleData = data.frame(
  condition = factor(rep(c("ctrl_minus", "ctrl_plus",
                           "ko_minus", "ko_plus"), 3)),
  batch = factor(rep(1:6, each = 2))))
#>     condition batch
#> 1  ctrl_minus     1
#> 2   ctrl_plus     1
#> 3    ko_minus     2
#> 4     ko_plus     2
#> 5  ctrl_minus     3
#> 6   ctrl_plus     3
#> 7    ko_minus     4
#> 8     ko_plus     4
#> 9  ctrl_minus     5
#> 10  ctrl_plus     5
#> 11   ko_minus     6
#> 12    ko_plus     6
vd <- VisualizeDesign(sampleData = sampleData,
                      designFormula = ~ 0 + batch + condition,
                      textSizeFitted = 4, lineWidthFitted = 20,
                      dropCols = "conditionko_minus")
cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

## El intercepto da un grupo de referencia que se debe eliminar.
## El planteamiento de la pregunta biológica es fundamental.

## ----EMM_example1_interactive, eval = FALSE-------------------
# ## Usaremos shiny otra ves
app <- ExploreModelMatrix::ExploreModelMatrix(
     sampleData = sampleData,
     designFormula = ~ genotype + treatment
 )
if (interactive()) shiny::runApp(app)


## ----download_SRP045638---------------------------------------
library("recount3")

human_projects <- available_projects()

rse_gene_SRP045638 <- create_rse(
  subset(
    human_projects,
    project == "SRP045638" & project_type == "data_sources"
  )
)
assay(rse_gene_SRP045638, "counts") <- compute_read_counts(rse_gene_SRP045638)


## ----describe_issue-------------------------------------------
rse_gene_SRP045638$sra.sample_attributes[1:3]

## dev_stage;;Fetal se encuentra solamente en el primer output. Hay que eliminarla.

## ----solve_issue----------------------------------------------
rse_gene_SRP045638$sra.sample_attributes <- gsub(
  "dev_stage;;Fetal\\|",
  "",
  rse_gene_SRP045638$sra.sample_attributes
  )
rse_gene_SRP045638$sra.sample_attributes[1:3]

## Limpieza de datos.

## ----attributes-----------------------------------------------
rse_gene_SRP045638 <- expand_sra_attributes(rse_gene_SRP045638)

colData(rse_gene_SRP045638)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP045638)))
]


## ----re_cast--------------------------------------------------
## Pasar de character a numeric o factor
rse_gene_SRP045638$sra_attribute.age <- as.numeric(rse_gene_SRP045638$sra_attribute.age)
rse_gene_SRP045638$sra_attribute.disease <- factor(tolower(rse_gene_SRP045638$sra_attribute.disease))
rse_gene_SRP045638$sra_attribute.RIN <- as.numeric(rse_gene_SRP045638$sra_attribute.RIN)
rse_gene_SRP045638$sra_attribute.sex <- factor(rse_gene_SRP045638$sra_attribute.sex)

## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP045638)[
  ,
  grepl("^sra_attribute.[age|disease|RIN|sex]", colnames(colData(rse_gene_SRP045638)))
]))


## ----new_variables--------------------------------------------
## Encontraremos diferencias entre muestra prenatalas vs postnatales
rse_gene_SRP045638$prenatal <- factor(ifelse(rse_gene_SRP045638$sra_attribute.age < 0, "prenatal", "postnatal"))
table(rse_gene_SRP045638$prenatal)

## http://rna.recount.bio/docs/quality-check-fields.html
rse_gene_SRP045638$assigned_gene_prop <- rse_gene_SRP045638$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP045638$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP045638$assigned_gene_prop)
with(colData(rse_gene_SRP045638), plot(assigned_gene_prop, sra_attribute.RIN))

## Hm... veamos si hay una diferencia entre los grupos
with(colData(rse_gene_SRP045638), tapply(assigned_gene_prop, prenatal, summary))


## ----filter_rse-----------------------------------------------
## Guardemos nuestro objeto entero por si luego cambiamos de opinión
rse_gene_SRP045638_unfiltered <- rse_gene_SRP045638

## Eliminemos a muestras malas
hist(rse_gene_SRP045638$assigned_gene_prop)
table(rse_gene_SRP045638$assigned_gene_prop < 0.3)
rse_gene_SRP045638 <- rse_gene_SRP045638[, rse_gene_SRP045638$assigned_gene_prop > 0.3]

## Calculemos los niveles medios de expresión de los genes en nuestras
## muestras.
## Ojo: en un análisis real probablemente haríamos esto con los RPKMs o CPMs
## en vez de las cuentas.
## En realidad usariamos:
# edgeR::filterByExpr() https://bioconductor.org/packages/edgeR/ https://rdrr.io/bioc/edgeR/man/filterByExpr.html
# genefilter::genefilter() https://bioconductor.org/packages/genefilter/ https://rdrr.io/bioc/genefilter/man/genefilter.html
# jaffelab::expression_cutoff() http://research.libd.org/jaffelab/reference/expression_cutoff.html
#
gene_means <- rowMeans(assay(rse_gene_SRP045638, "counts"))
summary(gene_means)

## Eliminamos genes
rse_gene_SRP045638 <- rse_gene_SRP045638[gene_means > 0.1, ]

## Dimensiones finales
dim(rse_gene_SRP045638)

## Porcentaje de genes que retuvimos
round(nrow(rse_gene_SRP045638) / nrow(rse_gene_SRP045638_unfiltered) * 100, 2)


## ----normalize------------------------------------------------
library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene_SRP045638, "counts"),
  genes = rowData(rse_gene_SRP045638)
)
dge <- calcNormFactors(dge)


## ----explore_gene_prop_by_age---------------------------------
library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP045638)), aes(y = assigned_gene_prop, x = prenatal)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Age Group")

## Para hacer el análisis de expresión, hay que definir el modelo estadístico.


## ----statiscal_model------------------------------------------
mod <- model.matrix(~ prenatal + sra_attribute.RIN + sra_attribute.sex + assigned_gene_prop,
                    data = colData(rse_gene_SRP045638)
)
colnames(mod)


## ----run_limma------------------------------------------------
library("limma") ## Se le pasa el objeto de dg list
vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP045638),
  sort.by = "none"
)
dim(de_results)
head


## Genes diferencialmente expresados entre pre y post natal con FDR < 5%
table(de_results$adj.P.Val < 0.05)

## Visualicemos los resultados estadísticos
plotMA(eb_results, coef = 2)

## Un valor positivo en el log-fold es un gen con mayor expresión en el grupo pre contra el post.
## prenatalprenatal -> variable y valor de contraste.

volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)
de_results[de_results$gene_name %in% c("ZSCAN2", "VASH2", "KIAA0922"), ]


## ----pheatmap-------------------------------------------------
## Extraer valores de los genes de interés
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

## Creemos una tabla con información de las muestras
## y con nombres de columnas más amigables
df <- as.data.frame(colData(rse_gene_SRP045638)[, c("prenatal", "sra_attribute.RIN", "sra_attribute.sex")])
colnames(df) <- c("AgeGroup", "RIN", "Sex")

## Hagamos un heatmap
library("pheatmap")
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)


## ----plot_mds-------------------------------------------------
## Para colores
library("RColorBrewer")

## Conviertiendo los grupos de edad a colores
col.group <- df$AgeGroup
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

## MDS por grupos de edad
plotMDS(vGene$E, labels = df$AgeGroup, col = col.group)

## Conviertiendo los valores de Sex a colores
col.sex <- df$Sex
levels(col.sex) <- brewer.pal(nlevels(col.sex), "Dark2")
col.sex <- as.character(col.sex)

## MDS por sexo
plotMDS(vGene$E, labels = df$Sex, col = col.sex)


## ----respuesta, out.height="1100px"---------------------------
## Tenemos que usar gene_id y gene_name
rowRanges(rse_gene_SRP045638)

## Alternativamente, podriamos haber usado de_results
head(de_results, n = 3)

## Es la misma información
identical(rowRanges(rse_gene_SRP045638)$gene_id, de_results$gene_id)
identical(rowRanges(rse_gene_SRP045638)$gene_name, de_results$gene_name)

## Guardemos los IDs de nuestros 50 genes
nombres_originales <- rownames(exprs_heatmap)

## Con match() podemos encontrar cual es cual
rownames(exprs_heatmap) <- rowRanges(rse_gene_SRP045638)$gene_name[
  match(rownames(exprs_heatmap), rowRanges(rse_gene_SRP045638)$gene_id)
]

## Vean que tambien podriamos haber usado rank()
identical(
  which(rank(de_results$adj.P.Val) <= 50),
  match(nombres_originales, rowRanges(rse_gene_SRP045638)$gene_id)
)

## Esta es otra solución
identical(
  de_results$gene_name[rank(de_results$adj.P.Val) <= 50],
  rownames(exprs_heatmap)
)

## Por último podemos cambiar el valor de show_rownames de FALSE a TRUE
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df
)

## Guardar la imagen en un PDF largo para poder ver los nombres de los genes
pdf("pheatmap_con_nombres.pdf", height = 14, useDingbats = FALSE)
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df
)
dev.off()


## ----"centered_and_scaled", out.height="1100px"---------------
## Versión con centering y scaling en los renglones (los genes)
pheatmap::pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df,
  scale = "row"
)


## ----"complexheatmap", out.height="1100px"--------------------
## Misma versión pero ahora con ComplexHeatmap en vez del paquete pheatmap
ComplexHeatmap::pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df,
  scale = "row"
)
