## ----'start', message=FALSE-----------------------------------
## Load recount3 R package
library("recount3")

## recount3 tiene muchas dependiencias y carga el paquete summarized
## experiment, por ello arroja varios mensajes como output

## ----'quick_example'------------------------------------------
## Revisemos todos los proyectos con datos de humano en recount3
human_projects <- available_projects()

## class(human_projects)
## dim(human_projects)
## head(human_projects[order(human_projects$n_samples), ])


## Encuentra tu proyecto de interés. Aquí usaremos
## SRP009615 de ejemplo
proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)

## class(proj_info)
## dim(proj_info)

## Crea un objeto de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)
## Explora el objeto RSE
rse_gene_SRP009615

## Tenemos 63856 genes en este objeto.
## Tenemos 12 muestras.
## Tenemos 175 variables con información de las muestras.
## Tenemos 10 variables de información de los genes.
## Tenemos 8 variables de información del objeto de metadata.

## ----"interactive_display", eval = FALSE----------------------
# ## Explora los proyectos disponibles de forma interactiva
# proj_info_interactive <- interactiveDisplayBase::display(human_projects)
# ## Selecciona un solo renglón en la tabla y da click en "send".
#
# ## Aquí verificamos que solo seleccionaste un solo renglón.
# stopifnot(nrow(proj_info_interactive) == 1)
# ## Crea el objeto RSE
# rse_gene_interactive <- create_rse(proj_info_interactive)


## ----"tranform_counts"----------------------------------------
## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)


## ----"expand_attributes"--------------------------------------
## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)
colData(rse_gene_SRP009615)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]

## ----"Ejercicio de clase"--------------------------------------

iSEE::iSEE(rse_gene_SRP009615)


## The following list of commands will generate the plots created in iSEE
## Copy them into a script or an R session containing your SingleCellExperiment.
## All commands below refer to your SingleCellExperiment object as `se`.

se <- rse_gene_SRP009615
colormap <- ExperimentColorMap()
se <- iSEE::cleanDataset(se)
colormap <- synchronizeAssays(colormap, se)
all_contents <- list()

################################################################################
# Defining brushes
################################################################################

all_active <- list()
all_active[['FeatureAssayPlot1']] <- list()
all_active[['ColumnDataPlot1']] <- list()
all_active[['RowDataPlot1']] <- list()
all_active[['SampleAssayPlot1']] <- list()

################################################################################
## Row data table 1
################################################################################


tab <- as.data.frame(rowData(se));
tab <- tab[,c("source", "type", "bp_length", "gene_id", "gene_type", "gene_name",
              "level", "havana_gene", "tag", "rowRanges_start", "rowRanges_end",
              "rowRanges_seqnames", "rowRanges_strand"),drop=FALSE]

# Saving data for transmission
all_contents[['RowDataTable1']] <- tab

################################################################################
## Column data plot 1
################################################################################


plot.data <- data.frame(Y=colData(se)[, "rail_id"], row.names=colnames(se));
plot.data$X <- factor(character(ncol(se)))

plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y,
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Avoid visual biases from default ordering by shuffling the points
set.seed(12);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot() +
  geom_violin(aes(x=X, y=Y, group=GroupBy), alpha=0.2, data=plot.data, scale='width', width=0.8) +
  geom_point(aes(y=Y, x=jitteredX), alpha=1, plot.data, color='#000000', size=1) +
  labs(x="", y="rail_id", title="rail_id ") +
  coord_cartesian(ylim=range(plot.data$Y, na.rm=TRUE), expand=TRUE) +
  scale_x_discrete(drop=FALSE) +
  theme_bw() +
  theme(legend.position='bottom', legend.text=element_text(size=9),
        legend.title=element_text(size=11), legend.box='vertical',
        axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['ColumnDataPlot1']] <- plot.data

################################################################################
## Row data plot 1
################################################################################


plot.data <- data.frame(Y=rowData(se)[, "source"], row.names=rownames(se));
plot.data$X <- factor(character(nrow(se)))

set.seed(100);
j.out <- iSEE:::jitterSquarePoints(plot.data$X, plot.data$Y);
summary.data <- j.out$summary;
plot.data$jitteredX <- j.out$X;
plot.data$jitteredY <- j.out$Y;

# Avoid visual biases from default ordering by shuffling the points
set.seed(63856);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot(plot.data) +
  geom_tile(aes(x=X, y=Y, height=2*YWidth, width=2*XWidth, group=interaction(X, Y)),
            summary.data, color='black', alpha=0, size=0.5) +
  geom_point(aes(x=jitteredX, y=jitteredY), alpha=1, plot.data, color='#000000', size=1) +
  labs(x="", y="source", title="source ") +
  scale_x_discrete(drop=FALSE) +
  scale_y_discrete(drop=FALSE) +
  theme_bw() +
  theme(legend.position='bottom', legend.text=element_text(size=9),
        legend.title=element_text(size=11), legend.box='vertical',
        axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['RowDataPlot1']] <- plot.data

################################################################################
## Sample assay plot 1
################################################################################


plot.data <- data.frame(Y=assay(se, "raw_counts")[,"SRR387777"], row.names=rownames(se));
plot.data$X <- factor(character(nrow(se)));

plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y,
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Avoid visual biases from default ordering by shuffling the points
set.seed(63856);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot() +
  geom_violin(aes(x=X, y=Y, group=GroupBy), alpha=0.2, data=plot.data, scale='width', width=0.8) +
  geom_point(aes(y=Y, x=jitteredX), alpha=1, plot.data, color='#000000', size=1) +
  labs(x="", y="SRR387777 (raw_counts)", title="SRR387777") +
  coord_cartesian(ylim=range(plot.data$Y, na.rm=TRUE), expand=TRUE) +
  scale_x_discrete(drop=FALSE) +
  theme_bw() +
  theme(legend.position='bottom', legend.text=element_text(size=9),
        legend.title=element_text(size=11), legend.box='vertical',
        axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['SampleAssayPlot1']] <- plot.data

################################################################################
## Column data table 1
################################################################################

################################################################################
## Complex heatmap 1
################################################################################


.chosen.rows <- "ENSG00000278704.1";
.chosen.columns <- colnames(se);
plot.data <- assay(se, "raw_counts")[.chosen.rows, .chosen.columns, drop=FALSE]
plot.data <- as.matrix(plot.data);

.assay_colors <- assayColorMap(colormap, "raw_counts", discrete=FALSE)(21L)
.assay_colors <- circlize::colorRamp2(breaks = seq(0, 1, length.out = 21L), colors = .assay_colors)

# Keep all samples to compute the full range of continuous annotations
.column_data <- colData(se)[, character(0), drop=FALSE]
.column_data[["Selected points"]] <- iSEE::multiSelectionToFactor(list(), colnames(se))

.column_col <- list()

.column_col[["Selected points"]] <- iSEE::columnSelectionColorMap(colormap, levels(.column_data[["Selected points"]]))

.column_data <- .column_data[colnames(plot.data), , drop=FALSE]
.column_data <- as.data.frame(.column_data, optional=TRUE)
.column_annot_order <- order(.column_data[["Selected points"]])
.column_data <- .column_data[.column_annot_order, , drop=FALSE]
plot.data <- plot.data[, .column_annot_order, drop=FALSE]
.column_annot <- ComplexHeatmap::columnAnnotation(df=.column_data, col=.column_col, annotation_legend_param=list(direction="horizontal", nrow=10))

hm <- ComplexHeatmap::Heatmap(matrix=plot.data, col=.assay_colors,
                              top_annotation=.column_annot, cluster_rows=FALSE, cluster_columns=FALSE,
                              name="raw_counts", show_row_names=TRUE, show_column_names=FALSE,
                              row_names_gp=grid::gpar(fontsize=10),
                              column_names_gp=grid::gpar(fontsize=10),
                              heatmap_legend_param=list(direction="horizontal"))

ComplexHeatmap::draw(hm, heatmap_legend_side="bottom", annotation_legend_side="bottom")

################################################################################
## Feature assay plot 1
################################################################################


plot.data <- data.frame(Y=assay(se, "counts")["ENSG00000168314.17", ], row.names=colnames(se))
plot.data$X <- colData(se)[, "sra_attribute.shRNA_expression"];

plot.data[["X"]] <- factor(plot.data[["X"]]);

plot.data$ColorBy <- colData(se)[, "sra_attribute.treatment"];
plot.data[["ColorBy"]] <- factor(plot.data[["ColorBy"]]);

plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y,
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Avoid visual biases from default ordering by shuffling the points
set.seed(12);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot() +
  geom_violin(aes(x=X, y=Y, group=GroupBy), alpha=0.2, data=plot.data, scale='width', width=0.8) +
  geom_point(aes(y=Y, color=ColorBy, x=jitteredX), alpha=1, plot.data, size=1) +
  labs(x="sra_attribute.shRNA_expression", y="ENSG00000168314.17 (counts)", color="sra_attribute.treatment", title="ENSG00000168314.17 vs sra_attribute.shRNA_expression") +
  coord_cartesian(ylim=range(plot.data$Y, na.rm=TRUE), expand=TRUE) +
  scale_color_manual(values=colDataColorMap(colormap, "sra_attribute.treatment", discrete=TRUE)(2), na.value='grey50', drop=FALSE) +
  scale_fill_manual(values=colDataColorMap(colormap, "sra_attribute.treatment", discrete=TRUE)(2), na.value='grey50', drop=FALSE) +
  scale_x_discrete(drop=FALSE) +
  theme_bw() +
  theme(legend.position='bottom', legend.text=element_text(size=9),
        legend.title=element_text(size=11), legend.box='vertical',
        axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['FeatureAssayPlot1']] <- plot.data

################################################################################
## To guarantee the reproducibility of your code, you should also
## record the output of sessionInfo()
sessionInfo()
