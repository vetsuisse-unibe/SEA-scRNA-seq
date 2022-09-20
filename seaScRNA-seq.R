library(Seurat)
library(clustree)
library(celldex)
library(SingleR)
library(cowplot)
library(patchwork)
library(dplyr)
library(Matrix)
library(gdata)
library(kableExtra)
library(ggplot2)
require("biomaRt")
library(RColorBrewer)
library(stringr)
library(viridis)
library(BiocParallel)
library(scDblFinder)
library(scater)
library(SingleCellExperiment)
library(magrittr)
library(scFeatureFilter)
library(glmGamPoi)

sea_sc2021 <- readRDS("initial_seurat_object.rds")

# Quality control and filtering 
# Add mitochondrial reads to the Seurat object's meta.data
mito.genes<-c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB","MT-tRNA-Phe",  
              "MT-s-rRNA","MT-tRNA-Val","MT-l-rRNA","MT-tRNA-Leu","MT-tRNA-Ile","MT-tRNA-Gln","MT-tRNA-Met","MT-tRNA-Trp","MT-tRNA-Ala",
              "MT-tRNA-Cys","MT-tRNA-Tyr","MT-tRNA-Ser","MT-tRNA-Asp","MT-tRNA-Lys","MT-tRNA-Gly","MT-tRNA-Arg","MT-tRNA-His","MT-tRNA-Ser.1",
              "MT-tRNA-Leu.1","MT-tRNA-Glu","MT-tRNA-Thr","MT-tRNA-Pro")#36/37 MT genes because tRNA-Asn is absent from this dataset

sea_sc2021[["percent.mito"]] <- PercentageFeatureSet(sea_sc2021, features=mito.genes)
summary(sea_sc2021$percent.mito)

sea_sc2021<- PercentageFeatureSet(sea_sc2021, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")

# Initial gene metrics 
kable(do.call("cbind", tapply(sea_sc2021$nFeature_RNA,
                              Idents(sea_sc2021),quantile,probs=seq(0,1,0.05))),
      caption = "5% Quantiles of Genes/Cell by Sample") %>% kable_styling()

kable(round(do.call("cbind", tapply(sea_sc2021$percent.mito, Idents(sea_sc2021),quantile,probs=seq(0,1,0.05))), digits = 3),
      caption = "5% Quantiles of Percent Mitochondria by Sample") %>% kable_styling()

# Sample wise Quality metrics 
VlnPlot(
  sea_sc2021,
  features = c("nFeature_RNA", "nCount_RNA","percent.mito"),
  ncol = 1, pt.size = 0.3)

RidgePlot(sea_sc2021, features=c("nFeature_RNA","nCount_RNA", "percent.mito"), ncol = 2)
plot(sort(Matrix::rowSums(GetAssayData(sea_sc2021) >= 3), decreasing = TRUE) , xlab="gene rank", ylab="number of cells", main="Cells per genes (reads/gene >= 3 )")

FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "percent.mito", shuffle = TRUE,cols=brewer.pal(n = 11, name = "Paired")) + geom_vline(xintercept = c(1000,12000)) + geom_hline(yintercept = 8)
FeatureScatter(sea_sc2021, feature1 = "nFeature_RNA", feature2 = "percent.mito", shuffle = TRUE,cols=brewer.pal(n = 11, name = "Paired")) + geom_vline(xintercept = 700) + geom_hline(yintercept = 8)
FeatureScatter(sea_sc2021, "nCount_RNA", "nFeature_RNA",pt.size = 0.5, shuffle = TRUE,cols=brewer.pal(n = 11, name = "Paired"))  + geom_vline(xintercept = c(1000,12000)) + geom_hline(yintercept = 700)

# Trying to find the best thresholds to filter for nFeature_RNA
(VlnPlot(sea_sc2021, features = "nFeature_RNA", pt.size = 0) + 
    geom_hline(yintercept=200, linetype="dashed", color = "red") +
    geom_hline(yintercept=8000, linetype="dashed", color = "red") +
    labs(title= NULL, x=NULL, y="Feature count")) &
  (scale_fill_discrete(limits = c("30", "33", "44","46","51","52","56","68","89","90","91"),
                       labels = c("Horse 30", "Horse 33", "Horse 44","Horse 46","Horse 51","Horse 52","Horse 56","Horse 68","Horse 89","Horse 90","Horse 91")))
FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept=200, linetype="dashed", color = "red") +
  geom_hline(yintercept=8000, linetype="dashed", color = "red") + 
  labs(title= NULL, x="RNA count", y="Feature count")

#200 and 8000 look good
#Creating nice figure to justify the choice of filters
Plot_FilterVln <- ((VlnPlot(sea_sc2021, features = "nFeature_RNA", pt.size = 0) + 
                      geom_hline(yintercept=200, linetype="dashed", color = "red") +
                      geom_hline(yintercept=8000, linetype="dashed", color = "red") +
                      labs(title= NULL, x=NULL, y="Feature count")) &
                     (VlnPlot(sea_sc2021, features = "percent.mito", pt.size = 0) + 
                        geom_hline(yintercept=15, linetype="dashed", color = "red") +
                        labs(title= NULL, x=NULL, y="Mitochondrial reads %")) +
                     plot_layout(ncol = 1)) &
  (theme(axis.title = element_text(size=10), axis.text.x = element_blank(), 
         axis.ticks.x = element_blank(), axis.text.y = element_text(size=8))) &
  scale_fill_discrete(limits = c("30", "33", "44","46","51","52","56","68","89","90","91"),
                      labels = c("Horse 30", "Horse 33", "Horse 44","Horse 46","Horse 51","Horse 52","Horse 56","Horse 68","Horse 89","Horse 90","Horse 91"))

Plot_FilterFeat <- ((FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
                       geom_hline(yintercept=200, linetype="dashed", color = "red") +
                       geom_hline(yintercept=8000, linetype="dashed", color = "red") + 
                       labs(title= NULL, x="RNA count", y="Feature count")) +
                      (FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "percent.mito") + 
                         geom_hline(yintercept=15, linetype="dashed", color = "red") 
                       + labs(title= NULL, x="RNA count", y="Mitochondrial reads %")) +
                      plot_layout(ncol = 1)) &
  theme(axis.title = element_text(size=10), axis.text = element_text(size=8)) &
  NoLegend()


(Plot_FilterFeat | Plot_FilterVln) + plot_annotation(tag_levels = 'A')

# QC metrics grouped by disease state
VlnPlot(sea_sc2021, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "disease_state",
        pt.size = 0.001, combine = T) & theme(legend.position = 'none',
                                              axis.title.x = element_blank())
VlnPlot(sea_sc2021, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "disease_state",
        pt.size = 0, combine = T) & theme(legend.position = 'none', 
                                          axis.title.x = element_blank())

FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "disease_state")
FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "disease_state")
FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "percent.ribo", group.by = "disease_state")
FeatureScatter(sea_sc2021, feature1 = "nFeature_RNA", feature2 = "percent.ribo", group.by = "disease_state")

# Find Doublets 
# find Doublets 
### Convert object into singlecellexperiment
sea.sce <- as.SingleCellExperiment(sea_sc2021)
#sea.sce <- scDblFinder(sea.sce, samples="sample_id", clusters=NULL, BPPARAM=MulticoreParam(10))
sea.sce <- scDblFinder(sea.sce, samples="sample_id", clusters=NULL)
## Convert sce object back to seurat
sea_seurat <- as.Seurat(sea.sce, counts = "counts", data = "logcounts")

head(sea_seurat@meta.data) %>%
  kable(format = "html", col.names = colnames(head(sea_seurat@meta.data))) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")

# Count the number of potential doublets found
table(sea_seurat$scDblFinder.class)
prop.table(table(sea_seurat$scDblFinder.class))
#12% doublets found with scDblFinder

# Remove all the potential doublets found 
sea_singlet <- subset(sea_seurat, subset = scDblFinder.class  == "singlet")
dim(sea_seurat)
## 32,835 genes and 85,961 cells.
dim(sea_singlet)
## 32,835 genes and 75,616 cells (10,345 cells = 12% filtered)

#Subset Data
#removing potential doublets and subsetting based on nFeature and percent.mt
sea_singlet_flt <- subset(sea_singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 9000
                          & percent.mito < 15 )
dim(sea_singlet_flt)
## 32,835 genes and 60,415 cells (15,201 additional cells filtered)
## 15,201/75,616 = 20% of the singlets have been filtered
## 15,201/85,961 = 17.6% of all cells have been additionally filtered

## In total, 1 - (60,415/85,961) = 30% cells filtered.

# Update Metadata with MT and Rb genes information
sea_singlet_flt[["percent.mito"]] <- PercentageFeatureSet(sea_singlet_flt, features=mito.genes)
summary(sea_singlet_flt$percent.mito)

sea_singlet_flt<- PercentageFeatureSet(sea_singlet_flt, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")
summary(sea_singlet_flt$percent.ribo)

VlnPlot(sea_singlet_flt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"), ncol = 4, group.by = "disease_state",
        pt.size = 0.001, combine = T) & theme(legend.position = 'none',
                                              axis.title.x = element_blank())

VlnPlot(sea_singlet_flt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"), ncol = 4, group.by = "disease_state",
        pt.size = 0, combine = T) & theme(legend.position = 'none', 
                                          axis.title.x = element_blank())

FeatureScatter(sea_singlet_flt, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "disease_state")
FeatureScatter(sea_singlet_flt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "disease_state") 

# QUALITY CONTROL USING SCATER

sc22 <- sea_singlet_flt

# Calculate a cell cycle score for each cell
convertHumanGeneList <- function(x){
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  horse = useEnsembl("ensembl", dataset = "ecaballus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("external_gene_name"), martL = horse, uniqueRows=T)
  humanx <- unique(genes[, 2])
  print(head(humanx))
  return(humanx)
}
m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

sc22 <- CellCycleScoring(sc22, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
head(sc22)
table(sc22@meta.data$Phase) %>% kable(caption = "Number of Cells in each Cell Cycle Stage", col.names = c("Stage", "Count"), align = "c") %>% kable_styling()
table(Idents(sc22))
Idents(sc22) <- "orig.ident"

# converting Seurat object to a SCE object
cts <- Seurat::GetAssayData(sc22, slot = "counts")

sc22_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = cts),
  colData = sc22@meta.data,
  rowData = rownames(sc22)
)
class(sc22)
class(sc22_sce)
sc22_sce

# add the Mt and RP genes info the SCE object
ribo_genes <- rownames(sc22)[grep(pattern = "^RP[S|L]", rownames(sc22), perl = T)] #should not add before conversion to SCE?
sc22_sce <- scuttle::addPerCellQC(sc22_sce,
                                  subsets=list(mito_genes=which(rownames(sc22_sce) %in% mito.genes),
 
                                                                                             ribo_genes=which(rownames(sc22_sce) %in% ribo_genes)))

SingleCellExperiment::colData(sc22_sce)

scater::plotColData(sc22_sce, x = "sum", y="detected")
scater::plotColData(sc22_sce, x = "detected", y="subsets_mito_genes_percent")
scater::plotColData(sc22_sce, x = "subsets_mito_genes_percent", y="subsets_ribo_genes_percent")
detected_genes <- rowSums(counts(sc22_sce)) > 0
table(detected_genes)
sc22_sce <- sc22_sce[detected_genes,]

plotColData(sc22_sce, x="orig.ident", y="sum",other_fields="disease_state") + 
  facet_wrap(~disease_state, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  ggtitle("Total count")

plotColData(sc22_sce, x="orig.ident", y="detected", other_fields="disease_state") + 
  facet_wrap(~disease_state, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  ggtitle("Detected features")

plotColData(sc22_sce, x="orig.ident", y="subsets_mito_genes_percent", other_fields="disease_state") + 
  facet_wrap(~disease_state, nrow=1, scales = "free_x") +
  ggtitle("Mito percent")

colData(sc22_sce) %>% 
  as.data.frame() %>% 
  arrange(subsets_mito_genes_percent) %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(colour = subsets_mito_genes_percent > 10)) + 
  facet_wrap(vars(disease_state))

#Highly expressed genes
# On the gene level, we can look at a plot that shows the top (by default 50) 
# most-expressed genes. Each row in the plot corresponds to a gene; 
# each bar corresponds to the expression of a gene in a single cell; 
# the circle indicates the median expression of each gene, with which genes are sorted. 
# We expect to see the "usual suspects", i.e., mitochondrial genes, ribosomal protein)
scater::plotHighestExprs(sc22_sce, exprs_values = "counts", colour_cells_by="detected")

# Normalization, so that the statistics reflect changes in relative expression.
sc22_sce <- scater::logNormCounts(sc22_sce)  # alternative to Seurat's normalization here using scater

#PCA plot 
# Before normalization 
#Without log-transformation or normalization, PCA plot fails to separate the datasets by replicate or 
#individual. We mostly see the effects of sequencing depth - samples (cells) with lots of expression, 
#and particularly highly expressed genes, dominate the PCs:
sc22_sce <- runPCA(sc22_sce, exprs_values = "counts")
dim(reducedDim(sc22_sce, "PCA"))
plotPCA(sc22_sce, colour_by = "disease_id", size_by = "detected", shape_by = "disease_state")

# After normalization 
#With log-transformation, we equalize the large difference between strongly and weakly expressed genes
# and  see cells still form groups by sequencing depth and not by individual or disease_state.
#log-transformation it reduces the variance on the first principal component. 
 
sc22_sce <- runPCA(sc22_sce, exprs_values = "logcounts")
dim(reducedDim(sc22_sce, "PCA"))
plotPCA(sc22_sce, colour_by = "disease_id", size_by = "detected", shape_by = "disease_state")
plotPCA(sc22_sce, colour_by = "disease_id", size_by = "sum", shape_by = "disease_state")

#TSNE 
# Again the plot shows that groups are formed by sequencing depth. Perlexity upto 250 still shows the same
set.seed(123456)
sc22_sce <- runTSNE(sc22_sce, exprs_values = "logcounts", perplexity = 130)
plotTSNE(sc22_sce, colour_by = "disease_id", size_by = "detected", shape_by = "disease_state")

#
getExplanatoryPCs(sc22_sce,variables = "sum")
plotExplanatoryPCs(sc22_sce,variables = "sum")
plotExplanatoryPCs(sc22_sce,variables = "breed")
plotExplanatoryPCs(sc22_sce,variables = "subsets_mito_genes_percent")
plotExplanatoryPCs(sc22_sce,variables = "subsets_ribo_genes_percent")


# look for explanatory variables
vars <- scater::getVarianceExplained(sc22_sce,
                                     variables = c("sum", "detected","disease_state","orig.ident", "Phase",
                                                   "subsets_mito_genes_percent","subsets_ribo_genes_percent"))
head(vars)
scater::plotExplanatoryVariables(vars)


# scTransform on sce object 
sc22_sparse <- as(counts(sc22_sce), "dgCMatrix")
sctnorm_data <- sctransform::vst(umi = sc22_sparse, min_cells = 1,
                                 cell_attr = as.data.frame(colData(sc22_sce)))
#Calculating cell attributes from input UMI matrix: log_umi
#Variance stabilizing transformation of count matrix of size 24665 by 60415
#Model formula is y ~ log_umi
#Get Negative Binomial regression parameters per gene
#Using 2000 genes, 60415 cells
#|===============================================================================================| 100%
#There are 3 estimated thetas smaller than 1e-07 - will be set to 1e-07
#Found 88 outliers - those will be ignored in fitting/regularization step

#Second step: Get residuals using fitted parameters for 24665 genes
#|===============================================================================================| 100%
#Calculating gene attributes
#Wall clock passed: Time difference of 16.33139 mins

dim(sctnorm_data$y)
dim(sc22_sce)
sctnorm_data$model_str
sctransform::plot_model_pars(sctnorm_data)

# Variance against mean 

ggplot(sctnorm_data$gene_attr,
       aes(log10(gmean), residual_variance)) +
  geom_point(alpha=0.3, shape=16) +
  geom_density_2d(size = 0.3)

# genes with large variance 
dd <- sctnorm_data$gene_attr %>%
  arrange(-residual_variance) %>%
  slice_head(n = 22) %>%
  mutate(across(where(is.numeric), round, 2))

dd %>% tibble::rownames_to_column("ID") %>%
  left_join(as.data.frame(rowData(sc22_sce))[,c("X")], "ID") %>%
  DT::datatable(rownames = FALSE)

# scTransform on seurat object 
sea_singlet_flt_sct<- SCTransform(sea_singlet_flt, method = "glmGamPoi", variable.features.n = 2000)
DefaultAssay(sea_singlet_flt_sct) 
#Calculating cell attributes from input UMI matrix: log_umi
#Variance stabilizing transformation of count matrix of size 21436 by 60415
#Model formula is y ~ log_umi
#Get Negative Binomial regression parameters per gene
#Using 2000 genes, 5000 cells
#|===============================================================================================| 100%
#Found 148 outliers - those will be ignored in fitting/regularization step

#Second step: Get residuals using fitted parameters for 21436 genes
#|===============================================================================================| 100%
#Computing corrected count matrix for 21436 genes
#|===============================================================================================| 100%
#Calculating gene attributes
#Wall clock passed: Time difference of 7.411177 mins
#Determine variable features
#Place corrected count matrix in counts slot
#Centering data matrix
#|===============================================================================================| 100%
#Set default assay to SCT

#dimension reduction
sea_singlet_flt_sct <- RunPCA(sea_singlet_flt_sct, verbose = FALSE)
sea_singlet_flt_sct <- RunUMAP(sea_singlet_flt_sct, dims = 1:30, verbose = FALSE)

DimPlot(pbmc, label = TRUE) + NoLegend()
