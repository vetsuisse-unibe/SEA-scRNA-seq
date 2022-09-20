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
library(scran)

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


# cell and gene sparsity 
sc22_sce<- addPerFeatureQC(sc22_sce)
rowData(sc22_sce)
colData(sc22_sce)$cell_sparsity <- 1 - (colData(sc22_sce)$detected / nrow(sc22_sce))
rowData(sc22_sce)$gene_sparsity <- (100 - rowData(sc22_sce)$detected) / 100

#the cell sparsity: for each cell, the proportion of genes that are not detected
#The cell sparsity plot shows that most cells have between 85% and 99% 0’s, which is typical.
hist(sc22_sce$cell_sparsity, breaks=50, col="grey80", xlab="Cell sparsity", main="")

#the gene sparsity: for each gene, the proportion of cells in which it is not detected
#The gene sparsity plot shows that a large number of genes are almost never detected, which is also regularly observed.
hist(rowData(sc22_sce)$gene_sparsity, breaks=50, col="grey80", xlab="Gene sparsity", main="")

lib.sf <- librarySizeFactors(sc22_sce)
summary(lib.sf)
dd <- data.frame("log10libSf"=log10(lib.sf))
ggplot(dd, aes(x=log10libSf)) + geom_histogram(bins=50)


# Normalization 
#Deconvolution (Ref:https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7)
set.seed(100) # clusters with PCA from irlba with approximation
clust <- quickCluster(sc22_sce) # slow with all cells.
table(clust)

#Compute size factors
sc22_sce <- computePooledFactors(sc22_sce,
                                 clusters = clust,
                                 min.mean = 0.1)
deconv.sf <- sizeFactors(sc22_sce)
summary(deconv.sf)

#Plot deconvolution size factors against library size factors:
sc22_sce <- addPerFeatureQC(sc22_sce) # PATCH

colData(sc22_sce)$cell_sparsity <- 1 - (colData(sc22_sce)$detected / nrow(sc22_sce))
rowData(sc22_sce)$gene_sparsity <- (100 - rowData(sc22_sce)$detected) / 100

deconvDf <- data.frame(lib.sf, deconv.sf,
                       "source_name" = sc22_sce$disease_state,
                       "sum" = sc22_sce$sum,
                       "mito_content" = sc22_sce$subsets_mito_genes_percent,
                       "cell_sparsity" = sc22_sce$cell_sparsity)

sp <- ggplot(deconvDf, aes(x=lib.sf, y=deconv.sf, col=source_name)) +
  geom_point()

sp

sp <- ggplot(deconvDf, aes(x=lib.sf, y=deconv.sf, col=cell_sparsity)) +
  geom_point()
sp

# log normalization using the size factors
#For each cell, raw counts for genes are divided by the size factor for that cell and 
#log-transformed so downstream analyses focus on genes with strong relative differences.

sc22_sce <- logNormCounts(sc22_sce)
print(assayNames(sc22_sce))


#ScTransform 
print(class(counts))

print(dim(counts))
## name columns (cells) with barcodes
colnames(counts) <- colData(sc22_sce)$Barcode

# gene attributes:
# prepare a data frame named e.g. 'gene_attr' to keep gene attributes, inc:
gene_attr <- data.frame(mean = rowMeans(counts), 
                        detection_rate = rowMeans(counts > 0),
                        var = rowVars(counts))
gene_attr$log_mean <- log10(gene_attr$mean)
gene_attr$log_var <- log10(gene_attr$var)
# name rows of the 'gene_attr' data frame:
rownames(gene_attr) <- rownames(counts)

# cell attributes:
cell_attr <- data.frame(n_umi = colSums(counts),
                        n_gene = colSums(counts > 0))
rownames(cell_attr) <- colnames(counts)

#Gene attributes 
dim(gene_attr)
head(gene_attr)

#Cell attributes
dim(cell_attr)
head(cell_attr)


#Mean-variance relationship
#For the genes, on the log10 scale we can see that up to a mean UMI count of 0 the variance 
#follows the line through the origin with slop one, i.e. variance and mean are roughly equal as 
#expected under a Poisson model. However, genes with a higher average UMI count show 
#overdispersion compared to Poisson.
ggplot(gene_attr, aes(log_mean, log_var)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_density_2d(size = 0.3) +
  geom_abline(intercept = 0, slope = 1, color='red')

# add the expected detection rate under Poisson model
x = seq(from = -3, to = 2, length.out = 1000)
poisson_model <- data.frame(log_mean = x,
                            detection_rate = 1 - dpois(0, lambda = 10^x))
ggplot(gene_attr, aes(log_mean, detection_rate)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_line(data=poisson_model, color='red') +
  theme_gray(base_size = 8)

#The plot below show the relationship between the to cell attributes computed: library size (n_umi) and number of genes detected (n_gene)
ggplot(cell_attr, aes(n_umi, n_gene)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_density_2d(size = 0.3)

future::plan(strategy = 'multicore', workers = 7)
options(future.globals.maxSize = 10 * 1024 ^ 3)


set.seed(44) # for reproducibility
sc22_sparse <- as(counts(sc22_sce), "dgCMatrix")
vst_out <- sctransform::vst(
  sc22_sparse, # A matrix of UMI counts with genes as rows and cells as columns
  latent_var = c('log_umi'), # The independent variables to regress out as a character vector
  return_gene_attr = TRUE, # Make cell attributes part of the output
  return_cell_attr = TRUE, # Calculate gene attributes and make part of output
  verbosity = 0 # An integer specifying what to show (0: nothing, 1: messages, 2: + progress bar)
)

sctransform::plot_model_pars(
  vst_out, # The output of a vst run
  verbosity = 1 # Messages only, no progress bar
)

rowData(sc22_sce) %>%
  as.data.frame %>%
  dplyr::filter(X %in% c('CD3E', 'SNX10'))

geneId <- rowData(sc22_sce) %>%
  as.data.frame %>%
  dplyr::filter(X %in% c('CD3E', 'SNX10')) %>%
  pull("X")

sctransform::plot_model(
  vst_out, # The output of a vst run
  counts, # UMI count matrix
  geneId, # Vector of genes to plot
  plot_residual = TRUE
)

#The distribution of residual mean is cetered around 0:
ggplot(vst_out$gene_attr, aes(residual_mean)) +
  geom_histogram(binwidth=0.01)

#The distribution of residual variance is centered around 1:
ggplot(vst_out$gene_attr, aes(residual_variance)) +
  geom_histogram(binwidth=0.1) +
  geom_vline(xintercept=1, color='red') +
  xlim(0, 10)
#The following plot of the residual variance against the mean: after transformation 
#there is no relationship between gene mean and variance.

ggplot(vst_out$gene_attr,
       aes(log10(gmean), residual_variance)) +
  geom_point(alpha=0.3, shape=16) +
  geom_density_2d(size = 0.3)

#Check genes with large residual variance. 
#These genes would be markers of expected cell populations. 
#Note how they represent a great range of mean UMI and detection rate values.

dd <- vst_out$gene_attr %>%
  arrange(-residual_variance) %>%
  slice_head(n = 10) %>%
  mutate(across(where(is.numeric), round, 2))
symbol<-rowData(sc22_sce)[,c("X")]
rowData(sc22_sce)$ID <-symbol
dd %>% tibble::rownames_to_column("ID") %>%
  left_join(as.data.frame(rowData(sc22_sce))[,c("ID", "X")], "ID") %>%
  DT::datatable(rownames = FALSE)


# genes that are expressed in fewer than 5 cells are not used and not returned
# so to add vst_out$y as an assay we need to ditch the missing genes completely.
# https://github.com/ChristophH/sctransform/issues/27
geneOverlap <- rownames(sc22_sce) %in% rownames(vst_out$y)
if(!all(geneOverlap))
{
  table(rownames(sc22_sce) %in% rownames(vst_out$y))
  tmpInd <- which(rownames(sc22_sce) %in% rownames(vst_out$y))
  sc22_sce <- sc22_sce[tmpInd,]
  assayNames(sc22_sce)
}
