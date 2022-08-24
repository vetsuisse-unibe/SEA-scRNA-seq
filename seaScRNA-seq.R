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
sea.sce <- scDblFinder(sea.sce, samples="sample_id", clusters=NULL, BPPARAM=MulticoreParam(10))
## Convert sce object back to seurat
sea_seurat <- as.Seurat(sea.sce, counts = "counts", data = "logcounts")

head(sea_seurat@meta.data) %>%
  kable(format = "html", col.names = colnames(head(sea_seurat@meta.data))) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")

# Count the number of potential doublets found
table(sea_seurat$scDblFinder.class)
prop.table(table(sea_seurat$scDblFinder.class))
#11% doublets found with scDblFinder
