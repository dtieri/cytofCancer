
library(readxl)
library(flowWorkspace)
library(CytoML)
library(SingleCellExperiment)
library(diffcyt)
library(lme4)
library(CATALYST)
#library("cytofWorkflow")



##Import metadata and marker pannel inf

md <- read_excel("samMetadata.xlsx")
#head(data.frame(md))

panel <- read_excel("David_B cell Panel_CyTOF.xlsx")
#head(data.frame(panel))



##import flowjo file to gating set

wsfile <- list.files("/home/david/Projects/cytofCancer/flojoFiles", pattern="Full patient compiled.wsp", full = TRUE)

ws <- open_flowjo_xml(wsfile)
#flowjo_ws_close(ws)
#ws
#fj_ws_get_samples(ws)

gs <- flowjo_to_gatingset(ws,name=1, path = "data"); #import the first group

#gs_get_pop_paths(gs, path = "full")
#plot(gs)

#plotGate(gs, "B cells")
#plotGate(gs[21], "B cells")

#gs_pop_get_count_fast(gs)

#fs <- gs_pop_get_data(gs,y = "CD45", inverse.transform = FALSE)
#fs <- gs_pop_get_data(gs,y = "CD45/B cells", inverse.transform = FALSE)
fs <- gs_pop_get_data(gs,y = "CD45/T cells", inverse.transform = FALSE)
#fs <- gs_pop_get_data(gs,y = "CD45+/T cells/CD4", inverse.transform = FALSE)
#fs <- gs_pop_get_data(gs,y = "CD45+/T cells/CD8", inverse.transform = FALSE)
#fs <- gs_pop_get_data(gs,y = "CD45+/T cells/gd T cells", inverse.transform = FALSE)
#fs



#Import gated flowset into a SingleCellExperiment so that we can make the plots from Cytof workflow

# specify levels for conditions & sample IDs to assure desired ordering
md$condition <- factor(md$condition, levels = c("Normal", "Tumor", "PBMC"))
md$sample_id <- factor(md$sample_id, levels = md$sample_id[order(md$condition)])

#md

# construct SingleCellExperiment
sce <- prepData(fs, panel,md, md_cols = list(file="file_name",id = "sample_id", factors = c("condition", "patient_id")), features = panel$fcs_colname)
#sce



#Diagnostic plots

plotFolder<-"plotsBCells"
dir.create(plotFolder)

p <- plotExprs(sce, color_by = "condition")
pdf(paste0(plotFolder,"/density.pdf"),width=10,height=10)
p
dev.off()

rm(p)

q=plotCounts(sce, color_by = "condition")
q

pdf(paste0(plotFolder,"/cellCounts.pdf"),width=9,height=9)
q
dev.off()


sce2<-sce[,which(sce$condition=="Tumor")]
sce2$sample_id<-droplevels(sce2$sample_id)
sce2$condition<-droplevels(sce2$condition)
sce2$patient_id<-droplevels(sce2$patient_id)
metadata(sce2)$experiment_info<-metadata(sce2)$experiment_info[which(metadata(sce2)$experiment_info$condition=="Tumor"),]

r<-plotMDS(sce2, color_by = "patient_id")

pdf(paste0(plotFolder,"/mds2.pdf"),width=9,height=9)
r
dev.off()

s<-plotExprHeatmap(sce, bin_anno = TRUE, row_anno = TRUE)
pdf(paste0(plotFolder,"/heatmap.pdf"),width=10,height=10)
s
dev.off()

s2<-plotExprHeatmap(sce2, bin_anno = TRUE, row_anno = TRUE)
pdf(paste0(plotFolder,"/heatmap2.pdf"),width=10,height=10)
s2
dev.off()


dim(sce[,which(sce$condition=="Tumor")])

dim(sce)

t<-plotNRS(sce, features = type_markers(sce), color_by = "condition")
t
pdf(paste0(plotFolder,"/nrs.pdf"),width=14,height=14)
t
dev.off()




#Clustering for cell type identification
library(ConsensusClusterPlus)

kclustering<-"meta10"

set.seed(1234)
sce <- cluster(sce, features = type_markers(sce),
    xdim = 10, ydim = 10, maxK = 20, seed = 1234)

sce2c<-cluster(sce2, features = type_markers(sce2),
    xdim = 10, ydim = 10, maxK = 20, seed = 1234)

type_markers(sce2)

metadata(sce2c)

pc<-plotClusterHeatmap(sce2c,
    hm2 = NULL, k = kclustering, m = NULL,
    cluster_anno = TRUE, draw_freqs = TRUE)
pc
pdf(paste0(plotFolder,"/lineageCluster.pdf"),width=10,height=10)
pc
dev.off()

sessionInfo()



#qc<-plotClusterExprs(sce, k = kclustering, features = "type")
#qc
#pdf(paste0(plotFolder,"/lineageIntensities.pdf"),width=14,height=14)
#qc
#dev.off()


#rc<-plotClusterHeatmap(sce, hm2 = "CD19", k = kclustering, draw_freqs = TRUE)
#rc
#pdf(paste0(plotFolder,"/functionalCluster.pdf"),width=10,height=10)
#rc
#dev.off()




#UMAP (t-sne)
library(cowplot)
library(ggplot2)
set.seed(1234)

sce <- runDR(sce, dr = "TSNE", cells = 500, features = "type")
sce <- runDR(sce, dr = "UMAP", cells = 1e3, features = "type")

p0<-plotDR(sce, "UMAP", color_by = "CD4")
p0
pdf(paste0(plotFolder,"/umap.pdf"),width=8,height=8)
p0
dev.off()

p1 <- plotDR(sce, "TSNE", color_by = "meta20") +
    theme(legend.position = "none")
p2 <- plotDR(sce, "UMAP", color_by = "meta20")
lgd <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")

p3<-plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))
p3
pdf(paste0(plotFolder,"/tsneANDumap.pdf"),width=8,height=8)
p3
dev.off()

p4<-plotDR(sce, "UMAP", color_by = "meta20") + facet_wrap("sample_id") +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
p4
pdf(paste0(plotFolder,"/umapSample.pdf"),width=8,height=8)
p4
dev.off()

p5<-plotDR(sce, "UMAP", color_by = "meta20") + facet_wrap("condition") +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
p5
pdf(paste0(plotFolder,"/umapCondition.pdf"),width=8,height=8)
p5
dev.off()

p6<-plotCodes(sce, k = "meta20")
p6
pdf(paste0(plotFolder,"/tsneANDpcaCodes.pdf"),width=8,height=8)
p6
dev.off()

p7<-plotClusterHeatmap(sce,
    hm2 = "CD19", k = "som100", m = "meta20",
    cluster_anno = FALSE, draw_freqs = TRUE)
p7
pdf(paste0(plotFolder,"/clusterHeatmap.pdf"),width=14,height=14)
p7
dev.off()



#Clustering only Tumor Samples
library(ConsensusClusterPlus)

maxK<-20
kclustering<-paste0("meta",maxK)

sce2<-sce[,which(sce$condition=="Tumor")]
sce2<-sce2[,-which(sce2$sample_id=="962_Tumor")]
sce2$sample_id<-droplevels(sce2$sample_id)
sce2$condition<-droplevels(sce2$condition)
sce2$patient_id<-droplevels(sce2$patient_id)
metadata(sce2)$experiment_info<-metadata(sce2)$experiment_info[which(metadata(sce2)$experiment_info$condition=="Tumor"),]
metadata(sce2)$experiment_info<-metadata(sce2)$experiment_info[-which(metadata(sce2)$experiment_info$sample_id=="962_Tumor"),]


n_cells(sce2)



set.seed(1234)
sce2 <- cluster(sce2, features = type_markers(sce2),
    xdim = 10, ydim = 10, maxK = maxK, seed = 1234)

pc<-plotClusterHeatmap(sce2,
    hm2 = NULL, k = "meta15", m = NULL,
    cluster_anno = TRUE, draw_freqs = TRUE)
pc
pdf(paste0(plotFolder,"/lineageCluster2.pdf"),width=10,height=10)
pc
dev.off()

#qc<-plotClusterExprs(sce, k = kclustering, features = "type")
#qc
#pdf(paste0(plotFolder,"/lineageIntensities.pdf"),width=14,height=14)
#qc
#dev.off()


#rc<-plotClusterHeatmap(sce, hm2 = "CD19", k = kclustering, draw_freqs = TRUE)
#rc
#pdf(paste0(plotFolder,"/functionalCluster.pdf"),width=10,height=10)
#rc
#dev.off()





#UMAP (t-sne) for just tumor samples
library(cowplot)
library(ggplot2)
set.seed(1234)

sce2 <- runDR(sce2, dr = "TSNE", cells = 500, features = "type")
sce2 <- runDR(sce2, dr = "UMAP", cells = 1e4, features = "type")

sce2

which(is.na(reducedDim(sce2)[,1]))

length(which(is.na(reducedDim(sce2)[,1])))

length(reducedDim(sce2)[,2])

p4<-plotDR(sce2, "UMAP", color_by = "meta10") + geom_point(position=position_jitter(h=0.1, w=0.1),
             shape = 21, alpha = 0.5, size = 0.1)  + facet_wrap("sample_id") +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
p4
pdf(paste0(plotFolder,"/umapSample2b.pdf"),width=8,height=8)
p4
dev.off()

dim(sce2)

p4<-plotDR(sce2, "TSNE", color_by = kclustering) + facet_wrap("sample_id") +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
p4
pdf(paste0(plotFolder,"/TSNESample3.pdf"),width=8,height=8)
p4
dev.off()

plotAbundances(sce2, k = "meta10", by = "cluster_id", shape = "patient_id")

metadata(sce2)




#correlations

#bAbundances<-plotAbundances(sce2, k = "meta15", by = "sample_id")
#bAbundances<-bAbundances$data[,c(1,2,3)]
#bAbundances$cluster_id<-paste0(bAbundances$cluster_id,"_B")

tAbundances<-plotAbundances(sce2, k = "meta15", by = "sample_id")
tAbundances<-tAbundances$data[,c(1,2,3)]
tAbundances$cluster_id<-paste0(tAbundances$cluster_id,"_T")


abundances<-rbind(bAbundances,tAbundances)

abundances

abundancePerCluster<-list()
for (ii in unique(abundances$sample_id)) {
    temp<-subset(abundances,sample_id==ii)
    abundancePerCluster[[ii]]<-temp$freq
}


abundancedf<-data.frame(do.call("cbind",abundancePerCluster))
rownames(abundancedf)<-temp$cluster_id

abundancedf

write.table(abundancedf,"abundancedf.csv",sep="\t",row.names = TRUE,col.names=NA)

abundancedf

as.matrix(t(abundancedf))

library(corrplot)
rst = cor(t(abundancedf))
rs = cor(abundancedf)
pdf("BT_Correlations_Tumor.pdf")
corrplot(rs, type = "upper", order = "alphabet", tl.col = "black", tl.srt = 45, tl.cex=.5)
dev.off()

library(gplots)
library(factoextra)
pdf("Heatmap.LN.Marker.Correlations_Cancer.Only.pdf")
heatmap.2(log2(as.matrix(t(abundancedf))+0.2), col=redgreen(75), trace="none", margins=c(10,10))
dev.off()

abundancedf

plotAbundances(sce2, k = "meta10", by = "cluster_id", shape = "patient_id")

### FINAL PCA for paper
###------------### All cells
### PCA, biplots, etc.

res.pca <- prcomp(t(abundancedf), scale = TRUE)
fviz_pca_ind(res.pca)

pdf(paste0("plots/",folder,"/",ISOTYPE,"_","All.Markers_Scree.pdf"))
#pdf(paste0("plots/",ISOTYPE,"/pcaJ/All.Markers_Scree.pdf"))
fviz_eig(res.pca)
dev.off()



fviz_pca_ind(res.pca)

pdf(paste0("plots/",folder,"/",ISOTYPE,"_","All.Markers_PC1.2_without.labels.ellipse.pdf"))
#pdf(paste0("plots/",ISOTYPE,"/pcaJ/All.Markers_PC1.2_without.labels.ellipse.pdf"))
fviz_pca_ind(res.pca, geom = "point",
           col.ind = rownames(t(abundancedf)),
           axes = c(1, 2),
            pointsize = 4,
            invisible = "quali",
            #repel = TRUE     # Avoid text overlapping
            ggtheme=theme(axis.text=element_text(size=30), axis.title=element_text(size=25, face="bold")))
dev.off()




#Merging Clusters

#DO NOT USE NOW!!!!!!!! 5.Dec.2019

merging_table1 <- "PBMC8_cluster_merging1.xlsx"
download.file(file.path(url, merging_table1),
    destfile = merging_table1, mode = "wb")
merging_table1 <- read_excel(merging_table1)
data.frame(merging_table1)

# convert to factor with merged clusters in desired order
merging_table1$new_cluster <- factor(merging_table1$new_cluster,
    levels = c("B-cells IgM+", "B-cells IgM-", "CD4 T-cells",
        "CD8 T-cells", "DC", "NK cells", "monocytes", "surface-"))

sce <- mergeClusters(sce, k = kclustering,
    table = merging_table1, id = "merging1")


mp1<-plotDR(sce, "UMAP", color_by = "merging1")
mp1
pdf(paste0(plotFolder,"/mergedclustersUMAP.pdf"),width=7,height=7)
mp1
dev.off()

mp2<-plotClusterHeatmap(sce, k = kclustering, m = "merging1")
mp2
pdf(paste0(plotFolder,"/mergedclustersHeatmapWithOldClusters.pdf"),width=7,height=7)
mp2
dev.off()

mp3<-plotClusterHeatmap(sce, k = "merging1")
mp3
pdf(paste0(plotFolder,"/mergedclustersHeatmap.pdf"),width=7,height=7)
mp3
dev.off()




#Differential analysis

library(diffcyt)

class(sce)

#dp1<-plotAbundances(sce, k = "merging1", by = "sample_id")
dp1<-plotAbundances(sce, k = kclustering, by = "sample_id")
dp1
pdf(paste0(plotFolder,"/differentialAbundanceSample.pdf"),width=7,height=7)
dp1
dev.off()

dp2<-plotAbundances(sce, k = kclustering, by = "cluster_id", shape = "patient_id")
dp2
pdf(paste0(plotFolder,"/differentialAbundanceClusters.pdf"),width=7,height=7)
dp2
dev.off()



library(SingleCellExperiment)
library(diffcyt)

library(lme4)
colData(sce)

sce

metadata
metadata(sce)

ei <- metadata(sce)$experiment_info
(da_formula1 <- createFormula(ei,
    cols_fixed = "condition",
    cols_random = "sample_id"))

?createContrast

ei

da_formula1
contrast <- createContrast(c(0, 1))

(da_formula2 <- createFormula(ei,
    cols_fixed = c("condition","response"),
    cols_random = "sample_id"))

da_formula2

contrast <- createContrast(c(0, 1))

contrast

da_res1 <- diffcyt(sce,
    formula = da_formula1, contrast = contrast,
    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
    clustering_to_use = kclustering, verbose = FALSE)
da_res2 <- diffcyt(sce,
    formula = da_formula2, contrast = contrast,
    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
    clustering_to_use = kclustering,verbose = FALSE)

rowData(da_res2$res)

da_res1$res

names(da_res1)
rowData(da_res1$res)
table(rowData(da_res1$res)$p_adj < .1)
