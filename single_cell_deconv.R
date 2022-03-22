# upload datasets:
moffitt_df      <- read.table(...)
collisson_df    <- read.table(...)
bailey_df       <- read.table(...)
icgc_array_df   <- read.table(...)
badea_df        <- read.table(...)
tcga_df         <- read.table(...)
yang_df         <- read.table(...)
sandhu_df       <- read.table(...)

# remove the 1000 genes less variables from Sandhu, ICGC and TCGA data since input too big for deconv function
sandhu_df$var <- rowVars(as.matrix(sandhu_df[,c(-1)]))
sandhu_df <- sandhu_df[order(sandhu_df$var),]
dim(sandhu_df)
sandhu_df <- sandhu_df[1000:dim(sandhu_df)[1], colnames(sandhu_df) != "var"]

icgc_array_df$var <- rowVars(as.matrix(icgc_array_df[,c(-1)]))
icgc_array_df <- icgc_array_df[order(icgc_array_df$var),]
icgc_array_df <- icgc_array_df[1000:dim(icgc_array_df)[1], colnames(icgc_array_df) != "var"]

tcga_df <- tcga_df[rownames(tcga_df) %in% rownames(icgc_array_df),]

library(Seurat)
library(omnideconv)
# get Tosti et al. human pancreas single-cell data
sc_file <- readRDS(file = ".../data/single-cells/adult_pancreas.rds")
DimPlot(sc_file, reduction = "umap")
sc_file[["percent.mt"]] <- PercentageFeatureSet(sc_file, pattern = "^MT-")
VlnPlot(sc_file, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# We filter cells that have unique feature counts over 2,500 or less than 200
sc_filter <- subset(sc_file, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
rna = sc_filter@assays$RNA@count

# build a signature matrix
rna_moffitt = rna[rownames(rna) %in% intersect(rownames(rna), rownames(moffitt_df)), ]
rna_collisson = rna[rownames(rna) %in% intersect(rownames(rna), rownames(collisson_df)), ]
rna_bailey = rna[rownames(rna) %in% intersect(rownames(rna), rownames(bailey_df)), ]
rna_icgc_array = rna[rownames(rna) %in% intersect(rownames(rna), rownames(icgc_array_df)), ]
rna_badea = rna[rownames(rna) %in% intersect(rownames(rna), rownames(badea_df)), ]
rna_tcga = rna[rownames(rna) %in% intersect(rownames(rna), rownames(tcga_df)), ]
rna_yang = rna[rownames(rna) %in% intersect(rownames(rna), rownames(yang_df)), ]
rna_sandhu = rna[rownames(rna) %in% intersect(rownames(rna), rownames(sandhu_df)), ]

sig_mat_moffitt    = omnideconv::build_model(rna_moffitt, cell_type_annotations = as.vector(sc_filter$Cluster),   method="cibersortx")
sig_mat_collisson  = omnideconv::build_model(rna_collisson, cell_type_annotations = as.vector(sc_filter$Cluster), method="cibersortx")
sig_mat_bailey     = omnideconv::build_model(rna_bailey, cell_type_annotations = as.vector(sc_filter$Cluster),    method="cibersortx")
sig_mat_icgc_array = omnideconv::build_model(rna_icgc_array, cell_type_annotations = as.vector(sc_filter$Cluster),method="cibersortx")
sig_mat_badea      = omnideconv::build_model(rna_badea, cell_type_annotations = as.vector(sc_filter$Cluster),     method="cibersortx")
sig_mat_tcga       = omnideconv::build_model(rna_tcga, cell_type_annotations = as.vector(sc_filter$Cluster),      method="cibersortx")
sig_mat_yang       = omnideconv::build_model(rna_yang, cell_type_annotations = as.vector(sc_filter$Cluster),      method="cibersortx")
sig_mat_sandhu     = omnideconv::build_model(rna_sandhu, cell_type_annotations = as.vector(sc_filter$Cluster),    method="cibersortx")

deconv_moffitt    = omnideconv::deconvolute(bulk_gene_expression=moffitt_df, signature=sig_mat_moffitt, method="cibersortx")
deconv_collisson  = omnideconv::deconvolute(bulk_gene_expression=collisson_df, signature=sig_mat_collisson, method="cibersortx")
deconv_bailey     = omnideconv::deconvolute(bulk_gene_expression=bailey_df, signature=sig_mat_bailey, method="cibersortx")
deconv_icgc_array = omnideconv::deconvolute(bulk_gene_expression=icgc_array_df, signature=sig_mat_icgc_array, method="cibersortx")
deconv_badea      = omnideconv::deconvolute(bulk_gene_expression=badea_df, signature=sig_mat_badea, method="cibersortx")
deconv_tcga       = omnideconv::deconvolute(bulk_gene_expression=tcga_df, signature=sig_mat_tcga, method="cibersortx")
deconv_yang       = omnideconv::deconvolute(bulk_gene_expression=yang_df, signature=sig_mat_yang, method="cibersortx")
deconv_sandhu     = omnideconv::deconvolute(bulk_gene_expression=sandhu_df, signature=sig_mat_sandhu, method="cibersortx")

write.csv(deconv_moffitt, ...)
write.csv(deconv_collisson, ...)
write.csv(deconv_bailey, ...)
write.csv(deconv_badea, ...)
write.csv(deconv_yang, ...)
write.csv(deconv_sandhu, ...)
write.csv(deconv_icgc_array, ...)
write.csv(deconv_tcga, ...)
