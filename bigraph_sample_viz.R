#Calculation
#source("/Users/dokada/Dropbox/analysis/2022.5/bigraph_sample_viz.R") #Run in 2024.5
out_path <- "/Users/dokada/Desktop/work/bigraph_sample_viz/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#load library
set.seed(520)
library(TabulaMurisSenisData)
library(SingleCellExperiment)
library(multigraph)

#Create graph
raw_data <- TabulaMurisSenisFACS(tissues = "All", processedCounts = TRUE)
count_mat <- as.matrix(logcounts(raw_data$All))
sample_annot <- colData(raw_data$All)
gene_annot <- rowData(raw_data$All)
tissue_all <- sample_annot$tissue
celltype_all <- sample_annot$cell_ontology_class
tc_all <- paste0(tissue_all, ".", celltype_all)
unq_tc <- sort(unique(tc_all))
tissue <- sapply(strsplit(unq_tc, "\\."), function(x){x[1]})
celltype <- sapply(strsplit(unq_tc, "\\."), function(x){x[2]})
dat <- data.frame(tissue, celltype)
ts <- unique(dat$tissue)
ct <- unique(dat$celltype)
m <- matrix(0, nrow=length(ts), ncol=length(ct))
rownames(m) <- ts
colnames(m) <- ct
for(i in 1:nrow(dat)){
    m[dat[i, "tissue"], dat[i, "celltype"]] <- 1
}

#Sub-dataset
selected_ts <- c('BAT', 'GAT', 'MAT', 'SCAT')
selected_ct <- c('B cell', 'NK cell', 'T cell', 'endothelial cell', 'epithelial cell', 'mesenchymal stem cell of adipose', 'myeloid cell')
m_sub <- m[selected_ts, selected_ct]
png(paste0(out_path, "mgraph_sub1.png"))
bmgraph(m_sub, showLbs=T, fsize=10, lwd=7) 
dev.off()

#Sub-dataset
selected_ts <- c('Brain_Non-Myeloid', 'Diaphragm', 'Limb_Muscle', 'MAT', 'Trachea')
selected_ct <- c('T cell', 'endothelial cell', 'macrophage')
m_sub <- m[selected_ts, selected_ct]
png(paste0(out_path, "mgraph_sub2.png"))
bmgraph(m_sub, showLbs=T, fsize=10, lwd=7) 
dev.off()