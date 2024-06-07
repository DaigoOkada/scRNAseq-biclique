#source("/Users/dokada/Dropbox/analysis/2022.5//bigraph_realdata_ts0520.R") ＃Run in 2024.5.27
out_path <- "/Users/dokada/Desktop/work/bigraph_realdata_ts0520/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

library(TabulaMurisSenisData)
library(WebGestaltR)
library(SingleCellExperiment)
library(effsize)
dbs <- c("geneontology_Biological_Process_noRedundant") #for GSEA


#Comn¥mpare two tissue (FACS dataset)
set.seed(520)
tmp_ds <-  "FACS"
raw_data <- TabulaMurisSenisFACS(tissues = "All", processedCounts = TRUE)
count_mat <- as.matrix(logcounts(raw_data$All))
sample_annot <- colData(raw_data$All)
gene_annot <- rowData(raw_data$All)
#Association analysis
id_all <- as.character(sample_annot$mouse.id)
if(tmp_ds == "FACS"){
    age_all <-  sapply(id_all, function(chr){as.numeric(strsplit(chr,"_")[[1]][1])})
}else{
    age_all <- sapply(id_all, function(chr){as.numeric(strsplit(chr,"-")[[1]][1])})
}
sex_all <- sample_annot$sex
tissue_all <- sample_annot$tissue
celltype_all <- sample_annot$cell_ontology_class
stcia_all <- paste0(sex_all, ".", tissue_all, ".", celltype_all, ".", id_all,".", age_all)
stc_all <- paste0(sex_all, ".", tissue_all, ".", celltype_all)
stca_all <- paste0(sex_all, ".", tissue_all, ".", celltype_all, ".", age_all)
tci_all <- paste0(tissue_all, ".", celltype_all, ".", id_all)
selected_id <- c('3_10_M', '3_11_M', '3_8_M', '3_9_M')
selected_tissue <- c('GAT', 'SCAT')
selected_ct <- c('endothelial cell', 'mesenchymal stem cell of adipose', 'myeloid cell')
select_stc_idx <- which((id_all %in% selected_id) & (tissue_all %in% selected_tissue) & (celltype_all %in% selected_ct))
coefs_mat <- NULL
for(i in 1:nrow(count_mat)){
    y <- count_mat[i,select_stc_idx]
    dat2 <- data.frame(y=y, id=id_all[select_stc_idx], tissue=tissue_all[select_stc_idx], celltype=celltype_all[select_stc_idx])
    model <- glm(y ~ id + tissue + celltype, data = dat2)
    coefs <- summary(model)$coefficients["tissueSCAT", c("Estimate", "Pr(>|t|)")]
    coefs_mat <- rbind(coefs_mat, coefs)
    cat(i, "\n")
}
rownames(coefs_mat) <- rownames(count_mat)
png(paste0(out_path, "facs.qqplot_tissueSCAT.png"))
par(mar = c(5, 5, 5, 2)) ##bottom, left, top, right
xlab = "Expected -log10(P value)"
ylab = "Observed -log10(P value)"
main <- "GAT vs SCAT"
qqman::qq(coefs_mat[,"Pr(>|t|)"], xlab="", ylab="", main="", cex=1, cex.axis=2, pch=19)
mtext(xlab, side=1, line=3, cex=2)
mtext(ylab, side=2, line=3, cex=2)
mtext(main, side=3, line=2, cex=2, adj=0)
dev.off()

#ORA analysis
adjp <- p.adjust(coefs_mat[,"Pr(>|t|)"], method="bonferroni")
gene_set <- names(adjp)[which(adjp < 0.05)] #829 genes
bg_geneset <- names(adjp)
dbs <- c("geneontology_Biological_Process_noRedundant")
enrichResult <- WebGestaltR(enrichMethod="ORA", organism="mmusculus", enrichDatabase=dbs, interestGene=gene_set, interestGeneType="genesymbol", 
referenceGene=bg_geneset, referenceGeneType="genesymbol", isOutput=FALSE, fdrThr=0.05)
write.csv(enrichResult, paste0(out_path, tmp_ds , ".tis.ORA_result.csv"), row.names=FALSE)

#Barplot visualization
if(!is.null(dim(enrichResult))){
    top_num <- min(10, nrow(enrichResult))
    tmp <- enrichResult[1:top_num, ]
    en <- tmp$enrichmentRatio
    names(en) <-tmp$description
    en <- sort(en)
    left <- 60
    png(paste0(out_path, "GATSCAT.ora.png"), width=960, height=960)
    par(mar = c(9, left, 2, 2)) ##bottom, left, top, right
    main = paste0("GAT vs SCAT")
    barplot(en, beside = TRUE, las = 1, horiz=T, col = "blue", names.arg = names(en), cex.axis=3, cex.lab=3, cex.names=3)
    dev.off()
}

#out_put
coefs_mat <- cbind(coefs_mat, adjp)
write.csv(coefs_mat, paste0(out_path, "FACS.coefs_mat.csv"), row.names=FALSE)


#Dropletで実施
tmp_ds <-  "Droplet"
raw_data <- TabulaMurisSenisDroplet(tissues = "All", processedCounts = TRUE)
count_mat <- as.matrix(logcounts(raw_data$All))
sample_annot <- colData(raw_data$All)
gene_annot <- rowData(raw_data$All)
id_all <- as.character(sample_annot$mouse.id)
if(tmp_ds == "FACS"){
    age_all <-  sapply(id_all, function(chr){as.numeric(strsplit(chr,"_")[[1]][1])})
}else{
    age_all <- sapply(id_all, function(chr){as.numeric(strsplit(chr,"-")[[1]][1])})
}
sex_all <- sample_annot$sex
tissue_all <- sample_annot$tissue
celltype_all <- sample_annot$cell_ontology_class
stcia_all <- paste0(sex_all, ".", tissue_all, ".", celltype_all, ".", id_all,".", age_all)
stc_all <- paste0(sex_all, ".", tissue_all, ".", celltype_all)
stca_all <- paste0(sex_all, ".", tissue_all, ".", celltype_all, ".", age_all)
tci_all <- paste0(tissue_all, ".", celltype_all, ".", id_all)
selected_id <- c('3-F-56', '3-F-57')
selected_tissue <- c('Limb_Muscle', 'Mammary_Gland')
selected_ct <- c('B cell', 'T cell', 'endothelial cell', 'macrophage')
select_stc_idx <- which(id_all %in% selected_id & tissue_all %in% selected_tissue & celltype_all %in% selected_ct)
coefs_mat <- NULL
for(i in 1:nrow(count_mat)){
    y <- count_mat[i,select_stc_idx]
    dat2 <- data.frame(y=y, id=id_all[select_stc_idx], tissue=tissue_all[select_stc_idx], celltype=celltype_all[select_stc_idx])
    model <- glm(y ~ id + tissue + celltype, data = dat2)
    coefs <- summary(model)$coefficients["tissueMammary_Gland", c("Estimate", "Pr(>|t|)")]
    coefs_mat <- rbind(coefs_mat, coefs)
    cat(i, "\n")
}
rownames(coefs_mat) <- rownames(count_mat)
png(paste0(out_path, "Droplet.qqplot_tissueMG.png"))
par(mar = c(5, 5, 5, 2)) ##bottom, left, top, right
xlab = "Expected -log10(P value)"
ylab = "Observed -log10(P value)"
main <- paste0("Limb_Muscle vs Mammary_Gland")
qqman::qq(coefs_mat[,"Pr(>|t|)"], xlab="", ylab="", main="", cex=1, cex.axis=2, pch=19)
mtext(xlab, side=1, line=3, cex=2)
mtext(ylab, side=2, line=3, cex=2)
mtext(main, side=3, line=2, cex=2, adj=0)
dev.off()

#ORA analysis
adjp <- p.adjust(coefs_mat[,"Pr(>|t|)"], method="bonferroni")
gene_set <- names(adjp)[which(adjp < 0.05)] #1199 genes
bg_geneset <- names(adjp)
dbs <- c("geneontology_Biological_Process_noRedundant")
enrichResult <- WebGestaltR(enrichMethod="ORA", organism="mmusculus", enrichDatabase=dbs, interestGene=gene_set, interestGeneType="genesymbol", 
referenceGene=bg_geneset, referenceGeneType="genesymbol", isOutput=FALSE, fdrThr=0.05)
write.csv(enrichResult, paste0(out_path, tmp_ds , ".tis.ORA_result.csv"), row.names=FALSE)

#Barplot visualization
if(!is.null(dim(enrichResult))){
    top_num <- min(10, nrow(enrichResult))
    tmp <- enrichResult[1:top_num, ]
    en <- tmp$enrichmentRatio
    names(en) <-tmp$description
    en <- sort(en)
    left <- 60
    png(paste0(out_path, "Limb_Muscle_Mammary_Gland.ora.png"), width=960, height=960)
    par(mar = c(9, left, 2, 2)) ##bottom, left, top, right
    main = paste0("Limb_Muscle vs Mammary_Gland")
    barplot(en, beside = TRUE, las = 1, horiz=T, col = "blue", names.arg = names(en), cex.axis=3, cex.lab=3, cex.names=3)
    dev.off()
}


#out_put
coefs_mat <- cbind(coefs_mat, adjp)
write.csv(coefs_mat, paste0(out_path, "Droplet.coefs_mat.csv"), row.names=FALSE)
