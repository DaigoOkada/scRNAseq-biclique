#Calculation
#source("/Users/dokada/Dropbox/analysis/2022.5//bigraph_realdata_ag0520.R")　 #Run in 2024.5.22
out_path <- "/Users/dokada/Desktop/work/bigraph_realdata_ag0520/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#Prepare data
library(TabulaMurisSenisData)
library(WebGestaltR)
library(SingleCellExperiment)
library(effsize)
set.seed(520)
dbs <- c("geneontology_Biological_Process_noRedundant") #for GSEA

#BAT vs SCAT
tmp_ds <- "FACS"
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
selected_sex <- c('female', 'male')
selected_tissue <- c('BAT', 'SCAT')
selected_ct <- c('B cell', 'endothelial cell', 'myeloid cell')
res_tissue_list <- list()
for(t in 1:2){
    tmp_tissue <- selected_tissue[t]
    select_stc_idx <- which(sex_all %in% selected_sex & tissue_all == tmp_tissue & celltype_all %in% selected_ct & age_all != 1)
    age_vec <- ifelse(age_all[select_stc_idx] == 3, "Young", "Old")
    res_mat <- NULL
    for(i in 1:nrow(count_mat)){
        y <- count_mat[i,select_stc_idx]
        dat2 <- data.frame(y=y, sex=sex_all[select_stc_idx], celltype=celltype_all[select_stc_idx], age=age_vec)
        model <- glm(y ~ sex +  celltype + age, data = dat2)
        res_tmp <- summary(model)$coefficients["ageYoung",]
        res_mat <- cbind(res_mat, res_tmp)
        cat(i, "\n")
    }
    adjp <- p.adjust(res_mat["Pr(>|t|)",],method="bonferroni")
    res_mat <- rbind(res_mat, adjp)
    colnames(res_mat) <- rownames(count_mat)
    res_tissue_list[[t]] <- res_mat
}
t1 <- res_tissue_list[[1]]
t2 <- res_tissue_list[[2]]
t1_adjp <- p.adjust(t1["Pr(>|t|)",], method="BH")
t2_adjp <- p.adjust(t2["Pr(>|t|)",], method="BH")
pn_genes_idx <- which(t1_adjp < 0.05 & t2_adjp < 0.05 & t1["Estimate",] > 0 & t2["Estimate",] < 0)
np_genes_idx <- which(t1_adjp < 0.05 & t2_adjp < 0.05 & t1["Estimate",] < 0 & t2["Estimate",] > 0)

#Output
output <- cbind(t1["Estimate",], t1["Pr(>|t|)",], t1_adjp, t2["Estimate",], t2["Pr(>|t|)",], t2_adjp)
colnames(output) <- c("Coefficient_Young(BAT)", "P-value_Young(BAT)", "BH_adjP_Young(BAT)", "Coefficient_Young(SCAT)", "P-value_Young(SCAT)", "BH_adjP_Young(SCAT)")
write.csv(output, paste0(out_path, "BATSCAT_aging_result.csv"), row.names=TRUE)

#Plots
cols <- rep("black", nrow(count_mat))
cols[pn_genes_idx] <- "red"
cols[np_genes_idx] <- "blue"
write.table(table(cols), paste0(out_path, "BAT_SCAT_aging.txt"), row.names=TRUE, col.names=TRUE)
png(paste0(out_path, "BAT_SCAT_aging.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
xlab = "BAT"
ylab = "SCAT"
main <- paste0("Coefficient of Young")
plot(t1["Estimate",], t2["Estimate",], col=cols, xlab="", ylab="", main="", cex=3, cex.axis=4, pch=19)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4, adj=0)
dev.off()


#ORA analysis
gene_set <- rownames(count_mat)[c(pn_genes_idx, np_genes_idx)]
bg_geneset <- rownames(count_mat)
dbs <- c("geneontology_Biological_Process_noRedundant")
enrichResult <- WebGestaltR(enrichMethod="ORA", organism="mmusculus", enrichDatabase=dbs, interestGene=gene_set, interestGeneType="genesymbol", 
referenceGene=bg_geneset, referenceGeneType="genesymbol", isOutput=FALSE, fdrThr=0.05)
write.csv(enrichResult, paste0(out_path, "facs.BATSCAT.ORA_aging.csv"), row.names=FALSE)

#Barplot visualization
if(!is.null(dim(enrichResult))){
    tmp <- enrichResult
    en <- tmp$enrichmentRatio
    names(en) <-tmp$description
    en <- sort(en)
    left <- 50
    png(paste0(out_path, "BATSCAT.ora.png"), width=960, height=960)
    par(mar = c(9, left, 2, 2)) ##bottom, left, top, right
    main = paste0("BAT vs SCAT")
    barplot(en, beside = TRUE, las = 1, horiz=T, col = "blue", names.arg = names(en), cex.axis=3, cex.lab=3, cex.names=3)
    dev.off()
}


#MAT vs SCAT
selected_sex <- c('female', 'male')
selected_tissue <- c('MAT', 'SCAT')
selected_ct <- c('B cell', 'endothelial cell', 'mesenchymal stem cell of adipose')
res_tissue_list <- list()
for(t in 1:2){
    tmp_tissue <- selected_tissue[t]
    select_stc_idx <- which(sex_all %in% selected_sex & tissue_all == tmp_tissue & celltype_all %in% selected_ct & age_all != 1)
    age_vec <- ifelse(age_all[select_stc_idx] == 3, "Young", "Old")
    res_mat <- NULL
    for(i in 1:nrow(count_mat)){
        y <- count_mat[i,select_stc_idx]
        dat2 <- data.frame(y=y, sex=sex_all[select_stc_idx], celltype=celltype_all[select_stc_idx], age=age_vec)
        model <- glm(y ~ sex +  celltype + age, data = dat2)
        res_tmp <- summary(model)$coefficients["ageYoung",]
        res_mat <- cbind(res_mat, res_tmp)
        cat(i, "\n")
    }
    adjp <- p.adjust(res_mat["Pr(>|t|)",],method="bonferroni")
    res_mat <- rbind(res_mat, adjp)
    colnames(res_mat) <- rownames(count_mat)
    res_tissue_list[[t]] <- res_mat
}
t1 <- res_tissue_list[[1]]
t2 <- res_tissue_list[[2]]
t1_adjp <- p.adjust(t1["Pr(>|t|)",], method="BH")
t2_adjp <- p.adjust(t2["Pr(>|t|)",], method="BH")
pn_genes_idx <- which(t1_adjp < 0.05 & t2_adjp < 0.05 & t1["Estimate",] > 0 & t2["Estimate",] < 0)
np_genes_idx <- which(t1_adjp < 0.05 & t2_adjp < 0.05 & t1["Estimate",] < 0 & t2["Estimate",] > 0)

#Output
output <- cbind(t1["Estimate",], t1["Pr(>|t|)",], t1_adjp, t2["Estimate",], t2["Pr(>|t|)",], t2_adjp)
colnames(output) <- c("Coefficient_Young(MAT)", "P-value_Young(MAT)", "BH_adjP_Young(MAT)", "Coefficient_Young(SCAT)", "P-value_Young(SCAT)", "BH_adjP_Young(SCAT)")
write.csv(output, paste0(out_path, "MATSCAT_aging_result.csv"), row.names=TRUE)

#Plots
cols <- rep("black", nrow(count_mat))
cols[pn_genes_idx] <- "red"
cols[np_genes_idx] <- "blue"
write.table(table(cols), paste0(out_path, "MAT_SCAT_aging.txt"), row.names=TRUE, col.names=TRUE)
png(paste0(out_path, "MAT_SCAT_aging.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
xlab = "MAT"
ylab = "SCAT"
main <- paste0("Coefficient of Young")
plot(t1["Estimate",], t2["Estimate",], col=cols, xlab="", ylab="", main="", cex=3, cex.axis=4, pch=19)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4, adj=0)
dev.off()

#ORA analysis
gene_set <- rownames(count_mat)[c(pn_genes_idx, np_genes_idx)]
bg_geneset <- rownames(count_mat)
dbs <- c("geneontology_Biological_Process_noRedundant")
enrichResult <- WebGestaltR(enrichMethod="ORA", organism="mmusculus", enrichDatabase=dbs, interestGene=gene_set, interestGeneType="genesymbol", 
referenceGene=bg_geneset, referenceGeneType="genesymbol", isOutput=FALSE, fdrThr=0.05)
write.csv(enrichResult, paste0(out_path, "facs.MATSCAT.ORA_aging.csv"), row.names=FALSE)

#Barplot visualization
if(!is.null(dim(enrichResult))){
    tmp <- enrichResult
    en <- tmp$enrichmentRatio
    names(en) <-tmp$description
    en <- sort(en)
    left <- 50
    png(paste0(out_path, "MATSCAT.ora.png"), width=960, height=960)
    par(mar = c(9, left, 2, 2)) ##bottom, left, top, right
    main = paste0("MAT vs SCAT")
    barplot(en, beside = TRUE, las = 1, horiz=T, col = "blue", names.arg = names(en), cex.axis=3, cex.lab=3, cex.names=3)
    dev.off()
}

