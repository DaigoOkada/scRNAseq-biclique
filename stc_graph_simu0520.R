#Calculation
#source("/Users/dokada/Dropbox/analysis/2022.5/stc_graph_simu0520.R") #Run in 2024.5
out_path <- "/Users/dokada/Desktop/work/stc_graph_simu0520/"
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
library(qqman)

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
edge_file <- paste0(out_path, "celltissue.csv")
write.csv(dat, edge_file, row.names=FALSE)
out_file <- paste0(out_path, "selected_celltissue.txt")
command <- paste0("python /Users/dokada/Desktop/work/haraguchi_lab/code_0510/main.py ", edge_file, " > ", out_file)
system(command)

#Visualizaiton
ts <- unique(dat$tissue)
ct <- unique(dat$celltype)
m <- matrix(0, nrow=length(ts), ncol=length(ct))
rownames(m) <- ts
colnames(m) <- ct
for(i in 1:nrow(dat)){
    m[dat[i, "tissue"], dat[i, "celltype"]] <- 1
}
colnames(m) <- rep(NA, ncol(m))
png(paste0(out_path, "mgraph.png"), width=960, height=960)
bmgraph(m, showLbs=T, fsize=20) 
dev.off()

#Simulation
set.seed(520)
pats <- expand.grid(c(0,1), c(0,1))
colnames(pats) <- c("tissue", "celltype")
pats_names <- paste0(pats$tissue, pats$celltype)
n_pats <- nrow(pats)
n_genes <- 10000
unq_tissue <- sort(unique(dat$tissue))
unq_ct <- sort(unique(dat$celltype))
n_tissue <- length(unq_tissue)
n_ct <- length(unq_ct)
n_each <- 10
for(i in 1:n_pats){
    pvals_mat <- NULL
    pvals_mat_sel <- NULL
    tmp_pats_name <- pats_names[i]
    for(j in 1:n_genes){
        eff_tissue <- pats[i, "tissue"] * rnorm(n_tissue, 0, 1)
        names(eff_tissue) <- unq_tissue
        eff_ct <- pats[i, "celltype"] * rnorm(n_ct, 0, 1)
        names(eff_ct) <- unq_ct
        y_each <- NULL
        tissue_each <- NULL
        celltype_each <- NULL
        for(k in 1:nrow(dat)){
            y_tmp <- rep(eff_tissue[dat$tissue[k]] + eff_ct[dat$celltype[k]], n_each)
            y_each <- c(y_each, y_tmp)
            tissue_each <- c(tissue_each, rep(dat$tissue[k], n_each))
            celltype_each <- c(celltype_each, rep(dat$celltype[k], n_each))
        }
        y_each <- y_each + rnorm(n_each*nrow(dat), 0, 0.01)
        dat2 <- data.frame(tissue_each, celltype_each, y_each)
        colnames(dat2) <- c("tissue", "celltype", "y")

        #Apply GLM for  full dataset
        model <- glm(y ~ tissue + celltype, data = dat2)
        pvals <- summary(model)$coefficients[,"Pr(>|t|)"]
        pvals_mat <- cbind(pvals_mat, pvals)

        #Apply GLM for  selected dataset
        selected_tissue <- c('BAT', 'GAT', 'MAT', 'SCAT')
        selected_ct <- c('B cell', 'NK cell', 'T cell', 'endothelial cell', 'epithelial cell', 'mesenchymal stem cell of adipose', 'myeloid cell')
        select_stc_idx <- which(dat2$tissue %in% selected_tissue & dat2$celltype %in% selected_ct)
        dat3 <- dat2[select_stc_idx,]
        if(length(selected_tissue)*length(selected_ct)*n_each != nrow(dat3)){
            stop("Error")
        }
        model <- glm(y ~ tissue + celltype, data = dat3)
        pvals <- summary(model)$coefficients[,"Pr(>|t|)"]
        pvals_mat_sel <- cbind(pvals_mat_sel, pvals)
        cat(i, j, "\n")
    }

    #全体を使った場合とサブレコードを遣った場合のQQ plot
    all_sel_can <- c("all", "sub")
    tmp2_can <- c("tissue", "celltype")
    command <- paste0("cd ", out_path, " ; montage -tile 2x2  -geometry +0+0 ")
    for(j2 in 1:length(all_sel_can)){
        for(j in 1:length(tmp2_can)){
            tmp2 <- tmp2_can[j]
            tmp_as <- all_sel_can[j2]
            if(tmp_as == "all"){
                pvals_mat_sub <- pvals_mat[grepl(paste0("^", tmp2), rownames(pvals_mat)),,drop=F]
            }else{
                pvals_mat_sub <- pvals_mat_sel[grepl(paste0("^", tmp2), rownames(pvals_mat_sel)),,drop=F]
            }
            png(paste0(out_path, tmp_pats_name, ".", tmp_as, ".", tmp2, "QQplot.png"))
            par(mar = c(5, 5, 5, 2)) ##bottom, left, top, right
            xlab = "Expected -log10(P value)"
            ylab = "Observed -log10(P value)"
            main <- paste0(tmp2, ".", tmp_as)
            qqman::qq(as.numeric(pvals_mat_sub), xlab="", ylab="", main="", cex=1, cex.axis=2, pch=19)
            mtext(xlab, side=1, line=3, cex=2)
            mtext(ylab, side=2, line=3, cex=2)
            mtext(main, side=3, line=2, cex=2, adj=0)
            dev.off()
            command <- paste0(command, paste0(out_path, tmp_pats_name, ".", tmp_as, ".", tmp2, "QQplot.png"), " ")
        }
    }
    command <- paste0(command, tmp_pats_name, "qqsimu.png")
    system(command)
}
