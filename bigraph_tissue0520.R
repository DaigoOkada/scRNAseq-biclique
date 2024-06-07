#source("/Users/dokada/Dropbox/analysis/2022.5//bigraph_tissue0520.R")　＃Run in 2-２４.５.21
out_path <- "/Users/dokada/Desktop/work/bigraph_tissue0520/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#analysis
library(TabulaMurisSenisData)
library(WebGestaltR)
library(SingleCellExperiment)
library(effsize)
set.seed(520)
dataset <- c("FACS", "Droplet")
for(ds in 1:2){

    #set parameters
    rm(list=setdiff(ls(), c("out_path", "dataset", "an_type", "ds", "dbs")))
    gc()
    tmp_ds <- dataset[ds]

    #load data
    if(tmp_ds == "FACS"){
        raw_data <- TabulaMurisSenisFACS(tissues = "All", processedCounts = TRUE)
    }else{
        raw_data <- TabulaMurisSenisDroplet(tissues = "All", processedCounts = TRUE)
    }
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

    tab <- table(tci_all[age_all == 3])
    used_stc <- names(tab)[tab > 10]
    sex <- sapply(strsplit(used_stc, "\\."), function(x){x[1]})
    tissue <- sapply(strsplit(used_stc, "\\."), function(x){x[2]})
    celltype <- sapply(strsplit(used_stc, "\\."), function(x){x[3]})
    id <- sapply(strsplit(used_stc, "\\."), function(x){x[4]})
    dat <- data.frame(id, tissue, celltype)
    edge_file <- paste0(out_path, "celltissue_young.csv")
    write.csv(dat, edge_file, row.names=FALSE)
    #Graph analysis (Haraghchi lab)
    out_file <- paste0(out_path, tmp_ds, "selected_celltissue.txt")
    command <- paste0("python /Users/dokada/Desktop/work/haraguchi_lab/code_0510/main.py ", edge_file, " > ", out_file)
    system(command)

    #TCI 
    tab <- table(tci_all[age_all == 3])
    used_stc <- names(tab)[tab > 10]
    tissue <- sapply(strsplit(used_stc, "\\."), function(x){x[1]})
    celltype <- sapply(strsplit(used_stc, "\\."), function(x){x[2]})
    id <- sapply(strsplit(used_stc, "\\."), function(x){x[3]})
    dat <- data.frame(id, tissue, celltype)
    edge_file <- paste0(out_path, "celltissue_young.csv")
    write.csv(dat, edge_file, row.names=FALSE)
    #Graph analysis (Haraghchi lab)
    out_file <- paste0(out_path, tmp_ds, "selected_celltissue.txt")
    command <- paste0("python /Users/dokada/Desktop/work/haraguchi_lab/code_0510/main.py ", edge_file, " > ", out_file)
    system(command)
}

