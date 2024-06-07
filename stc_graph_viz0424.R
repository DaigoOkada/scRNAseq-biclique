#Calculation
#source("/Users/dokada/Dropbox/analysis/2022.5/stc_graph_viz0424.R") #Run in 2024.5.4
out_path <- "/Users/dokada/Desktop/work/stc_graph_viz0424/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#load edge list
library(igraph)
dat <- read.csv("/Users/dokada/Desktop/work/stc_graph_simu0424/celltissue.csv", header=TRUE)
ts <- unique(dat$tissue)
ct <- unique(dat$celltype)
node <- c(ts, ct)
tc <- c(rep("tisssue", length(ts)), rep("celltype", length(ct)))
node_type <- data.frame(node, tc)
write.csv(node_type, "/Users/dokada/Desktop/work/stc_graph_simu0424/node_type.csv", row.names=FALSE)


#Vizualization
library(multigraph)
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