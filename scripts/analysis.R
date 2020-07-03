require(data.table);require(tidyverse);require(xlsx);require(stringdist);require(Matrix);require(deMULTIplex)
matrix_dir = '~/Documents/NMDA10x/'
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2





bclist=c('GGAGAAGA
CCACAATG
TGAGACCT
GCACACGC
AGAGAGAG
TCACAGCA
GAAAAGGG
CGAGATTC
GTAGCACT
CGACCAGC
TTAGCCAG
GGACCCCA
CCAACCGG
TGACCGAT
GCAACGCC
CAATCGGT
ATAGCGTC
GAATCTCG
CTAGCTGA
AGACCTTG
GAAGGAAG
CTACGACA
AGAAGAGG
GTACGCAT
CGAAGCCC
ACATGCGT
TAAGGCTC
GGAAGGAA
CCATGGCG
AAAGGGGA
TTACGGTG
GCATGTAC')
bclistParse=str_split(bclist, pattern = "\n", simplify = T) %>% t() %>% as.character()
write(bclistParse, file = '~/Documents/NMDA10x/multiSeq/data/multiSeqBarcodes_1_to_32.txt')
cellIDs=gsub('-1','', barcode.names$V1)#[sample(10000)]
readTable = fread('~/Documents/NMDA10x/test.csv')

bar.table <- as.data.frame(matrix(0L, ncol=length(bclistParse)+2, nrow=length(cellIDs)))
colnames(bar.table) <- c(paste("Bar", 1:length(bclistParse),sep=""), "nUMI", "nUMI_total")
rownames(bar.table) <- cellIDs
for (cell in cellIDs){
  r1.ind <- which(readTable$Cell == cell)
  umis <- readTable$Umi[r1.ind]
  tags <- readTable$Sample[r1.ind]
  bar.table[cell, "nUMI_total"] <- length(umis)
  tag.dists <- stringdistmatrix(a=tags, b=bclistParse)
  
  tag.hds <- apply(tag.dists, 1, min)
  tag.ind <- which(tag.hds <= 1)
  if (length(tag.ind) == 0) { next } ## Skip cellIDs with 0 aligned read (bug fix)
  tag.dists <- tag.dists[tag.ind, ]
  umis <- umis[tag.ind]
  umi.ind <- which(duplicated(umis) == FALSE)
  bar.table[cell,"nUMI"] <- length(umi.ind)
  for (tag in 1:length(tag.dists)){
    bar.table[cell,tag] <- length(which(tag.dists[tag] <= 1))
  }
}




readTable

tag.dists=stringdistmatrix(a = , b = bclistParse)
tag.dists
tag.hds <- apply(tag.dists, 1, min)
tag.ind <- which(tag.hds <= 1)
tag.dists <- tag.dists[tag.ind, ]


bar.table2 = fread('~/Documents/NMDA10x/multiSeq/outs/barTable.csv')
bar.table2
## Note: Exclude columns 97:98 (assuming 96 barcodes were used) which provide total barcode UMI counts for each cell. 
barTSNE <- function(barTable) {
  require(Rtsne)
  
  ## Normalize barcode count matrix
  n_BC <- ncol(barTable)
  barTable.n <- as.data.frame(log2(barTable))
  for (i in 1:n_BC) {
    ind <- which(is.finite(barTable.n[,i]) == FALSE)
    barTable.n[ind,i] <- 0
    barTable.n[,i] <- barTable.n[,i]-mean(barTable.n[,i])
  }
  
  ## Run tSNE
  tsne.res <- Rtsne(barTable.n, dims=2, initial_dims=n_BC, verbose=TRUE, check_duplicates=FALSE, max_iter=2500)
  tsne.embedding <- as.data.frame(tsne.res$Y)
  colnames(tsne.embedding) <- c("TSNE1","TSNE2")
  tsne.embedding[,3:(n_BC+2)] <- barTable.n
  rownames(tsne.embedding) <- rownames(barTable)
  
  return(tsne.embedding)
}

bar.tsne <- barTSNE(bar.table2[,2:33]) 


pdf("~/Documents/NMDA10x/multiSeq/outs/bc.check.pdf")
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}
dev.off()


bar.table.full <- bar.table2[,2:35]
rownames(bar.table.full) = bar.table2$cellID
good.bars <- paste("Bar",1:32,sep="")  # NOTE: In this hypothetical example, barcodes 91-96 were not detected
bar.table <- bar.table2[, 2:33]  # Remove missing bars and summary columns
rownames(bar.table) = bar.table2$cellID
names(bar.table) = good.bars
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))


round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")] %>% as.numeric()
bar.table <- bar.table[-which(rownames(bar.table) %in% bar.table2$cellID[neg.cells]), ]

## Round 2 -----------------------------------------------------------------------------------------------------
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])

## Repeat until all no negative cells remain (usually 3 rounds)...
final.calls <- c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls) <- c(names(round2.calls),neg.cells)

reclass.cells <- findReclassCells(bar.table, as.numeric(names(final.calls)[which(final.calls=="Negative")]))
reclass.res <- rescueCells(barTable = bar.table, final.calls, reclass.cells)
