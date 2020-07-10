require(data.table);require(tidyverse);require(xlsx);require(stringdist);require(Matrix);require(deMULTIplex)
matrix_dir = '~/Documents/NMDA10x/multiSeq/data/'
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
#readTable = fread('~/Documents/NMDA10x/test.csv')

readTable <- MULTIseq.preProcess(R1 = 'R1.fastq', R2 = 'R2.fastq', cellIDs = cellIDs, cell=c(1,16), umi=c(17,28), tag=c(1,8))
bar.table <- MULTIseq.align(readTable, cellIDs, ref = bclistParse)
setDT(bar.table, keep.rownames = T)
processed.bartable = fread('~/Documents/NMDA10x/multiSeq/outs/barTableCompare.csv')












## Note: Exclude columns 97:98 (assuming 96 barcodes were used) which provide total barcode UMI counts for each cell. 
barTSNE <- function(barTable) {
  require(Rtsne)
  
  ## Normalize barcode count matrix
  ptm <- proc.time()
  n_BC <- ncol(barTable)
  barTable.n <- as.data.frame(log2(barTable))
  ptm <- proc.time()
  for (i in 1:n_BC) {
    ind <- which(is.finite(barTable.n[,i]) == FALSE)
    barTable.n[ind,i] <- 0
    barTable.n[,i] <- barTable.n[,i]-mean(barTable.n[,i])
  }
  proc.time() - ptm
  
  ## Run tSNE
  tsne.res <- Rtsne(barTable.n, dims=2, initial_dims=n_BC, verbose=TRUE, check_duplicates=FALSE, max_iter=2500)
  tsne.embedding <- as.data.frame(tsne.res$Y)
  colnames(tsne.embedding) <- c("TSNE1","TSNE2")
  tsne.embedding[,3:(n_BC+2)] <- barTable.n
  rownames(tsne.embedding) <- rownames(barTable)
  
  return(tsne.embedding)
}

bar.tsne <- barTSNE(bar.table2[,2:33]) 


pdf("~/Documents/NMDA10x/multiSeq/outs/bcAllStrict.check.pdf")
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}
dev.off()


## Round 1 -----------------------------------------------------------------------------------------------------
## Perform Quantile Sweep
bar.table.full <- bar.table[, 1:32]#as.data.frame(row.names = bar.table2$cellID, bar.table)
good.bars <- paste("Bar",1:32,sep="")  # NOTE: In this hypothetical example, barcodes 91-96 were not detected
bar.table <- bar.table.full[, good.bars]  # Remove missing bars and summary columns
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

## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

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



## Perform semi-supervised negative cell reclassification
reclass.cells <- findReclassCells(bar.table.full, neg.cells = unique(names(final.calls)[which(final.calls=="Negative")]))
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
  geom_point() + xlim(c(nrow(reclass.res)-1,1)) + 
  ylim(c(0,1.05)) +
  geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_vline(xintercept = 16, color="blue",lty=2)

## Finalize negative cell rescue results
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 20) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]
### check align
cellIDS = c(bar.table2$cellID[1:30])
bar.ref = bclistParse
readTable = fread('~/Documents/NMDA10x/multiSeq/data/Multi-seq1_S89_L003_ReadTable.csv.gz')
readTable = as.data.frame(readTable)
bar.table.de <- MULTIseq.align(readTable, cellIDS, bar.ref)


barTable = fread('~/Documents/NMDA10x/multiSeq/outs/barTable.csv')
barTable
vec.bc1 = log2(barTable$TGAGACCT)#[,1]
vec.bc1[vec.bc1 == -Inf] =0
vec.bc1 = scale(vec.bc1, center = T, scale = F)[,1]
hist(vec.bc1)
ggplot(data=NULL, aes((vec.bc1)))+geom_density()
require(KernSmooth)
model <- tryCatch( { approxfun(bkde(vec.bc1, kernel="normal")) },
                   error=function(e) { print(paste0("No threshold found for","...")) } )
x <-  seq(from=quantile(vec.bc1,0.001), to=quantile(vec.bc1,0.999), length.out=100)

bkde(vec.bc1, kernel="normal")
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

extrema <- localMaxima(model(x))
low.extreme <- extrema[which.max(model(x)[extrema])]
high.extreme <- max(extrema)
thresh <- quantile(c(x[high.extreme], x[low.extreme]), 0.70)
cell_i <- which(vec.bc1 >= thresh)
out=data.frame(counts= barTable$TGAGACCT, norm.counts=vec.bc1, class='neg')
out$class[cell_i] = 'Pos'
ggplot(out, aes(norm.counts, fill=factor(class)))+geom_density()

require(reshape2)
call.list = barTable_sweep.list
res <- as.data.frame(matrix(0L, nrow=length(call.list), ncol=4))
colnames(res) <- c("q","pDoublet","pNegative","pSinglet")
q.range <- unlist(strsplit(names(call.list), split="q="))
res$q <- as.numeric(q.range[grep("0", q.range)])
nCell <- length(call.list[[1]])

for (i in 1:nrow(res)) {
  temp <- table(call.list[[i]])
  if ( "Doublet" %in% names(temp) == TRUE ) { res$pDoublet[i] <- temp[which(names(temp) == "Doublet")] }
  if ( "Negative" %in% names(temp) == TRUE ) { res$pNegative[i] <- temp[which(names(temp) == "Negative")] }
  res$pSinglet[i] <- sum(temp[which(names(temp) != "Doublet" & names(temp) != "Negative")])
}

res <- melt(res, id.vars="q")
res[,4] <- res$value/nCell
colnames(res)[2:4] <- c("Subset","nCells","Proportion")

extrema <- res$q[localMaxima(res$Proportion[which(res$Subset == "pSinglet")])]





bar.tsne=barTSNE(barTable = barTable[, 2:33])
bar.tsne = as.data.frame(bar.tsne)
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


setDT(bar.tsne)
ggplot(bar.tsne, aes(TSNE1, TSNE2, color=GGAGAAGA))+geom_point(alpha=0.2)+
  scale_color_gradient(low = "black", high = "red")
bar.tsne$CELL.ID = barTable$cellID
bar.tsne.melt=melt(bar.tsne, id.vars = c("TSNE1","TSNE2", "CELL.ID"))
require(tidyverse)
ggplot(bar.tsne.melt, aes(TSNE1, TSNE2, color=value))+geom_point(size=0.1)+
  scale_color_viridis_c()+facet_wrap(~variable, nrow = 8, ncol = 4)+theme_gdocs(base_size = 8)
ggsave('~/Documents/NMDA10x/multiSeq/outs/bc.check.png', device = 'png', width = 8, height = 10, dpi = 500)


