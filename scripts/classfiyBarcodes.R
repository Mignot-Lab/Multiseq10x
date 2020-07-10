require(data.table);require(tidyverse);require(Matrix);require(deMULTIplex)
args = commandArgs(trailingOnly=TRUE)
args[1] = '~/Documents/NMDA10x/multiSeq/outs/barTable.csv'
args[2] = '~/Documents/NMDA10x/multiSeq/data/barcodes.tsv.gz'
args[3] = '~/Documents/NMDA10x/multiSeq/outs/test1'
barTableMain = read.csv(args[1])
rownames(barTableMain) = barTableMain$cellID
cell.barCodes = read.delim(args[2], header = F,stringsAsFactors = F)
outFile = args[3]

## Tsne Check 

## Visualize barcode space
bar.tsne <- barTSNE(barTableMain[,2:33]) 
## Note: Exclude columns 97:98 (assuming 96 barcodes were used) which provide total barcode UMI counts for each cell. 

pdf(paste0(outFile, "_bc.check.pdf"))
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
barTable.full <- barTableMain[, 2:33]#as.data.frame(row.names = barTable2$cellID, barTable)
good.bars <- paste("Bar",1:32,sep="")  # NOTE: In this hypothetical example, barcodes 91-96 were not detected
names(barTable.full) = good.bars
barTable <- barTable.full[, good.bars]  # Remove missing bars and summary columns
barTable_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  barTable_sweep.list[[n]] <- classifyCells(barTable, q=q)
  names(barTable_sweep.list)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=barTable_sweep.list)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
ggsave(filename = paste0(outFile, '.quantile.png'), height = 6, width = 8, dpi = 300, device = "png")

## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(barTable, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
barTable <- barTable[-which(rownames(barTable) %in% neg.cells), ]


## Round 2 -----------------------------------------------------------------------------------------------------
barTable_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  barTable_sweep.list[[n]] <- classifyCells(barTable, q=q)
  names(barTable_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results2 <- findThresh(call.list=barTable_sweep.list)
round2.calls <- classifyCells(barTable, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])

## Repeat until all no negative cells remain (usually 3 rounds)...
final.calls <- c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls) <- c(names(round2.calls),neg.cells)



## Perform semi-supervised negative cell reclassification
reclass.cells <- findReclassCells(barTable.full, neg.cells = unique(names(final.calls)[which(final.calls=="Negative")]))
reclass.res <- rescueCells(barTable.full, final.calls, reclass.cells)
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
  geom_point() + xlim(c(nrow(reclass.res)-1,1)) + 
  ylim(c(0,1.05)) +
  geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_vline(xintercept = 16, color="blue",lty=2)
ggsave(filename = paste0(outFile, '.classStab.png'), height = 6, width = 8, dpi = 300, device = "png")

## Finalize negative cell rescue results
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 16) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]
sampleDF=data.frame(good.bars, bc=names(barTableMain)[2:33])
setDT(sampleDF, key="good.bars")
outData=data.frame(barcode.num=final.calls.rescued, cell.Ids = names(final.calls.rescued))
setDT(outData, key="barcode.num")
outData[sampleDF, barcodeSample:=bc, by=.EACHI]

