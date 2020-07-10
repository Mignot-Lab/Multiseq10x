# Multiseq10x
A pre-processing pipeline companion to https://github.com/chris-mcginnis-ucsf/MULTI-seq. This collection of scripts is meant to replace 2 pre-processing functions in the R package deMultiplex.  
>```makeReadtable.py``` is analogous to ```MULTIseq.preProcess``` from deMultiplex but atleast 10x faster and has minimal RAM footprint.  
>```makeBarMatrix.py``` is analogous to ```MULTIseq.align``` from deMultiplex but atleast 10x faster and uses hamming distance only if there is no sample barcode match.  

## Preprocessing multi-seq FASTQ
The `makeReadtable.py` requires 4 arguments.  
1. The `-C` should be provided a text file of cell id barcodes with single cell barcode per line.  
2. The `-R1` is the multi-seq fraction Fastq pair 1. (can be compressed in gz format).  
3. The `-R2` is the multi-seq fraction Fastq pair 2. (can be compressed in gz format).  
4. The `-O` is the output filename for the read table.  
A typical command would be
usage `makeReadtable.py -C cellIds.txt -R1 R1.fastq.gz -R2 R2.fastq.gz -O readTable.csv`. 