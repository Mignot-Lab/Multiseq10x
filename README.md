# Multiseq10x
A pre-processing pipeline companion to https://github.com/chris-mcginnis-ucsf/MULTI-seq. This collection of scripts is meant to replace 2 pre-processing functions in the R package deMultiplex.
>```makeReadtable.py``` is analogous to ```MULTIseq.preProcess``` from deMultiplex but atleast 10x faster and has minimal RAM footprint.
>```makeBarMatrix.py``` is analogous to ```MULTIseq.align``` from deMultiplex but atleast 10x faster and uses hamming distance only if there is no sample barcode match.

