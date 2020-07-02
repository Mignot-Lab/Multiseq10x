require(data.table);require(tidyverse);require(xlsx)
install.packages('stringdist')
require(stringdist)
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
bclistParse
readTable = fread('~/Documents/NMDA10x/test.csv')
readTable

tag.dists=stringdistmatrix(a = , b = bclistParse)
tag.dists
tag.hds <- apply(tag.dists, 1, min)
tag.ind <- which(tag.hds <= 1)
tag.dists <- tag.dists[tag.ind, ]


