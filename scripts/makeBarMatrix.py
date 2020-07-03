import gzip
from collections import defaultdict
import distance

def tagDist(tag, bcList):
    '''returns 1 mismatch tag closest sample barcode
    example tag = 'GGAGAAGG'
    tagDist(tag, bcList)
    '''
    if tag:
        tagDist=[distance.hamming(tag, i) for i in bcList]
        if min(tagDist) == 1:
            simTagind = tagDist.index(min(tagDist))
            outTag = bcList[simTagind]
            return outTag
        else:
            return None


cellIds = gzip.open('data/barcodes.tsv.gz', 'rb')
bcList = [i.strip() for i in open('data/multiSeqBarcodes_1_to_32.txt')]
readTable = gzip.open('data/Multi-seq1_S89_L003_ReadTable.csv.gz')

cellSet = defaultdict(str)
for cell in cellIds:
    cellParse = cell.decode().strip().split('-')[0]
    cellSet[cellParse]= ''
cellIds.close()

barDict = defaultdict(lambda:defaultdict(int))
umiDict = defaultdict(set)
umiCounterDict = defaultdict(int)
umiTotalDict = defaultdict(int)

for n, line in enumerate(readTable):
    if n > 0:
        cell, umi, tag = line.decode().strip().split(',')
        umiDict[cell].add('')
        umiTotalDict[cell] += 1
        if cell in cellSet:
            if tag in bcList:
                if umi not in umiDict.get(cell):
                    barDict[cell][tag] += 1
                    umiDict[cell].add(umi)
                    umiCounterDict[cell] += 1
            else:
                #continue
                outTag=tagDist(tag, bcList)
                if outTag is not None:
                    if umi not in umiDict.get(cell):
                        barDict[cell][tag] += 1
                        umiDict[cell].add(umi)
                        umiCounterDict[cell] += 1
        if n % 10000 == 0:
            print('PROCESSED {} READS '.format(n))
            
outFile = open('outs/barTable.csv', 'w')
outFile.write('cellID,'+','.join(bcList)+',nUMI,nUMITotal'+'\n')
for k, v in barDict.items():
    outString =k
    for bc in bcList:
        if bc in v:
            outString += ','+str(v.get(bc))
        else:
            outString += ','+'0'
    if k in umiCounterDict:
        outString += ','+str(umiCounterDict.get(k))
    else:
        outString += ','+'0'
    if k in umiTotalDict:
        outString += ','+str(umiTotalDict.get(k))
    else:
        outString += ','+'0'
    outFile.write(outString +'\n')
outFile.close()




