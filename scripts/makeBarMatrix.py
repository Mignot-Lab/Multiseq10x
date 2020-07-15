import gzip
from collections import defaultdict
import argparse
from multiprocessing import Pool
import json

def gzipHandle(fileName):
    if '.gz' in fileName:
        fileOut = gzip.open(fileName, 'rt')
    else:
        fileOut = open(fileName, 'rt')
    return fileOut

def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

def lineCounter(fileName):
    f = gzipHandle(fileName)
    lineCount=sum(bl.count("\n") for bl in blocks(f))
    print('LINE COUNT AT {}'.format(lineCounter))
    return lineCount

def bucketReadTable(readTable, buckets):
    lineC = lineCounter(readTable)
    payLoad = int(lineC/buckets)
    readTableHandle = gzipHandle(readTable)
    readList = []
    chunkId = 1
    outName ='{}.temp_{}'.format(readTable.replace('.gz', ''), chunkId)
    chunkFile = open(outName, 'w')
    for n, line in enumerate(readTableHandle):
        if n > 0:
            chunkFile.write(line) # keep writing 
            if n % payLoad == 0: ## if payload is reached spawn new file
                chunkFile.close()
                print(' BUCKETED {} CHUNKS'.format(chunkId))
                chunkId += 1
                outName = '{}.temp_{}'.format(readTable.replace('.gz', ''), chunkId)
                chunkFile = open(outName, 'w')
            if outName not in readList: # we will need this list to trigger multiprocessing run so return this iterator
                readList.append(outName)
    chunkFile.close()
    return readList

def hamming(seq1, seq2):
	"""
    adapted from Distance package by https://raw.githubusercontent.com/doukremt/distance/master/distance/_simpledists.py
    Compute the Hamming distance between the two sequences `seq1` and `seq2`.
	The Hamming distance is the number of differing items in two ordered
	sequences of the same length. If the sequences submitted do not have the
	same length, an error will be raised.
	"""
	L = len(seq1)
	if L != len(seq2):
		raise ValueError("expected two strings of the same length")
	if L == 0:
		return 0  # equal
	dist = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
	return dist

def tagDist(tag, bcList):
    '''returns 1 mismatch tag closest sample barcode
    example tag = 'GGAGAAGG'
    tagDist(tag, bcList)
    '''
    if tag:
        tagDist=[hamming(tag, i) for i in bcList]
        if min(tagDist) == 1:
            simTagind = tagDist.index(min(tagDist))
            outTag = bcList[simTagind]
            return outTag
        else:
            return None
    else:
        raise ValueError('NO VALID TAG WAS FOUND')


def processCellIds(cellIdFile):
    '''takes in a cell ranger output such as barcodes.tsv.gz
    '''
    cellIds = gzipHandle(cellIdFile)
    cellSet = defaultdict(str)
    for n,cell in enumerate(cellIds):
        cellParse = cell.strip().split('-')[0]
        cellSet[cellParse]= ''
    print('PROCESSED {} CELLS '.format(n))
    cellIds.close()
    return cellSet


def processReadtable(readTable, cellSet, bcList):
    '''this function quantifies the frequency of cellids and sample barcode combinations
    if the sample barcode does not match a hamming distance is called to check 
    if there is 1 nucleotide mismatch and if found will be added to respective count
    '''
    barDict = defaultdict(lambda:defaultdict(int))
    umiDict = defaultdict(set)
    umiCounterDict = defaultdict(int)
    umiTotalDict = defaultdict(int)
    readTableHandle = gzipHandle(readTable)
    print(bcList)
    for n, line in enumerate(readTableHandle):
        if n > 0:
            cell, umi, tag = line.strip().split(',')
            #print(tag)
            umiDict[cell].add('')
            umiTotalDict[cell] += 1
            if cell in cellSet:
                if tag in bcList:
                    if umi not in umiDict.get(cell):
                        barDict[cell][tag] += 1
                        umiDict[cell].add(umi)
                        umiCounterDict[cell] += 1
                else:
                    outTag=tagDist(tag, bcList)
                    if outTag is not None:
                        if umi not in umiDict.get(cell):
                            barDict[cell][tag] += 1
                            umiDict[cell].add(umi)
                            umiCounterDict[cell] += 1
            if n % 10000 == 0:
                print('PROCESSED {} READS '.format(n))
    jsonOut = readTable.replace('csv', 'json')
    with open(jsonOut, 'w') as outDict:
        json.dump(barDict, outDict)
    #return[barDict, umiCounterDict, umiTotalDict]
            
def writeBartable(outFile, bcList, barDict, umiCounterDict, umiTotalDict):
    outFile = open(outFile, 'w')
    outFile.write('cellID,'+','.join(bcList)+',nUMI,nUMI_Total'+'\n')
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


cellIds = 'Multi-seq7_S95_L003_cell.barcode'
cellSet = processCellIds(cellIds)
bcList = [i.strip() for i in open('multiSeqBarcodes_1_to_32.txt')]
readTable = 'test_ReadTable.csv.gz'
readList=bucketReadTable(readTable, buckets=50)
parallelArgs=[(i, cellSet, bcList) for i in readList][:5]

barDict, umiCounterDict, umiTotalDict = processReadtable(readList[0], cellSet, bcList)
with Pool() as pool:
    pool.starmap(processReadtable, parallelArgs)
def main():
    parser = argparse.ArgumentParser(description='A script to make CELL.ID vs Barcode count when a ReadTable is input ')
    parser.add_argument('-C', help='Cell Ids list', required=True)
    parser.add_argument('-R', help='ReadTable as generated by makeReadTable.py', required=True)
    parser.add_argument('-B', help='Sample Barcode file', required=True)
    parser.add_argument('-O', help='Output file name', required=True)
    args=parser.parse_args()
    cellIds = args.C
    bcList = [i.strip() for i in open(args.B)]
    readTable = args.R
    outFile = args.O
    cellSet = processCellIds(cellIds)
    barDict, umiCounterDict, umiTotalDict = processReadtable(readTable, cellSet, bcList)
    writeBartable(outFile, bcList, barDict, umiCounterDict, umiTotalDict)

if __name__ == "__main__":main()
# cellIds = gzip.open('data/barcodes.tsv.gz', 'rb')
# bcList = [i.strip() for i in open('data/multiSeqBarcodes_1_to_32.txt')]
# readTable = gzip.open('data/testCompare.csv.gz')
