import os
import glob
from collections import defaultdict
import gzip
import re

readTableList=sorted(glob.glob('/labs/mignot/10xNMDA/191231_A00351_0308_BH2LW7DSXY_10x/YS-JY-16/filtered_readtables/*ReadTable.csv.gz'))
filteredBarcodeLoc = '/labs/mignot/10xNMDA/191231_A00351_0308_BH2LW7DSXY_10x/cDNA_Aggr_NMDA/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'

def processCellbarcode(barCodefile, readTableList):
    '''takes in processed cell ranger barcode raw file 
    from processed folder /outs/raw_feature_bc_matrix/barcodes.tsv.gz and spits out a set of barcodes
    '''
    cellidFile = defaultdict(set)
    lineN = 0
    for sampleId in gzip.open(barCodefile, 'rb'):
        lineN += 1
        if lineN % 100000 == 0:
            print('PROCESSED LINES {} FROM {}'.format(lineN, barCodefile))
        parseLine=sampleId.decode().strip()
        s=re.sub('-.*', '', parseLine)
        fileId = int(parseLine.split('-')[1])-1
        fileName = readTableList[fileId]
        fileName=fileName.replace('ReadTable.csv.gz', 'cell.barcode')
        cellidFile[fileId].add(fileName)
        if os.path.exists(fileName):
            outFile = open(fileName, 'a')
        else:
            outFile = open(fileName, 'w')
        outFile.write(s+'\n')
    outFile.close()
    return cellidFile

cellidFile=processCellbarcode(filteredBarcodeLoc, readTableList)

def writeCommands(cellidFile, readTableList, sampleBarcodeFile):
    os.chdir('/labs/mignot/10xNMDA/191231_A00351_0308_BH2LW7DSXY_10x/YS-JY-16/filtered_readtables')
    for k, v in cellidFile.items():
        readTableCurrent = readTableList[k]
        cellidFileCurrent = v.pop()
        outFile = readTableCurrent.replace('ReadTable.csv.gz', 'barTable.csv')
        commandFile = open('slurmCommands_{}.sh'.format(k), 'w')
        commandFile.write('#!/bin/bash -l'+'\n')
        commandFile.write('#SBATCH --job-name=BARTABLE_{}'.format(k)+'\n')
        commandFile.write('#SBATCH --mem-per-cpu=16000'+'\n')
        commandFile.write('#SBATCH --time=120:00:00'+'\n')
        commandFile.write('#SBATCH --account=mignot'+'\n')
        commandFile.write('conda activate sseq'+'\n')
        commandFile.write('python /oak/stanford/scg/lab_mignot/10xNMDA/191231_A00351_0308_BH2LW7DSXY_10x/Multiseq10x/scripts/makeBarMatrix.py -C {} -R {} -B {} -O {}'.format(cellidFileCurrent, readTableCurrent, sampleBarcodeFile, outFile))
        commandFile.close()


writeCommands(cellidFile = cellidFile, readTableList=readTableList, sampleBarcodeFile='/labs/mignot/10xNMDA/191231_A00351_0308_BH2LW7DSXY_10x/Multiseq10x/data/multiSeqBarcodes_1_to_32.txt')