## a script to make a UMI count matrix from fastq
import gzip
import re
import os
import argparse
from multiprocessing import Pool

def processCellbarcode(barCodefile):
    '''takes in processed cell ranger barcode raw file 
    from processed folder /outs/raw_feature_bc_matrix/barcodes.tsv.gz and spits out a set of barcodes
    '''
    cellIDs = set()
    lineBuff = 0
    lineN = 0
    for sampleId in gzip.open(barCodefile, 'rb'):
        lineN += 1
        if lineN == 100000:
            lineBuff += 100000
            lineN =0
            print('PROCESSED LINES {} FROM {}'.format(lineBuff, barCodefile))
        s=re.sub('-.*', '', sampleId.decode().strip())
        cellIDs.add(s)
    return cellIDs

def readTable(R1path, R2path, outReadTable, cellIDs):
    cell, umi, tag= [0, 16], [16, 28], [0, 8]
    cell_start, cell_end = cell
    umi_start, umi_end = umi
    tag_start, tag_end = tag
    R1File = gzip.open(R1path, 'rb')
    R2File = gzip.open(R2path, 'rb')
    readTrue = 0
    outFile = open(outReadTable, 'w')
    outFile.write("Cell,Umi,Sample\n")
    for lineN, (r1, r2) in enumerate(zip(R1File, R2File)):
        if lineN % 1000000 == 0:
            print('{} PROCESSED LINES {}'.format(R1path, lineN))
        if lineN % 4 == 1: # every four lines
            fastaLineR1 = r1.strip()
            fastaLineR2 = r2.strip()
            if fastaLineR1 and fastaLineR2:
                cellNuc=fastaLineR1[cell_start:cell_end].decode()
                if cellNuc in cellIDs:
                    umiNuc=fastaLineR1[umi_start:umi_end].decode()
                    sampleNuc=fastaLineR2[tag_start:tag_end].decode()
                    outFile.write("{},{},{}\n".format(cellNuc,umiNuc,sampleNuc))
                    readTrue += 1
                    # print('CELL is {} \n UMI is {} \n & TAG is {}'.format(cellNuc, umiNuc, sampleNuc))
        lineN += 1
    print('PROCESSED {} FASTA RECORDS AND WROTE {} RECORDS TO {} readTable '.format(lineN//4, readTrue, outReadTable))
    outFile.close()

def preprocess(paths):
    R1path, R2path = paths
    R1path, R2path = os.path.join(multiseq_path, R1path), os.path.join(multiseq_path, R2path)
    outReadTable = os.path.join(multiseq_path, R1path[:-16]+'_ReadTable.csv')
    readTable(R1path, R2path, outReadTable, cellIDs)

if __name__ == "__main__":

    multiseq_path = "/labs/mignot/10xNMDA/191231_A00351_0308_BH2LW7DSXY_10x/YS-JY-16"
    cellIDs = processCellbarcode("/labs/mignot/10xNMDA/191231_A00351_0308_BH2LW7DSXY_10x/10xRUN_cDNA1_3/outs/raw_feature_bc_matrix/barcodes.tsv.gz")
    print("There are {} cellIDs".format(len(cellIDs)))

    # Get all eligible R1 and R2 paths
    r1_paths = []
    r2_paths = []
    for path in os.listdir(multiseq_path):
        if path.startswith("Multi-seq"):
            if path.endswith("R1_001.fastq.gz"):
                r1_paths.append(path)
            elif path.endswith("R2_001.fastq.gz"):
                r2_paths.append(path)
    r1_paths.sort()
    r2_paths.sort()
    
    with Pool(32) as p:
        print(p.map(preprocess, zip(r1_paths, r2_paths)))