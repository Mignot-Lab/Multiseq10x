## a script to make a UMI count matrix from fastq
import gzip
import re
import os
import argparse
from makeBarMatrix import gzipHandle

def processCellbarcode(barCodefile):
    '''takes in processed cell ranger barcode raw file 
    from processed folder /outs/raw_feature_bc_matrix/barcodes.tsv.gz and spits out a set of barcodes
    '''
    cellIDs = set()
    lineN = 0
    for sampleId in gzipHandle(barCodefile):
        lineN += 1
        if lineN % 100000 == 0:
            print('PROCESSED LINES {} FROM {}'.format(lineN, barCodefile))
        s=re.sub('-.*', '', sampleId.decode().strip())
        cellIDs.add(s)
    return cellIDs

def readTable(R1path, R2path, outReadTable, cellIDs):
    cell, umi, tag= [0, 16], [16, 28], [0, 8]
    cell_start, cell_end = cell
    umi_start, umi_end = umi
    tag_start, tag_end = tag
    R1File = gzipHandle(R1path)
    R2File = gzipHandle(R2path)
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


def main():
    parser = argparse.ArgumentParser(description='Reads multiseq FASTQ R1 R2 ')
    parser.add_argument('-C', help='Cell Ids list', required=True)
    parser.add_argument('-R1', help='FASTQ R1', required=True)
    parser.add_argument('-R2', help='FASTQ R2', required=True)
    parser.add_argument('-O', help='Output file name', required=True)
    args=parser.parse_args()
    barCodefile = args.C
    R1path = args.R1 
    R2path = args.R2
    outReadTable = args.O
    cellIDs=processCellbarcode(barCodefile)
    readTable(R1path, R2path, outReadTable, cellIDs)

if __name__ == "__main__":main()

