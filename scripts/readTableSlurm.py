import os 
multiseq_path = "/labs/mignot/10xNMDA/191231_A00351_0308_BH2LW7DSXY_10x/YS-JY-16"
commandpaths = []
for path in os.listdir(multiseq_path):
    if path.startswith("Multi-seq"):
        seqFile = '_'.join(path.split('_')[:3])
        outFileName = multiseq_path+'/filtered_readtables/'+seqFile+'_ReadTable.csv'
        makeBarcode = multiseq_path+'/filtered_readtables/'+seqFile+'_cell.barcode'
        makeR1 = multiseq_path+'/'+seqFile+"_R1_001.fastq.gz"
        makeR2 = multiseq_path+'/'+seqFile+"_R2_001.fastq.gz"
        out= tuple([makeBarcode, makeR1, makeR2, outFileName])
        if out not in commandpaths:
            commandpaths.append(out)


def writeCommandsReadTable(commandpaths):
    os.chdir('/labs/mignot/10xNMDA/191231_A00351_0308_BH2LW7DSXY_10x/YS-JY-16/filtered_readtables')
    for n, item in enumerate(commandpaths):
        C, R1, R2, O = item
        commandFile = open('slurmReadTable_{}.sh'.format(n), 'w')
        commandFile.write('#!/bin/bash -l'+'\n')
        commandFile.write('#SBATCH --job-name=ReadTable_{}'.format(n)+'\n')
        commandFile.write('#SBATCH --mem-per-cpu=16000'+'\n')
        commandFile.write('#SBATCH --time=120:00:00'+'\n')
        commandFile.write('#SBATCH --account=mignot'+'\n')
        commandFile.write('conda activate sseq'+'\n')
        commandFile.write('python /oak/stanford/scg/lab_mignot/10xNMDA/191231_A00351_0308_BH2LW7DSXY_10x/Multiseq10x/scripts/makeReadtable.py -C {} -R1 {} -R2 {} -O {}'.format(C, R1, R2, O))
        commandFile.close()

writeCommandsReadTable(commandpaths)