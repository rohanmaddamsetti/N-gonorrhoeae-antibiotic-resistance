#!/usr/bin/env python

## NOTE! first, source activate evcouplings_env to get python3, AND
## module load seq/breseq/0.31.0

import subprocess
import os

'''
run breseq on training and validation data.
use masked FA1090 reference genome.
run breseq in consensus as well as polymorphism modes.
we will use consensus calls, and substitute rRNA polymorphism calls
to generate the alignment and the rRNA SNP frequency files.
'''

def get_training_sra_ids():
    sra_ids = []
    sra_data_fh = open("training-data.csv")
    for l in sra_data_fh:
        l = l.strip()
        if l.startswith('WGS_ID'): ## skip header
            continue
        ldata = l.split(',')
        sra_id = ldata[22]
        sra_ids.append(sra_id)
    return sra_ids

def get_validation_sra_ids():
    sra_ids = []
    sra_data_fh = open("validation-data.csv")
    for l in sra_data_fh:
        l = l.strip()
        ## skip header
        if not l.startswith('UK') and not l.startswith('Canada'):
            continue
        ldata = l.split(',')
        sra_id = ldata[2]
        sra_ids.append(sra_id)
    return sra_ids

def get_outdir(mode,dataset):
    if mode == 'consensus' and dataset == 'training':
        outdir = './breseq-training/'
    elif mode == 'consensus' and dataset == 'validation':
        outdir = './breseq-validation/'
    elif mode == 'polymorphism' and dataset == 'training':
        outdir = './breseq-poly-training/'
    elif mode == 'polymorphism' and dataset == 'validation':
        outdir = './breseq-poly-validation/'
    else:
        print("ERROR: unknown mode")
        quit()
    return outdir

def run_breseq_on_data(mode,dataset,sra_id):
    
    refgenome = 'masked-gonorrhoeae-FA1090.gbk'
    outdir = get_outdir(mode,dataset)
    if dataset == 'training':
        datadir = './training-data/'
    elif dataset == 'validation':
        datadir = './validation-data/'
    else:
        print('ERROR: unknown dataset')
        quit()

    if mode not in ['consensus', 'polymorphism']:
        print('ERROR: unknown mode')
        quit()

    bsub_args = ['bsub', '-o', '/dev/null', '-q', 'medium', '-W', '120:0', 'breseq']

    fulloutdir = os.path.join(outdir,sra_id)
    datafile1 = os.path.join(datadir,sra_id+'_pass_1.fastq.gz')
    datafile2 = os.path.join(datadir,sra_id+'_pass_2.fastq.gz')

    breseq_args =  ['-r', './masked-gonorrhoeae-FA1090.gbk', '-o', outdir+sra_id, datadir+sra_id+'_pass_1.fastq.gz', datadir+sra_id+'_pass_2.fastq.gz']

    if mode == 'polymorphism':
        ''' run in polymorphism mode for rRNA SNP calls '''
        breseq_args = ['-p'] + breseq_args

    breseq_args = bsub_args + breseq_args

    output_done = os.path.join(outdir,sra_id,'output/output.done')    
    if not os.path.exists(output_done):
        print(' '.join(breseq_args))
        ##subprocess.run(breseq_args)
        ## for checking out failed/unfinished runs.
        return sra_id
    else:
        return ''

def delete_unfinished_runs(mode, dataset, unfinished_runs):
    breseq_outdir = get_outdir(mode, dataset)
    for x in unfinished_runs:
        delete_args = ['rm', '-rf', os.path.join(breseq_outdir,x)]
        subprocess.run(delete_args)
    ##print(unfinished_runs)
    print(len(unfinished_runs))
        
def main():
    masked_FA1090 = 'masked-gonorrhoeae-FA1090.gbk'
    ## check out unfinished/failed runs.
    my_unfinished_consensus_training_runs = []
    my_unfinished_polymorphism_training_runs = []
    my_unfinished_consensus_validation_runs = []
    my_unfinished_polymorphism_validation_runs = []

    training_sra_ids = get_training_sra_ids()

    for sra_id in training_sra_ids:
        consensus_training_stat = run_breseq_on_data('consensus','training',sra_id)
        if len(consensus_training_stat):
            my_unfinished_consensus_training_runs.append(consensus_training_stat)

        polymorphism_training_stat = run_breseq_on_data('polymorphism','training',sra_id)
        if len(polymorphism_training_stat):
            my_unfinished_polymorphism_training_runs.append(polymorphism_training_stat)

    validation_sra_ids = get_validation_sra_ids()

    for sra_id in validation_sra_ids:
        consensus_validation_stat = run_breseq_on_data('consensus','validation',sra_id)
        if len(consensus_validation_stat):
            my_unfinished_consensus_validation_runs.append(consensus_validation_stat)

        polymorphism_validation_stat = run_breseq_on_data('polymorphism','validation',sra_id)
        if len(polymorphism_validation_stat):
            my_unfinished_polymorphism_validation_runs.append(polymorphism_validation_stat)


    ## These lines should be commented out unless re-running
    ## failed or unfinished runs.
    delete_unfinished_runs('consensus','training',my_unfinished_consensus_training_runs)
    delete_unfinished_runs('polymorphism','training',my_unfinished_polymorphism_training_runs)
    delete_unfinished_runs('consensus','validation',my_unfinished_consensus_validation_runs)
    delete_unfinished_runs('polymorphism','validation',my_unfinished_polymorphism_validation_runs)

main()
