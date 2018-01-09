#!/usr/bin/env python

## NOTE! first, source activate evcouplings_env to get python3, AND
## module load seq/breseq/0.30.0

import subprocess
import os

## first filter all mutations in gd files except for SNPs. 
##run gdtools APPLY -f FASTA

def get_sra_ids(sra_data_f):
    if sra_data_f not in ["./validation-data.csv","./training-data.csv"]:
        print("Error: file name does not match")
        quit()
    sra_ids = []
    sra_data_fh = open(sra_data_f)
    if sra_data_f == "./validation-data.csv":
        for l in sra_data_fh:
            l = l.strip()
            if not l.startswith('UK') and not l.startswith('Canada'):
                continue
            ldata = l.split(',')
            sra_id = ldata[2]
            sra_ids.append(sra_id)
    elif sra_data_f == "./training-data.csv":
        sra_data_fh = open(sra_data_f)
        for l in sra_data_fh:
            l = l.strip()
            if l.startswith('WGS_ID'): ## skip header
                continue
            ldata = l.split(',')
            sra_id = ldata[22]
            sra_ids.append(sra_id)
    return sra_ids

def filter_for_snps(sra_data_f, breseq_outdir,filtered_outdir):
    sra_ids = get_sra_ids(sra_data_f)
    for sra_id in sra_ids:
        gdloc = os.path.join(breseq_outdir,sra_id, 'output/output.gd')
        filtered_gdout = os.path.join(filtered_outdir,sra_id+'_filtered.gd')
        with open(filtered_gdout, 'w') as filter_handle:
            with open(gdloc) as gd_handle:
                for line in gd_handle:
                    if line.startswith('#') or line.startswith('SNP'):
                        filter_handle.write(line)

def gd_to_fasta(sra_data_f, filtered_outdir,fasta_outdir):
    sra_ids = get_sra_ids(sra_data_f)
    for sra_id in sra_ids:
        filtered_gd = os.path.join(filtered_outdir,sra_id+'_filtered.gd')
        fastaout = os.path.join(fasta_outdir,sra_id+'.fasta')
        breseq_args = ['bsub', '-o', '/dev/null', '-q', 'short', '-W', '12:0', 'gdtools','APPLY', '-o', fastaout, '-f', 'FASTA', '-r', 'gonorrhoeae-FA1090.gb.txt', filtered_gd]
        print(' '.join(breseq_args))
        subprocess.run(breseq_args)

def main():
    filter_for_snps('./validation-data.csv', './breseq-validation','./filtered-validation-gds')
    filter_for_snps('./training-data.csv', './breseq-training','./filtered-training-gds')
    gd_to_fasta('./validation-data.csv', './filtered-validation-gds', './fasta-validation')
    gd_to_fasta('./training-data.csv', './filtered-training-gds', './fasta-training')

main()
