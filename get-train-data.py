#!/usr/bin/env python

## NOTE! first, source activate evcouplings_env to get python3, AND
## module load seq/sratoolkit/2.8.1

import subprocess

## get read data for genomes in Grad et al. (2016) Journal of Infectious Diseases.
sra_data_fh = open("training-data.csv")

fsize_dict = {}
for l in sra_data_fh:
    l = l.strip()
    if l.startswith('WGS_ID'): ## skip header
        continue
    ldata = l.split(',')
    sra_id = ldata[22]
    mbytes = int(ldata[19])
    fsize_dict[sra_id] = mbytes
    ##print(sra_id)
    ## NOTE: send output to /dev/null so that I don't get 1100 emails about jobs.
    ## write errors to file called train-dump-fails.txt.
    sra_args = ['bsub', '-o', '/dev/null', '-e', './train-dump-fails.txt', '-q', 'short', '-W', '12:0', 'fastq-dump', '--outdir', './training-data', '--gzip', '--skip-technical', '--readids', '--read-filter', 'pass', '--dumpbase', '--split-files', '--clip', sra_id]
    print(' '.join(sra_args))
    subprocess.run(sra_args)
