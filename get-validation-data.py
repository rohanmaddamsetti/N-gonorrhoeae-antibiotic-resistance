#!/usr/bin/env python

## NOTE! first, source activate evcouplings_env to get python3, AND
## module load seq/sratoolkit/2.8.1

import subprocess

## get read data for Brighton genomes
sra_data_fh = open("validation-data.csv")
for l in sra_data_fh:
    l = l.strip()
    ## skip header
    if not l.startswith('UK') and not l.startswith('Canada'):
        continue
    ldata = l.split(',')
    sra_id = ldata[2]
    sra_args = ['bsub', '-q', 'short', '-W', '12:0', 'fastq-dump', '--gzip', '--outdir', './validation-data',
                '--skip-technical', '--readids', '--read-filter', 'pass', '--dumpbase', '--split-files', '--clip', sra_id]
    print(' '.join(sra_args))
    subprocess.run(sra_args)
