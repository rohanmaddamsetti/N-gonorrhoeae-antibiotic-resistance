#!/usr/bin/env python

''' annotate-rRNA-SNPs.py by Rohan Maddamsetti

- bsub to O2 because CPU usage is high.

 - Get locations of the first (the reference) 
 5S, 16S, 23S rRNAs in FA1090.

 For training or validation data:

 Open rRNA-SNP-validation-frequencies.csv for writing.
 Or rRNA-SNP-training-frequencies.csv.
 For each strain/genome:
  -- read in the SNP-called FASTA file and keep in memory.
  -- get the 5S, 16S, 23S rRNAs for this genome.
  -- filter for rRNA SNPs and keep in memory.
  -- for each rRNA:
       -- for each SNP:
            -- find the position of the genome to modify,
               keeping in mind strand and 0-based indexing.
            -- write a line to rRNA-SNP-frequencies.csv
               with the following attributes:
               strain id (different for training and validation data),
               rRNA type, rRNA position,
               GENOME position, reference state, evolved state, frequency.
            -- modify that position in the genome sequence.
  -- write edited_genome to file. 
 close csv file for writing.
'''

import os
from Bio import SeqIO

def get_sra_ids(sra_data_f):
    if sra_data_f not in ["./validation-data.csv","./training-data.csv"]:
        print("Error: file name does not match")
        quit()
    sra_ids = []
    sra_data_fh = open(sra_data_f)
    if sra_data_f == "./validation-data.csv":
        for l in sra_data_fh:
            l = l.strip()
            if not l.startswith('Brighton') and not l.startswith('Canada'):
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

def get_rRNA_annotation(gbkf):
    ''' returns a list of SeqFeatures '''
    rec = SeqIO.read(gbkf,"genbank")
    locus_tags = ['NGO_r01', 'NGO_r02', 'NGO_r03']
    rRNAlist = [feat for feat in rec.features if feat.type == 'rRNA']
    rRNAlist = [feat for feat in rRNAlist if feat.qualifiers['locus_tag'][0] in locus_tags]
    rRNA_5S = [x for x in rRNAlist if '5S' in x.qualifiers['product'][0]][0]
    rRNA_16S = [x for x in rRNAlist if '16S' in x.qualifiers['product'][0]][0]
    rRNA_23S = [x for x in rRNAlist if '23S' in x.qualifiers['product'][0]][0]
    rRNAdict = {'5S':rRNA_5S,'16S':rRNA_16S, '23S':rRNA_23S}
    return rRNAdict

def get_rRNA_pos(SNP, rRNA_ref_dict):
    '''returns a 1-indexed rRNA position for printing to file. '''
    my_id, my_rRNAtype, my_ref_pos, my_evolstate, my_freq = SNP
    SNP_genome_pos = int(my_ref_pos) - 1
    rRNA_annotation = rRNA_ref_dict[my_rRNAtype]
    if rRNA_annotation.strand == 1:
        rRNA_pos = SNP_genome_pos - int(rRNA_annotation.location.start) + 1
    elif rRNA_annotation.strand == -1: ## closed interval on 'start', open interval on 'end'
        rRNA_pos = int(rRNA_annotation.location.end) - SNP_genome_pos
    assert rRNA_pos > 0
    return str(rRNA_pos)

def get_rRNA_SNPs(sra_id,rRNA_dict,training=True):
    ''' 
    snps are represented as tuples like:
    (ERR191730, 16S, 12345674, T, 0.7654).
    (sra_id, rRNAtype, ref_pos, evolstate, freq    
    '''
    if training == True:
        my_dataset = 'training'
    else:
        my_dataset = 'validation'
    my_datadir = 'breseq-poly-'+my_dataset
    my_gd = os.path.join('.',my_datadir,sra_id,"output","output.gd")
    my_snps = []
    with open(my_gd) as gd_fh:
        for l in gd_fh:
            l = l.strip()
            if l.startswith('SNP'):
                ldata = l.split()
                ref_pos = ldata[4]
                zerorefpos = int(ref_pos) - 1
                if rRNA_dict['5S'].location.start <= zerorefpos and zerorefpos < rRNA_dict['5S'].location.end:
                    rRNAtype = '5S'
                elif rRNA_dict['16S'].location.start <= zerorefpos and zerorefpos < rRNA_dict['16S'].location.end:
                    rRNAtype = '16S'
                elif rRNA_dict['23S'].location.start <= zerorefpos and zerorefpos < rRNA_dict['23S'].location.end:
                    rRNAtype = '23S'
                else:
                    continue
                evolstate = ldata[5]
                freq = ldata[6].split('=')[-1]
                snp = (sra_id, rRNAtype, ref_pos, evolstate, freq)
                my_snps.append(snp)
    return my_snps


def write_edited_fasta_files(is_training=True):
    gbkf = "./gonorrhoeae-FA1090.gb.txt"
    rRNA_ref_dict = get_rRNA_annotation(gbkf)
    csvheader = "ID,rRNA_type,rRNA_position,FA1090_position,consensus_call,polymorphism_call,frequency\n"   

    if is_training:
        freq_file = "./rRNA-SNP-training-frequencies.csv"
        fasta_aln = "./training-with-rRNAs.fasta"
        datafile = "./training-data.csv"
        fasta_dir = "./fasta-training"
    else:
        freq_file = "./rRNA-SNP-validation-frequencies.csv"
        fasta_aln = "./validation-with-rRNAs.fasta"
        datafile = "./validation-data.csv"
        fasta_dir = "./fasta-validation"

    with open(freq_file,"w") as csvhandle, open(fasta_aln,"w") as edited_genomes_h:
        csvhandle.write(csvheader)
        sra_ids = get_sra_ids(datafile)
        for sraid in sra_ids:
            genomef = os.path.join(fasta_dir,sraid+".fasta")
            genome_fasta = SeqIO.read(genomef,"fasta")
            genome_fasta_list = list(genome_fasta.seq)
            rRNA_SNPs = get_rRNA_SNPs(sraid,rRNA_ref_dict,is_training)            
            for SNP in rRNA_SNPs:
                my_id, my_rRNAtype, my_ref_pos, my_evolstate, my_freq = SNP
                assert my_id == sraid
                zero_ref_pos = int(my_ref_pos) - 1
                my_rRNA_pos = get_rRNA_pos(SNP, rRNA_ref_dict)
                ref_state = genome_fasta_list[zero_ref_pos]
                ''' STRINGS ARE IMMUTABLE: operate on a list. '''
                genome_fasta_list[zero_ref_pos] = my_evolstate
                snpfreqstring = ','.join([training_id,my_rRNAtype,my_rRNA_pos,my_ref_pos,ref_state,my_evolstate,my_freq])
                csvhandle.write(snpfreqstring+"\n")
            my_genome_header = '>'+sraid+"\n"
            edited_genomes_h.write(my_genome_header)
            edited_genomes_h.write(''.join(genome_fasta_list)+"\n")

def main(): 
    write_edited_fasta_files(is_training=True)
    write_edited_fasta_files(is_training=False)

main()

