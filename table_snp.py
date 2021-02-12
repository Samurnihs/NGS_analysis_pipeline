import time
import os
import numpy as np
import pandas as pd
import datetime
import Bio
from Bio import SeqIO
from functools import partial
import multiprocessing as mp
from functools import partial
import argparse

def seq_in_fasta(se, fasta_p): # finding at least 1 subsequence in fasta
    
    fasta = list(SeqIO.parse(fasta_p, "fasta"))
    
    num = 0
    n_pos = 0
    for i in range(len(fasta)):
        num+=1
        c = fasta[i].seq.count(se)
        if c == 1:
            return 1
    return 0


def get_loc_seq(ref, pos, alt): # getting local subsequence of length 10 with mutation
    if pos <= 10:
        return alt + ref[pos:pos+19]
    elif pos >= (len(ref) - 10):
        return ref[pos-20:pos-1] + alt
    else:
        return ref[pos-10:pos-1] + alt + ref[pos:pos+10]


def check_fastas(fastas_fold, ref, pos, alt): # searching particular mutation in particular position in local database 
    fastas = os.listdir(fastas_fold)
    subseq = get_loc_seq(ref, pos, alt)
        
    
    for path in list(map(lambda x: os.path.join(fastas_fold, x), fastas)):
        if seq_in_fasta(subseq, path) == 1:
            return 1 
    return 0


def pd_snps(snps): # getting pandas DataFrame with snps
    df = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
    lines = snps.readlines()
    for line in lines:
        if line[0] != '#':
            attr = line.split('\t')[1:6]
            if float(attr[-1]) >= 100: # only high quality (>=100)
                df.loc[len(df)] = line.split('\t')[0:9]
    return df


# filtration  of all positions with several bases
def filter_mutations(file_in_name, positions, refer, fastas_fold, read_thresh=50, freq_thresh=1):

    columns = ['File', 'Position', 'Reads','Reference', 'Alternative', 'A', 'T', 'G', 'C', \
               'Max_base', 'GISAIDloc'] # head of dataframe
    df = pd.DataFrame(columns=columns) # empty dataframe

    fin = open(file_in_name, 'r') # opening file obtained after samtools mpileup
    lines = fin.readlines() # reading lines of the file
    
    for line in lines: # counting different bases for every position 
        sample = line.split('\t')
        if sample[1] in list(positions['POS']):
            
            bases = sample[4].lower() # all letters to lower register
            num_A = bases.count('a')
            num_T = bases.count('t')
            num_G = bases.count('g')
            num_C = bases.count('c')
            num_bases = num_A + num_T + num_G + num_C
            
            ref = str(list(positions[positions['POS'] == sample[1]]['REF'])[0])
            alt = str(list(positions[positions['POS'] == sample[1]]['ALT'])[0])
            
            reference = list(SeqIO.parse(refer, "fasta"))[0].seq
            
            gisaidloc = check_fastas(fastas_fold, reference, int(sample[1]), alt)
            
            order = sorted((num_A, num_T, num_G, num_C))
            max_b, second_b, third_b = order[-1], order[-2], order[-3] # getting 3 highest scores for bases

            # we are interested only in high-reliability data with coverage more than read_thresh reads
            if (num_bases >= read_thresh): 
                max_freq = max_b/num_bases # getting maximum frequency
                if (max_freq <= freq_thresh):
                    df.loc[len(df)] = [file_in_name.split('/')[-1].split('.')[0], sample[1], \
                                    num_bases, ref, alt, num_A/num_bases, num_T/num_bases, \
                                       num_G/num_bases, num_C/num_bases, max_b/num_bases, gisaidloc]

    fin.close()
    
    print('{} Successfully processed!'.format(file_in_name))
    
    #fout.close()
    df.sort_values(by=['Position', 'File'])
    #df.drop(labels=['Max'], axis=1)
    return df


def proc_all(filename, ref, fastas_fold): # using previous function in context of files
    vcf = open(filename)
    st = os.path.join(os.getcwd(), 'mpileups', filename.replace('.flt.vcf', '.txt').split('/')[-1])
    print(st)
    df = filter_mutations(st, pd_snps(vcf), \
                          ref, fastas_fold)
    return df


def procpar(filenames, ref, fastas_fold): # parallel adaptation for previous function
    t1 = time.time()
    pool = mp.Pool(mp.cpu_count())
    procfix = partial(proc_all, ref=ref, fastas_fold=fastas_fold)
    print(str(time.time() - t1))
    return pd.concat(pool.map(procfix, filenames))


if __name__ =='__main__':
    parser = argparse.ArgumentParser(description='Getting information about mutations.')
    parser.add_argument('ref', type=str, help='Path to reference fasta file.')
    parser.add_argument('vcfs', type=str, nargs='+', help='Path to vcf files.')
    parser.add_argument('-o', type=str, help='Path output folder.', default='csv_results')
    parser.add_argument('-db', type=str, help='Fasta database to search mutations.', default='fastas')
    args = parser.parse_args()


    df = procpar(args.vcfs, args.ref, args.db) # getting result as DataFrame
    df.to_csv(os.path.join(args.o, ('sars_result' + \
                  str(datetime.datetime.now()).replace('-', '_').replace(':', \
                      '_').replace(' ', '_')[0:19] + '.csv'))) # saving result as csv
