import time
import os
import numpy as np
import pandas as pd
import datetime
import Bio
from Bio import SeqIO
from functools import partial
import multiprocessing as mp

def counter_fasta(se, fasta_p):
    
    fasta = list(SeqIO.parse(fasta_p, "fasta"))
    
    num = 0
    n_pos = 0
    for i in range(len(fasta)):
        num+=1
        c = fasta[i].seq.count(se)
        if c > 0:
            n_pos+=1
            #print(i, fasta[i].seq.find(se))
    return n_pos


def get_loc_seq(ref, pos, alt):
    if pos <= 10:
        #print(ef[pos-1:pos+19], alt + ref[pos:pos+19])
        return alt + ref[pos:pos+19]
    elif pos >= (len(ref) - 10):
        #print(ref[pos-20:pos], ref[pos-20:pos-1] + alt)
        return ref[pos-20:pos-1] + alt
    else:
        #print(ref[pos-10:pos+10], ref[pos-10:pos-1] + alt + ref[pos:pos+10])
        return ref[pos-10:pos-1] + alt + ref[pos:pos+10]


def check_fastas(fastas_fold, ref, pos, alt):
    #scores = list()
    fastas = os.listdir(fastas_fold)
    subs = get_loc_seq(ref, pos, alt)
        
    #p = counter_fasta(subs, list(map(lambda x: os.path.join(fastas_fold, x), fastas)))
        #print(p)
    #    scores.append(p)
    countpar = partial(counter_fasta, subs)    
        
    pool = mp.Pool(mp.cpu_count())
    scores = list(pool.map(countpar, list(map(lambda x: os.path.join(fastas_fold, x), fastas))))
    
    if max(scores) > 0:
        return 1
    else:
        return 0


def pd_snps(snps):
    df = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
    lines = snps.readlines()
    for line in lines:
        if line[0] != '#':
            attr = line.split('\t')[1:6]
            if float(attr[-1]) >= 100:
                #print(line.split('\t')[0:9])
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


def proc_all(filename, vcf_dir, txt_dir, ref, fastas_fold):
    vcf = open(os.path.join(vcf_dir, filename))
    df = filter_mutations(os.path.join(txt_dir, filename.replace('.flt.vcf', '.txt')), pd_snps(vcf), \
                          ref, fastas_fold)
    return df


t1 = time.time()

fnames = os.listdir('snp_results')

pd.concat(map(lambda x: proc_all(x, 'snp_results', 'txt_results', 'sars-ref.fasta', 'fastas'), \
              fnames)).to_csv(os.path.join('csv_results', ('sars_result' + \
                  str(datetime.datetime.now()).replace('-', '_').replace(':', \
                      '_').replace(' ', '_')[0:19] + '.csv')))   

print(str(time.time() - t1))

