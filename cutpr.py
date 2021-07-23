import os
import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp
import pandas as pd
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from distutils.core import setup
from Cython.Build import cythonize
import c2
import argparse


def comp(sym): #returns complementary base
    pairs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    sym = sym.upper()
    if sym in pairs.keys():
        return pairs[sym]
    else:
        return sym.upper()


def rev_comp(string): # returns reverse compliment for sequence
    return ''.join(list(map(comp, string)))[:: -1]


def cutprimers_fq(file_p, primers):
    
    f_out = open(file_p.replace('.fastq', '_cut.fastq'), 'w')
    file = open(file_p, 'r')
    lines = file.readlines()
    
    thr_len = 50
    thr_mis = 0.5
    x_b = 'zzzzzzzzzzzzzzz'
    
    rc_primers = list()
    for i in range(len(primers)):
        rec = primers[i]
        rec.seq = rev_comp(primers[i].seq)
        rc_primers.append(rec)
    
    for i in range(len(lines)//4):
        n = lines[4*i].strip()
        s = lines[4*i+1].strip()
        q = lines[4*i+3].strip()
    
        if len(s) >= thr_len:
            for j in range(len(primers)):
                
                q1 = c2.mismatch(str(primers[j].seq).encode(), (x_b+s[:35]).encode())
                q2 = c2.mismatch(str(rc_primers[j].seq).encode(), (s[-35:]+x_b).encode())

                if (q1[1] >= thr_mis):
                    s = s[(q1[0] + len(primers[j].seq) - 5) :]
                    q = q[(q1[0] + len(primers[j].seq) - 5):]

                if (q2[1] >= thr_mis):
                    s = s[: (q2[0]-35)]
                    q = q[: (q2[0]-35)]

            f_out.write(n + '\n' + s + '\n+\n' + q + '\n')
        else:
            f_out.write(n + '\n' + '' + '\n+\n' + '' + '\n')
    file.close()
    f_out.close()
    print('{} is processed!'.format(file_p))
    return file_p.replace('.fastq', '_cut.fastq')


if __name__ =='__main__':
    parser = argparse.ArgumentParser( #parsing arguments
    description='Cuts off all primers.')
    parser.add_argument('inputs', type=str, nargs='+', help='Path to fastq files.')
    args = parser.parse_args()

    prm = list(SeqIO.parse("s_primers.fasta", "fasta"))
    
    
    def cut_fq(file_p, primers=prm):
        return cutprimers_fq(file_p, primers)
    
    
    pool = mp.Pool(32)
    pool.map(cut_fq, args.inputs)




