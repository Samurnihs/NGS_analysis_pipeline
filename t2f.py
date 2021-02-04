from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import multiprocessing as mp
import argparse
import datetime
import os


def letter(A, T, G, C):
    d = {A: 'A', T: 'T', G: 'G', C: 'C'}
    mel = max(A, T, G, C)
    return d[mel]


def txt2fa(fname):
    s = ['N']*30000
    
    ch = ''
    dis = ''
    
    f = open(fname, 'r')
    lines = f.readlines()
    
    for line in lines: 
        sample = line.split('\t')
        
        if ch == '':
            ch = sample[0]
        
        if dis == '':
            dis = fname.split('/')[-1].replace('.txt', '')
            
        bases = sample[4].upper()
        
        s[int(sample[1])-1] = letter(bases.count('A'), bases.count('T'), bases.count('G'), bases.count('C'))
        
        f.close()
        
    return SeqRecord(Seq("".join(s).lstrip('N').rstrip('N')), id=ch, description=dis)


parser = argparse.ArgumentParser(description='Mpileup txt file to fasta')
parser.add_argument('inputs', type=str, nargs='+', help='Path to bam files.')
parser.add_argument('-o', type=str, help='Path output folder.', default='out_fastas')
args = parser.parse_args()

files = list(map(str, args.inputs))


if __name__ =='__main__':
    pool = mp.Pool(mp.cpu_count())
    records = list(pool.map(txt2fa, files))
    SeqIO.write(records, os.path.join(args.o, ('result'+str(datetime.datetime.now()).replace('-', '_').replace(':', '_').replace(' ', '_')[0:19]+'.fasta')), "fasta")

