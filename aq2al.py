import argparse 
import fq2snp as f2s # import pipelines
import m2f
import table_snp as ts 
import os
import time
import pandas
import datetime

parser = argparse.ArgumentParser(description='Complete pipeline from fastq files (paired reads) to fasta alignment.')
parser.add_argument('inputs', type=str, nargs='+', help='Path to fq files.')
parser.add_argument('ref', type=str, help='Path to reference fasta file.')
args = parser.parse_args()

f2s.fqbamsnp(args.ref, args.inputs, 'mpileups', 'snps') # obtaining mpileups and snps

mpileups=list(map(lambda x: os.path.join('mpileups', x), os.listdir('mpileups')))

m2f.f2fa_par(mpileups, 'generated_fastas') # obtaining fasta files from pileups

out_f = os.listdir('generated_fastas') # list of fasta files in directory generated_fastas

for path in list(map(lambda x: os.path.join('generated_fastas', x), out_f)): # making alignments
    os.system('mafft --6merpair --thread -4 --addfragments {} {} > {}'.format(path, args.ref, \
        path.replace('generated_fastas', 'alignments').replace('.fasta', '_al.fasta')))

snps = list(map(lambda x: os.path.join('snps', x), os.listdir('snps')))
df = ts.procpar(snps, args.ref, 'fastas') # table snps
df.to_csv(os.path.join('csv_results', ('sars_result' + \
                str(datetime.datetime.now()).replace('-', '_').replace(':', \
                    '_').replace(' ', '_')[0:19] + '.csv'))) # saving as csv
