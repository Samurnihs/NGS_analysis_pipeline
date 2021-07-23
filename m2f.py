from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import multiprocessing as mp
import argparse
import datetime
import os


def find_pluses(string): # returns all indices of + signs in string
    indices = list()
    for i in range(len(string)):
        if string[i] == '+':
            indices.append(i)
    return indices


def find_minuses(string): # returns all indices of - signs in string
    indices = list()
    for i in range(len(string)):
        if string[i] == '-':
            indices.append(i)
    return indices


def parse_num(string): # parses number in the beginning of string
    i = 1
    while string[0:i].isdigit():
        i+=1
    
    if i > 1:
        return int(string[0:i-1])
    else:
        return -1

def find_ins(string):
    pluses = find_pluses(string)
    nums_pos = list(map(lambda x: x+1, pluses))
    poses_lens = list(zip(nums_pos, list(map(lambda x: parse_num(string[x:x+4]), nums_pos))))
    insertions = list(map(lambda x: string[x[0]+len(str(x[1])):x[0]+len(str(x[1]))+x[1]], poses_lens))
    return insertions


def letter(A, T, G, C, N, emp):
    d = {A: 'A', T: 'T', G: 'G', C: 'C', N: 'N', emp: '-'}
    mel = max(A, T, G, C, N, emp)
    if (mel > 0):
        if mel/(A + T + G + C + N + emp) >= 0.25:
            return d[mel]
        else:
            return 'N'
    else:
        return '-'


def insertion(inses):
    in_con = list()
    if len(inses) > 0:
        for i in range(min(list(map(len, inses)))):
            in_dict = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0,  'other': 0}
            for st in inses:
                if i < (len(st)-1):
                    ch = st[i].upper()
                    if ch in in_dict.keys():
                        in_dict[ch]+=1
            in_con.append(letter(in_dict['A'], in_dict['T'], in_dict['G'], in_dict['C'], in_dict['N'], in_dict['other']))
    return ''.join(in_con).replace('-', '')




def find_dels(string):
    minuses = find_minuses(string)
    #print(minuses)
    lens = list(map(lambda x: parse_num(string[x+1:x+5]), minuses))
    #print(lens)
    #poses_lens = list(zip(nums_pos, list(map(lambda x: parse_num(st[x:x+4]), nums_pos))))
    #deletions = list(map(lambda x: string[x[0]+len(str(x[1])):x[0]+len(str(x[1]))+x[1]], poses_lens))
    out_st = string
    for i in range(1, len(minuses)+1):
        #print((minuses[-i], minuses[-i]+len(str(lens[-i]))+lens[-i]))
        out_st = out_st[0:minuses[-i]] + out_st[minuses[-i]+len(str(lens[-i]))+lens[-i]+1:]
    return out_st


def find_inses(string):
    pluses = find_pluses(string)
    #print(minuses)
    lens = list(map(lambda x: parse_num(string[x+1:x+5]), pluses))
    #print(lens)
    #poses_lens = list(zip(nums_pos, list(map(lambda x: parse_num(st[x:x+4]), nums_pos))))
    #deletions = list(map(lambda x: string[x[0]+len(str(x[1])):x[0]+len(str(x[1]))+x[1]], poses_lens))
    out_st = string
    for i in range(1, len(pluses)+1):
        #print((minuses[-i], minuses[-i]+len(str(lens[-i]))+lens[-i]))
        out_st = out_st[0:pluses[-i]] + out_st[pluses[-i]+len(str(lens[-i]))+lens[-i]+1:]
    return out_st


def txt2fa(fname):
    s = ['N']*30000
    
    ch = ''
    dis = ''
    
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    
    for line in lines: 
        sample = line.split('\t')
        
        if ch == '':
            ch = sample[0]
        
        if dis == '':
            dis = fname.split('/')[-1].replace('.txt', '')
        
        bases = sample[4].upper()
        #bases = find_dels(bases)
        inses = find_ins(bases)
        l_inses = len(inses)
        insert = find_inses(insertion(inses))
        #print(sample[1] + ' ' + insert)
        bases = find_inses(find_dels(bases))
        syms = [bases.count('A'), bases.count('T'), bases.count('G'), bases.count('C'), bases.count('*')]
        if l_inses > 0:
            if l_inses/(l_inses + sum(syms)) >= 0.25:
                ins = insert
            else:
                ins = ''
        else:
            ins = ''
        
        s[int(sample[1])-1] = letter(bases.count('A'), bases.count('T'), bases.count('G'), bases.count('C'), bases.count('N'), bases.count('*')) + ins
        #print(sample[1] + ' ' + s[int(sample[1])-1], len(bases), l_inses)
        
        
    return SeqRecord(Seq("".join(s).strip('N')), id=ch, description=dis)
    #return SeqRecord(Seq("".join(s)), id=ch, description=dis)


def f2fa_par(inputs):#, output, fname=('result'+str(datetime.datetime.now()).replace('-', \
    #    '_').replace(':', '_').replace(' ', '_')[0:19]+'.fasta')):
    pool = mp.Pool(mp.cpu_count())
    records = list(pool.map(txt2fa, inputs))
    #file_path = os.path.join(os.getcwd(), output, fname)
    #SeqIO.write(records, file_path, "fasta")
    return records #file_path


if __name__ =='__main__':

    parser = argparse.ArgumentParser(description='Mpileup txt file to fasta')
    parser.add_argument('inputs', type=str, nargs='+', help='Path to mpileup txt files.')
    parser.add_argument('-o', type=str, help='Path output file.', default=os.path.join('generated_fastas', \
        ('result'+str(datetime.datetime.now()).replace('-', '_').replace(':', '_').replace(' ', '_')[0:19]+'.fasta')))
    args = parser.parse_args()
    
    SeqIO.write(f2fa_par(args.inputs), args.o, "fasta")# f2fa_par(args.inputs, args.o)

