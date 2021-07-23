import Bio
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio import AlignIO
import multiprocessing as mp
import numpy as np
import pandas as pd
import argparse

# list of all mutations and deletions
mutations = [(23063, 'a', 't', 'N501Y'), # N501Y
             (23271, 'c', 'a', 'A570D'), # A570D
             (22813, 'g', 't', 'K417N'), # K417N 
             (23012, 'g', 'a', 'E484K'), # E484K 
             (22812, 'a', 'c', 'K417T'), # K417T 
             (21614, 'c', 't', 'L18F'), # L18F
             (21801, 'a', 'c', 'D80A'), # D80A
             (22206, 'a', 'g', 'D215G'), # D215G
             (23664, 'c', 't', 'A701V'), # A701V
             (21621, 'c', 'a', 'T20N'), # T20N
             (21638, 'c', 't', 'P26S'), # P26S
             (21974, 'g', 't', 'D138Y'), # D138Y
             (22132, 'g', 't', 'R190S'), # R190S
             (23525, 'c', 't', 'H655Y'), # H655Y
             (24642, 'c', 't', 'T1027I'), # T1027I
             (23042,'t', 'c', 'S494P'), # S494P
             (22992, 'g', 'a', 'S477N'), # S477N
             (23012, 'g', 'c', 'E484Q'), # E484Q
             (22917, 't', 'g', 'L452R'), # L452R
             (21618, 'c', 'g', 'T19R'), # T19R
             (22995, 'c', 'a', 'T478K'), # T478K
             (23604, 'c', 'g', 'P681R') # P681R
             ]

dels = [(21765, 21770, 'del HV69-70'), # del HV69-70
        (21991, 21993, 'del Y144'), # del Y144
        (22028, 22036, 'del ER156-158'), # del ER156-158
        (21967, 21993, 'del CY136-144') # del CY136-144
        ] 


inses = [(23598, 'aggggatagcac', 'INS23598')] # The only unsertion can be detected now


def comm(x, y): # returns 1 if variables are equal
    if x == y:
        return 1
    else:
        return 0


def r_co_st(a, b): # returns number of matches between 2 strings
    r = 0
    l = min(len(a), len(b))
    for i in range(min(len(a), len(b))):
        r+=comm(a[i], b[i])
    return r / l


# function that finds out presence of mutations and deletions
def find_muts_dels_single_seq(seq, refPoses, muts, dels, inses):
    mlist = dict()
    bases = ['a', 't', 'g', 'c']
    for m in muts: # going through all mutations 
        if seq[refPoses[m[0] - 1]] == m[2]:
            mlist[m[-1]] = ['ДА']
        elif (seq[refPoses[m[0] - 1]] == m[1]):
            mlist[m[-1]] = ['НЕТ']
        elif (seq[refPoses[m[0] - 1]] in bases):
            mlist[m[-1]] = ['иное']
        else:
            mlist[m[-1]] = ['не определено']

    for d in dels: # going through all deletions 
        isempty = True # emptiness of all positions of deletion
        for s in seq[refPoses[(d[0]-1)]:refPoses[d[1]]]:
            if s in bases: # thereis no emptiness if there is at least one base in place of deletion
                isempty = False
                mlist[d[-1]] = ['НЕТ']
                break
        if isempty: # if there is empty positions before or after deletions, then this place is just wasn't sequenced correctly,
            if (seq[refPoses[d[0]-2]] not in bases) or (seq[refPoses[d[1]+1]] not in bases): # so we cannot determine presence of deletion
                mlist[d[-1]] = ['не определено'] 
            else:
                mlist[d[-1]] = ['ДА']

    for ins in inses: # going through all insertions 
        if r_co_st(seq[refPoses[ins[0] - 1]:refPoses[ins[0] - 1]+len(ins[1])], ins[1]) >= 0.5:
            mlist[ins[-1]] = ['ДА']
        else:
            mlist[ins[-1]] = ['НЕТ']

    return mlist # returns information about mutations and deletions


def get_string(rec, refPoses, muts, dels, inses): # returns one-row-DataFrame for one SeqRecord
    cols = ['Sequence', 'N501Y', 'A570D', 'K417N', 'E484K', 'K417T', 
        'del HV69-70', 'del Y144', 
        'L18F', 'D80A', 'D215G', 'A701V', 'T20N', 'P26S', 'D138Y', 'R190S', 'H655Y', 'T1027I', 'S494P', 'S477N']
    df = pd.DataFrame()#columns=cols)
    df = pd.concat([df, pd.DataFrame({**{'Sequence': rec.description.split('_S')[0].split(' ')[-1].split('-')[0]}, **find_muts_dels_single_seq(rec.seq, refPoses, muts, dels, inses)})])
    return df


def get_df(aln, refPoses, muts=mutations, dels=dels, inses=inses): # returns pandas.DataFrame with results for all sequences
    df = pd.concat(map(lambda x: get_string(x, refPoses, muts, dels, inses), aln))
    return df


def no_ins_ind(ref, index):
    a = ref[:index].count('-')
    b = index
    while ref[b] == '-':
        b+=1
    return a+b


def ref_poses(sec_rec):
    poses = dict()
    for i in range(len(sec_rec.seq)):
        poses[i] = no_ins_ind(str(sec_rec.seq), i)
    return poses



if __name__ =='__main__':
    parser = argparse.ArgumentParser(description='Turns alignment to table of stamm-important mutations.')
    parser.add_argument('input', type=str, help='Path to alignment file.')
    parser.add_argument('output', type=str, help='Path to output excel file.')
    args = parser.parse_args()

    alignments = list(SeqIO.parse(args.input, "fasta"))

    refPoses = ref_poses(alignments[0])

    get_df(alignments[1:], refPoses).to_excel(args.output)
