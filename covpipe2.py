import argparse
import mil
import m2f
import os
import cutpr as cut
import shutil
import multiprocessing as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import mimuts4 as mut 
import offreps2 as rep


def proc_s(seqs, prefix):
    recs = list()
    zero = True
    ind_fa = prefix+'_ind_fa'
    os.mkdir(ind_fa)
    for record in seqs:
        record.seq = Seq(str(record.seq[21608:25272]).strip('-n').upper().replace('-', ''))
        if len(record.id.split('_')) > 1:
            record.description = record.description.split('_')[1].split('-')[0].split()[-1]
            record.id = record.description
            if zero:
                zero = False
            else:
                recs.append(record)
                SeqIO.write(record, os.path.join(ind_fa, record.id+'.fasta'), "fasta")
    return recs


def fq2m(pth_files, prefix):
    analyze_p = prefix+'_analyze' # names of folders
    mpileups_p = prefix+'_mpileups'
    vcfs_p = prefix+'_vcfs'
    
    for name in [analyze_p, mpileups_p, vcfs_p]: # making new directories
        if not os.path.exists(os.path.join(os.getcwd(), name)):
            os.mkdir(name)
    
    for i in range(len(pth_files)):
        shutil.copy(pth_files[i], os.path.join(os.getcwd(), analyze_p))
    #map(lambda x: shutil.copy(x, os.path.join(os.getcwd(), analyze_p)), pth_files) # making copies of files
    
    os.system('gunzip {}/*'.format(os.path.join(os.getcwd(), analyze_p)))

    fastqs = sorted(list(map(lambda x: os.path.join(analyze_p, x), \
        filter(lambda x: x.endswith('.fastq') or x.endswith('.fastq.gz'), os.listdir(analyze_p)))))
    
    pool = mp.Pool(64)
    pool.map(cut_fq, fastqs)
    
    fastqs = list(map(lambda x: x.replace('.fastq', '_cut.fastq'), fastqs))
    
    mil.fqbamsnp(os.path.join('ref', 'sars-ref.fasta'), fastqs, \
        os.path.join(os.getcwd(), mpileups_p), os.path.join(os.getcwd(), vcfs_p)) # processing fastq to mpileups
    
    fname = prefix+'_fastas.fasta'

    SeqIO.write(m2f.f2fa_par(list(map(lambda x: os.path.join(mpileups_p, x), os.listdir(mpileups_p)))), os.path.join("generated_fastas", fname), "fasta") # obtaining 
    # previous fasta for alignment
    os.system('mafft --6merpair --thread -32 --addfragments generated_fastas/{} ref/sars-ref.fasta > {}'.format(\
        fname, os.path.join("alns", prefix+'_aln.fasta'))) # alignment

    alignments = list(SeqIO.parse(os.path.join("alns", prefix+'_aln.fasta'), "fasta")) # getting alignment

    refPoses = mut.ref_poses(alignments[0]) # getting dictionary of all positions

    mutations = mut.get_df(alignments[1:], refPoses) # getting table of mutations
    mutations_strains = rep.get_table(mutations).to_excel(os.path.join('tables', prefix+'_report.xlsx'), index=False) # getting previous table with strain column and saving it


    SeqIO.write(proc_s(alignments, prefix), prefix+'.fasta', "fasta")





if __name__ =='__main__':

    parser = argparse.ArgumentParser(description='Pipeline.')
    parser.add_argument('inputs', type=str, nargs='+', help='Path to fastq files.')
    parser.add_argument('prefix', type=str, help='Prefix.')
    args = parser.parse_args()
    
    
    prm = list(SeqIO.parse("s_primers.fasta", "fasta"))
    def cut_fq(file_p, primers=prm):
        return cut.cutprimers_fq(file_p, primers)
        
    fq2m(args.inputs, args.prefix)
    
