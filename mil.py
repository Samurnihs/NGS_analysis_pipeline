import os
import argparse

def comm(string): # useful debugging function
    os.system(string)

def single_assay_processing(ref, file1, file2, results, snps): #ref=sfor.fasta minlen 50
    pref1 = file1.split('.')[0] # making prefixes
    pref2 = file2.split('.')[0]
    pre = pref1[: pref1.rfind('R') - 1] # the most common prefix

    clean1 = pref1 + '_qc.fastq.gz' # cleaning files
    clean2 = pref2 + '_qc.fastq.gz'
    comm('bbduk.sh in1={} in2={} out1={} out2={} maq=22 literal=AGATCGGAAGAG,TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG,GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG ref=sfor.fasta ktrim=l k=31 mink=4 hdist=1 tpe tbo minlen=80'.format(file1, file2, clean1, clean2))
    
    #comm('rm {}'.format(file1)) # we already do not need original files
    #comm('rm {}'.format(file2))

    sam = pre + '.sam' # generating .sam
    comm('bwa mem {} {} {} > {}'.format(ref, clean1, clean2, sam))

    bam = pre + '.bam' # generating binary
    comm('samtools view -S -b {} > {}'.format(sam, bam)) 
    
    comm('rm {}'.format(sam)) # sam removal    

    s1 = pre + '_sorted.bam' # collating
    comm('samtools sort {} > {}'.format(bam, s1)) # sorting file +

    comm('rm {}'.format(bam)) # collate removal

    comm('samtools index {}'.format(s1)) # indexing file +

    groups = pre + '_groups.bam' # adding groups
    comm('java -jar /export/home/kotov/picard/picard.jar AddOrReplaceReadGroups I={} O={}  RGLB=lib1  RGPL=illumina  RGPU=unit1 RGSM=ILoveYou'.format(s1, groups)) 
    # +
    
    comm('samtools index {}'.format(groups)) # indexing file +

    vcf = os.path.join(snps, (pre.split('/')[-1] + '.vcf'))
    comm('gatk HaplotypeCaller -R {} -I {} -O {}'.format(ref, groups, vcf)) # making snp-calling +
    
    output = os.path.join(results, os.path.split(pre)[1]) + '.txt' # obtaining all variants of base for every position
    comm('samtools mpileup {} > {}'.format(groups, output)) # + 


def analyze_list(file_list, ref, results, snps): # analyzing list of paired files
    for i in range(len(file_list)//2):
        single_assay_processing(ref, file_list[i * 2], file_list[i * 2 + 1], results, snps)

def analyze_folder(folder, ref, results, snps): # analyzing all files in folder
    file_list = list(map(lambda x: os.path.join(folder, x), os.listdir(folder)))
    analyze_list(file_list, ref, results, snps)

def analyze_folders_inside(path, ref, results, snps): # analyzing all subfolders in folder
    fol_list = list(map(lambda x: os.path.join(path, x), os.listdir(path)))
    for fol in fol_list:
        analyze_folder(fol, ref, results, snps)


def fqbamsnp(ref, inputs, txt_output, snp_output): # complete analysis
    analyze_list(inputs, ref, txt_output, snp_output)

if __name__ =='__main__':
    parser = argparse.ArgumentParser( #parsing arguments
    description='Process pair-reads NGS data from fastq to bam. Includes QC, adapters and duplicates removal.')
    parser.add_argument('ref', type=str, help='Path to reference fasta file. Supporting files shold be at the same place.')
    parser.add_argument('inputs', type=str, nargs='+', help='Path to fastq files.')
    parser.add_argument('-ot', type=str, help='Path to output folder for mpileup.', default=os.path.join(os.getcwd(), 'mpileups'))
    parser.add_argument('-os', type=str, help='Path to output folder for snps.', default=os.path.join(os.getcwd(), 'snps'))
    args = parser.parse_args()

    print(args.ot)

    fqbamsnp(args.ref, sorted(args.inputs), args.ot, args.os)
