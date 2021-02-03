import os

def comm(string):
    os.system(string)

def single_assay_processing(ref, file1, file2, results, snps):
    pref1 = file1.split('.')[0] # making prefixes
    pref2 = file2.split('.')[0]
    pre = pref1[: pref1.rfind('R') - 1] # the most common prefix

    clean1 = pref1 + '_qc.fastq.gz' # cleaning files
    clean2 = pref2 + '_qc.fastq.gz'
    comm('bbduk.sh in1={} in2={} out1={} out2={} ftl=15 ftm=15 maq=30 literal=AGATCGGAAGAG,TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG,GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG ktrim=r k=31 mink=4 hdist=1 tpe tbo'.format(file1, \
                file2, clean1, clean2))
    
    #comm('rm {}'.format(file1)) # we already do not need original files
    #comm('rm {}'.format(file2))
    
    #bam = pre + '.bam' # generating binary
    #comm('bbmap.sh in1={} in2={} out={} ref={}'.format(clean1, clean2, bam, ref))

    sam = pre + '.sam' # generating .sam
    comm('bwa mem {} {} {} > {}'.format(ref, clean1, clean2, sam))

    bam = pre + '.bam' # generating binary
    comm('samtools view -S -b {} > {}'.format(sam, bam))
    
    comm('rm {}'.format(sam)) # sam removal

    coll = pre + '_collate.bam' # collating
    comm('samtools collate -o {} {}'.format(coll, bam))

    fixmate = pre + '_fixmate.bam' # fixmate
    comm('samtools fixmate -m {} {}'.format(coll, fixmate))
    
    comm('rm {}'.format(coll)) # collate removal

    positionsort = pre + '_positionsort.bam' # position sotrt
    comm('samtools sort -o {} {}'.format(positionsort, fixmate))
    
    comm('rm {}'.format(fixmate)) # fixmate removal

    markdup = pre + '_markdup.bam' # duplicates removal
    comm('samtools markdup -r {} {}'.format(positionsort, markdup))
    
    comm('rm {}'.format(positionsort)) # positionsort removal

    comm('samtools index {}'.format(markdup)) # indexing file

    sort = pre + '_markdup.sorted.bam' #sorting markedup-file
    comm('samtools sort {} > {}'.format(markdup, sort))
    
    comm('rm {}'.format(markdup)) # markdup removal
    comm('samtools index {}'.format(sort)) # indexing sorted markdup file

    output = os.path.join(results, os.path.split(pre)[1]) + '.txt' # obtaining all variants of base for every position
    comm('samtools mpileup {} > {}'.format(sort, output))

    comm("bcftools mpileup -Ou -f {} {} | bcftools call -Ou -mv | bcftools filter -s LowQual -e '%QUAL<20 || DP>100' > {}.flt.vcf".format(ref, \
		 sort, os.path.join(snps, pre.split('/')[-1])))

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


comm('samtools faidx sars-ref.fasta')
#os.system('bwa index xx0620/sars-ref.fasta')
#single_assay_processing('xx0620/sars-ref.fasta', 'xx0620/Cov_reseq_D186A1184_S1_R1_001.fastq', 'xx0620/Cov_reseq_D186A1184_S1_R1_001.fastq', 'results')
#analyze_folder('xx0620', 'sars-ref.fasta', 'results')
#analyze_folder('covid_files/17_11_2020', 'sars-ref.fasta', 'txt_results', 'snp_results')
analyze_folders_inside('covid_files', 'sars-ref.fasta', 'txt_results', 'snp_results')

