import os
import pandas as pd
import numpy as np
import pysam
import argparse
import multiprocessing as mp

def get_cov(bam, reg, contig='NC_045512.2'): #returns the characteristic coverage for region
    samfile = pysam.AlignmentFile(bam, "rb") #parsing bam file
    # counts coverages for A C G T for each position 
    covs = samfile.count_coverage(contig, start=reg[0], stop=reg[1])
    # summary coverage
    sum_covs = np.array(covs[0]) + np.array(covs[1]) + np.array(covs[2]) + np.array(covs[3])
    # getting median coverage as characteristic
    median_coverage = int(np.median(sum_covs))
    return median_coverage


def coverfiles(files, regions): # calculates coverage 
    cols = ['FILE'] + list(map(lambda x: 'COV'+str(x+1), range(len(regions)))) 
    cov = pd.DataFrame(columns=cols) # preparing pandas DataFrame for results
    
    for file in files:
        quan = [0]*len(regions)
        for i in range(len(regions)): # obtains coverage for all regions in file
            quan[i] = get_cov(file, regions[i])
        cov.loc[len(cov)] = [file] + quan
        print('{} is finished!'.format(file))
    
    cov['FILE'] = cov['FILE'].apply(lambda x: os.path.split(x)[1].split('_')[0].split('-')[0])
    cov = cov.sort_values(by=['FILE'])
    return cov


if __name__ =='__main__':

    parser = argparse.ArgumentParser(description='Obtains coverage for amplicons in bam files from SARS-CoV-2 NGS-sequencing.')
    parser.add_argument('inputs', type=str, nargs='+', help='Path to bam files.')
    parser.add_argument('output', type=str, help='Path output excel file (.xlsx).')
    parser.add_argument('pattern', type=str, help='Amplicons pattern.', choices= ['five', 'sgen'])
    args = parser.parse_args()
    
    patterns = {'five': [(21598, 21837), 
			  (21899, 22217), 
			  (22811, 23063), 
			  (23271, 23438), 
			  (23518, 23743)],
    		  'sgen': [(21609, 21689), 
			   (21747, 21857), 
			   (21913, 22074), 
			   (22129, 22248), 
			   (22317, 22468),
			   (22529, 22682),
			   (22756, 22850),
			   (22903, 22977),
			   (23136, 23194),
			   (23260, 23372),
			   (23468, 23608),
			   (23663, 23744),
			   (23881, 23935),
			   (23997, 24155),
			   (24221, 24400),
			   (24450, 24537),
			   (24625, 24733),
			   (24821, 24926),
			   (24988, 25097),
			   (25173, 25272)
			   ]}
    
    files = list(filter(lambda x: x.endswith('.bam'), args.inputs))
    
    def cov(files, reg=patterns[args.pattern]):
        return coverfiles(files, reg)
    
    pool = mp.Pool(32)
    df = pd.concat(pool.map(cov, files))
    df.to_excel(args.output)

