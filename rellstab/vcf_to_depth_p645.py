import gzip
import numpy as np
import scipy.stats

def vcf_to_allele_depth(vcf_file, out_file):
    if vcf_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    with opener(vcf_file) as INFILE:
        with open(out_file, 'w') as OUTFILE:
            header_lines = 0
            for line in INFILE:
                if line.startswith('##'):
                    header_lines += 1
                elif line.startswith('#CHROM'):
                    print 'skipped {} header lines'.format(header_lines)
                    header = line.strip().strip('#').split('\t')
                    inds = header[9:]
                    print 'found {} individuals'.format(len(inds))
                else:
                    tabs = line.split('\t')
                    contig = tabs[0]
                    pos = tabs[1]
                    locus_ID = tabs[2]
                    genotypes = tabs[9:]
                    depth_a_of_ind = dict()
                    depth_b_of_ind = dict()
                    for gen_idx, gen in enumerate(genotypes):
                        #if gen.split(':')[0] in ['1/0', '0/1']: # if het  # orig
                        if gen.split(':')[0] in ['1/0', '0/1'] and gen.split(':')[1] not in ['.']: # if het   # for p316
                            # depth_a_of_ind[inds[gen_idx]] = int(gen.split(':')[2].split(',')[0]) # orig where data is in GT:DP:AD:GL
                            # depth_b_of_ind[inds[gen_idx]] = int(gen.split(':')[2].split(',')[1]) # orig
			    # print 'gen {}'.format(gen)
			    # print int(gen.split(':')[1].split(',')[0]) 
                            depth_a_of_ind[inds[gen_idx]] = int(gen.split(':')[1].split(',')[0])   # for p316 where data is GT:AD:DP:GQ:PL
                            depth_b_of_ind[inds[gen_idx]] = int(gen.split(':')[1].split(',')[1])   # for p316 
                    sum_a = sum(depth_a_of_ind.values())
                    sum_b = sum(depth_b_of_ind.values())
                    num_hets = len(depth_b_of_ind.values())
                    num_samples = len(genotypes)
                    if sum_a+sum_b > 0:
                        OUTFILE.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(contig, pos, locus_ID, sum_a, sum_b, np.float(sum_a)/(sum_a+sum_b), num_hets, num_samples))
                        
                        
def vcf_to_allele_depth_by_ind(vcf_file, out_file):
    if vcf_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    with opener(vcf_file) as INFILE:
        with open(out_file, 'w') as OUTFILE:
            header_lines = 0
            for line in INFILE:
                if line.startswith('##'):
                    header_lines += 1
                elif line.startswith('#CHROM'):
                    print 'skipped {} header lines'.format(header_lines)
                    header = line.strip().strip('#').split('\t')
                    inds = header[9:]
                    print 'found {} individuals'.format(len(inds))
                else:
                    tabs = line.split('\t')
                    contig = tabs[0]
                    pos = tabs[1]
                    locus_ID = tabs[2]
                    genotypes = tabs[9:]
                    depth_a_of_ind = dict()
                    depth_b_of_ind = dict()
                    for gen_idx, gen in enumerate(genotypes):
                        if gen.split(':')[0] in ['1/0', '0/1']: # if het
                            depth_a = int(gen.split(':')[2].split(',')[0])
                            depth_b = int(gen.split(':')[2].split(',')[1])
                            std = scipy.stats.binom(n = depth_a + depth_b, p = .5).std()
                            z = ((depth_a + depth_b)/2. - depth_a)/ -std          
                            OUTFILE.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                    contig, pos, locus_ID, inds[gen_idx], depth_a, depth_b, std, z))
