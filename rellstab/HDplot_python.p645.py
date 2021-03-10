import scipy.stats
from scipy.stats import gaussian_kde
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import statsmodels
import os.path
import sys
#from IPython.core.pylabtools import figsize

#%matplotlib inline




#LOAD VCF FILE AND EXTRACT SEQUENCE COUNTS FROM HETEROZYGOUS INDIVIDUALS
from vcf_to_depth_p645  import vcf_to_allele_depth
#vcf_file = ('Bpendula_Oulu_v2.vcf')
vcf_file = sys.argv[1]
filebase=os.path.splitext(os.path.basename(vcf_file))[0]
#print "", filebase
#depths_file = ('Bpendula_Oulu_v2.depths')
depths_file = (filebase + '.depths')
vcf_to_allele_depth(vcf_file=vcf_file, out_file = depths_file)

#LOAD COUNT DATA INTO DATAFRAME
depths = pd.read_csv(depths_file, sep = '\t', header = None)
depths.columns = ['contig', 'pos', 'locus_ID', 'depth_a' , 'depth_b', 'ratio', 'num_hets', 'num_samples']
depths.head()

#SUM READ COUNTS PER LOCUS
depths['total_depth'] = depths['depth_a'] + depths['depth_b']
depths['depth_per_het'] = depths['total_depth']/[np.float(xx) for xx in depths['num_hets']]
depths.head()

#CALCULATE HETEROZYGOSITY
depths['hetPerc']=depths['num_hets']/depths['num_samples']
depths.head()

#CALCULATE EXPECTED STANDARD DEVIATION BASED ON BINOMIAL DISTRIBUTION
depths['std'] = scipy.stats.binom(n = depths['total_depth'], p = .5).std()
depths.head()

#CALCULATE Z-SCORE BASED ON STANDARD DEVIATION
depths['z'] = -(depths['total_depth']/2. - depths['depth_a'])/ depths['std']
depths.head()

#WRITE OUTPUT FILE CONTAINING DEPTH AND BIAS INFORMATION
depthsBiasFile = filebase + '.depthsBias'
depths.to_csv(depthsBiasFile, sep="\t")

##FRACTION OF 'A' ALLELES
#sum_a = sum(depths['depth_a'])
#sum_b = sum(depths['depth_b'])
#frac_a = np.float(sum_a)/(sum_a + sum_b)
#print sum_a, sum_b, frac_a

#PLOTS
figsize(8,4)

#PLOT READ-RATIO DEVIATION (D) AGAINST HETEROZYGOSITY
figsize(6,6)
plt.ylim(-80, 80)
plt.scatter(depths['hetPerc'], depths['z'], alpha = .1, s=15, color= 'green')
plt.xlabel('H')
plt.ylabel('D')
plt.title('HDplot Results')
#plt.show()
plt.savefig('HDplot_results.pdf')
plt.clf()

# plot read-ratio deviation (D) with density 
figsize(6,6)
x=depths['hetPerc']
y=depths['z']
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
plt.ylim(-80, 80)
plt.scatter(x, y, c=z, s=15, alpha = 0.1, edgecolor='')
plt.xlabel('H')
plt.ylabel('D')
plt.title('HDplot Results; Density')
plt.savefig('HDplot_results_density.pdf')
plt.clf()

#PLOT RATIO AGAINST HETEROZYGOSITY
figsize(6,6)
plt.scatter(depths['hetPerc'], depths['ratio'], alpha = .1, s=15, color= 'green')
plt.xlabel('H')
plt.ylabel('Allele Ratio')
plt.title('Allele Ratio vs. Heterozygosity for All loci')
#plt.show()
plt.savefig('H_vs_AlleleRatio_results.pdf')
plt.clf()

# plot ratio with density
figsize(6,6)
x=depths['hetPerc']
y=depths['ratio']
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
plt.scatter(x, y, c=z, s=15, alpha = 0.1, edgecolor='')
plt.xlabel('H')
plt.ylabel('Allele Ratio')
plt.title('Allele Ratio vs. Heterozygosity for All loci; Density')
plt.savefig('H_vs_AlleleRatio_results_density.pdf')
plt.clf()


