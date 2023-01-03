"""
Visualization of clones and samples classification performances.
"""

# Code
import pickle
import re
import os
import sys
import gc
from Cellula.plotting._plotting import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from MI_TO.preprocessing import *
from MI_TO.diagnostic_plots import sturges
from MI_TO.heatmaps_plots import *
from MI_TO.utils import *
from MI_TO.diagnostic_plots import *
#matplotlib.use('macOSX')


##


# Set paths
path_main = sys.argv[1]
overwrite = True if sys.argv[2] == 'over' else False
#path_main = '/Users/IEO5505/Desktop/MI_TO/'
sample_names = ['MDA', 'PDX', 'AML']

path_clones = path_main + '/results_and_plots/clones_classification/'
path_samples = path_main + '/results_and_plots/samples_classification/'
path_results = path_main + '/results_and_plots/classification_performance/'

# Read reports
clones = pd.read_excel(path_clones + 'report_classification_clones.xlsx', index_col=0)
samples = pd.read_excel(path_samples + 'report_classification_samples.xlsx', index_col=0)

# Re-format analysis
clones['analysis'] += '_' + clones['model']
samples['analysis'] += '_' + samples['model']

############## 
# For each sample (3x) clones, what are the top 3 analyses (median f1 score across clones)? 
# Intersection among selected SNVs??

# Save top3 for easy quering
top_3 = {}
for sample in clones['sample'].unique():
    top_3[sample] = clones.query('sample == @sample').groupby(['analysis']).agg(
        {'f1':np.median}).sort_values(
        'f1', ascending=False).index[:3].to_list()

# Load top3 variants for each sample clones, and visualize their intersection (i.e., J.I.), by sample
top3_sample_variants = {}
for sample in sample_names:
        var_dict = {}
        for x in os.listdir(path_results + f'top_3/{sample}/'):
            if x.endswith('.xlsx'):
                n = '_'.join(x.split('.')[0].split('_')[2:-1])
                df_ = pd.read_excel(path_results + f'top_3/{sample}/{x}', index_col=0)
                var_dict[n] = df_.index.unique().to_list()
        top3_sample_variants[sample] = var_dict
##############


##


############## 
# For each sample top3 analysis on the clone task, what are the AF profiles of the variants selected?
# Which relatinship can we visualize among clone cells, using:
# 1) hclustering of cell x var AFM 
# 2) hclustering of a cell x cell similarity matrix?
path_data = path_main + 'data/'
path_distances = path_main + 'results_and_plots/distances/'

# Here we go
for sample in sample_names:

    if not os.path.exists(path_results + 'top_3'):
        os.mkdir(path_results + 'top_3')
    os.chdir(path_results + 'top_3')
    
    if not os.path.exists(sample):
        os.mkdir(sample)
    os.chdir(sample)
            
    # For all top3 analysis of that sample...:
    for analysis in top3_sample_variants[sample]:
        a_ = analysis.split('_')[:-1]
        filtering = a_[0]  
        min_cell_number = int(a_[1])
        min_cov_treshold = int(a_[2])
 
        # Control vars...
        print(analysis)

        # Get info!
        if not os.path.exists(path_results + f'top_3/{sample}/{analysis}/'):
            os.mkdir(path_results + f'top_3/{sample}/{analysis}/')
            os.chdir(path_results + f'top_3/{sample}/{analysis}/')
            for x in os.listdir(path_distances):
                if bool(re.search(f'{sample}_{"_".join(analysis.split("_")[:-1])}_', x)):

                    # Control dists
                    print(x)

                    a_ = x.split('_')[:-1]
                    metric = a_[-1]
                    with_nans = 'w/i nans' if a_[-2] == 'yes' else 'w/o nans'
        else:
            print(f'Analysis {analysis} hclusts have been already computed...')
##############
    
