"""
Heatmaps.
"""

import re
import os
import sys
import scanpy as sc
import anndata
from scipy.stats import zscore
from Cellula._utils import Timer
from Cellula.plotting._plotting import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *
from MI_TO.preprocessing import *
from MI_TO.heatmaps_plots import *
matplotlib.use('MacOSX')


##


# Read data
path_main = sys.argv[1]
sample = sys.argv[2]

path_main = '/Users/IEO5505/Desktop/MI_TO/'
sample = 'MDA'

path_results = path_main + 'results_and_plots/hclust/'
path_dists = path_main + f'results_and_plots/distances/'
path_top3 = path_main + 'results_and_plots/top3_analysis/'

# Set sample colors
samples = ['MDA', 'AML', 'PDX']
sample_colors = {'MDA':'#DA5700', 'AML':'#0074DA', 'PDX':'#0F9221'}

# Read AFM(s)
if sample != 'all':

    if sample in samples:
        afm = read_one_sample(path_main, sample=sample)
        options = [ x.split('_')[2:5] for x in os.listdir(f'{path_top3}/{sample}') ]
        options_ = [ x.split('_')[1:4] for x in os.listdir(path_dists) ]

        with PdfPages(path_top3 + f'/{sample}.pdf') as pdf:

            for opt in options_:

                filtering, min_cell_number, min_cov_treshold = opt

                if filtering != 'density':
                    afm = filter_cells_coverage(afm, mean_coverage=int(min_cov_treshold))

                    if int(min_cell_number) > 0:
                        cell_counts = afm.obs.groupby('GBC').size()
                        clones_to_retain = cell_counts[cell_counts>min_cell_number].index 
                        cells_to_retain = afm.obs.query('GBC in @clones_to_retain').index
                        afm = afm[cells_to_retain, :].copy()

                    if filtering == 'miller2022':
                        a = filter_miller2022(afm)
                        a = nans_as_zeros(a)
                        clone_colors = create_palette(a.obs, 'GBC', palette='dark')
                        cell_anno_clones = [ clone_colors[clone] for clone in a.obs['GBC'] ]
                        # Viz cells x vars
                        g = cells_vars_heatmap(a, cell_anno=cell_anno_clones, anno_colors=clone_colors, 
                            heat_label='AF', legend_label='Clone', figsize=(11, 8), title=f'{sample} clones'
                        )
                        pdf.savefig() 

                    for x in os.listdir(path_dists):

                        if bool(re.search('_'.join([sample] + opt), x)):
                            print(x)
                            metric = x.split('_')[5].capitalize()
                            D = sc.read(path_dists + x)
                            assert (D.obs_names == a.obs_names).all()

                            g = cell_cell_dists_heatmap(D, 
                                cell_anno=cell_anno_clones, anno_colors=clone_colors, 
                                heat_label=f'Similarity (1-{metric})', legend_label='Clone', figsize=(11, 6.5), 
                                title=f'{sample} clones'
                            )
                            pdf.savefig() 
            plt.close()






            



# else:
#     afm = read_all_samples(path_main, sample_list=samples)
# 
# # Retrieve top analyses, for each sample "clones" task
# if sample != 'all':
# else:
#     print('Find top... Or plot all!')


##


################################################################

# Cell_vars 

# One sample

# Filter, convert nans, create clone-colors
a = filter_velten2021(afm)
a = nans_as_zeros(a)
clone_colors = create_palette(afm.obs, 'GBC', palette='dark')
cell_anno_clones = [ clone_colors[clone] for clone in afm.obs['GBC'] ]

# Viz 
g = cells_vars_heatmap(a, cell_anno=cell_anno_clones, anno_colors=clone_colors, 
    heat_label='AF', legend_label='Clone', figsize=(11, 8), title=f'{sample} clones'
)
g.savefig(path_results + 'cell_vars/single_sample.pdf')


##


# All samples

# Filter, convert nans, create col-colors
a = filter_velten2021(afm)
a = nans_as_zeros(a)
cell_anno_samples = [ sample_colors[sample] for sample in afm.obs['sample'] ]

# Viz 
g = cells_vars_heatmap(a, cell_anno=cell_anno_samples, anno_colors=sample_colors, 
    heat_label='AF', legend_label='Sample', figsize=(11, 8), title='All samples'
)
g.savefig(path_results + 'cell_vars/all_samples.pdf')


##


################################################################

# Dists

# Read one: example
D = sc.read(path_dists + 'MDA_miller2022_0_30_False_cosine_no_kernel.h5ad')

## NBB: to fix this
cell_anno_clones = [ clone_colors[clone] for clone in a.loc[D.obs_names]['GBC'] ]

# Viz 
g = cell_cell_dists_heatmap(D, cell_anno=cell_anno_clones, anno_colors=sample_colors, 
    heat_label='Similarity', legend_label='Clone', figsize=(11, 6.5), 
    title=f'{sample} clones'
)

g.savefig(path_results + 'dists/single_sample.pdf')




