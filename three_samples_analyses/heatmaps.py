"""
Heatmaps.
"""

import os
import scanpy as sc
import anndata
from scipy.stats import zscore
from Cellula._utils import Timer
from Cellula.plotting._plotting import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *
from MI_TO.preprocessing import *
matplotlib.use('MacOSX')


##


# Read data
path_main = '/Users/IEO5505/Desktop/MI_TO/'
path_results = path_main + 'results_and_plots/three_samples/exploratory/'

ORIG = {}
AFMs = {}

samples = ['MDA_clones', 'AML_clones', 'PDX']
for x in samples:
    orig = sc.read(path_main + f'data/AFMs/{x}_afm.h5ad')
    CBC_GBC = pd.read_csv(path_main + f'data/CBC_GBC_cells/CBC_GBC_{x}.csv', index_col=0)
    afm, variants = format_matrix(orig, CBC_GBC)
    afm.obs = afm.obs.assign(sample=x)
    AFMs[x] = afm
    ORIG[x] = orig

# Set colors
sample_colors = {'MDA_clones':'#DA5700', 'AML_clones':'#0074DA', 'PDX':'#0F9221'}


##


################################################################

# Heatmaps for MT-SNVs clusters/clones vizualization
afm = AFMs['MDA_clones']

# Filter, convert nans, create col-colors
a = filter_miller2022(afm)
a = nans_as_zeros(a)
clone_colors = create_palette(afm.obs, 'GBC', palette='dark')
cell_anno_clones = [ clone_colors[clone] for clone in afm.obs['GBC'] ]

# Viz 
g = sns.clustermap(pd.DataFrame(data=a.X, columns=a.var_names), 
    cmap='viridis', yticklabels=False, xticklabels=True, 
    dendrogram_ratio=(.3, .04), figsize=(11, 8), row_cluster=True, col_cluster=True, 
    annot=False, 
    cbar_kws={'use_gridspec' : False, 'orientation' : 'vertical', 'label' : 'AFM'}, 
    colors_ratio=0.05, row_colors=cell_anno_clones
)

g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=5)
g.fig.subplots_adjust(right=0.7)
g.ax_col_dendrogram.set_visible(False) 
g.fig.suptitle('MDA_clones', x=0.47)
g.ax_cbar.set_position((0.82, 0.2, 0.02, 0.25))

handles = create_handles(clone_colors.keys(), colors=clone_colors.values())
g.fig.legend(handles, clone_colors.keys(), loc='lower center', 
    bbox_to_anchor=(0.825, 0.5), ncol=1, frameon=False, title='Clones'
)
g.savefig('/Users/IEO5505/Desktop/a.pdf')


##


# Samples
del ORIG

ORIG = {}
samples = ['MDA_clones', 'AML_clones', 'PDX']
for x in samples:
    orig = sc.read(path_main + f'data/AFMs/{x}_afm.h5ad')
    orig.obs = orig.obs.assign(sample=x)
    ORIG[x] = orig
    meta_vars = orig.var

orig = anndata.concat(ORIG.values(), axis=0)
orig.var = meta_vars

del ORIG

afm, variants = format_matrix(orig, no_clones=True)

# Filter, convert nans, create col-colors
a = filter_ludwig2019(afm)
a = nans_as_zeros(a)
cell_anno_samples = [ sample_colors[sample] for sample in afm.obs['sample'] ]

# Viz 
g = sns.clustermap(pd.DataFrame(data=a.X, columns=a.var_names), 
    cmap='viridis', yticklabels=False, xticklabels=True, 
    dendrogram_ratio=(.3, .04), figsize=(11, 8), row_cluster=True, col_cluster=True, 
    annot=False, 
    cbar_kws={'use_gridspec' : False, 'orientation' : 'vertical', 'label' : 'AFM'}, 
    colors_ratio=0.05, row_colors=cell_anno_samples
)

g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=5)
g.fig.subplots_adjust(right=0.7)
g.ax_col_dendrogram.set_visible(False) 
g.fig.suptitle('All samples', x=0.47)
g.ax_cbar.set_position((0.82, 0.2, 0.02, 0.25))

handles = create_handles(sample_colors.keys(), colors=sample_colors.values())
g.fig.legend(handles, sample_colors.keys(), loc='lower center', 
    bbox_to_anchor=(0.825, 0.5), ncol=1, frameon=False, title='Samples'
)

plt.show()
g.savefig('/Users/IEO5505/Desktop/samples.pdf')






