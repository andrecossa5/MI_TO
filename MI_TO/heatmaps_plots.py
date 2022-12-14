"""
Utils and plotting functions to visualize (clustered and annotated) cells x vars AFM matrices
or cells x cells distances/affinity matrices.
"""

import os
import scanpy as sc
import anndata
from scipy.stats import zscore
from Cellula._utils import Timer
from Cellula.plotting._plotting import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *


##


# Cells x vars AFMs
def cells_vars_heatmap(afm, covariate='GBC', palette_anno='dark', cmap='mako',
    var_names=True, cell_names=False, cluster_cells=True, cluster_vars=True,
    dendrogram_ratio=(.3, .04), colors_ratio=0.05, figsize=(11, 8), var_names_size=5,
    cut_from_right=0.7, halign_title=0.47, position_cbar=(0.82, 0.2, 0.02, 0.25), 
    title_legend='Clones', loc_legend='lower center', bbox_legend=(0.825, 0.5)
    ):
    """
    Given a (filtered) cells x vars AFM, produce its (clustered, or ordered) 
    annotated heatmap visualization.
    """
    # Create annot colors
    colors = create_palette(afm.obs, covariate, palette=palette_anno)
    cells_anno = [ colors[x] for x in afm.obs[covariate] ]

    # Viz 
    g = sns.clustermap(pd.DataFrame(data=afm.X, columns=afm.var_names), 
        cmap=cmap, yticklabels=cell_names, xticklabels=var_names, 
        dendrogram_ratio=dendrogram_ratio, figsize=figsize, row_cluster=cluster_cells, col_cluster=cluster_vars, 
        annot=False, 
        cbar_kws={'use_gridspec' : False, 'orientation' : 'vertical', 'label' : 'AFM'}, 
        colors_ratio=colors_ratio, row_colors=cells_anno
    )

    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=var_names_size)
    g.fig.subplots_adjust(right=cut_from_right)
    g.ax_col_dendrogram.set_visible(False) 
    g.fig.suptitle('MDA_clones', x=halign_title)
    g.ax_cbar.set_position(position_cbar)

    handles = create_handles(colors.keys(), colors=colors.values())
    g.fig.legend(handles, colors.keys(), loc=loc_legend, 
        bbox_to_anchor=bbox_legend, ncol=1, frameon=False, title=title_legend
    )

    return g