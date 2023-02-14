""""
Miscellaneous utilities.
"""

import pickle
import numpy as np
import pandas as pd
import scanpy as sc

from .preprocessing import read_one_sample


##


def summary_stats_vars(afm, variants=None):
    """
    Calculate the most important summary stats for a bunch of variants, collected for
    a set of cells.
    """
    if variants is not None:
        test = afm.var_names.isin(variants)
        density = (~np.isnan(afm[:, test].X)).sum(axis=0) / afm.shape[0]
        median_vafs = np.nanmedian(afm[:, test].X, axis=0)
        median_coverage_var = np.nanmedian(afm[:, test].layers['coverage'], axis=0)
        fr_positives = np.sum(afm[:, test].X > 0, axis=0) / afm.shape[0]
        var_names = afm.var_names[test]
    else:
        density = (~np.isnan(afm.X)).sum(axis=0) / afm.shape[0]
        median_vafs = np.nanmedian(afm.X, axis=0)
        median_coverage_var = np.nanmedian(afm.layers['coverage'], axis=0)
        fr_positives = np.sum(afm.X > 0, axis=0) / afm.shape[0]
        var_names = afm.var_names

    df = pd.DataFrame(
        {   
            'density' : density,
            'median_coverage' : median_coverage_var,
            'median_AF' : median_vafs,
            'fr_positives' : fr_positives
        }, index=var_names
    )

    return df


##


def prep_things_for_umap(top_runs_per_sample, i, solutions, connectivities, path_main=None):
    """
    Utility used in leiden performance viz.
    """
    # Get top solutions
    d_run = top_runs_per_sample.iloc[i, :].to_dict()

    # Prepare ingredients for embs calculations
    s = d_run['sample']
    a = '_'.join(d_run['analysis'].split('_')[1:])

    path_ = path_main + f'results_and_plots/classification_performance/top_3/{s}/{a}/cell_x_var_hclust.pickle'

    with open(path_, 'rb') as f:
        d_cell_x_var = pickle.load(f)

    cells = d_cell_x_var['cells']
    variants = d_cell_x_var['vars']

    afm = read_one_sample(path_main, sample=s)
    X = afm[cells, variants].X.copy()

    conn_name = f'{d_run["analysis"]}_{d_run["with_nans"]}_{d_run["metric"]}_None'
    leiden_pickle_name = f'{d_run["analysis"]}_{d_run["with_nans"]}_{d_run["metric"]}_None|{d_run["k"]}|{d_run["res"]}'

    labels, true_clones, ARI = solutions[s][leiden_pickle_name]
    conn = connectivities[s][conn_name]

    return X, conn, cells, true_clones, labels, ARI, d_run
