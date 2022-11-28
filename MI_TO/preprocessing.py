"""
Module to preprocess AFMs: reformat original AFM; filter variants/cells.
"""

from timeit import repeat
import numpy as np
import pandas as pd
import scanpy as sc
import anndata


##


def format_matrix(A, cbc_gbc_df):
    """
    Create a full cell x variant AFM from the original maegatk output. 
    Add lentiviral clones' labels to resulting .obs.
    """

    # Add labels to .obs, add A fw counts to layers
    A.obs = A.obs.join(cbc_gbc_df)
    A.obs['GBC'] = pd.Categorical(A.obs['GBC'])
    A.layers['A_counts_fw'] = A.X

    # Create a coverage layer for each feature (each of the possible base per position)
    coverage = np.repeat(A.layers['coverage'].toarray(), repeats=4, axis=1)

    # Create lists for all possible position-base and position-reference combos 
    all_features = [ f'{x}_{base}' for x in range(16569) for base in ['A', 'C', 'T', 'G'] ]
    ref_features = A.var.reset_index().loc[:, ['index', 'refAllele']].agg('_'.join, axis=1).tolist()

    # Initialize AFM
    afm = anndata.AnnData(
        X=np.zeros((len(A.obs_names), len(all_features))),
        obs=A.obs
    )
    afm.var_names = all_features
    afm.layers['coverage'] = coverage

    # For each position and cell, compute each base AF
    A_x = np.array((A.layers['A_counts_fw'] + A.layers['A_counts_rev']) / A.layers['coverage'])
    C_x = np.array((A.layers['C_counts_fw'] + A.layers['C_counts_rev']) / A.layers['coverage'])
    T_x = np.array((A.layers['T_counts_fw'] + A.layers['T_counts_rev']) / A.layers['coverage'])
    G_x = np.array((A.layers['G_counts_fw'] + A.layers['G_counts_rev']) / A.layers['coverage'])

    # Fill afm and store variants names
    variants = []
    j = 0
    for i in range(afm.shape[1]):
        pos_features = all_features[j:j+4]
        for var in pos_features:
            if var.endswith('A'):
                afm[:, var].X = A_x[:, i]
            elif var.endswith('C'):
                afm[:, var].X = C_x[:, i]
            elif var.endswith('T'):
                afm[:, var].X = T_x[:, i]
            elif var.endswith('G'):
                afm[:, var].X = G_x[:, i]
            else:
                print('N found. I am not adding anything')
            if var not in ref_features:
                variants.append(var)
        j += 4
    
    # Subset afm for variants only
    afm_variants = afm[:, variants].copy()
    

    return afm_variants, variants


##


def nans_as_zeros(afm):
    """
    Fill nans with zeros. Technical holes are considered as biologically meaningul, absent quantities.
    """
    X = afm.X
    X_copy = X.copy()
    X_copy[np.isnan(X_copy)] = 0
    afm.X = X_copy

    return afm


##


def filter_CV(afm, n=1000):
    """
    Filter variants based on their coefficient of variation.
    """
    X = afm.X
    CV = X.mean(axis=0) / X.std(axis=0)
    idx_to_filter = np.argsort(CV)[::-1][:n]
    afm = afm[:, idx_to_filter].copy()
    
    return afm


##


def filter_tresholds(afm, median_coverage=100, median_af=0.01, min_perc_cells=0.1):
    """
    Filter variants based on their coefficient of variation.
    """
    X = afm.X
    cov = afm.layers['coverage']

    test_median_af = np.median(X, axis=0) > median_af
    test_min_perc = (np.sum(X>0, axis=0) / X.shape[0]) > min_perc_cells
    test_median_coverage = np.median(cov, axis=0) > median_coverage

    idx_to_filter = test_median_af & test_min_perc & test_median_coverage
    afm = afm[:, idx_to_filter].copy()

    return afm

