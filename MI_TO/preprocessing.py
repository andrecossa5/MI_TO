"""
Module to preprocess AFMs: reformat original AFM; filter variants/cells.
"""

import sys
import gc
import re
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import anndata


##


def create_one_base_tables(A, base):
    '''
    Create a full cell x variant AFM from the original maegatk output, and one of the 4 bases,
    create the allelic frequency df for that base, a cell x {pos}_{base} table.
    '''
    X = np.array((A.layers[f'{base}_counts_fw'] + A.layers[f'{base}_counts_rev']) / A.layers['coverage'])
    q = A.layers[f'{base}_qual_fw'] + A.layers[f'{base}_qual_rev']
    m = np.where(A.layers[f'{base}_qual_fw'].toarray() > 0, 1, 0) + np.where(A.layers[f'{base}_qual_rev'].toarray() > 0, 1, 0)
    qual = q.toarray() / m
    df_x = pd.DataFrame(data=X, index=A.obs_names, columns=[ f'{pos}_{base}' for pos in range(A.shape[1]) ])
    df_qual = pd.DataFrame(data=qual, index=A.obs_names, columns=[ f'{pos}_{base}' for pos in range(A.shape[1]) ])
    gc.collect()

    return df_x, df_qual


##


def format_matrix(A, cbc_gbc_df=None, no_clones=False):
    """
    Create a full cell x variant AFM from the original maegatk output. 
    Add lentiviral clones' labels to resulting .obs.
    """

    # Add labels to .obs
    if not no_clones:
        A.obs = A.obs.join(cbc_gbc_df)
        A.obs['GBC'] = pd.Categorical(A.obs['GBC'])
    # add A fw counts to layers
    A.layers['A_counts_fw'] = A.X

    # Create lists for all possible position-base and position-reference combos 
    all_features = [ f'{x}_{base}' for x in range(16569) for base in ['A', 'C', 'T', 'G'] ]
    ref_features = A.var.reset_index().loc[:, ['index', 'refAllele']].agg('_'.join, axis=1).tolist()

    # For each position and cell, compute each base AF and quality tables
    A_x, A_qual = create_one_base_tables(A, 'A')
    C_x, C_qual = create_one_base_tables(A, 'C')
    T_x, T_qual = create_one_base_tables(A, 'T')
    G_x, G_qual = create_one_base_tables(A, 'G')

    # Create a base quality layer, at each site
    quality = np.zeros(A.shape)
    n_times = np.zeros(A.shape)
    for k in A.layers:
        if bool(re.search('qual', k)):
            quality += A.layers[k].toarray()
            r, c = np.nonzero(A.layers[k].toarray())
            n_times[r, c] += 1
    quality = quality / n_times

    # Initialize AFM, its uns and layers slots
    afm = anndata.AnnData(
        X=np.zeros((len(A.obs_names), len(all_features))),
        obs=A.obs
    )
    # Add per position coverage and quality
    afm.var_names = all_features
    afm.uns['per_position_coverage'] = A.layers['coverage'].toarray()
    afm.uns['per_position_quality'] = quality # nans matained, NB
    afm.layers['quality'] = np.zeros((len(A.obs_names), len(all_features)))

    # Fill afm and store variants names
    variants = []
    for i, x in enumerate(afm.var_names):
        if x.endswith('A'):
            afm[:, x].X = A_x.loc[:, x].values
            afm.layers['quality'][:, i] = A_qual.loc[:, x].values
        elif x.endswith('C'):
            afm[:, x].X = C_x.loc[:, x].values
            afm.layers['quality'][:, i] = C_qual.loc[:, x].values
        elif x.endswith('T'):
            afm[:, x].X = T_x.loc[:, x].values
            afm.layers['quality'][:, i] = T_qual.loc[:, x].values
        elif x.endswith('G'):
            afm[:, x].X = G_x.loc[:, x].values
            afm.layers['quality'][:, i] = G_qual.loc[:, x].values
        else:
            print('N found. I am not adding anything')
        if x not in ref_features:
            variants.append(x)
    
    # Subset afm for variants only
    afm_variants = afm[:, variants].copy()
    gc.collect()
    
    return afm_variants, variants


##


def nans_as_zeros(afm):
    """
    Fill nans with zeros. Technical holes are considered as biologically meaningul zeroes.
    """
    X_copy = afm.X.copy()
    X_copy[np.isnan(X_copy)] = 0
    afm.X = X_copy

    return afm


##


def filter_cells_coverage(afm, mean_coverage=100):
    """
    Simple filter to subset an AFM only for cells with at least n mean site coverage. 
    """
    test_cells = np.mean(afm.uns['per_position_coverage'], axis=1) > mean_coverage 
    filtered = afm[test_cells, :].copy()
    return filtered


##


def filter_CV(afm, n=100):
    """
    Filter variants based on their coefficient of variation.
    """
    CV = np.nanmean(afm.X, axis=0) / np.nanvar(afm.X, axis=0)
    variants = np.argsort(CV)[::-1][:n]
    filtered = afm[:, variants].copy()
    return filtered


##


def filter_ludwig2019(afm, mean_AF=0.5, mean_qual=0.2):
    """
    Filter variants based on fixed tresholds adopted in Ludwig et al., 2019, 
    in the xperiment without ATAC-seq reference.
    """
    test_vars_het = np.nanmean(afm.X, axis=0) > mean_AF # highly heteroplasmic variants
    test_vars_qual = np.nanmean(afm.layers['quality'], axis=0) > mean_qual # high quality vars
    test_vars = test_vars_het & test_vars_qual
    filtered = afm[:, test_vars].copy()
    return filtered


##


def filter_velten2021(afm, mean_AF=0.1, min_cell_perc=0.2):
    """
    Filter variants based on fixed tresholds adopted in Ludwig et al., 2021.
    """
    test_vars_considering_site = []
    test_sites = np.sum(afm.uns['per_position_coverage'] > 5, axis=0) > 20
    for x in afm.var_names:
        i = int(x.split('_')[0])
        if test_sites[i]:
            test_vars_considering_site.append(True)
        else:
            test_vars_considering_site.append(False)
    test_vars_considering_site = np.array(test_vars_considering_site)
    test_vars_het = np.nanmean(afm.X, axis=0) > mean_AF
    test_vars_exp = (np.sum(afm.X > 0, axis=0) / afm.shape[0]) > min_cell_perc
    test_vars = test_vars_considering_site & test_vars_het & test_vars_exp
    filtered = afm[:, test_vars].copy()

    return filtered

##


def filter_miller2022(afm, mean_coverage=100, mean_qual=0.3, perc_1=0.01, perc_99=0.1):
    """
    Filter variants based on adaptive adopted in Miller et al., 2022.
    """
    test_vars_considering_site = []
    test_sites = np.mean(afm.uns['per_position_coverage'], axis=0) > mean_coverage
    for x in afm.var_names:
        i = int(x.split('_')[0])
        if test_sites[i]:
            test_vars_considering_site.append(True)
        else:
            test_vars_considering_site.append(False)
    test_vars_considering_site = np.array(test_vars_considering_site)
    test_vars_qual = np.nanmean(afm.layers['quality'], axis=0) > mean_qual
    test_vars_het = (np.percentile(afm.X, q=1, axis=0) < perc_1) & (np.percentile(afm.X, q=99, axis=0) > perc_99)
    test_vars = test_vars_considering_site & test_vars_qual & test_vars_het
    filtered = afm[:, test_vars].copy()

    return filtered


##


def filter_density(afm, density=0.5, steps=np.Inf):
    """
    Jointly filter cells and variants based on the iterative filtering algorithm adopted by Moravec et al., 2022.
    """
    # Get AF matrix, convert into a df
    logger = logging.getLogger("my_logger")
    X_bool = np.where(~np.isnan(afm.X), 1, 0)

    # Check initial density not already above the target one
    d0 = X_bool.sum() / X_bool.size
    if d0 >= density:
        logger.info(f'Density is already more than the desired target: {d0}')
        return afm
        
    else:
        print(f'Initial density: {d0}')

    # Iteratively remove lowest density rows/cols, until desired density is reached
    densities = []
    i = 0
    while i < steps:

        print(f'Step {i}:')
        rowsums = X_bool.sum(axis=1)
        colsums = X_bool.sum(axis=0)
        d = X_bool.sum() / X_bool.size
        densities.append(d)
        print(f'Density: {d}')

        if d >= density or (len(rowsums) == 0 or len(colsums) == 0):
            break

        rowmin = rowsums.min()
        colmin = colsums.min()
        if rowmin <= colmin:
            lowest_density_rows = np.where(rowsums == rowmin)[0]
            X_bool = X_bool[ [ i for i in range(X_bool.shape[0]) if i not in lowest_density_rows ], :]
            afm = afm[ [ i for i in range(afm.shape[0]) if i not in lowest_density_rows ], :].copy()
        else:
            lowest_density_cols = np.where(colsums == colmin)[0]
            X_bool = X_bool[:, [ i for i in range(X_bool.shape[1]) if i not in lowest_density_cols ] ]
            afm = afm[:, [ i for i in range(afm.shape[1]) if i not in lowest_density_cols ] ].copy()
        i += 1
    
    return afm