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
from mquad.mquad import *
from scipy.sparse import coo_matrix
from pegasus.tools.hvf_selection import fit_loess


##
    

def create_one_base_tables(A, base, full=False):
    '''
    Create a full cell x variant AFM from the original maegatk output, and one of the 4 bases,
    create the allelic frequency df for that base, a cell x {pos}_{base} table.
    '''
    cov = A.layers[f'{base}_counts_fw'].toarray() + A.layers[f'{base}_counts_rev'].toarray()
    X = cov / A.layers['coverage'].toarray()
    q = A.layers[f'{base}_qual_fw'] + A.layers[f'{base}_qual_rev']
    m = np.where(A.layers[f'{base}_qual_fw'].toarray() > 0, 1, 0) + np.where(A.layers[f'{base}_qual_rev'].toarray() > 0, 1, 0)
    qual = q.toarray() / m

    df_cov = pd.DataFrame(data=cov, index=A.obs_names, columns=[ f'{pos}_{base}' for pos in range(A.shape[1]) ])
    df_x = pd.DataFrame(data=X, index=A.obs_names, columns=[ f'{pos}_{base}' for pos in range(A.shape[1]) ])
    df_qual = pd.DataFrame(data=qual, index=A.obs_names, columns=[ f'{pos}_{base}' for pos in range(A.shape[1]) ])
    gc.collect()

    if not full:
        test = (A.var['refAllele'] != base).values
        return df_cov.loc[:, test], df_x.loc[:, test], df_qual.loc[:, test]
    else:
        return df_cov, df_x, df_qual


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

    # For each position and cell, compute each base AF and quality tables
    A_cov, A_x, A_qual = create_one_base_tables(A, 'A')
    C_cov, C_x, C_qual = create_one_base_tables(A, 'C')
    T_cov, T_x, T_qual = create_one_base_tables(A, 'T')
    G_cov, G_x, G_qual = create_one_base_tables(A, 'G')

    # Concat all of them in three complete coverage, AF and quality matrices, for each variant from the ref
    cov = pd.concat([A_cov, C_cov, T_cov, G_cov], axis=1)
    X = pd.concat([A_x, C_x, T_x, G_x], axis=1)
    qual = pd.concat([A_cov, C_cov, T_cov, G_cov], axis=1)
    assert (cov.shape[1] % 3 == 1) and (X.shape[1] % 3 == 1) and (qual.shape[1] % 3 == 1) # Check dimensions

    # Reorder columns...
    ref_allele = [ f'{i}_{base}' for i, base in enumerate(A.var['refAllele']) ]
    variants = []
    for i in range(A.shape[1]):
        variants += [ f'{i}_{base}' for base in ['A', 'C', 'T', 'G'] if f'{i}_{base}' not in ref_allele ] 
    cov = cov.loc[:, variants]
    X = X.loc[:, variants]
    qual = qual.loc[:, variants]

    # Create the per position quality matrix
    quality = np.zeros(A.shape)
    n_times = np.zeros(A.shape)
    for k in A.layers:
        if bool(re.search('qual', k)):
            quality += A.layers[k].toarray()
            r, c = np.nonzero(A.layers[k].toarray())
            n_times[r, c] += 1
    quality = quality / n_times

    # Create AnnData with variants and sites matrices
    afm = anndata.AnnData(X=X, obs=A.obs)
    afm.layers['coverage'] = cov
    afm.layers['quality'] = qual

    # Per site slots, in 'uns'. Each matrix is a ncells x nsites matrix
    afm.uns['per_position_coverage'] = A.layers['coverage'].toarray()
    afm.uns['per_position_quality'] = quality

    gc.collect()
    
    return afm


##


def read_one_sample(path_main, sample=None):
    """
    Read and format one sample AFM.
    """
    orig = sc.read(path_main + f'data/AFMs/{sample}_afm.h5ad')
    CBC_GBC = pd.read_csv(path_main + f'data/CBC_GBC_cells/CBC_GBC_{sample}.csv', index_col=0)
    afm = format_matrix(orig, CBC_GBC)
    afm.obs = afm.obs.assign(sample=sample)

    return afm


##


def read_all_samples(path_main, sample_list=None):
    """
    Read and format all samples AFMs. 
    """
    ORIG = {}
    for sample in sample_list:
        orig = sc.read(path_main + f'data/AFMs/{sample}_afm.h5ad')
        orig.obs = orig.obs.assign(sample=sample)
        ORIG[sample] = orig
        meta_vars = orig.var
    orig = anndata.concat(ORIG.values(), axis=0)
    orig.var = meta_vars
    del ORIG
    afm = format_matrix(orig, no_clones=True)

    return afm


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
    filtered.uns['per_position_coverage'] = filtered.uns['per_position_coverage'][test_cells, :]

    return filtered


##


def filter_minimal(afm):
    """
    Minimal filter on variants, applied before any other one.
    """
    # Site covered by at least (median) 10 UMIs
    test_vars_considering_site = []
    test_sites = np.median(afm.uns['per_position_coverage'], axis=0) > 10 
    for x in afm.var_names:
        i = int(x.split('_')[0])
        if test_sites[i]:
            test_vars_considering_site.append(True)
        else:
            test_vars_considering_site.append(False)
    test_vars_considering_site = np.array(test_vars_considering_site)
    
    # Variants seen by at least 5 UMIs and with an AFM of at least 0.01 in at least 5 cells
    test_vars_coverage = np.sum(afm.layers['coverage'] > 5, axis=0) > 5
    test_vars_AF = np.sum(afm.X > 0.01, axis=0) > 5

    # Final test
    test = test_vars_considering_site & test_vars_coverage & test_vars_AF
    filtered = afm[:, test].copy()

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


def filter_seurat(afm, nbins=50, n=2000, log=True):
    """
    Filter in a procedure akin to scanpy/seurat.
    """
    # Calc stats
    X = afm.X
    X[np.isnan(X)] = 0
    mean = X.mean(axis=0)
    var = X.var(axis=0)
    dispersion = np.full(X.shape[1], np.nan)
    idx_valid = (mean > 0.0) & (var > 0.0)
    dispersion[idx_valid] = var[idx_valid] / mean[idx_valid] 

    # Optional, put them in the log space
    if log:
        mean = np.log1p(mean) 
        dispersion = np.log(dispersion)
    
    # Bin dispersion values, subtract from each value the bin dispersion mean and divide by the bin std
    df = pd.DataFrame({"log_dispersion": dispersion, "bin": pd.cut(mean, bins=nbins)})
    log_disp_groups = df.groupby("bin")["log_dispersion"]
    log_disp_mean = log_disp_groups.mean()
    log_disp_std = log_disp_groups.std(ddof=1)
    log_disp_zscore = (
        df["log_dispersion"].values - log_disp_mean.loc[df["bin"]].values
    ) / log_disp_std.loc[df["bin"]].values
    log_disp_zscore[np.isnan(log_disp_zscore)] = 0.0

    # Rank, order and slice first n. Subset AFM
    hvf_rank = np.full(X.shape[1], -1, dtype=int)
    ords = np.argsort(log_disp_zscore)[::-1]
    hvf_rank[ords[:n]] = range(n)
    select = np.where(hvf_rank != -1)[0]
    filtered = afm[:, select].copy()
    
    return filtered


##


def filter_pegasus(afm, span=0.02, n=2000):
    """
    Filter in a procedure akin to scanpy/seurat.
    """
    # Get means and vars
    X = afm.X
    X[np.isnan(X)] = 0
    mean = X.mean(axis=0)
    var = X.var(axis=0)

    # Fit var to the mean with loess regression, readjusting the span
    span_value = span
    while True:
        lobj = fit_loess(mean, var, span=span_value, degree=2)
        if lobj is not None:
            break
        span_value += 0.01

    # Create two ranks
    rank1 = np.zeros(mean.size, dtype=int)
    rank2 = np.zeros(mean.size, dtype=int)
    delta = var - lobj.outputs.fitted_values # obs variance - fitted one
    fc = var / lobj.outputs.fitted_values # obs variance / fitted one
    rank1[np.argsort(delta)[::-1]] = range(mean.size) # Rank in desc order
    rank2[np.argsort(fc)[::-1]] = range(mean.size)

    # Rank according to the sum of the two ranks, and filter AFM.
    hvf_rank = rank1 + rank2
    hvf_index = np.zeros(mean.size, dtype=bool)
    hvf_index[np.argsort(hvf_rank)[:n]] = True
    filtered = afm[:, hvf_index].copy()

    return filtered


##


def filter_Mquad(afm, nproc=8, minDP=10, minAD=1, minCell=2, path_=None):
    """
    Filter variants using the Mquad method.
    """
    # Get DP and AD counts matrices
    DP = afm.uns['per_position_coverage'].T
    DP = coo_matrix(DP)

    x = afm.var_names.to_list()
    ad_vars = []
    pdio = 3106
    i = 0
    while i+3 <= afm.shape[1]:
        if i != pdio*3:
            vars_ = x[i:i+3]
            cum_coverage = afm.layers['coverage'][:,i:i+3].sum(axis=0)
            i += 3
        elif i == pdio*3:
            vars_ = x[i:i+4]
            cum_coverage = afm.layers['coverage'][:, i:i+4].sum(axis=0)
            i += 4
        ad_vars.append(vars_[np.argmax(cum_coverage)])
    AD = afm[:, ad_vars].layers['coverage'].T
    AD = coo_matrix(AD)

    # Select variants
    M = Mquad(AD=AD, DP=DP)
    df = M.fit_deltaBIC(out_dir=path_, nproc=nproc, minDP=minDP, beta_mode=False)
    best_ad, best_dp = M.selectInformativeVariants(
        min_cells=minCell, out_dir=path_, tenx_cutoff=None, export_heatmap=False, export_mtx=False
    )
    selected_idx = M.final_df.index.to_list()
    selected_vars = [ ad_vars[i] for i in selected_idx ]
    # Subset matrix
    filtered = afm[:, selected_vars].copy()

    return filtered


##


def filter_DADApy(afm, ):
    """
    Filter using DADApy.
    """

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


##


def filter_cells_and_vars(afm, filtering=None, min_cell_number=None, min_cov_treshold=None, 
    nproc=8, path_=None, n=2000):
    """
    Filter cells and vars from an afm.
    """ 
    if filtering in [
        'CV', 'ludwig2019', 'velten2021', 'miller2022', 'seurat', 'pegasus', 'MQuad', 'DADApy'
    ]:

        # Cells
        a_cells = filter_cells_coverage(afm, mean_coverage=min_cov_treshold)
        if min_cell_number > 0:
            cell_counts = a_cells.obs.groupby('GBC').size()
            clones_to_retain = cell_counts[cell_counts>min_cell_number].index 
            test = a_cells.obs['GBC'].isin(clones_to_retain)
            a_cells = a_cells[a_cells.obs_names[test], :].copy()
            a_cells.uns['per_position_coverage'] = a_cells.uns['per_position_coverage'][test.astype(bool).values, :]
       
        # Variants
        a_cells = filter_minimal(a_cells)
        if filtering == 'CV':
            a = filter_CV(a_cells, n=100)
        elif filtering == 'ludwig2019':
            a = filter_ludwig2019(a_cells, mean_AF=0.5, mean_qual=0.2)
        elif filtering == 'velten2021':
            a = filter_velten2021(a_cells, mean_AF=0.1, min_cell_perc=0.2)
        elif filtering == 'miller2022':
            a = filter_miller2022(a_cells, mean_coverage=100, mean_qual=0.3, perc_1=0.01, perc_99=0.1)
        elif filtering == 'seurat':
            a = filter_seurat(a_cells, nbins=50, n=n, log=True)
        elif filtering == 'pegasus':
            a = filter_pegasus(a_cells, span=0.02, n=n)
        elif filtering == 'MQuad':
            a = filter_Mquad(a_cells, nproc=nproc, minDP=10, minAD=1, minCell=2, path_=path_)
        elif filtering == 'DADApy':
            a = filter_DADApy(a_cells,...)

    elif filtering == 'density':
        a = filter_density(afm, density=0.5, steps=np.Inf)
        if min_cell_number > 0:
            cell_counts = a.obs.groupby('GBC').size()
            clones_to_retain = cell_counts[cell_counts>min_cell_number].index 
            cells_to_retain = a.obs.query('GBC in @clones_to_retain').index
            a = a[cells_to_retain, :].copy()

    return a_cells, a