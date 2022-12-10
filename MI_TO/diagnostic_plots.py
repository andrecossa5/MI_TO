"""
Utils and plotting functions to visualize and inspect SNVs from a MAESTER experiment and maegatk output.
"""

import os
import scanpy as sc
import anndata
from itertools import product
from Cellula._utils import Timer
from Cellula.plotting._plotting import *
from Cellula.plotting._plotting_base import *
from MI_TO.preprocessing import *


##


# Cell level diagnostics
def sturges(x):   
    return round(1 + 2 * 3.322 * np.log(len(x))) 


##


def cell_depth_dist(afm, ax=None, color='b', title=None):
    """
    Cell depth distribution plot.
    """
    hist(afm.obs, 'depth', n=sturges(afm.obs['depth']), ax=ax, c=color)

    if title is None:
        t = 'MASTER library cell depth (nUMIs)'
    else:
        t = title
    format_ax(afm.obs, ax=ax, title=t, xlabel='total nUMIs', ylabel='n cells')

    median = np.median(afm.obs['depth'])
    r = (afm.obs['depth'].min(), afm.obs['depth'].max())
    ax.text(0.5, 0.9, f'Median: {median:.2f}', transform=ax.transAxes)
    ax.text(0.5, 0.85, f'Min-max: {r[0]:.2f}-{r[1]:.2f}', transform=ax.transAxes)
    ax.text(0.5, 0.8, f'Total cells: {afm.shape[0]}', transform=ax.transAxes)

    return ax


##


def cell_n_sites_covered_dist(afm, ax=None, color='b', title=None):
    """
    n covered positions/total position. Cell distribution.
    """

    df_ = pd.DataFrame(
        np.sum(afm.uns['per_position_coverage'] > 0, axis=1), 
        columns=['n']
    )

    hist(df_, 'n', n=sturges(df_['n']), ax=ax, c=color)

    if title is None:
        t = 'n covered sites (total=16569), per cell'
    else:
        t = title
    format_ax(df_, ax=ax, title=t, xlabel='n sites covered', ylabel='n cells')

    r = (df_['n'].min(), df_['n'].max())
    ax.text(0.05, 0.9, f'Median: {round(np.median(df_["n"]))}', transform=ax.transAxes)
    ax.text(0.05, 0.85, f'Min-max: {r[0]}-{r[1]}', transform=ax.transAxes)
    ax.text(0.05, 0.8, f'Total cells: {afm.shape[0]}', transform=ax.transAxes)

    return ax


##


def cell_n_vars_detected_dist(afm, ax=None, color='b', title=None):
    """
    n covered variants. Cell distribution.
    """
    df_ = pd.DataFrame(np.sum(afm.X > 0, axis=1), columns=['n'])

    hist(df_, 'n', n=sturges(df_['n']), ax=ax, c=color)

    if title is None:
        t = 'n detected variants (total=16569*3), per cell'
    else:
        t = title
    format_ax(df_, ax=ax, title=t, xlabel='n variants detected', ylabel='n cells')

    r = (df_['n'].min(), df_['n'].max())
    ax.text(0.6, 0.9, f'Median: {round(np.median(df_["n"]))}', transform=ax.transAxes)
    ax.text(0.6, 0.85, f'Min-max: {r[0]}-{r[1]}', transform=ax.transAxes)
    ax.text(0.6, 0.8, f'Total cells: {afm.shape[0]}', transform=ax.transAxes)

    return ax


##


def cell_median_site_quality_dist(afm, ax=None, color='b', title=None):
    """
    Median base quality per cell distribution.
    """
    df_ = pd.Series(np.nanmedian(afm.uns['per_position_quality'], axis=1)).to_frame(
        ).rename(columns={0:'qual'})

    hist(df_, 'qual', n=sturges(df_['qual']), ax=ax, c=color)

    if title is None:
        t = 'Median base quality, per cell'
    else:
        t = title
    format_ax(df_, ax=ax, title=t, xlabel='Phred score', ylabel='n cells')

    median = np.median(df_['qual'])
    r = (df_['qual'].min(), df_['qual'].max())
    ax.text(0.6, 0.9, f'Median: {median:.2f}', transform=ax.transAxes)
    ax.text(0.6, 0.85, f'Min-max: {r[0]:.2f}-{r[1]:.2f}', transform=ax.transAxes)
    ax.text(0.6, 0.8, f'Total cells: {afm.shape[0]}', transform=ax.transAxes)

    return ax


##


# Site level diagnostics
def site_median_coverage_dist(afm, ax=None, color='b', title=None):
    """
    Median site coverage across cells, distribution.
    """
    df_ = pd.Series(np.median(afm.uns['per_position_coverage'], axis=0)).to_frame(
        ).rename(columns={0:'cov'})

    hist(df_, 'cov', n=sturges(df_['cov']), ax=ax, c=color)

    if title is None:
        t = 'Median (across cells) MT sites (total=16569) coverage (nUMIs)'
    else:
        t = title
    format_ax(df_, ax=ax, title=t, xlabel='total nUMIs', ylabel='n sites')

    ax.set_xlim((-50, 800))
    median = round(np.median(df_['cov']))
    r = (df_['cov'].min(), df_['cov'].max())
    ax.text(0.5, 0.9, f'Median: {median}', transform=ax.transAxes)
    ax.text(0.5, 0.85, f'Min-max: {r[0]}-{r[1]}', transform=ax.transAxes)

    return ax


##


def site_median_quality_dist(afm, ax=None, color='b', title=None):
    """
    Median site quality across cells, distribution.
    """
    df_ = pd.Series(np.nanmedian(afm.uns['per_position_quality'], axis=0)).to_frame(
        ).rename(columns={0:'qual'})

    hist(df_, 'qual', n=sturges(df_['qual']), ax=ax, c=color)

    if title is None:
        t = 'MASTER library median quality, per site'
    else:
        t = title
    format_ax(df_, ax=ax, title=t, xlabel='Phred score', ylabel='n sites')

    median = np.median(df_['qual'])
    r = (df_['qual'].min(), df_['qual'].max())
    ax.text(0.08, 0.9, f'Median: {median:.2f}', transform=ax.transAxes)
    ax.text(0.08, 0.85, f'Min-max: {r[0]:.2f}-{r[1]:.2f}', transform=ax.transAxes)

    return ax


##


# Variant level diagnostics
def vars_n_positive_dist(afm, ax=None, color='b', title=None):
    """
    Percentage of positive cells per variant, distribution.
    """
    df_ = pd.DataFrame(
        np.sum(afm.X > 0, axis=0), 
        columns=['n']
    )

    hist(df_, 'n', n=sturges(df_['n']), ax=ax, c=color)

    if title is None:
        t = 'n positive cells, per variant'
    else:
        t = title
    format_ax(df_, ax=ax, title=t, xlabel='n cells+', ylabel='n variants')
    
    #ax.set_ylim((0, 1.5))
    r = (df_['n'].min(), df_['n'].max())
    ax.text(0.6, 0.9, f'Median: {round(np.median(df_["n"]))}', transform=ax.transAxes)
    ax.text(0.6, 0.85, f'Min-max: {r[0]}-{r[1]}', transform=ax.transAxes)
    ax.text(0.6, 0.80, f'total cells: {afm.shape[0]}', transform=ax.transAxes)

    return ax


#


def vars_AF_dist(afm, ax=None, color='b', title=None):
    """
    Ranked AF distributions (VG-like).
    """
    to_plot = afm.X.copy()
    to_plot[np.isnan(to_plot)] = 0

    for i in range(to_plot.shape[1]):
        x = to_plot[:, i]
        x = np.sort(x)
        ax.plot(x, '--', color=color, linewidth=0.5)

    if title is None:
        t = 'Ranked AF distributions, per variant'
    else:
        t = title
    format_ax(pd.DataFrame(x), ax=ax, title=t, xlabel='Cell rank', ylabel='AF')

    return ax


##


def strand_concordances(orig, variants):
    """
    Utils for strand concordances calculation.
    """
    strand_conc = {}
    for base in ['A', 'C', 'T', 'G']:
        fw = orig.layers[f'{base}_counts_fw'].toarray()
        rev = orig.layers[f'{base}_counts_rev'].toarray()
        corr = []
        for i in range(fw.shape[1]):
            x = fw[:, i]
            y = rev[:, i]
            corr.append(np.corrcoef(x, y)[0,1])
        corr = np.array(corr)
        corr[np.isnan(corr)] = 0
        strand_conc[base] = corr

    strand_conc = pd.DataFrame(
        strand_conc).reset_index().rename(
        columns={'index':'MT_site'}).melt(
        id_vars='MT_site', var_name='base', value_name='corr'
    )
    strand_conc['MT_site'] = strand_conc['MT_site'].astype(str)
    strand_conc['var'] = strand_conc.loc[:, ['MT_site', 'base']].apply('_'.join, axis=1)
    strand_conc = strand_conc.drop(['MT_site', 'base'], axis=1).set_index('var')
    strand_conc = strand_conc.loc[variants, :]

    return strand_conc


##


def vars_strand_conc_dist(orig, variants, ax=None, color='b', title=None):
    """
    Variant strand concordances distribution.
    """
    orig.layers['A_counts_fw'] = orig.X

    # Calculation of strand concordances, at all {site}_{base} combo.
    strand_conc = strand_concordances(orig, variants)

    # Viz
    hist(strand_conc, 'corr', n=sturges(strand_conc['corr']), ax=ax, c=color)

    if title is None:
        t = 'Strand concordances distribution'
    else:
        t = title
    format_ax(strand_conc, ax=ax, title=t, xlabel="Pearson's r", ylabel='n variants') 

    median = np.median(strand_conc['corr'])
    r = (strand_conc['corr'].min(), strand_conc['corr'].max())
    ax.text(0.6, 0.9, f'Median: {median:.2f}', transform=ax.transAxes)
    ax.text(0.6, 0.85, f'Min-max: {r[0]:.2f}-{r[1]:.2f}', transform=ax.transAxes)

    return ax


##


def AF_mean_strand_conc_corr(orig, variants, afm, ax=None, color='b', title=None):
    """
    Mean AF/strand concordance relationship.
    """
    strand_conc = strand_concordances(orig, variants)

    df_ = strand_conc.assign(mean=np.nanmean(afm.X, axis=0))

    scatter(df_, 'mean', 'corr', c=color, ax=ax)

    if title is None:
        t = 'AF mean-strand concordance trend, per variant'
    else:
        t = title
    format_ax(df_, ax=ax, title=t, xlabel='Mean', ylabel='pho')

    test = df_['mean'] < 0.6
    x = df_['mean'][test]
    y = df_['corr'][test]
    fitted_coefs = np.polyfit(x, y, 1)
    y_hat = np.poly1d(fitted_coefs)(x)
    ax.plot(x, y_hat, linestyle='--', color='black')
    corr = np.corrcoef(x, y)[0,1]
    ax.text(0.7, 0.9, f"Pearson's r: {corr:.2f}", transform=ax.transAxes)

    return ax


##


def AF_mean_var_corr(afm, ax=None, color='b', title=None):
    """
    Mean AF/strand concordance relationship.
    """
    df_ = pd.DataFrame(
        data = np.stack((np.nanmean(afm.X, axis=0), np.nanvar(afm.X, axis=0)), axis=1),
        columns = ['mean', 'variance']
    )

    scatter(df_, 'mean', 'variance', c='b', ax=ax)

    if title is None:
        t = 'AF mean-variance trend'
    else:
        t = title
    format_ax(df_, ax=ax, title=t, xlabel='Mean', ylabel='Variance')

    test = df_['mean'] < 0.6
    x = df_['mean'][test]
    y = df_['variance'][test]
    fitted_coefs = np.polyfit(x, y, 1)
    y_hat = np.poly1d(fitted_coefs)(x)
    ax.plot(x, y_hat, linestyle='--', color='black')
    corr = np.corrcoef(x, y)[0,1]
    ax.text(0.7, 0.9, f"Pearson's r: {corr:.2f}", transform=ax.transAxes)

    return ax


##


def positive_events_by_var_type(afm, orig, ax=None, color=None, title=None):
    """
    % of +events over total + events, stratified by variant type.
    """
    # Compute
    bases = ['A', 'C', 'T', 'G']
    var_types = [ '>'.join(x) for x in product(bases, bases) if x[0] != x[1] ]

    var_type = []
    for x in afm.var_names:
        idx = x.split('_')[0]
        mt_base = x.split('_')[1]
        ref_base = orig.var.loc[idx, 'refAllele']
        t = '>'.join([ref_base, mt_base])
        var_type.append(t)

    afm.var['var_type'] = var_type

    n_positive = {}
    total_positive_events = np.sum(afm.X > 0)
    for x in afm.var['var_type'].unique():
        if not x.startswith('N'):
            test = afm.var['var_type'] == x
            n = np.sum(afm[:, test].X.toarray() > 0) / total_positive_events
            n_positive[x] = n

    # Viz 
    df_ = pd.Series(n_positive).to_frame().rename(
        columns={0:'freq'}).sort_values(by='freq', ascending=False)
    df_['freq'] = df_['freq'].map(lambda x: round(x*100, 2))

    bar(df_, 'freq', ax=ax, s=0.75, c=color, annot_size=8)

    if title is None:
        t = '\n% of + events over total + events'
    else:
        t = title
    format_ax(df_, ax=ax, title=t, xticks=df_.index, ylabel='\n%')

    return ax