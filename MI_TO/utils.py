""""
Miscellaneous utilities.
"""

import numpy as np
import pandas as pd
import scanpy as sc


##


def summary_stats_vars(afm, variants=None):
    """
    Calculate the most important summary stats for a bunch of variants, collected for
    a set of cells.
    """
    if variants is not None:
        test = afm.var_names.isin(variants)
        median_vafs = np.nanmedian(afm[:, test].X, axis=0)
        median_coverage_var = np.nanmedian(afm[:, test].layers['coverage'], axis=0)
        fr_positives = np.sum(afm[:, test].X > 0, axis=0) / afm.shape[0]
        var_names = afm.var_names[test]
    else:
        median_vafs = np.nanmedian(afm.X, axis=0)
        median_coverage_var = np.nanmedian(afm.layers['coverage'], axis=0)
        fr_positives = np.sum(afm.X > 0, axis=0) / afm.shape[0]
        var_names = afm.var_names

    df = pd.DataFrame(
        {
            'median_coverage' : median_coverage_var,
            'median_AF' : median_vafs,
            'fr_positives' : fr_positives
        }, index=var_names
    )

    return df