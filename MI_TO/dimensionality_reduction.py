"""
Module for dimensionality reduction.
"""

import gc
from joblib import cpu_count
import numpy as np
import pandas as pd
import scanpy as sc
from umap.umap_ import simplicial_set_embedding, find_ab_params


##


def umap_from_X_conn(X, conn, ncomps=2, metric='euclidean'):
    """
    Wrapper around umap.umap_.simplicial_set_embedding() to create a umap embedding of the 
    feature matrix X using a precomputed fuzzy graph.
    """
    a, b = find_ab_params(1.0, 0.5)
    X_umap, _ = simplicial_set_embedding(
        X, conn, ncomps, 1.0, a, b, 1.0, 5, 200, 'spectral', 
        random_state=np.random.RandomState(0), metric=metric, metric_kwds={},
        densmap=None, densmap_kwds=None, output_dens=None
    )

    return X_umap


##