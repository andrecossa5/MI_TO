"""
Module to create custom distance function among cell AF profiles.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from numba import njit
from scipy.spatial.distance import pdist

from Cellula._utils import Timer


##


@njit(fastmath=True)
def sqeuclidean_(X):
    """
    Squared euclidean distance.
    """
    # result = np.sum((x - y) ** 2)
    return D


##


@njit(fastmath=True)
def euclidean(X):
    """
    Standard euclidean distance.
    """
    #result = np.sum((x - y) ** 2)
    return D


##

