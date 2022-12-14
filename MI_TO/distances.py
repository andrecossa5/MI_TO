"""
Module to create custom distance function among cell AF profiles.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from numba import njit

from Cellula._utils import Timer


##


def sqeuclidean_(x, y):
    """
    Squared euclidean distance.
    """
    result = np.sum((x - y) ** 2)
    return result


##


def euclidean(x, y):
    """
    Standard euclidean distance.
    """
    result = np.sum((x - y) ** 2)
    return np.sqrt(result)


##

