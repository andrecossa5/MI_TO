"""
Module to create custom distance function among cell AF profiles.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from numba import njit

from Cellula._utils import Timer

##


@njit(fastmath=True)
def sqeuclidean_nb(x, y):
    """
    Squared euclidean distance.
    """
    result = 0.0
    for i in range(len(x)):
        result += (x[i] - y[i]) ** 2
    return result


##


def sqeuclidean_np(x, y):
    """
    Squared euclidean distance.
    """
    result = np.sum((x - y) ** 2)
    return result


##


@njit(fastmath=True)
def euclidean_nb(x, y):
    """
    Standard euclidean distance.
    """
    result = 0.0
    for i in range(len(x)):
        result += (x[i] - y[i]) ** 2
    return np.sqrt(result)


##


def euclidean_np(x, y):
    """
    Standard euclidean distance.
    """
    result = np.sum((x - y) ** 2)
    return np.sqrt(result)


##


def test_functions(X, fun1, fun2):
    """
    Test numpy vs numba performance. No parallelization here. Nested for loop.
    """
    t = Timer()

    n = X.shape[0]
    D = np.zeros((n, n))

    # Numba
    t.start()
    for i in range(n):
        for j in range(n):
            if j > i:
                x = X[i, :]
                y = X[j, :]
                test_covered = (~np.isnan(x)) & (~np.isnan(y)) # Here only not nans
                x = x[test_covered]
                y = y[test_covered]
                D[i, j] = D[j, i] = fun1(x, y)

    print(f'Elapsed time numba: {t.stop()/60} min.')

    # Numpy
    t.start()
    for i in range(n):
        for j in range(n):
            if j > i:
                x = X[i, :]
                y = X[j, :]
                test_covered = (~np.isnan(x)) & (~np.isnan(y)) # Here only not nans
                x = x[test_covered]
                y = y[test_covered]
                D[i, j] = D[j, i] = fun2(x, y)

    print(f'Elapsed time numpy: {t.stop()/60} min.')


##


def pairwise_dist_serial(X, fun):
    """
    Returns a pairwise distance matrix from a feature matrix X.
    """

    t = Timer()
    t.start()

    n = X.shape[0]
    D = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if j>i:
                x = X[i, :]
                y = X[j, :]
                test_covered = (~np.isnan(x)) & (~np.isnan(y))
                x = x[test_covered]
                y = y[test_covered]
                D[i, j] = D[j, i] = fun(x, y)

    print(f'Elapsed time: {t.stop()/60} min.')

    return D


##
