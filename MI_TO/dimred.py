"""
Dimensionality reduction utils to compress (pre-filtered) AFMs.
"""

from joblib import cpu_count
import numpy as np
import pandas as pd
from umap.umap_ import UMAP
from sklearn.metrics import pairwise_distances
from scipy.linalg import eigh
from sklearn.preprocessing import StandardScaler

from Cellula.preprocessing._pp import my_PCA


##


def find_diffusion_matrix(X=None, alpha=5):
    """
    Function to find the diffusion matrix P.
    """

    dists = pairwise_distances(X)
    K = np.exp(-dists**2 / alpha)
    
    r = np.sum(K, axis=0)
    Di = np.diag(1/r)
    P = np.matmul(Di, K)
    
    D_right = np.diag((r)**0.5)
    D_left = np.diag((r)**-0.5)
    P_prime = np.matmul(D_right, np.matmul(P,D_left))

    return P_prime, P, Di, K, D_left


##


def find_diffusion_map(P_prime, D_left, n_eign=3):
    """
    Function to find the diffusion coordinates in the diffusion space.
    """   
    eigenValues, eigenVectors = eigh(P_prime)
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    diffusion_coordinates = np.matmul(D_left, eigenVectors)
    
    return diffusion_coordinates[:,:n_eign]


##


def reduce_dimensions(afm, method='PCA', metric='euclidean', n_comps=30, alpha=0.01, sqrt=False, scale=True):
    """
    Util to create dimension-reduced representation of the input SNVs AFM.
    """
    # Get AFM np.array
    X = afm.X

    # Sqrt and scale, optionally
    if sqrt:
        X = np.sqrt(X)
    if scale:
        scaler = StandardScaler()
        X = scaler.fit_transform(X)
    
    # Reduce
    if method == 'PCA':
        PCA = my_PCA()
        PCA.calculate_PCA(X, n_components=n_comps)
        X_reduced = PCA.embs
        feature_names = [ f'PC{i}' for i in range(1, X_reduced.shape[1]+1)]

    elif method == 'UMAP':
        umap = UMAP(n_components=n_comps, metric=metric, random_state=1234)
        X_reduced = umap.fit_transform(X)
        feature_names = [ f'UMAP{i}' for i in range(1, X_reduced.shape[1]+1)]

    elif method == 'diffmap':
        P_prime, P, Di, K, D_left = find_diffusion_matrix(X, alpha=alpha)
        X_reduced = find_diffusion_map(P_prime, D_left, n_eign=n_comps)
        feature_names = [ f'Diff{i}' for i in range(1, X_reduced.shape[1]+1)]

    return X_reduced, feature_names

