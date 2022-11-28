import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from numba import njit
from joblib import Parallel, delayed, cpu_count
import leidenalg
import igraph as ig


from Cellula._utils import Timer
from Cellula.preprocessing._neighbors import kNN_graph
from Cellula.preprocessing._metrics import custom_ARI




# Read data
path_data = '/Users/IEO5505/Desktop/MI_TO/data/'
GBC_labels = pd.read_csv(path_data + 'CBC_GBC_cells/CBC_GBC_PDX.csv', index_col=0)
orig = sc.read(path_data + '/AFMs/PDX_afm.h5ad')
orig.obs = orig.obs.join(GBC_labels)
afm.obs['GBC'] = pd.Categorical(afm.obs['GBC'])
orig.layers['A_counts_fw'] = orig.X

# afm = anndata.AnnData(X=np.zeros(afm.shape), obs=afm.obs, var=afm.var)

all_features = [ f'{x}_{base}' for x in range(16569) for base in ['A', 'C', 'T', 'G'] ]
ref_features = orig.var.reset_index().loc[:, ['index', 'refAllele']].agg('_'.join, axis=1).tolist()
cells = orig.obs_names

afm = anndata.AnnData(X=np.zeros((len(cells), len(all_features))), obs=orig.obs)
afm.var_names = all_features

A_x = np.array((orig.layers['A_counts_fw'] + orig.layers['A_counts_rev']) / orig.layers['coverage'])
C_x = np.array((orig.layers['C_counts_fw'] + orig.layers['C_counts_rev']) / orig.layers['coverage'])
T_x = np.array((orig.layers['T_counts_fw'] + orig.layers['T_counts_rev']) / orig.layers['coverage'])
G_x = np.array((orig.layers['G_counts_fw'] + orig.layers['G_counts_rev']) / orig.layers['coverage'])

variants = []
j = 0
for i in range(A_x.shape[1]):
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
            print('N found. not adding anything')
        if var not in ref_features:
            variants.append(var)
    j += 4

X = afm[:, variants].X.toarray()
X.shape

# np.median(np.isnan(X).sum(axis=0) / X.shape[0])

np.sum(~np.isnan(X[0, :])) / X.shape[1]
np.sum(~np.isnan(X[11, :])) / X.shape[1]

x = X[0, :]
y = X[11, :]

test_covered = (~np.isnan(x)) & (~np.isnan(y))
x = x[test_covered]
y = y[test_covered]


@njit(fastmath=True)
def euclidean_numba(x, y):
    """
    Standard euclidean distance.
    """
    result = 0.0
    for i in range(len(x)):
        result += (x[i] - y[i]) ** 2
    return result

def euclidean_numpy(x, y):
    """
    Standard euclidean distance.
    """
    result = np.sum((x - y) ** 2)
    return result


def test_functions(X):

    t = Timer()

    n = X.shape[0]
    D_numba = np.zeros((n, n))
    D_numpy = np.zeros((n, n))

    # Numba
    t.start()
    for i in range(n):
        for j in range(n):
            x = X[i, :]
            y = X[j, :]
            test_covered = (~np.isnan(x)) & (~np.isnan(y))
            x = x[test_covered]
            y = y[test_covered]
            D_numba[i, j] = euclidean_numba(x, y)
    print(f'Elapsed time numba: {t.stop()/60} min.')

    # Numpy
    t.start()
    for i in range(n):
        for j in range(n):
            x = X[i, :]
            y = X[j, :]
            test_covered = (~np.isnan(x)) & (~np.isnan(y))
            x = x[test_covered]
            y = y[test_covered]
            D_numpy[i, j] = euclidean_numpy(x, y)
    print(f'Elapsed time numpy: {t.stop()/60} min.')

    return D_numba, D_numpy


# test_functions(X[:50, :])


def pairwise_dist_euclidean(X):

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
                d = euclidean_numpy(x, y)
                D[i, j] = D[j, i] = d
    print(f'Elapsed time numpy: {t.stop()/60} min.')

    return D


D = pairwise_dist_euclidean(X)

k = 15
for i in range(D.shape[0]):
    NN_idx = np.argsort(D[i, :])[:k]
    others_idx = [ x for x in range(D.shape[0]) if x not in NN_idx ]
    D[i, others_idx] = 0

# With nans as zeros, and UMAP connectivities as edges weights
X_copy = X.copy()
X_copy[np.isnan(X_copy)] = 0

knn_indices_with_zeros, distances, connectivities = kNN_graph(X_copy, k=15, n_components=X_copy.shape[1])





def compute_leiden(A, res=0.5):
    """
    Compute leiden clustering, at some resolution.
    """
    g = sc._utils.get_igraph_from_adjacency(A, directed=True)
    part = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=res,
        seed=1234
    )
    labels = np.array(part.membership)
    return labels


for r in np.linspace(0.1, 1.0, 10):
    np.unique(compute_leiden(connectivities, res=r))

for r in np.linspace(0.1, 1.0, 10):
    np.unique(compute_leiden(D, res=r))

labels_11_solution = compute_leiden(connectivities, res=0.5)
labels_7_solution = compute_leiden(connectivities, res=0.4)
labels_GBCs = afm.obs['GBC'].cat.codes.values



for r in np.linspace(0.1, 10, 10):
    sol = compute_leiden(connectivities, res=r)
    custom_ARI(afm.obs['GBC'], sol)






































