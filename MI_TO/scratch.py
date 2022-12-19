"""
Scratch for AFM construction, filtering, distance matrices computations, and vizualization
"""

from Cellula._utils import Timer
from Cellula.plotting._plotting import *
from Cellula.dist_features._dist_features import *
from MI_TO.preprocessing import *
from MI_TO.kNN import *
from MI_TO.distances import *
from MI_TO.heatmaps_plots import *
from MI_TO.spectral_clustering import *
matplotlib.use('MacOSX')

# Read data
path_main = '/Users/IEO5505/Desktop/MI_TO/'
path_data = path_main + 'data/'
path_results = path_main + 'results_and_plots/viz_clones/'
sample = 'MDA'
filtering = 'miller2022'
min_cov_treshold = 30
min_cell_number = 10

# Read data
orig = sc.read(path_data + f'/AFMs/{sample}_afm.h5ad')
CBC_GBC = pd.read_csv(path_data + f'CBC_GBC_cells/CBC_GBC_{sample}.csv', index_col=0)

# Format variants AFM
afm, variants = format_matrix(orig, CBC_GBC)
ncells0 = afm.shape[0]
n_all_clones = len(afm.obs['GBC'].unique())

# Filter cells and vars
if filtering in ['CV', 'ludwig2019', 'velten2021', 'miller2022']:
    # Cells
    a_cells = filter_cells_coverage(afm, mean_coverage=min_cov_treshold) 
    if min_cell_number > 0:
        cell_counts = a_cells.obs.groupby('GBC').size()
        clones_to_retain = cell_counts[cell_counts>min_cell_number].index 
        cells_to_retain = a_cells.obs.query('GBC in @clones_to_retain').index
        a_cells = a_cells[cells_to_retain, :].copy()
    # Variants
    if filtering == 'CV':
        a = filter_CV(a_cells, n=50)
    elif filtering == 'ludwig2019':
        a = filter_ludwig2019(a_cells, mean_AF=0.5, mean_qual=0.2)
    elif filtering == 'velten2021':
        a = filter_velten2021(a_cells, mean_AF=0.1, min_cell_perc=0.2)
    elif filtering == 'miller2022':
        a = filter_miller2022(a_cells, mean_coverage=100, mean_qual=0.3, perc_1=0.001, perc_99=0.15)
        
elif filtering == 'density':
    afm = filter_density(afm, density=0.5, steps=np.Inf)
    if min_cell_number > 0:
        cell_counts = afm.obs.groupby('GBC').size()
        clones_to_retain = cell_counts[cell_counts>min_cell_number].index 
        cells_to_retain = afm.obs.query('GBC in @clones_to_retain').index
        a = afm[cells_to_retain, :].copy()




# Nans robust primitives


# def euclidean_python_nans(x, y):
#     result = 0.0
#     for i in range(x.shape[0]):
#         if (~np.isnan(x[i])) and (~np.isnan(y[i])):
#             result += (x[i] - y[i]) ** 2
#     return np.sqrt(result)

# @njit(fastmath=True, parallel=True)
# def euclidean_pynn_nans(x, y):
#     result = 0.0
#     for i in prange(x.shape[0]):
#         if (~np.isnan(x[i])) and (~np.isnan(y[i])):
#             result += (x[i] - y[i]) ** 2
#     return np.sqrt(result)

# def test_primitive(f, size=10000, n_times=1000, nans=False):
#     """
#     Test a primitive function making a d(x,y) calculations n_times among two size sized 
#     vectors, with or without a 0.5 ratio on nans.
#     """
#     np.random.seed(1234)
#     x = np.random.random(size)
#     np.random.seed(1453)
#     y = np.random.random(size)
# 
#     t = Timer()
#     t.start()
# 
#     for _ in range(n_times):
#         d = f(x, y)
# 
#     print(f'Size: {size}; n_times: {n_times}; Execution: {t.stop():.3f} s')



def euclidean_nans(x, y):
    ix = np.where(~np.isnan(x))[0]
    iy = np.where(~np.isnan(y))[0]
    idx = list(set(ix) & set(iy))
    x_ = x[idx]
    y_ = y[idx]
    return euclidean_std(x_, y_)

def sqeuclidean_nans(x, y):
    ix = np.where(~np.isnan(x))[0]
    iy = np.where(~np.isnan(y))[0]
    idx = list(set(ix) & set(iy))
    x_ = x[idx]
    y_ = y[idx]
    return sqeuclidean(x_, y_)

def correlation_nans(x, y):
    ix = np.where(~np.isnan(x))[0]
    iy = np.where(~np.isnan(y))[0]
    idx = list(set(ix) & set(iy))
    x_ = x[idx]
    y_ = y[idx]
    return correlation(x_, y_)

def cosine_nans(x, y):
    ix = np.where(~np.isnan(x))[0]
    iy = np.where(~np.isnan(y))[0]
    idx = list(set(ix) & set(iy))
    x_ = x[idx]
    y_ = y[idx]
    return cosine(x_, y_)


##


x = np.linspace(0, 1, 100)
x[np.random.randint(0, 100, size=10)] = np.nan
y = np.linspace(2, 3, 100)
y[np.random.randint(0, 100, size=10)] = np.nan

cosine_nans(x, y)


##


def pair_d(X, **kwargs):
    """
    Function for calculating pairwise distances within the row vectors of a matrix X.
    """
    # Get kwargs
    try:
        metric = kwargs['metric']
    except:
        metric = 'euclidean'
    try:
        ncores = kwargs['ncores']
    except:
        ncores = 8
    try:
        nans = kwargs['nans']
    except:
        nans = False

    print(f'pair_d arguments: metric={metric}, ncores={ncores}, nans={nans}')

    # Compute D
    if not nans:
        D = pairwise_distances(X, metric=metric, n_jobs=ncores)
    else:
        print(f'Custom, nan-robust {metric} metric will be used here...')
        if metric == 'euclidean':
            D = pairwise_distances(X, metric=euclidean_nans, n_jobs=ncores, force_all_finite=False)
        elif metric == 'sqeuclidean':
            D = pairwise_distances(X, metric=sqeuclidean_nans, n_jobs=ncores, force_all_finite=False)
        elif metric == 'correlation':
            D = pairwise_distances(X, metric=correlation_nans, n_jobs=ncores, force_all_finite=False)
        elif metric == 'cosine':
            D = pairwise_distances(X, metric=cosine_nans, n_jobs=ncores, force_all_finite=False)

    return D

##

def test_matrix(f, rows=1000, cols=100, **kwargs):
    """
    Test a primitive function making a d(x,y) calculations n_times among two size sized 
    vectors, with or without a 0.5 ratio on nans.
    """
    np.random.seed(1234)
    X = np.random.rand(rows, cols)

    if kwargs['nans']:
        np.random.seed(1234)
        rows_idx = np.random.randint(0, rows, size=rows//2)
        np.random.seed(1234)
        cols_idx = np.random.randint(0, cols, size=cols//2)
        X[np.ix_(rows_idx, cols_idx)] = np.nan

    t = Timer()
    t.start()

    D = f(X, **kwargs)

    print(f'Size: {rows * cols}; nrows {rows}; ncols {cols}; Execution: {t.stop():.3f} s')


##


test_matrix(pair_d, rows=10000, cols=1000, metric='euclidean', nans=False)
test_matrix(pair_d, rows=10000, cols=1000, metric='euclidean', nans=True)





























































































# Filtering: CV, tresholds
afm = nans_as_zeros(afm)
afm = filter_CV(afm, n=1000)
afm = filter_tresholds(afm, median_coverage=100, median_af=0.01, min_perc_cells=0.01)
afm

# kNN computation
X = afm.X
g = kNN_graph(X, k=15, n_components=X.shape[1])
conn = g['connectivities']

# Clustering and testing
solutions = []
for r in np.linspace(0.2, 0.8, 20):
    sol = leiden_clustering(conn, res=r)
    solutions.append(list(sol))
    print(len(np.unique(sol)))
    print(custom_ARI(afm.obs['GBC'].cat.codes.values, sol))

# Create df solutions
solutions = pd.DataFrame(
    data=np.array(solutions).T, 
    columns=[ f'res_{round(x, 2)}' for x in np.linspace(0.2, 0.8, 20) ],
    index = afm.obs_names,
    dtype='category'
)
solutions['GBC'] = afm.obs['GBC']

# Viz 
fig, ax = plt.subplots()
ji_cells_one_couple(solutions, 'res_0.res_0.71', 'GBC', ax=ax)
plt.show()










































