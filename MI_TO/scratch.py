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

# Read data
afm = read_one_sample(path_main, 'MDA')
a_cells, a = filter_cells_and_vars(afm, filtering='miller2022', min_cell_number=50, min_cov_treshold=50)




from pegasus.tools.hvf_selection import select_hvf_pegasus


def select_hvf_pegasus(
    data: Union[MultimodalData, UnimodalData], batch: str, n_top: int = 2000, span: float = 0.02
    ) -> None:
    """ 
    Select highly variable features using the pegasus method
    """

    estimate_feature_statistics(data, batch)

    robust_idx = data.var["robust"].values
    hvf_index = np.zeros(afm.shape[1], dtype=bool)

    mean = data.var.loc[robust_idx, "mean"]
    var = data.var.loc[robust_idx, "var"]

    span_value = span
    while True:
        lobj = fit_loess(mean, var, span = span_value, degree = 2)
        if lobj is not None:
            break
        span_value += 0.01
    if span_value > span:
        logger.warning("Leoss span is adjusted from {:.2f} to {:.2f} to avoid fitting errors.".format(span, span_value))

    rank1 = np.zeros(hvf_index.size, dtype=int)
    rank2 = np.zeros(hvf_index.size, dtype=int)

    delta = var - lobj.outputs.fitted_values
    fc = var / lobj.outputs.fitted_values

    rank1[np.argsort(delta)[::-1]] = range(hvf_index.size)
    rank2[np.argsort(fc)[::-1]] = range(hvf_index.size)
    hvf_rank = rank1 + rank2

    hvf_index[np.argsort(hvf_rank)[:n_top]] = True

    data.var["hvf_loess"] = 0.0
    data.var.loc[robust_idx, "hvf_loess"] = lobj.outputs.fitted_values

    data.var["hvf_rank"] = -1
    data.var.loc[robust_idx, "hvf_rank"] = hvf_rank
    data.var["highly_variable_features"] = False
    data.var.loc[robust_idx, "highly_variable_features"] = hvf_index








































































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










































