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
a_cells, a = filter_cells_and_vars(
    afm, filtering='pegasus', min_cell_number=50, min_cov_treshold=50, n=5000
)

from dadapy.plot import plot_inf_imb_plane
from dadapy.metric_comparisons import MetricComparisons

d = MetricComparisons(a.X)
best_sets, best_imbs, all_imbs = d.greedy_feature_selection_full(n_coords=10, n_best=1, k=1)

len(best_sets)

best_sets[0]



a.X.shape