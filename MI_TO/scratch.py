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
filtering = 'ludwig2019'
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

a
a = nans_as_zeros(a)
ncells = a.shape[0]
n_clones_analyzed = len(a.obs['GBC'].unique())

# Viz cells x vars
g = cells_vars_heatmap(a, covariate='GBC', palette_anno='dark', cmap='mako',
    var_names=True, cell_names=False, cluster_cells=True, cluster_vars=True,
    dendrogram_ratio=(.3, .04), colors_ratio=0.05, figsize=(11, 8), var_names_size=5,
    cut_from_right=0.7, halign_title=0.47, position_cbar=(0.82, 0.2, 0.02, 0.25), 
    title_legend='Clones', loc_legend='lower center', bbox_legend=(0.825, 0.5)
)

# Viz cells x cells distances

g = cells_vars_heatmap(a, covariate='GBC', palette_anno='dark', cmap='mako',
    var_names=True, cell_names=False, cluster_cells=True, cluster_vars=True,
    dendrogram_ratio=(.3, .04), colors_ratio=0.05, figsize=(11, 8), var_names_size=5,
    cut_from_right=0.7, halign_title=0.47, position_cbar=(0.82, 0.2, 0.02, 0.25), 
    title_legend='Clones', loc_legend='lower center', bbox_legend=(0.825, 0.5)
)

plt.show()






















































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










































