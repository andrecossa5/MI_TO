"""
Scratch for AFM construction, filtering, distance matrices computations, and spectral clustering.
"""

from Cellula._utils import Timer
from Cellula.preprocessing._metrics import custom_ARI
from Cellula.plotting._plotting import *
from MI_TO.preprocessing import *
from MI_TO.kNN import *
from MI_TO.spectral_clustering import *
matplotlib.use('MacOSX')

# Read data
path_data = '/Users/IEO5505/Desktop/MI_TO/data/'
orig = sc.read(path_data + '/AFMs/MDA_clones_afm.h5ad')
CBC_GBC = pd.read_csv(path_data + 'CBC_GBC_cells/CBC_GBC_MDA.csv', index_col=0)

# Create variants AFM
afm, variants = format_matrix(orig, CBC_GBC)

# Filtering: CV, tresholds
afm = nans_as_zeros(afm)
# afm = filter_CV(afm, n=50)
afm = filter_tresholds(afm, median_coverage=100, median_af=0.005, min_perc_cells=0.001)
afm

# kNN computation
X = afm.X
knn, dist, conn = kNN_graph(X, k=50, n_components=X.shape[1])

# Clustering and testing
solutions = []
for r in np.linspace(0.1, 1.5, 10):
    sol = leiden_clustering(conn, res=r)
    solutions.append(list(sol))
    print(len(np.unique(sol)))
    print(custom_ARI(afm.obs['GBC'].cat.codes.values, sol))

# Create df solutions
solutions = pd.DataFrame(
    data=np.array(solutions).T, 
    columns=[ f'res_{round(x, 2)}' for x in np.linspace(0.1, 1.5, 10) ],
    index = afm.obs_names,
    dtype='category'
)
solutions['GBC'] = afm.obs['GBC']

# Viz 
fig, ax = plt.subplots()
ji_cells_one_couple(solutions, 'res_0.88', 'GBC', ax=ax)
plt.show()










































