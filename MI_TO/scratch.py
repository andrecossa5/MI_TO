"""
Scratch for AFM construction, filtering, distance matrices computations, and spectral clustering.
"""

from Cellula._utils import Timer
from Cellula.preprocessing._metrics import custom_ARI
from Cellula.plotting._plotting import *
from Cellula.ML._ML import *
from Cellula.dist_features._dist_features import *
from MI_TO.preprocessing import *
from MI_TO.kNN import *
from MI_TO.spectral_clustering import *
matplotlib.use('MacOSX')

# Read data
path_data = '/Users/IEO5505/Desktop/MI_TO/data/'
orig = sc.read(path_data + '/AFMs/MDA_clones_afm.h5ad')
CBC_GBC = pd.read_csv(path_data + 'CBC_GBC_cells/CBC_GBC_MDA_clones.csv', index_col=0)

# Create variants AFM
afm, variants = format_matrix(orig, CBC_GBC)

filter_CV(afm, mean_coverage=100, n=50)
filter_ludwig2019(afm, mean_coverage=100, mean_AF=0.5, mean_qual=0.2)
filter_velten2021(afm, mean_coverage=100, mean_AF=0.1, min_cell_perc=0.2)
filter_miller2022(afm, mean_coverage=100, mean_qual=0.3, perc_1=0.01, perc_99=0.1)
filter_density(afm, density=0.7, steps=np.Inf)

# Class
a = filter_miller2022(afm, mean_coverage=100, mean_qual=0.3, perc_1=0.01, perc_99=0.1)
a = nans_as_zeros(a)
X = a.X

feature_names = a.var_names
y = pd.Categorical(a.obs['GBC'])

if len(y.categories) > 2:
    Y = one_hot_from_labels(y)
    # Here we go
    DF = []
    for i in range(Y.shape[1]):
        comparison = f'{y.categories[i]}_vs_rest' 
        y_ = Y[:, i]
        df = classification(X, y_, feature_names, key='xgboost', GS=True, 
            score='f1', n_combos=10, cores_model=8, cores_GS=1)
        df = df.assign(comparison=comparison, feature_type='CV')          
        df = df.loc[:,
            ['feature_type', 'rank', 'evidence', 'evidence_type', 'effect_size', 'es_rescaled',
            'effect_type', 'comparison']
        ]
        DF.append(df)

df = pd.concat(DF, axis=0)
df['evidence'].describe()
























































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










































