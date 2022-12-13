"""
Summary of classification performances.
"""

# Code
import os
import re
from Cellula._utils import Timer, set_logger
from Cellula.plotting._plotting import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *
from Cellula.dist_features._dist_features import *
from MI_TO.preprocessing import *

# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/'
path_data = path_main + '/data/'
path_clones = path_main + '/results_and_plots/clones_classification/'
path_samples = path_main + '/results_and_plots/samples_classification/'

res = []
for x in os.listdir(path_clones):
    if bool(re.search('PDX', x)):
        df = pd.read_excel(path_clones + x, index_col=0)
        d = { k : v for k, v in zip(['sample', 'filter', 'classifier'], x.split('_')[1:-2]) }
        d['median_f1'] = np.median(df['evidence'])
        d['n_clones'] = df['comparison'].unique().size
        res.append(d)
    elif bool(re.search('MDA', x)) or bool(re.search('AML', x)):
        df = pd.read_excel(path_clones + x, index_col=0)
        l = x.split('_')[1:-2]
        d = {}
        d['sample'] = '_'.join(l[:2])
        d['filter'] = l[2]
        d['classifier'] = l[3]
        d['median_f1'] = np.median(df['evidence'])
        d['n_clones'] = df['comparison'].unique().size
        res.append(d)

res = pd.DataFrame(res).sort_values('median_f1', ascending=False)
res.to_excel(path_clones + 'report_classification_clones.xlsx')

res = []
for x in os.listdir(path_samples):
    df = pd.read_excel(path_samples + x, index_col=0)
    d = { k : v for k, v in zip(['filter', 'classifier'], x.split('_')[:2]) }
    d['n_samples'] = 3
    d['median_f1'] = np.median(df['evidence'])
    res.append(d)

res = pd.DataFrame(res).sort_values('median_f1', ascending=False)
res.to_excel(path_clones + 'report_classification_samples.xlsx')

    


