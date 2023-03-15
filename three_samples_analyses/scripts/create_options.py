#!/usr/bin/python

import pandas as pd
from itertools import product


samples = ['MDA']# , 'AML']#, 'PDX']
models = ['kNN']#, 'logit'] # 'SVM, 'xgboost']
filtering = ['miller2022']#, 'MQuad', 'velten2021', 'CV', 'pegasus', 'seurat', 'density']
min_cell_number = [0]#, 10, 50]

jobs = list(product(samples, models, filtering, min_cell_number)) 
pd.DataFrame(jobs, columns=['sample', 'model', 'filtering', 'min_cell_number']).to_csv('jobs.csv')




