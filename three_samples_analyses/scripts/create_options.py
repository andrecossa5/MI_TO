#!/usr/bin/python

import pandas as pd
from itertools import product

# Lists
samples = ['MDA']# , 'AML']#, 'PDX']
input_mode = ['less_stringent']#, 'more_stringent']
filtering = ['miller2022']#, 'MQuad', 'velten2021', 'CV', 'pegasus', 'seurat', 'density']
dimred = ['no_dimred']#, 'more_stringent']
models = [ 'logit'] #, 'kNN', 'SVM, 'xgboost']
min_cell_number = [0]#, 10, 50]

# Product, and write
jobs = list(product(samples, input_mode, filtering, dimred, models, min_cell_number)) 
pd.DataFrame(
    jobs, 
    columns=['sample', 'input_mode', 'filtering', 'dimred', 'model', 'min_cell_number']
).to_csv('jobs.csv')



