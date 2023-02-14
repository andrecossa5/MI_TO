#!/usr/bin/python

# Dimred computation script

########################################################################

# Parsing CLI args 

# Libraries
import sys 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='dimred',
    description=
        '''
        Calculate reduced-dimension representation(s) of the original data (potentially reduced by feature selection only).
        '''
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_main', 
    type=str,
    default='..',
    help='The path to the main project directory. Default: .. .'
)

# Filter
my_parser.add_argument(
    '--sample', 
    type=str,
    default='MDA',
    help='Sample to use. Default: MDA.'
)

# Model
my_parser.add_argument(
    '--methods', 
    type=str,
    default='all',
    help='Which dimensionality reduction method to use. Default: all.'
)

# n components
my_parser.add_argument(
    '--n_comps', 
    type=int,
    default=30,
    help='Number of reduced dimension to retain. Default: 30.'
)

# Filter
my_parser.add_argument(
    '--filtering', 
    action='store_true',
    help='Wheter to pre-filter MT-SNVs. Default: False.'
)

# skip
my_parser.add_argument(
    '--skip', 
    action='store_true',
    help='Skip analysis. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
sample = args.sample
if args.methods == 'all':
    methods = ['umap', 'pca', 'vae']
else:
    methods = args.methods.split(':')
n_comps = args.n_comps

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    from joblib import cpu_count
    from Cellula._utils import Timer, set_logger
    from MI_TO.preprocessing import *
    from MI_TO.dimred import *

    #-----------------------------------------------------------------#

    # Set other paths
    path_main = '/Users/IEO5505/Desktop/MI_TO/'
    sample = 'MDA'
    n_comps = 30

    path_data = path_main + '/data/'
    path_results = path_main + '/data/dimred/'
    path_runs = path_main + '/runs/'

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, f'logs_{sample}_dimred.txt')

########################################################################

# Main
def main():

    T = Timer()
    T.start()

    # Load data
    t = Timer()
    t.start()

    logger.info(f'Execute dimensionality reduction: --sample {sample} --filtering {args.filtering}')

    # Read data
    afm = read_one_sample(path_main, sample=sample)

    # Filter 'good quality' cells and variants
    if args.filtering:
        _, a = filter_cells_and_vars(
            afm, 
            filtering='pegasus',
            n=5000, 
            min_cell_number=0, 
            min_cov_treshold=50
        )
    else:
        a = filter_cells_coverage(afm, mean_coverage=50) 

    # Perform dimred
    a = nans_as_zeros(a) 
    X = a.X


    # Dimred
    #...


    from Cellula.preprocessing._pp import my_PCA


    # 
    PCA = my_PCA()
    PCA.calculate_PCA(X, n_components=n_comps)
    X_red = PCA.embs


    from umap.umap_ import UMAP

    umap = UMAP(n_components=n_comps, metric='euclidean', random_state=1234)
    umap.fit_transform(X)


    import scanpy as sc


    sc.tl.diffmap()




    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        main()

#######################################################################

