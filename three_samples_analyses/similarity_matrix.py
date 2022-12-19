#!/usr/bin/python

# Calculate cell-cell distance/affinity matrix script

########################################################################

# Parsing CLI args 

# Libraries
import sys 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='similarity_matrix',
    description=
        '''
        Calculates a cell-to-cell pairwise distance or affinity matrix for all (included) cells in a sample(s).
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

# Sample
my_parser.add_argument(
    '--sample', 
    type=str,
    default='MDA',
    help='Sample to use. Default: MDA. If None, this is done for all cells in all samples found in $path_main/data.'
)

# Filtering
my_parser.add_argument(
    '--filtering', 
    type=str,
    default='ludwig2019',
    help='Method to filter MT-SNVs. Default: ludwig2019.'
)

# metric
my_parser.add_argument(
    '--metric', 
    type=str,
    default='euclidean',
    help='Distance metric chosen. Default: euclidean.'
)

# ncores
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=8,
    help='ncores to use for model training. Default: 8.'
)

# min_cell_number
my_parser.add_argument(
    '--min_cell_number', 
    type=int,
    default=0,
    help='Include in the analysis only cells with membership in clones with >= min_cell_number. Default: 0.'
)

# min_cov_treshold
my_parser.add_argument(
    '--min_cov_treshold', 
    type=int,
    default=30,
    help='Include in the analysis only cells MAESTER sites mean coverage > min_cov_treshold. Default: 30.'
)

# Kernel
my_parser.add_argument(
    '--kernel', 
    type=str,
    default=None,
    help='Transform some distance metric with some kernel. Default: None'
)

# Skip
my_parser.add_argument(
    '--retain_nans', 
    action='store_true',
    help='Retain nans as missing values, and vot zeroes. Default" False',
)

# Skip
my_parser.add_argument(
    '--skip', 
    action='store_true',
    help='Skip analysis. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
sample = args.sample
filtering = args.filtering
metric = args.metric
ncores = args.ncores 
kernel = args.kernel if args.kernel is not None else 'no_kernel'
retain_nans = args.retain_nans
min_cell_number = args.min_cell_number
min_cov_treshold = args.min_cov_treshold

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    from Cellula._utils import Timer, set_logger
    from MI_TO.preprocessing import *
    from MI_TO.distances import *

    #-----------------------------------------------------------------#

    # Set other paths
    path_main = '/Users/IEO5505/Desktop/MI_TO/'
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/distances/'
    path_runs = path_main + '/runs/'

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, f'logs_{sample}_{filtering}_{min_cell_number}_{min_cov_treshold}_{metric}_{retain_nans}_{kernel}.txt')

########################################################################

# Main
def main():

    T = Timer()
    T.start()

    # Load data
    t = Timer()
    t.start()

    logger.info(f'Execute classification: --sample {sample} --filtering {filtering} --metric {metric} --kernel {kernel} --retain_nans {retain_nans} --min_cell_number {min_cell_number} --min_cov_treshold {min_cov_treshold}')

    # Read data
    orig = sc.read(path_data + f'/AFMs/{sample}_afm.h5ad')
    CBC_GBC = pd.read_csv(path_data + f'CBC_GBC_cells/CBC_GBC_{sample}.csv', index_col=0)

    # Format variants AFM
    afm, variants = format_matrix(orig, CBC_GBC)
    ncells0 = afm.shape[0]
    n_all_clones = len(afm.obs['GBC'].unique())

    # Filter 'good quality cells': 
    # 1: passing transcriptional QC (already done); 
    # 2: having a mean MEASTER site coverage >= min_cov_treshold;
    # 3: being from clones with >= min_cell_number cells;
    # Only on 'good quality cells', filter variants. If density method is used, cells and genes are filtered jointly, 
    # and no fixed treshold on MAESTER sites mean coverage is applied.

    if filtering in ['CV', 'ludwig2019', 'velten2021', 'miller2022']:

        # Cells
        afm = filter_cells_coverage(afm, mean_coverage=min_cov_treshold) 
        if min_cell_number > 0:
            cell_counts = afm.obs.groupby('GBC').size()
            clones_to_retain = cell_counts[cell_counts>min_cell_number].index 
            cells_to_retain = afm.obs.query('GBC in @clones_to_retain').index
            afm = afm[cells_to_retain, :].copy()

        # Variants
        if filtering == 'CV':
            a = filter_CV(afm, n=50)
        elif filtering == 'ludwig2019':
            a = filter_ludwig2019(afm, mean_AF=0.5, mean_qual=0.2)
        elif filtering == 'velten2021':
            a = filter_velten2021(afm, mean_AF=0.1, min_cell_perc=0.2)
        elif filtering == 'miller2022':
            a = filter_miller2022(afm, mean_coverage=100, mean_qual=0.3, perc_1=0.01, perc_99=0.1)

    elif filtering == 'density':
        a = filter_density(afm, density=0.5, steps=np.Inf)
        if min_cell_number > 0:
            cell_counts = afm.obs.groupby('GBC').size()
            clones_to_retain = cell_counts[cell_counts>min_cell_number].index 
            cells_to_retain = afm.obs.query('GBC in @clones_to_retain').index
            a = afm[cells_to_retain, :].copy()
    
    # Format afm for D computation
    if not retain_nans:
        a = nans_as_zeros(a)
        logger.info('Convert nans into zeros...')
    else:
        logger.info('nans mantained in the filtered feature matrix...')
    ncells = a.shape[0]
    n_clones_analyzed = len(a.obs['GBC'].unique())

    logger.info(f'Reading and formatting AFM: total {t.stop()} s.')
    logger.info(f'Total cells and clones in the original QCed sample (only transcriptional and perturb seq QC metrics): {ncells0}; {n_all_clones}.')
    logger.info(f'Total cells, clones and variants in final filtered sample: {ncells}; {n_clones_analyzed}, {a.shape[1]}.')

    # Calculate distance matrix
    t.start()
    logger.info('Begin pairwise distances calculations...')
    D = pair_d(a.X, metric=metric, ncores=ncores, nans=retain_nans)
    logger.info(f'Finished with distances: {t.stop()} s.')

    # Save as .h5ad adata object
    df_cells = pd.DataFrame({'cell' : a.obs_names}).set_index('cell')
    D = anndata.AnnData(X=D, obs=df_cells, var=df_cells)
    D.uns['distance_matrix'] = {
        'filtering' : filtering, 
        'retain_nans' : retain_nans, 
        'metric' : metric, 
        'kernel' : kernel, 
        'min_cell_number' : min_cell_number, 
        'min_cov_treshold' : min_cov_treshold
    }
    D.write(path_results + f'{sample}_{filtering}_{min_cell_number}_{min_cov_treshold}_{retain_nans}_{metric}_{kernel}.h5ad')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        main()

#######################################################################


