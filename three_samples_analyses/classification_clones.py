#!/usr/bin/python

# Clones classification script

########################################################################

# Parsing CLI args 

# Libraries
import sys 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='clones_classification',
    description=
        '''
        Systematically testing the ability of (filtered) MT-SNVs to distinguish ground truth clonal labels
        from lentiviral barcoding.
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

# Filter
my_parser.add_argument(
    '--filtering', 
    type=str,
    default='ludwig2019',
    help='Method to filter MT-SNVs. Default: ludwig2019.'
)

# Model
my_parser.add_argument(
    '--model', 
    type=str,
    default='xgboost',
    help='Classifier chosen. Default: xgboost.'
)

# ncombos
my_parser.add_argument(
    '--ncombos', 
    type=int,
    default=50,
    help='n combinations of hyperparameters to test. Default: 50.'
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

# Score
my_parser.add_argument(
    '--score', 
    type=str,
    default='f1',
    help='Classification performance scoring type. Default: f1-score.'
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
filtering = args.filtering
model = args.model
ncombos = args.ncombos
ncores = args.ncores 
score = args.score
min_cell_number = args.min_cell_number
min_cov_treshold = args.min_cov_treshold

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    from Cellula._utils import Timer, set_logger
    from Cellula.ML._ML import *
    from Cellula.dist_features._dist_features import *
    from MI_TO.preprocessing import *

    #-----------------------------------------------------------------#

    # Set other paths
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/clones_classification/'
    path_runs = path_main + '/runs/'

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, f'logs_{sample}_{filtering}_{min_cell_number}_{min_cov_treshold}_{model}_{score}.txt')

########################################################################

# Main
def main():

    T = Timer()
    T.start()

    # Load data
    t = Timer()
    t.start()

    logger.info(f'Execute classification: --sample {sample} --filtering {filtering} --model {model} --ncombos {ncombos} --score {score} --min_cell_number {min_cell_number} --min_cov_treshold {min_cov_treshold}')

    # Read data
    afm = read_one_sample(path_main, sample=sample)
    ncells0 = afm.shape[0]
    n_all_clones = len(afm.obs['GBC'].unique())

    # Filter 'good quality' cells and variants
    a_cells, a = filter_cells_and_vars(
        afm, 
        filtering=filtering, 
        min_cell_number=min_cell_number, 
        min_cov_treshold=min_cov_treshold, 
        nproc=ncores, 
        path_=path_results
    )
       
    # Format X and Y for classification
    a = nans_as_zeros(a) # For sklearn APIs compatibility
    ncells = a.shape[0]
    n_clones_analyzed = len(a.obs['GBC'].unique())
    X = a.X
    feature_names = a.var_names
    y = pd.Categorical(a.obs['GBC'])
    Y = one_hot_from_labels(y)

    logger.info(f'Reading and formatting AFM, X and y complete, total {t.stop()} s.')
    logger.info(f'Total cells and clones in the original QCed sample (only transcriptional and perturb seq QC metrics): {ncells0}; {n_all_clones}.')
    logger.info(f'Total cells, clones and variants in final filtered sample: {ncells}; {n_clones_analyzed}, {a.shape[1]}.')

    # Here we go
    DF = []
    for i in range(Y.shape[1]):

        t.start()
        comparison = f'{y.categories[i]}_vs_rest' 
        logger.info(f'Starting comparison {comparison}, {i+1}/{Y.shape[1]}...')

        y_ = Y[:, i]

        # Check numbers 
        if np.sum(y_) > min_cell_number:
            df = classification(X, y_, feature_names, key=model, GS=True, 
                score=score, n_combos=ncombos, cores_model=ncores, cores_GS=1)
            df = df.assign(comparison=comparison, feature_type=filtering)          
            df = df.loc[:,
                [
                    'feature_type', 'rank', 'evidence', 'evidence_type', 
                    'effect_size', 'es_rescaled', 'effect_type', 'comparison'
                ]
            ]
            DF.append(df)
            logger.info(f'Comparison {comparison} finished: {t.stop()} s.')
        else:
            logger.info(f'Clone {y.categories[i]} does not reach {min_cell_number} cells. This should not happen here... {t.stop()} s.')

    df = pd.concat(DF, axis=0)
    df['evidence'].describe()

    # Save results
    df.to_excel(path_results + f'clones_{sample}_{filtering}_{min_cell_number}_{min_cov_treshold}_{model}_{score}.xlsx')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        main()

#######################################################################

