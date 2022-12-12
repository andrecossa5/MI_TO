#!/usr/bin/python

# Samples classification script

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
        Systematically testing the ability of (filtered) MT-SNVs to distinguish sample memberships.
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

# Model
my_parser.add_argument(
    '--ncombos', 
    type=int,
    default=50,
    help='n combinations of hyperparameters to test. Default: 50.'
)

# Model
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=8,
    help='ncores to use for model training. Default: 8.'
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
filtering = args.filtering
model = args.model
ncombos = args.ncombos
ncores = args.ncores 
score = args.score

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
    path_results = path_main + '/results_and_plots/samples_classification/'
    path_runs = path_main + '/runs/'

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, f'logs_{filtering}_{model}_{ncombos}_{score}.txt')

########################################################################

# Main
def main():

    T = Timer()
    T.start()

    # Load data
    t = Timer()
    t.start()

    logger.info(f'Execute classification: --filtering {filtering} --model {model} --ncombos {ncombos} --score {score}')

    # Read and format data
    ORIG = {}
    samples = ['MDA_clones', 'AML_clones', 'PDX']
    for x in samples:
        orig = sc.read(path_data + f'AFMs/{x}_afm.h5ad')
        orig.obs = orig.obs.assign(sample=x)
        ORIG[x] = orig
        meta_vars = orig.var
    orig = anndata.concat(ORIG.values(), axis=0)
    orig.var = meta_vars
    del ORIG

    # Create variants AFM and filter variants
    afm, variants = format_matrix(orig, no_clones=True)

    if filtering == 'CV':
        a = filter_CV(afm, mean_coverage=100, n=50)
    elif filtering == 'ludwig2019':
        a = filter_ludwig2019(afm, mean_coverage=100, mean_AF=0.5, mean_qual=0.2)
    elif filtering == 'velten2021':
        a = filter_velten2021(afm, mean_coverage=100, mean_AF=0.1, min_cell_perc=0.2)
    elif filtering == 'miller2022':
        a = filter_miller2022(afm, mean_coverage=100, mean_qual=0.3, perc_1=0.01, perc_99=0.1)
    elif filtering == 'density':
        a = filter_density(afm, density=0.7, steps=np.Inf)
    
    # Format X and Y for classification
    a = nans_as_zeros(a)
    X = a.X
    feature_names = a.var_names
    y = pd.Categorical(a.obs['sample'])
    Y = one_hot_from_labels(y)

    logger.info(f'Reading and formatting AFM, X and y complete, total {t.stop()} s.')

    # Here we go
    DF = []
    for i in range(Y.shape[1]):

        t.start()
        comparison = f'{y.categories[i]}_vs_rest' 
        logger.info(f'Starting comparison {comparison}, {i+1}/{Y.shape[1]}...')

        y_ = Y[:, i]

        # Check numbers 
        if np.sum(y_) > 50:
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
            logger.info(f'Finished comparison {comparison}, {t.stop()} s.')
        
        else:
            print(f'Sample{y.categories[i]} does not reach 50 cells. Skip this one...')

    df = pd.concat(DF, axis=0)
    df['evidence'].describe()

    # Save results
    df.to_excel(path_results + f'{filtering}_{model}_{ncombos}_{score}.xlsx')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        main()

#######################################################################
