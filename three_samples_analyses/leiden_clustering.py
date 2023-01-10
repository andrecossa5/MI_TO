#!/usr/bin/python

# Leiden clustering script

########################################################################

# Parsing CLI args 

# Libraries
import sys 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='leiden_clustering',
    description=
    '''
    Takes each sample AFM and its top analysis distance matrices, subset for each cell (row)
    k neighbors, apply the UMAP kernel and partition each resulting kNN graph across some resolutions. 
    Save the output.
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
    '--analysis', 
    type=str,
    default=None,
    help='Analysis from which distance matrices have to be taken. Default: None.'
)

# Resolution range
my_parser.add_argument(
    '--range', 
    type=str,
    default=0.2:1,
    help='Resolution range for leiden clustering. Default: 0.2:1.'
)

# ncores
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=8,
    help='ncores to use for model training. Default: 8.'
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
analysis = args.analysis
res_range = args.range
ncores = args.ncores 

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import pickle
    from Cellula._utils import Timer, set_logger
    from MI_TO.preprocessing import *
    from MI_TO.kNN import *
    from MI_TO.spectral_clustering import *

    #-----------------------------------------------------------------#

    # Set other paths
    path_data = path_main + '/data/'
    path_distances = path_main + '/results_and_plots/distances/'
    path_clones = path_main + '/results_and_plots/clones_classification/'
    path_results = path_main + '/results_and_plots/leiden_clustering/'
    path_runs = path_main + '/runs/'

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, f'logs_{sample}_{analysis}.txt')

########################################################################

# Main
def main():

    T = Timer()
    T.start()

    # Load data
    t = Timer()

    logger.info(f'Execute leiden clustering: --sample {sample} --analysis {analysis}')

    # Read top3 dictionary
    with open(path_clones + 'top3.pkl', 'rb') as f:
        top3 = pickle.load(f)

    # Separate clones and analysis 
    samples = list(top3.keys())
    top_analysis = [ top3[k][0] for k in top3 ]
    
    # For each sample and its top analysis...
    samples_d = {}
    for i in range(len(samples)): 

        top = top_analysis[i]
        sample = samples[i]
        
        # Read in a dictionary the distance matrices of its top analysis
        d = {}
        DISTANCES = {
            x.split('.')[0] : \
            sc.read(path_distances + x) for x in os.listdir(path_distances) \
            if x.startswith(sample) and bool(re.search(top, x))
        }
    
        # For each distance...
        for key in DISTANCES:
            D = DISTANCES[key]
            logger.info(f'Go with distance {key}...')
            t.start()

            # Compute kNN graphs... 
            for k in [5, 15, 30, 50, 100]:
                kNN = kNN_graph(D.X, k=k, n_components=None)

                # And partition them with different settings...
                for res in np.linspace(res_range[0], res_range[1], 5):
                    labels = leiden_clustering(kNN['connectivities'], res=res)
                    d[f'{key}|{k}|{res}'] = labels
                    logger.info(f'Finished with distance matrix {key}: {t.stop()}')

        samples_d[sample] = d

    # Save dictionary
    with open(path_results + 'leiden_clustering.pkl', 'wb') as f:
        pickle.dump(samples_d, f)
                    
    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        main()

#######################################################################