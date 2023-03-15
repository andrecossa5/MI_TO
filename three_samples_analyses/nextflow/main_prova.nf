#!/usr/bin/env nextflow

//

// Paths
params.p = '/Users/IEO5505/Desktop/MI_TO/'
params.path_code = '/Users/IEO5505/Desktop/MI_TO/three_samples_analyses/scripts/'

//

// Create jobs options
process createOptions {

    output:
    path 'jobs.csv'

    script:
    """
    python $params.path_code/create_options.py > 'jobs.csv'
    """

}

// Run!
process runJobs {
  
    // Directives
    cpus 4
    memory '4 G' 
    // errorStrategy 'retry'
    // maxRetries 1
    maxForks 8

    input:
    val x

    output:
    path 'out.csv'
    
    script:
    """
    python $params.path_code/classification_clones.py \
        --p $params.p \
        --min_cov_treshold 50 \
        --ncombos 5 \
        --sample ${x[1]} \
        --model ${x[2]} \
        --filtering ${x[3]} \
        --min_cell_number ${x[4]} > out.csv
    """

}

// Workflow
workflow { 

    jobs = createOptions().splitCsv(skip:1)
    runJobs(jobs)
    view

}


