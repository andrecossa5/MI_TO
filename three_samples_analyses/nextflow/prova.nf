ch = Channel.fromPath('/Users/IEO5505/Desktop/cellula_example/data/default/full_embs.csv')

process foo {

    conda '/Users/IEO5505/mambaforge/envs/MI_TO'

    input:
    path x
    output:
    stdout

    //"""
    //cat $x | cut -d , -f 1,4-5 | head
    //"""

    """
    #!/usr/bin/python 
    import scanpy as sc
    print('Hello world!')
    """

}

workflow {
    ch | foo | view
}
