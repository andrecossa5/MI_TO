#!/usr/bin/python

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


##

# Paths
path_CBCs = sys.argv[1]
path_UMIs = sys.argv[2]
path_GBCs = sys.argv[3]
path_aligned_reads = sys.argv[4]

# Read components and aligned reads
CBCs = pd.read_csv(path_CBCs, index_col=0, sep='\t', header=None)
UMIs = pd.read_csv(path_UMIs, index_col=0, sep='\t', header=None)
GBCs = pd.read_csv(path_GBCs, index_col=0, sep='\t', header=None)
aligned_names = pd.read_csv(path_aligned_reads, index_col=0, sep='\t', header=None)

# Chek read names are identical
if not (CBCs.index == UMIs.index).all() and (UMIs.index == GBCs.index).all():
    raise ValueError('Read names of filtered CBCs, UMIs and GBCs tables are not identical. Something is wrong...')
# Merge in a single table
df = pd.concat([CBCs, UMIs, GBCs], axis=1)
df.columns = ['CBC', 'UMI', 'GBC']

# Filter for only reads aligned to the bulk reference, write all CBC-UMI-GBC table
df = df.loc[aligned_names.index,:]
df.to_csv('CBC_UMI_GBC_by_read.tsv.gz', sep='\t')

# Compute unique CBC-GBC combo table
df_read_counts = df.groupby(['CBC', 'GBC']).size().to_frame('read_counts')
df_umi_counts = df.reset_index().loc[:, ['CBC', 'GBC', 'UMI']].drop_duplicates().groupby(['CBC', 'GBC']).size().to_frame('umi_counts')
df_combos = df_read_counts.join(df_umi_counts).reset_index()
df_combos['coverage'] = df_combos['read_counts'] / df_combos['umi_counts']
df_combos.to_csv('CBC_GBC_combos.tsv.gz', sep='\t')

# Cell classification: from Adamson et al. 2016
test = (df_combos['read_counts'] > 30) & \
        (df_combos['umi_counts'] > 3) & \
        (df_combos['coverage'] > 1.5) # TOTUNE
df_combos['status'] = np.where(test, 1, 0)
df_combos_agg = df_combos.reset_index().loc[:, ['CBC', 'GBC', 'status']].groupby(['CBC', 'GBC']).agg({'status': 'sum'})
unique_combos = df_combos_agg.query('status == 1').index
unique_combos = [ [k, v] for k, v in zip(unique_combos.get_level_values('CBC'), unique_combos.get_level_values('GBC')) ]

# Cell assignment plot
fig, ax = plt.subplots()
x = np.log10(df_combos['read_counts'])
y = np.log10(df_combos['umi_counts'])
ax.plot(x[df_combos['status'] == 1], y[df_combos['status'] == 1], '.', label='assigned', color='blue')
ax.plot(x[df_combos['status'] == 0], y[df_combos['status'] == 0], '.', label='not-assigned', color='grey')
ax.set(title='CBC-GBC combination status', xlabel='log10_read_counts', ylabel='log10_umi_counts')
ax.legend()
fig.savefig('CBC_GBC_combo_status.png')


# Compute summary tables: cells 
df_cells = pd.merge(
    pd.DataFrame(unique_combos, columns=['CBC', 'GBC']), 
    df_combos.drop(columns='status'),
    on=['CBC', 'GBC']
).set_index('CBC')
df_cells.to_csv('cells_summary_table.csv')

# Compute summary tables: clones
df_clones = df_cells.reset_index().groupby('GBC').size().to_frame('n_cells')
df_clones['prevalence'] = df_clones['n_cells'] / df_clones['n_cells'].sum()
df_clones.to_csv('clones_summary_table.csv')

# Plots: correlation with bulk and distribution
fig, ax = plt.subplots()
x = np.log10(df_combos['read_counts'])
y = np.log10(df_combos['umi_counts'])
ax.plot(x[df_combos['status'] == 1], y[df_combos['status'] == 1], '.', label='assigned', color='blue')
ax.plot(x[df_combos['status'] == 0], y[df_combos['status'] == 0], '.', label='not-assigned', color='grey')
ax.set(title='CBC-GBC combination status', xlabel='log10_read_counts', ylabel='log10_umi_counts')
ax.legend()
fig.savefig('CBC_GBC_combo_status.png')