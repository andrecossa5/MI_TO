"""
Visualization of clones and samples classification performances.
"""

# Code
import pickle
import re
import os
import sys
from Cellula.plotting._plotting import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *
from MI_TO.preprocessing import *
from MI_TO.diagnostic_plots import sturges
from MI_TO.heatmaps_plots import *
from MI_TO.utils import *


##


# Set paths
path_main = sys.argv[1]
sample_names = ['MDA', 'PDX', 'AML']

path_clones = path_main + '/results_and_plots/clones_classification/'
path_samples = path_main + '/results_and_plots/samples_classification/'
path_results = path_main + '/results_and_plots/classification_performance/'

# Read reports
clones = pd.read_excel(path_clones + 'report_classification_clones.xlsx', index_col=0)
samples = pd.read_excel(path_samples + 'report_classification_samples.xlsx', index_col=0)

# Re-format analysis
clones['analysis'] += '_' + clones['model']
samples['analysis'] += '_' + samples['model']


##


############## Median f1 by comparision and task
fig, ax = plt.subplots(figsize=(8,8))

np.random.seed(123)

clones_agg = clones.groupby('comparison').agg({'f1' : np.median}).assign(task='Clones')
samples_agg = samples.groupby('comparison').agg({'f1' : np.median}).assign(task='Samples')
df_ = pd.concat([clones_agg, samples_agg], axis=0)
df_['Median f1'] = df_['f1']

strip(df_, 'task', 'Median f1', by='task', c={'Clones':'b', 'Samples':'r'}, a=1, l=None, s=7, ax=ax, with_stats=False, pairs=None)
ax.set(title='Median f1-score by comparison and task')

ax.text(0.42, 0.37, f'-Total analyses for clones (samples): {clones["analysis"].unique().size} ({samples["analysis"].unique().size})', transform=ax.transAxes)
ax.text(0.42, 0.34, f'-Total comparisons for clones (samples): {clones["comparison"].unique().size} ({samples["comparison"].unique().size})', transform=ax.transAxes)
n_clones = clones.loc[:, ['sample', 'comparison']].drop_duplicates().groupby('sample').count()['comparison']
ax.text(0.42, 0.31, f'-n clones by sample: MDA {n_clones["MDA"]}; AML {n_clones["AML"]}; PDX {n_clones["PDX"]};', transform=ax.transAxes)
ax.text(0.42, 0.22, f'-Median n of analyses a comparison occurred in,', transform=ax.transAxes)
ax.text(0.42, 0.19, f' for clones (samples): {np.median(clones.groupby("comparison").size())} ({np.median(samples.groupby("comparison").size())})', transform=ax.transAxes)
ax.text(0.42, 0.16, f'-Median f1-scores for clones (samples):', transform=ax.transAxes)
ax.text(0.42, 0.13, f' {np.median(clones["f1"]):.3f} ({clones["f1"].size} values across all analyses)', transform=ax.transAxes)
ax.text(0.42, 0.10, f' {np.median(samples["f1"]):.3f} ({samples["f1"].size} values across all analyses)', transform=ax.transAxes)

ax.text(1, np.median(samples.query('comparison == "PDX_vs_rest"')['f1']), 'PDX') 
ax.text(0.85, np.median(samples.query('comparison == "MDA_vs_rest"')['f1']), 'MDA')
ax.text(1.08, np.median(samples.query('comparison == "AML_vs_rest"')['f1']), 'AML')

top_3_clones = clones_agg.sort_values('f1', ascending=False).head(3)
ax.text(-0.4, np.median(clones.query('comparison == "CGCCGAACAGCTTCAGTG_vs_rest"')['f1']) + 0.02, 'CGCCGAACAGCTTCAGTG')
ax.text(-0.4, np.median(clones.query('comparison == "TCCCTGGAGTCTTCGAAC_vs_rest"')['f1']) + 0.04, 'TCCCTGGAGTCTTCGAAC')
ax.text(-0.4, np.median(clones.query('comparison == "CTCCTCCGCGGCGAAACG_vs_rest"')['f1']) - 0.06, 'CTCCTCCGCGGCGAAACG')

# Save
fig.savefig(path_results + 'median_f1_by_task.pdf')
############## 


##


############## f1 by clone and feat_type
feat_type_colors = create_palette(clones, 'feature_type', 'dark')

fig, ax = plt.subplots(figsize=(12, 5))

strip(clones, 'comparison', 'f1', by='feature_type', c=feat_type_colors, s=3, ax=ax)
box(clones, 'comparison', 'f1', c='grey', ax=ax)
format_ax(clones, ax, title='f1-scores by clone and variant selection method', rotx=90, xsize=5)
create_handles(feat_type_colors.keys(), marker='o', colors=None, size=10, width=0.5)
handles = create_handles(feat_type_colors.keys(), colors=feat_type_colors.values())
fig.legend(handles, feat_type_colors.keys(), loc='upper right', 
    bbox_to_anchor=(0.9, 0.9), ncol=2, frameon=False, title='Feature selection'
)

v = 0.8
np.sum([ np.sum(clones.query('comparison == @x')['f1'] > v) > 0 for x in clones['comparison'].unique() ])

v = 0.5
np.sum([ np.median(clones.query('comparison == @x')['f1']) > v for x in clones['comparison'].unique() ])

ax.text(0.2, 0.8, f'-n clones with more than one classification above 0.5 f1: 6', transform=ax.transAxes)
ax.text(0.2, 0.75, f'-n clones with more than one classification above 0.8 f1: 6', transform=ax.transAxes)
ax.text(0.2, 0.7, f'-n clones with median f1 > 0.5: 1', transform=ax.transAxes)

# Save
fig.tight_layout()
fig.savefig(path_results + 'clones_f1.pdf')
##############


##


############## f1 by sample and feat_type
fig, ax = plt.subplots(figsize=(8, 6.8))

strip(clones, 'sample', 'f1', by='feature_type', c=feat_type_colors, s=5, ax=ax)
box(clones, 'sample', 'f1', c='grey', ax=ax, s=0.3, a=0.001)
format_ax(clones, ax, title='Clones f1-scores by sample and variant selection method')
create_handles(feat_type_colors.keys(), marker='o', colors=None, size=10, width=0.5)
handles = create_handles(feat_type_colors.keys(), colors=feat_type_colors.values())
fig.legend(handles, feat_type_colors.keys(), loc='center', 
    bbox_to_anchor=(0.8, 0.65), ncol=1, frameon=False, title='Feature selection'
)
ax.text(0.7, 0.5, f'Mean f1 clones MDA: {clones.loc[clones["sample"]=="AML"]["f1"].mean():.3f}', transform=ax.transAxes)
ax.text(0.7, 0.47, f'Mean f1 clones AML: {clones.loc[clones["sample"]=="MDA"]["f1"].mean():.3f}', transform=ax.transAxes)
ax.text(0.7, 0.44, f'Mean f1 clones PDX: {clones.loc[clones["sample"]=="PDX"]["f1"].mean():.3f}', transform=ax.transAxes)

# Save
fig.tight_layout()
fig.savefig(path_results + 'clones_f1_by_sample.pdf')
##############


##



############## f1 by model and feat_type
fig, ax = plt.subplots(figsize=(8, 6.5))

strip(clones, 'model', 'f1', by='feature_type', c=feat_type_colors, s=5, ax=ax)
box(clones, 'model', 'f1', c='grey', ax=ax, s=0.3, a=0.001)
format_ax(clones, ax, title='Clones f1-scores by model and variant selection method')
create_handles(feat_type_colors.keys(), marker='o', colors=None, size=10, width=0.5)
handles = create_handles(feat_type_colors.keys(), colors=feat_type_colors.values())
fig.legend(handles, feat_type_colors.keys(), loc='center', 
    bbox_to_anchor=(0.5, 0.6), ncol=1, frameon=False, title='Feature selection'
)
ax.text(0.35, 0.45, f'Mean f1 clones xgboost: {clones.loc[clones["model"]=="xgboost"]["f1"].mean():.3f}', transform=ax.transAxes)
ax.text(0.35, 0.42, f'Mean f1 clones logit: {clones.loc[clones["model"]=="logit"]["f1"].mean():.3f}', transform=ax.transAxes)

# Save
fig.tight_layout()
fig.savefig(path_results + 'clones_f1_by_model.pdf')
##############


##


############## f1 by sample and feat_type
# Sizes
res = []
for x in os.listdir(path_main + 'data/CBC_GBC_cells'):
    if x.endswith('csv'):
        d = pd.read_csv(path_main + f'data/CBC_GBC_cells/{x}', index_col=0)
        res.append(d.assign(sample=x.split('_')[-1].split('.')[0]))
CBC_GBC = pd.concat(res, axis=0)
clones_sizes = CBC_GBC.groupby('GBC').size()

clones['GBC'] = clones['comparison'].map(lambda x: x.split('_')[0])

csizes = []
for x in clones['GBC']:
    csizes.append(clones_sizes[x])
clones['size'] = csizes

# Viz
fig, ax = plt.subplots(figsize=(6, 5))
scatter(clones, 'size', 'f1', c='black', s=10, ax=ax)
x = clones['size']
y = clones['f1']
fitted_coefs = np.polyfit(x, y, 1)
y_hat = np.poly1d(fitted_coefs)(x)
ax.plot(x, y_hat, linestyle='--', color='red')
corr = np.corrcoef(x, y)[0,1]
ax.text(0.6, 0.9, f"Pearson's r: {corr:.2f}", transform=ax.transAxes)
format_ax(clones, ax, title='f1-clone size correlation', xlabel='Clone size', ylabel='f1')

# Save
fig.tight_layout()
fig.savefig(path_results + 'clones_size_f1_corr.pdf')
##############



##


############## 
# For each sample (3x) clones, what are the top 3 analyses (median f1 score across clones)? 
# Intersection among selected SNVs??

# Save top3 for easy quering
top_3 = {}
for sample in clones['sample'].unique():
    top_3[sample] = clones.query('sample == @sample').groupby(['analysis']).agg(
        {'f1':np.median}).sort_values(
        'f1', ascending=False).index[:3].to_list()

# Load top3 variants for each sample clones, and visualize their intersection (i.e., J.I.), by sample
top3_sample_variants = {}
for sample in sample_names:
        var_dict = {}
        for x in os.listdir(path_results + f'top_3/{sample}/'):
            if x.endswith('.xlsx'):
                n = '_'.join(x.split('.')[0].split('_')[2:-1])
                df_ = pd.read_excel(path_results + f'top_3/{sample}/{x}', index_col=0)
                var_dict[n] = df_.index.unique().to_list()
        top3_sample_variants[sample] = var_dict

# Sample a
fig, axs = plt.subplots(1, 3, figsize=(10,5))

for k, sample in enumerate(top3_sample_variants):
    n_analysis = len(top3_sample_variants[sample].keys())
    JI = np.zeros((n_analysis, n_analysis))
    for i, l1 in enumerate(top3_sample_variants[sample]):
        for j, l2 in enumerate(top3_sample_variants[sample]):
            x = top3_sample_variants[sample][l1]
            y = top3_sample_variants[sample][l2]
            JI[i, j] = ji(x, y)
    JI = pd.DataFrame(data=JI, index=None, columns=top3_sample_variants[sample].keys())

    plot_heatmap(JI, palette='mako', ax=axs[k], title=sample, y_names=False,
        x_names_size=10, y_names_size=0, annot=True, annot_size=10, cb=True, label='JI variants'
    )

fig.tight_layout()
fig.savefig(path_results + 'overlap_selected_vars.pdf')
##############


##


############## 
# For each sample top3 analysis on the clone task, what are the AF profiles of the variants selected?
# Which relatinship can we visualize among clone cells, using:
# 1) hclustering of cell x var AFM 
# 2) hclustering of a cell x cell similarity matrix?

path_data = path_main + 'data/'
path_distances = path_main + 'results_and_plots/distances/'

# Here we go
for sample in sample_names:

    if not os.path.exists(path_results + 'top_3'):
        os.mkdir(path_results + 'top_3')
    os.chdir(path_results + 'top_3')
    
    if not os.path.exists(path_results + f'top_3/{sample}'):
        os.mkdir(sample)
    os.chdir(sample)

    # Read data
    orig = sc.read(path_data + f'/AFMs/{sample}_afm.h5ad')
    CBC_GBC = pd.read_csv(path_data + f'CBC_GBC_cells/CBC_GBC_{sample}.csv', index_col=0)

    # Format variants AFM
    afm, variants = format_matrix(orig, CBC_GBC)

    # For all top3 analysis of that sample...:
    for analysis in top3_sample_variants[sample]:
        print(analysis)
        a_ = analysis.split('_')[:-1]
        filtering = a_[0]  
        min_cell_number = int(a_[1])
        min_cov_treshold = int(a_[2])

        # Filter cells and vars
        if filtering in ['CV', 'ludwig2019', 'velten2021', 'miller2022']:
            # Cells
            a_cells = filter_cells_coverage(afm, mean_coverage=min_cov_treshold) 
            if min_cell_number > 0:
                cell_counts = a_cells.obs.groupby('GBC').size()
                clones_to_retain = cell_counts[cell_counts>min_cell_number].index 
                cells_to_retain = a_cells.obs.query('GBC in @clones_to_retain').index
                a_cells = a_cells[cells_to_retain, :].copy()
            # Variants
            if filtering == 'CV':
                a = filter_CV(a_cells, n=50)
            elif filtering == 'ludwig2019':
                a = filter_ludwig2019(a_cells, mean_AF=0.5, mean_qual=0.2)
            elif filtering == 'velten2021':
                a = filter_velten2021(a_cells, mean_AF=0.1, min_cell_perc=0.2)
            elif filtering == 'miller2022':
                a = filter_miller2022(a_cells, mean_coverage=100, mean_qual=0.3, perc_1=0.01, perc_99=0.1)

        elif filtering == 'density':
            a = filter_density(afm, density=0.5, steps=np.Inf)
            if min_cell_number > 0:
                cell_counts = a.obs.groupby('GBC').size()
                clones_to_retain = cell_counts[cell_counts>min_cell_number].index 
                cells_to_retain = a.obs.query('GBC in @clones_to_retain').index
                a = a[cells_to_retain, :].copy()

        # Control vars...
        print(analysis)
        print(a)
        assert all([ var in top3_sample_variants[sample][analysis] for var in a.var_names ])

        # 1-Viz selected variants properties
        fig, axs = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)

        # Set colors 
        colors = {'non-selected':'grey', 'selected':'red'}

        # To nans
        to_plot = a_cells.copy()
        to_plot.X[np.isnan(to_plot.X)] = 0

        # Vafs ditribution
        for i, var in enumerate(a_cells.var_names):
            x = to_plot.X[:, i]
            x = np.sort(x)
            if var in a.var_names:
                axs[0].plot(x, '--', color=colors['selected'], linewidth=0.5)
            else:
                axs[0].plot(x, '--', color=colors['non-selected'], linewidth=0.2)

        format_ax(pd.DataFrame(x), ax=axs[0], title='Ranked AFs', xlabel='Cell rank', ylabel='AF')

        # Vafs summary stats
        df_ = summary_stats_vars(to_plot, variants=None).drop('median_coverage', axis=1).reset_index(
            ).rename(columns={'index' : 'variant'}).assign(
            is_selected=lambda x: np.where(x['variant'].isin(a.var_names), 'selected', 'non-selected')).melt(
            id_vars=['variant', 'is_selected'], var_name='summary_stat')

        #strip(df_, 'summary_stat', 'value', by='is_selected', s=2, c=colors, ax=axs[1])
        violin(df_, 'summary_stat', 'value', by='is_selected', c=colors, ax=axs[1])
        format_ax(df_, ax=axs[1], title='Summary statistics', 
            xticks=df_['summary_stat'].unique(), xlabel='', ylabel='Value'
        )
        handles = create_handles(colors.keys(), marker='o', colors=colors.values(), size=10, width=0.5)
        axs[1].legend(handles, colors.keys(), title='Selection', loc='center left', 
            bbox_to_anchor=(1, 0.5), ncol=1, frameon=False
        )
        fig.suptitle(f'{sample}: analysis {analysis}')

        # Save
        fig.savefig(f'{analysis}_variants.pdf')
        
        # 2-Viz cell x var and cell x cell heatmaps
        with PdfPages(f'{sample}_{analysis}_heatmaps.pdf') as pdf:

            a = nans_as_zeros(a)
            clone_colors = create_palette(a.obs, 'GBC', palette='dark')
            cell_anno_clones = [ clone_colors[clone] for clone in a.obs['GBC'] ]

            # Viz 
            g = cells_vars_heatmap(a, cell_anno=cell_anno_clones, anno_colors=clone_colors, 
                heat_label='AF', legend_label='Clone', figsize=(11, 8), title=f'{sample}: {analysis}'
            )
            pdf.savefig() 

            # 3-Viz all cell x cell similarity matrices obtained from the filtered AFM one.
            for x in os.listdir(path_distances):
                if bool(re.search(f'{sample}_{"_".join(analysis.split("_")[:-1])}', x)):
                    a_ = x.split('_')[:-1]
                    metric = a_[-1]
                    with_nans = 'w/i nans' if a_[-2] == 'yes' else 'w/o nans'
                    D = sc.read(path_distances + x)

                    assert (a.obs_names == D.obs_names).all()
                    D.obs['GBC'] = a.obs['GBC']

                    # Draw clustered similarity matrix heatmap 
                    heat_title = f'{sample} clones: {filtering}_{min_cell_number}_{min_cov_treshold}, {metric} {with_nans}'
                    g = cell_cell_dists_heatmap(D, cell_anno=cell_anno_clones, anno_colors=clone_colors, 
                        heat_label='Similarity', legend_label='Clone', figsize=(11, 6.5), 
                        title=heat_title
                    )
                pdf.savefig() 

            plt.close()
############## 


##


############## 
# For each sample (3x) clones, what are the clones which are consistently predictable in the top analyses? 
top_clones = {}
for sample in clones['sample'].unique():
    top = top_3[sample]
    top_clones[sample] = clones.query('sample == @sample and analysis in @top').groupby(['comparison']).agg(
        {'f1':np.median}).sort_values(
        'f1', ascending=False).query('f1 > 0.5').index.to_list()
print(f'Top clones: {top_clones}')

# Load top3 analyses variants for each sample, and visualize their intersection (i.e., J.I.)
D = {}
for sample in sample_names:
    var_dict = {}
    for x in os.listdir(path_results + f'top_3/{sample}/'):
        if x.endswith('.xlsx'):
            n = '_'.join(x.split('.')[0].split('_')[2:-1])
            df_ = pd.read_excel(path_results + f'top_3/{sample}/{x}', index_col=0)

            df_ = df_.loc[df_['comparison'].str.contains('|'.join(top_clones[sample]))]
            var_dict[n] = df_.index.to_list()
    D[sample] = var_dict

# Sample a
fig, axs = plt.subplots(1, 3, figsize=(10,5))

for k, sample in enumerate(D):
    n_analysis = len(D[sample].keys())
    JI = np.zeros((n_analysis, n_analysis))
    for i, l1 in enumerate(D[sample]):
        for j, l2 in enumerate(D[sample]):
            x = D[sample][l1]
            y = D[sample][l2]
            JI[i, j] = ji(x, y)
    JI = pd.DataFrame(data=JI, index=None, columns=D[sample].keys())

    plot_heatmap(JI, palette='mako', ax=axs[k], title=sample, y_names=False,
        x_names_size=10, y_names_size=0, annot=True, annot_size=10, cb=True, label='JI variants'
    )

fig.tight_layout()
fig.savefig(path_results + 'overlap_selected_vars_only_top_clones.pdf')
################