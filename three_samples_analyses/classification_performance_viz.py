"""
Visualization of clones and samples classification performances.
"""

# Code
import pickle
import re
import os
import sys
import gc
from Cellula.plotting._plotting import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from MI_TO.preprocessing import *
from MI_TO.diagnostic_plots import sturges
from MI_TO.heatmaps_plots import *
from MI_TO.utils import *
from MI_TO.diagnostic_plots import *


##


# Set paths
path_main = sys.argv[1]
sample_names = sys.argv[2].split(':')

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

ax.text(0.42, 0.32, f'-Total analyses for clones (samples): {clones["analysis"].unique().size} ({samples["analysis"].unique().size})', transform=ax.transAxes)
ax.text(0.42, 0.29, f'-Total comparisons for clones (samples): {clones["comparison"].unique().size} ({samples["comparison"].unique().size})', transform=ax.transAxes)
n_clones = clones.loc[:, ['sample', 'comparison']].drop_duplicates().groupby('sample').count()['comparison']
ax.text(0.42, 0.26, f'-n clones by sample: MDA {n_clones["MDA"]}; AML {n_clones["AML"]}; PDX {n_clones["PDX"]};', transform=ax.transAxes)
ax.text(0.42, 0.23, f'-Median n of analyses a comparison occurred in,', transform=ax.transAxes)
ax.text(0.42, 0.20, f' for clones (samples): {np.median(clones.groupby("comparison").size())} ({np.median(samples.groupby("comparison").size())})', transform=ax.transAxes)
ax.text(0.42, 0.17, f'-Median f1-scores for clones (samples):', transform=ax.transAxes)
ax.text(0.42, 0.14, f' {np.median(clones["f1"]):.3f} ({clones["f1"].size} values across all analyses)', transform=ax.transAxes)
ax.text(0.42, 0.11, f' {np.median(samples["f1"]):.3f} ({samples["f1"].size} values across all analyses)', transform=ax.transAxes)

pdx_y = np.median(samples.query('comparison == "PDX_vs_rest"')['f1'])
ax.annotate('PDX', xy=(0.935, pdx_y-0.001), xytext=(0.7, 0.95), arrowprops={"arrowstyle":"->", "color":"black"})
aml_y = np.median(samples.query('comparison == "AML_vs_rest"')['f1'])
ax.annotate('AML', xy=(1.03, aml_y), xytext=(1.1, 0.6), arrowprops={"arrowstyle":"->", "color":"black"})
mda_y = np.median(samples.query('comparison == "MDA_vs_rest"')['f1'])
ax.annotate('MDA', xy=(0.99, mda_y), xytext=(1.15, 0.87), arrowprops={"arrowstyle":"->", "color":"black"})

top_3_clones = clones_agg.sort_values('f1', ascending=False).head(3)
clone_y = top_3_clones['f1'][0]
ax.annotate(top_3_clones.index[0].split('_')[0], xy=(-0.03, clone_y+0.01), xytext=(-0.25, 0.82), arrowprops={"arrowstyle":"->", "color":"black"})
clone_y = top_3_clones['f1'][1]
ax.annotate(top_3_clones.index[1].split('_')[0], xy=(-0.08, clone_y-0.01), xytext=(-0.45, 0.53), arrowprops={"arrowstyle":"->", "color":"black"})
clone_y = top_3_clones['f1'][2]
ax.annotate(top_3_clones.index[2].split('_')[0], xy=(0.09, clone_y), xytext=(0.1, 0.65), arrowprops={"arrowstyle":"->", "color":"black"})

# Save
fig.savefig(path_results + 'median_f1_by_task.pdf')
############## 


##


############## f1 by clone and feat_type
fig, ax = plt.subplots(figsize=(12, 5))

feat_type_colors = create_palette(clones, 'feature_type', 'Set1')
params = {   
            'showcaps' : True,
            'fliersize': 0,
            'boxprops' : {'edgecolor': 'black', 'linewidth': 0.3}, 
            'medianprops': {"color": "black", "linewidth": 1},
            'whiskerprops':{"color": "black", "linewidth": 1}
        }
box(clones, 'comparison', 'f1', c='#E9E7E7', ax=ax, params=params)
strip(clones, 'comparison', 'f1', by='feature_type', c=feat_type_colors, s=2, ax=ax)
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

ax.text(0.25, 0.8, f'-n clones with more than one classification above 0.5 f1: 6', transform=ax.transAxes)
ax.text(0.25, 0.75, f'-n clones with more than one classification above 0.8 f1: 6', transform=ax.transAxes)
ax.text(0.25, 0.7, f'-n clones with median f1 > 0.5: 1', transform=ax.transAxes)

# Save
fig.tight_layout()
fig.savefig(path_results + 'clones_f1.pdf')
##############


##


############## f1 by sample and feat_type
fig, ax = plt.subplots(figsize=(8, 6.8))

box(clones, 'sample', 'f1', ax=ax, s=0.5, c='#E9E7E7', params=params)
strip(clones, 'sample', 'f1', by='feature_type', c=feat_type_colors, s=3, ax=ax)
format_ax(clones, ax, title='Clones f1-scores by sample and variant selection method', ylabel='f1')
create_handles(feat_type_colors.keys(), marker='o', colors=None, size=10, width=0.75)
handles = create_handles(feat_type_colors.keys(), colors=feat_type_colors.values())
fig.subplots_adjust(right=0.75)
fig.legend(handles, feat_type_colors.keys(), loc='center', 
    bbox_to_anchor=(0.85, 0.6), ncol=1, frameon=False, title='Feature selection'
)
ax.text(1.05, 0.4, f'Mean MDA: {clones.loc[clones["sample"]=="AML"]["f1"].mean():.3f}', transform=ax.transAxes)
ax.text(1.05, 0.36, f'Mean AML: {clones.loc[clones["sample"]=="MDA"]["f1"].mean():.3f}', transform=ax.transAxes)
ax.text(1.05, 0.32, f'Mean PDX: {clones.loc[clones["sample"]=="PDX"]["f1"].mean():.3f}', transform=ax.transAxes)

# Save
fig.savefig(path_results + 'clones_f1_by_sample.pdf')
##############


##



############## f1 by model and feat_type
fig, ax = plt.subplots(figsize=(8, 6.5))

box(clones, 'model', 'f1', ax=ax, c='#E9E7E7', params=params)
strip(clones, 'model', 'f1', by='feature_type', c=feat_type_colors, s=3, ax=ax)
format_ax(clones, ax, title='Clones f1-scores by model and variant selection method', ylabel='f1')
create_handles(feat_type_colors.keys(), marker='o', colors=None, size=10, width=0.5)
handles = create_handles(feat_type_colors.keys(), colors=feat_type_colors.values())
fig.subplots_adjust(right=0.75)
fig.legend(handles, feat_type_colors.keys(), loc='center', 
    bbox_to_anchor=(0.85, 0.6), ncol=1, frameon=False, title='Feature selection'
)
ax.text(1.05, 0.39, f'Mean xgboost: {clones.loc[clones["model"]=="xgboost"]["f1"].mean():.3f}', transform=ax.transAxes)
ax.text(1.05, 0.36, f'Mean logit: {clones.loc[clones["model"]=="logit"]["f1"].mean():.3f}', transform=ax.transAxes)

# Save
fig.savefig(path_results + 'clones_f1_by_model.pdf')
##############


##


############## f1 by sample and feat_type
# Sizes
res = []
for x in os.listdir(path_main + '/data/CBC_GBC_cells/'):
    if x.endswith('csv'):
        d = pd.read_csv(path_main + f'/data/CBC_GBC_cells/{x}', index_col=0)
        res.append(d.assign(sample=x.split('_')[-1].split('.')[0]))
CBC_GBC = pd.concat(res, axis=0)
clones_sizes = CBC_GBC.groupby('GBC').size()

clones['GBC'] = clones['comparison'].map(lambda x: x.split('_')[0])

csizes = []
for x in clones['GBC']:
    csizes.append(clones_sizes[x])
clones['size'] = csizes

# Viz
fig, ax = plt.subplots(figsize=(6, 6))
scatter(clones, 'size', 'f1', c='#606060', s=3, ax=ax)
x = clones['size']
y = clones['f1']
fitted_coefs = np.polyfit(x, y, 1)
y_hat = np.poly1d(fitted_coefs)(x)
ax.plot(x, y_hat, linestyle='dotted', linewidth=2, color='r')
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

#with open(path_clones + 'top3.pkl', 'wb') as f:
#    pickle.dump(top_3, f)

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
    
    if not os.path.exists(sample):
        os.mkdir(sample)
    os.chdir(sample)

    # Read data and create colors
    print(sample)
    afm = read_one_sample(path_main, sample)
    clone_colors = create_palette(afm.obs, 'GBC', palette=sc.pl.palettes.default_20)
    gc.collect()
 
    # For all top3 analysis of that sample...:
    for analysis in top3_sample_variants[sample]:

        a_ = analysis.split('_')[:-1]
        filtering = a_[0]  
        min_cell_number = int(a_[1])
        min_cov_treshold = int(a_[2])
 
        # Filter cells and vars
        a_cells, a = filter_cells_and_vars(
            afm, 
            filtering=filtering,
            min_cell_number=min_cell_number,
            min_cov_treshold=min_cov_treshold,
            path_=path_results
        )
        gc.collect()
 
        # Control vars...
        print(analysis)
        print(a)
        assert all([ var in top3_sample_variants[sample][analysis] for var in a.var_names ])

        # Get info!
        if not os.path.exists(path_results + f'top_3/{sample}/{analysis}/'):
            os.mkdir(path_results + f'top_3/{sample}/{analysis}/')
            os.chdir(path_results + f'top_3/{sample}/{analysis}/')
    
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
                cell_anno_clones = [ clone_colors[clone] for clone in a.obs['GBC'] ]

                # Viz 
                g = cells_vars_heatmap(a, cell_anno=cell_anno_clones, anno_colors=clone_colors, 
                    heat_label='AF', legend_label='Clone', figsize=(11, 8), title=f'{sample}: {analysis}'
                )

                # Prep d for savings
                analysis_d = {}
                analysis_d['cells'] = a.obs_names.to_list()
                analysis_d['vars'] = a.var_names.to_list()
                analysis_d['dendrogram'] = g.dendrogram_row.dendrogram
                analysis_d['linkage'] = g.dendrogram_row.linkage

                with open('cell_x_var_hclust.pickle', 'wb') as f:
                    pickle.dump(analysis_d, f)

                pdf.savefig() 

                # 3-Viz all cell x cell similarity matrices obtained from the filtered AFM one.
                for x in os.listdir(path_distances):

                    if bool(re.search(f'{sample}_{"_".join(analysis.split("_")[:-1])}_', x)):

                        print(x)
                        a_ = x.split('_')[:-1]
                        metric = a_[-1]
                        with_nans = 'w/i nans' if a_[-2] == 'yes' else 'w/o nans'
                        D = sc.read(path_distances + x)
                        gc.collect()

                        if a.shape[0] == D.shape[0]:
                            print(a)
                            print(D)
                            assert (a.obs_names == D.obs_names).all()
                            D.obs['GBC'] = a.obs['GBC']

                            # Draw clustered similarity matrix heatmap 
                            heat_title = f'{sample} clones: {filtering}_{min_cell_number}_{min_cov_treshold}, {metric} {with_nans}'
                            g = cell_cell_dists_heatmap(D, 
                                cell_anno=cell_anno_clones, anno_colors=clone_colors, 
                                heat_label='Similarity', legend_label='Clone', figsize=(11, 6.5), 
                                title=heat_title
                            )

                            analysis_d = {}
                            analysis_d['dendrogram'] = g.dendrogram_row.dendrogram
                            analysis_d['linkage'] = g.dendrogram_row.linkage

                            with open(f'similarity_{"_".join(a_)}_hclust.pickle', 'wb') as f:
                                pickle.dump(analysis_d, f)

                            pdf.savefig()

                        else:
                            print(f'{x} not added...')
                plt.close()
        
        else:
            print(f'Analysis {analysis} hclusts have been already computed...')
##############
    

##


############## 
# For each sample top3 analysis on the clone task, what is the number of selected variants?
d = {}
for sample in  top3_sample_variants:
    n_vars = {}
    for analysis in top3_sample_variants[sample]:
        n_vars[analysis] = len(top3_sample_variants[sample][analysis])
    d[sample] = n_vars
df_ = pd.DataFrame(d).reset_index().rename(columns={'index':'analysis'}).melt(
    id_vars='analysis', var_name='sample', value_name='n_vars').dropna()

# Viz 
colors = {'MDA':'#DA5700', 'AML':'#0074DA', 'PDX':'#0F9221'}
fig, ax = plt.subplots(figsize=(6,7))
bar(df_, 'n_vars', x=None, by='sample', c=colors, ax=ax, s=0.75, annot_size=10)
format_ax(df_, ax=ax, xticks=df_['analysis'], rotx=90, 
    ylabel='n variants', title='n variants selected by the top 3 analyses of each sample'
)
handles = create_handles(colors.keys(), marker='o', colors=colors.values(), size=10, width=0.5)
ax.legend(handles, colors.keys(), title='Sample', loc='upper right', 
    bbox_to_anchor=(0.25, 0.95), ncol=1, frameon=False
)
fig.tight_layout()
fig.savefig(path_results + 'n_top3_selected_variants.pdf')
################


##


############## 
# For each sample (3x) clones, what are the clones that are consistently predictable in the top analyses? 
# What are their features? 
top_clones_d = {}
for sample in clones['sample'].unique():
    top_analyses = top_3[sample]
    top_clones = clones.query('sample == @sample and analysis in @top_analyses').groupby(['comparison']).agg(
        {'f1':np.median}).sort_values(
        'f1', ascending=False).query('f1 > 0.5').index.to_list()
    top_clones_stats = clones.query(
        'sample == @sample and analysis in @top_analyses and comparison in @top_clones'
        )
    top_clones_d[sample] = {'clones' : top_clones, 'stats' : top_clones_stats}

print(f'Top clones: {top_clones_d}')

# Here we go...
for sample in sample_names: 

    if len(top_clones_d[sample]['clones']) > 0:

        afm = read_one_sample(path_main, sample=sample)
        best_for_each_top_clone_df = top_clones_d[sample]['stats'].sort_values(
            'f1', ascending=False).groupby('comparison').head(1).loc[
                :, ['feature_type', 'min_cell_number', 'min_cov_treshold', 'comparison', 'model']
            ].assign(
                clone=lambda x: x['comparison'].map(lambda y: y.split('_')[0])
            ).set_index('clone').drop('comparison', axis=1)

        # For each top classified clone in that sample, and its top analysis...
        for i in range(best_for_each_top_clone_df.shape[0]):

            # Get top clone id, and its top classification analysis options
            topper = best_for_each_top_clone_df.index[i]
            filtering = best_for_each_top_clone_df['feature_type'][i]
            min_cell_number = best_for_each_top_clone_df['min_cell_number'][i]
            min_cov_treshold = best_for_each_top_clone_df['min_cov_treshold'][i]
            model = best_for_each_top_clone_df['model'][i]

            # Viz
            print(topper)
            fig = viz_clone_variants(
                afm, topper, 
                sample=sample, path=path_clones, 
                filtering=filtering, min_cell_number=min_cell_number, 
                min_cov_treshold=min_cov_treshold, model=model, 
                figsize=(10,10)
            )

            fig.savefig(path_results + f'top_3/{sample}/{topper}_features.pdf')
################


##


############## f1 by clone and feat_type

# top_3['MDA']
# top_clones_d['MDA']['clones']
# 
# clones.query('analysis in @top_3["MDA"]')




















fig, ax = plt.subplots(figsize=(12, 5))
params = {   
            'showcaps' : True,
            'fliersize': 0,
            'boxprops' : {'edgecolor': 'black', 'linewidth': 0.3}, 
            'medianprops': {"color": "black", "linewidth": 1},
            'whiskerprops':{"color": "black", "linewidth": 1}
        }
box(clones, 'comparison', 'f1', c='#E9E7E7', ax=ax, params=params)
strip(clones, 'comparison', 'f1', by='feature_type', c=feat_type_colors, s=2, ax=ax)
format_ax(clones, ax, title='f1-scores by clone and variant selection method', rotx=90, xsize=5)
create_handles(feat_type_colors.keys(), marker='o', colors=None, size=10, width=0.5)
handles = create_handles(feat_type_colors.keys(), colors=feat_type_colors.values())
fig.legend(handles, feat_type_colors.keys(), loc='upper right', 
    bbox_to_anchor=(0.9, 0.9), ncol=2, frameon=False, title='Feature selection'
)

# Save
fig.tight_layout()
fig.savefig(path_results + 'clones_f1_only_top3_per_sample.pdf')
##############



