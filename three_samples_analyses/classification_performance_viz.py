"""
Visualization of clones and samples classification performances.
"""

# Code
import pickle
import os
from Cellula.plotting._plotting import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *
from MI_TO.diagnostic_plots import sturges
matplotlib.use('MacOSX')

# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/'
path_clones = path_main + '/results_and_plots/clones_classification/'
path_samples = path_main + '/results_and_plots/samples_classification/'

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
fig.savefig(path_main + '/results_and_plots/classification_performance/median_f1_by_task.pdf')
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
fig.savefig(path_main + '/results_and_plots/classification_performance/clones_f1.pdf')
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
fig.savefig(path_main + '/results_and_plots/classification_performance/clones_f1_by_sample.pdf')
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
fig.savefig(path_main + '/results_and_plots/classification_performance/clones_f1_by_model.pdf')
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
fig.savefig(path_main + '/results_and_plots/classification_performance/clones_size_f1_corr.pdf')
##############



##


############## 
# For each sample (3x) clones, what are the top 3 analyses (median f1 score across clones)? 
# Intersection among selected SNVs?? --> Take out from cluster

# Save top3 for easy quering on the cluster
top_3 = {}
for sample in clones['sample'].unique():
    top_3[sample] = clones.query('sample == @sample').groupby(['analysis']).agg(
        {'f1':np.median}).sort_values(
        'f1', ascending=False).index[:3].to_list()
top_3

# Load top3 variants for each sample clones, and visualize their intersection (i.e., J.I.), by sample
D = {}
for sample in os.listdir(path_clones + 'top3_analysis/'):
    var_dict = {}
    for x in os.listdir(path_clones + f'top3_analysis/{sample}/'):
        n = '_'.join(x.split('.')[0].split('_')[2:-1])
        df_ = pd.read_excel(path_clones + f'top3_analysis/{sample}/{x}', index_col=0)
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


fig.savefig(path_main + '/results_and_plots/classification_performance/overalp_selected_vars.pdf')

##############
