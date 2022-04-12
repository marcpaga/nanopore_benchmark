import os

from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches


def readoutcomes(df, output_path):
    
    output_file = os.path.join(output_path, 'combined', 'readoutcomes.pdf')
    metrics = list()
    for c in df.columns:
        if c != 'model':
            metrics.append(c)
    models_in_df = sorted(np.unique(df['model']))[::-1]

    fig, axes = plt.subplots(1, len(metrics), figsize=(3*len(metrics), 1.5*len(models_in_df)), sharey=True)

    for metric, ax in zip(metrics, axes):
        for j, model in enumerate(models_in_df):
            plot_df = df[df['model'] == model]

            try:
                v = int(plot_df[metric].item())
            except ValueError:
                v = 0
            v_perc = v/plot_df[metrics].sum(1).item()

            ax.barh(
                y = j, 
                width= v, 
                color = '#ccc5b4', 
                edgecolor = 'black'
            )
            
            ax.annotate(
                text = str(v)+"\n({:.1%})".format(v_perc), 
                xy = (plot_df[metric] - np.max(df[metric])*0.25, j), 
                verticalalignment = 'center', 
                horizontalalignment='left',
                multialignment = 'center',
            )

            if j == 0:
                ax.set_title(metric)
            
    axes[0].set_yticks(ticks = np.arange(0, len(set(df['model'])), 1))
    axes[0].set_yticklabels(labels = models_in_df)

    fig.suptitle('Basecalled reads failure rates')
    fig.supxlabel('Number of reads')
    fig.tight_layout()
    plt.savefig(output_file) 


def eventrates(df, output_path):

    output_file = os.path.join(output_path, 'combined', 'eventrates.pdf')

    metrics = ['match_rate', 'mismatch_rate', 'insertion_rate', 'deletion_rate']
    metric_name_transformer = {
        'match_rate': 'Match',
        'mismatch_rate': 'Mismatch',
        'insertion_rate': 'Insertion',
        'deletion_rate': 'Deletion'
    }

    models_in_df = sorted(np.unique(df['model']).astype(str))[::-1]
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(10, 1.5*len(models_in_df)), sharey=True)
    axes = [ax1, ax2, ax3, ax4]
   
    for i, model in enumerate(models_in_df):
        subdf = df[df['model'] == model]

        for metric, ax in zip(metrics, axes):
            stats = {k:v for k, v in zip(subdf['stats'], subdf[metric])}
            bp = ax.bxp(
                [stats],
                positions = [i], 
                showfliers = False, 
                vert = False, 
                widths = [0.6], 
                patch_artist = True
            )

            if i == 0:
                ax.set_xlabel(metric_name_transformer[metric])
                
            for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
                plt.setp(bp[element], color='black',  linewidth=1.5)
            plt.setp(bp['boxes'], facecolor = '#ccc5b4')

    ax1.set_yticks(ticks = np.arange(0, len(set(df['model'])), 1))
    ax1.set_yticklabels(labels = models_in_df)
    fig.suptitle('Alignment event rates')
    fig.tight_layout()
    plt.savefig(output_file) 


def homopolymerrates(df, output_path):

    output_file = os.path.join(output_path, 'combined', 'homopolymerrates.pdf')

    metrics = [
        'total_homo_error_rate', 
        'A_homo_error_rate', 
        'C_homo_error_rate', 
        'G_homo_error_rate', 
        'T_homo_error_rate'
    ]
    metric_name_transformer = {
        'total_homo_error_rate': 'Total',
        'A_homo_error_rate': 'A',
        'C_homo_error_rate': 'C',
        'G_homo_error_rate': 'G',
        'T_homo_error_rate': 'T'
    }
    
    models_in_df = sorted(np.unique(df['model']).astype(str))[::-1]
    fig, axes = plt.subplots(1, 5, figsize=(12.5, 1.5*len(models_in_df)), sharey=True)

    for i, model in enumerate(models_in_df):
        subdf = df[df['model'] == model]

        for metric, ax in zip(metrics, axes):
            stats = {k:v for k, v in zip(subdf['stats'], subdf[metric])}
            bp = ax.bxp(
                [stats], 
                positions = [i], 
                showfliers = False, 
                vert = False, 
                widths = [0.6], 
                patch_artist=True
            )

            if i == 0:
                ax.set_xlabel(metric_name_transformer[metric])
                
            for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
                plt.setp(bp[element], color='black',  linewidth=1.5)
            plt.setp(bp['boxes'], facecolor = '#ccc5b4')

    axes[0].set_yticks(ticks = np.arange(0, len(set(df['model'])), 1))
    axes[0].set_yticklabels(labels = models_in_df)
    fig.suptitle('Homopolymer Error rates')
    fig.tight_layout()
    plt.savefig(output_file) 

def phredq(df, output_path):
    
    output_file = os.path.join(output_path, 'combined', 'phredq.pdf')
    models_in_df = sorted(np.unique(df['model']))[::-1]

    colorpalettetab10 = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    fig, ax = plt.subplots(1, 1, figsize=(5, 1*len(models_in_df)), sharey=True)
    metrics = ['phred_mean_correct', 'phred_mean_error']
    for i, model in enumerate(models_in_df):
        subdf = df[df['model'] == model]

        for j, metric in enumerate(metrics):
            stats = {k:v for k, v in zip(subdf['stats'], subdf[metric])}
            bp = ax.bxp([stats], positions = [i + 1 + (j-1)/len(models_in_df)], showfliers = False, vert = False, widths = [0.3], patch_artist=True)

            for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
                plt.setp(bp[element], color='black',  linewidth=1)
            for patch in bp['boxes']:
                patch.set(facecolor=colorpalettetab10[j])       

    ax.set_xlabel('PhredQ Scores')
    ax.set_yticks(ticks = np.arange(1, len(set(df['model'])) + 1, 1))
    ax.set_yticklabels(labels = models_in_df)
    handles, _ = plt.gca().get_legend_handles_labels()
    for c, n in zip(colorpalettetab10, ['Correct', 'Incorrect']):
        handles.extend([mpatches.Patch(facecolor=c, label=n, edgecolor = 'black')]) 
    plt.legend(handles=handles, edgecolor = 'black', bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)

    fig.tight_layout()
    plt.savefig(output_file) 

def signatures(df, output_path):
    
    sns.set_palette(sns.color_palette("tab10"))

    for model in np.unique(df['model']):
        output_file = os.path.join(output_path, model, 'signatures.pdf')

        signatures_df = df[df['model'] == model].copy()
        signatures_df.loc[signatures_df['Context'] == 'NNN', 'Context'] = 'NZN'
        bases = ['A', 'C', 'G', 'T']

        f, axs = plt.subplots(1, 5, sharey=False, gridspec_kw={'width_ratios': [1, 1, 1, 1, 0.32]}, figsize =(15, 3))
        for i, b in enumerate(bases + ['Grouped']):
            plot_df = signatures_df[signatures_df["Base"] == b].copy()
            plot_df = plot_df[plot_df["Error"].isin(["Missmatch_A", "Missmatch_C", "Missmatch_G", "Missmatch_T", "Deletion", "Insertion"])]
            plot_df = plot_df.sort_values(['Error', 'Context'])
            plot_df.loc[plot_df['Context'] == 'NZN', 'Context'] = 'NNN'
            sns.histplot(x = "Context", hue = "Error", weights='Rate', data = plot_df, multiple='stack', ax = axs[i], line_kws=dict(linewidth=0.5))
            
            axs[i].set_xlabel(None)
            for tick in axs[i].get_xticklabels():
                tick.set_rotation(90)
                tick.set_fontfamily('monospace')
            if i == 0:
                axs[i].set_ylabel("Fraction")
            else:
                axs[i].set_ylabel(None)
            if i != 4:
                axs[i].get_legend().remove()
            else:
                legend = axs[i].get_legend()
                legend.set_title(None)
                legend.texts[0].set_text("Deletion")
                legend.texts[1].set_text("Insertion")
                legend.texts[2].set_text("Mismatch: A")
                legend.texts[3].set_text("Mismatch: C")
                legend.texts[4].set_text("Mismatch: G")
                legend.texts[5].set_text("Mismatch: T")
                legend.set_bbox_to_anchor((1, 0.75))
                
        f.tight_layout()
        plt.savefig(output_file) 

def auc(df, output_path):

    def integrate(x, y):
        sm = 0
        for i in range(1, len(x)):
            h = x[i] - x[i-1]
            sm += h * (y[i-1] + y[i]) / 2

        return sm
    
    output_file = os.path.join(output_path, 'combined', 'auc.pdf')
    colorpalettetab10 = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    f, ax1 = plt.subplots(figsize=(8, 5))
    aucs = dict()
    for j, modelname in enumerate(np.sort(np.unique(df['model']))):

        model_df = df[df['model'] == modelname]
        aucs[modelname] = -integrate(model_df['fraction'], model_df['match_rate'])

        if j == 0:
            ax2 = ax1.twinx()
            ax1.set_xlabel('Fraction of reads')
            ax1.set_ylabel('Average match rate')
            ax2.set_ylabel('Minimum read-mean PhredQ score')

        ax1.plot(model_df['fraction'], model_df['match_rate'], color = colorpalettetab10[j])
        ax2.plot(model_df['fraction'], model_df['phred_mean'], color = colorpalettetab10[j], linestyle='--')        
        
    handles, _ = plt.gca().get_legend_handles_labels()
    for c, modelname in zip(colorpalettetab10, np.sort(np.unique(df['model']))):
        handles.extend([Line2D([0], [0], label= modelname + ' ('+ str(round(aucs[modelname], 3)) +')', color=c)]) 
    plt.legend(handles=handles, edgecolor = 'black', bbox_to_anchor=(1.2,0.5), loc="center left", borderaxespad=0)

    plt.tight_layout()
    plt.savefig(output_file) 

    for modelname in np.sort(np.unique(df['model'])):
        
        output_file = os.path.join(output_path, modelname, 'auc.pdf')
        _, ax1 = plt.subplots(figsize=(8, 5))
        model_df = df[df['model'] == modelname]
        auc = -integrate(model_df['fraction'], model_df['match_rate'])

        ax2 = ax1.twinx()
        ax1.set_xlabel('Fraction of reads')
        ax1.set_ylabel('Average match rate')
        ax2.set_ylabel('Minimum read-mean PhredQ score')

        ax1.plot(model_df['fraction'], model_df['match_rate'], color = 'black')
        ax2.plot(model_df['fraction'], model_df['phred_mean'], color = 'red', linestyle='--')        
        
        ax1.set_title('AUC: ' + str(round(auc, 3)))
        plt.savefig(output_file) 