import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import os


# Locations of the data folders
diff_folder = os.getenv('DIFF_FOLDER') or "../data/CookerOut/"
genotypes_folder = os.getenv('GENO_FOLDER') or "../data/CookerOut/"


def load_genotypes(cols, disease):
    case_path = os.path.join(genotypes_folder, f'{disease}_case.csv')
    control_path = os.path.join(genotypes_folder, f'{disease}_control.csv')
    case = pd.read_csv(case_path, usecols=[cols])
    control = pd.read_csv(control_path, usecols=[cols])
    return case, control

### GENOTYPE GRAPHING FUNCTIONS ###

def graph_genotypes(locus, disease):
    """
    Graph the genotypes of a single locus for a single disease.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    case, control = load_genotypes(locus, disease)
    graph_genotypes_helper(case, control, locus, ax)
    plt.show()


def graph_genotypes_helper(case, control, name, ax):
    ax.hist(case, bins=100, alpha=0.5, label='case', color='red')
    ax.hist(control, bins=100, alpha=0.5, label='control', color='blue')
    ax.set_title(f'{name} Distribution')
    ax.set_xlabel('Genotype')
    ax.set_ylabel('Frequency')
    ax.legend(loc='upper right')


def load_graphGenotypesLoci(loci, disease):
    num_loci = len(loci)
    num_cols = 4
    num_rows = math.ceil(num_loci / num_cols)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(30, num_rows*5)) # Adjust the figure size as per your requirement
    axs = axs.ravel() # flatten the array of axes

    for idx, loc in enumerate(loci):
        ax = axs[idx]
        graph_genotypes_helper(loc, disease, ax)

    plt.title(f'{disease} Distribution for {num_loci} Loci')
    plt.tight_layout()  # To prevent overlap of subplots
    plt.show()

def calculate_data(loci, case_df, control_df):
    data_dict = {}
    max_diff = 0
    
    for i, loc in enumerate(loci):
        # Extract data for current loci and drop null rows
        df = pd.concat([control_df[[loc]], case_df[[loc]]], axis=1)
        df.columns = ['control', 'case']
        df.dropna(inplace=True)

        # Melt the dataframe to long format for plotting
        df_melt = df.melt(value_name='Repeat Length', var_name='group')
        bins = np.arange(df_melt['Repeat Length'].min() - 0.5, df_melt['Repeat Length'].max() + 0.5)

        
        # Get histogram data for case and control groups
        hist_case, case_bins = np.histogram(df['case'], bins=bins)
        hist_control, _ = np.histogram(df['control'], bins=bins)

        # Calculate the difference in bin counts
        diff = hist_case - hist_control

        # create a convolution of the difference array with a gaussian kernel
        # to smooth out the line plot
        conv_diff = np.convolve(diff, np.ones(9)/9, mode='same')

        # if magnitude of max of conv_diff is greater than max_diff, update max_diff
        if np.abs(conv_diff).max() > max_diff:
            max_diff = np.abs(conv_diff).max()

        if len(conv_diff) <= 9: 
            conv_diff = diff
        
        bins = bins - 0.5
        bins = bins[1:]
        scatter_data = pd.DataFrame({
            'bin_midpoints': bins,
            'diff': conv_diff
        })
        scatter_data = scatter_data.loc[scatter_data['diff'] != 0]

        data_dict[loc] = (df_melt, scatter_data)

    return data_dict, max_diff

def plot_data(data_dict, max_diff, num_cols, disease, motifs, sample_case, sample_control):
    num_loci = len(data_dict.keys())
    num_rows = np.ceil(num_loci / num_cols).astype(int)
    
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(24, 3+num_rows*5))
    axs = axs.ravel()

    plt.subplots_adjust(wspace=0.3, hspace=0.5)


    for i, (loc, (df_melt, scatter_data)) in enumerate(data_dict.items()):
        sns.histplot(data=df_melt, x='Repeat Length', hue='group', element='step', common_norm=False, binwidth=1, ax=axs[i])
    
        ax2 = axs[i].twinx()  # Create a second axes that shares the same x-axis
        ax2.set_ylim(-max_diff * 2, max_diff * 2)

        sns.lineplot(data=scatter_data, x='bin_midpoints', y='diff', ax=ax2, color='purple')
        ax2.axhline(0, color='gray', linestyle='--')
        ax2.set_ylabel('Smoothed Case - Control Difference')
        ax2.set_xlim(0, scatter_data['bin_midpoints'].max() + 2)

        legend_elements = [Patch(facecolor='b', edgecolor='b', label='Control', alpha=0.3),
                            Patch(facecolor='orange', edgecolor='orange', label='Case', alpha=0.3),
                            Line2D([0], [0], color='purple', label='Smoothed Difference'),
                            Line2D([0], [0], color='gray', linestyle='--', label='Zero Difference')]

        if sample_case is not None:
            vals = sample_case.loc[:, loc].dropna().values
            for val in vals:
                axs[i].axvline(val, color='red', linestyle='--')
            legend_elements.append(Line2D([0], [0], color='red', linestyle='--', label='Sample Case Length'))
            

        if sample_control is not None:
            vals = sample_control.loc[:, loc].dropna().values
            for val in vals:
                axs[i].axvline(val, color='blue', linestyle='--')
            legend_elements.append(Line2D([0], [0], color='blue', linestyle='--', label='Sample Control Length'))

        axs[i].legend(handles=legend_elements, loc='upper right')
        axs[i].set_title(f'{loc}, {motifs["chr" + loc]}')
        # add vertical line at sample repeat length
        

        axs[i].set_xlabel('Repeat Length')

    plt.suptitle(f'Genotypes of {disease} Across Loci')
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.show()

def graphLociGenotypes(loci, disease, case_df=None, control_df=None, sample=None, num_cols=4):
    motifs = get_motifs(loci)

    case_df = case_df or load_genotypes(loci.append('sample_id'), disease)
    sample_cases = case_df[case_df['sample_id'] == sample] if sample else None
    sample_controls = control_df[control_df['sample_id'] == sample] if sample else None

    data_dict, max_diff = calculate_data(loci, case_df, control_df)
    plot_data(data_dict, max_diff, num_cols, disease, motifs, sample_cases, sample_controls)

# def _graph_multi_genotypes(title, subvalues, case_df, control_df, sample=None, num_cols=4):



def get_motifs(loci):
    df = pd.read_csv('data/other/locus_structures.csv')

    # add 'chr' to start of loci
    loci = ['chr' + loc for loc in loci]

    # return dict of loci in loci and their motifs
    df = df.loc[df['ReferenceRegion'].isin(loci)]
    return dict(zip(df['ReferenceRegion'], df['LocusStructure']))


def graphdiseasesGenotypes(locus, diseases):
    """
    Plots the genotypes of a single locus for each disease type.
    args:
        locus: the locus to plot
        diseases: a list of disease types to plot
    """
    num_diseases = len(diseases)
    num_cols = 4
    num_rows = np.ceil(num_diseases / num_cols).astype(int)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, num_rows*5))  # Adjust the figure size as per your requirement
    axs = axs.ravel()  # flatten the array of axes

    # Adjust padding between subplots
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    plt.subplots_adjust(top=0.95)  # adjust top margin

    for i, disease in enumerate(diseases):
        # Read in the case and control data for each disease
        case_df, control_df = load_genotypes(locus, disease)

        # Combine the case and control dataframes along columns
        df = pd.concat([case_df, control_df], axis=1)
        df.columns = ['case', 'control']
        df.dropna(inplace=True)

        # Melt the dataframe to long format for plotting
        df_melt = df.melt(value_name='Repeat Length', var_name='group')
        bins = np.arange(df_melt['Repeat Length'].min() - 0.5, df_melt['Repeat Length'].max() + 0.5)
        # Create the overlaid histograms
        sns.histplot(data=df_melt, x='Repeat Length', hue='group', element='step', common_norm=False, bins=bins, ax=axs[i])

        
        # Get histogram data for case and control groups
        hist_case, bins_case = np.histogram(df['case'], bins=bins)
        hist_control, bins_control = np.histogram(df['control'], bins=bins)

        # Calculate the difference in bin counts
        diff = hist_case - hist_control

        max_hist = np.max([hist_case.max(), hist_control.max()])
        ax2 = axs[i].twinx()  # Create a second axes that shares the same x-axis
        ax2.set_ylim(-max_hist/2, max_hist/2)  # Set the limits of the y-axis to the maximum difference

        # Prepare data for the smoothed line plot
        scatter_data = pd.DataFrame({
            'bin_midpoints': (bins_case[:-1] + bins_case[1:])/2,
            'diff': diff
        })
        scatter_data = scatter_data.loc[scatter_data['diff'] != 0]

        sns.regplot(data=scatter_data, x='bin_midpoints', y='diff', order=4, ax=ax2, color='purple', scatter=False, ci = None, truncate=True)
        ax2.set_ylabel('Difference in Bin Counts')
        ax2.set_xlim(0, scatter_data['bin_midpoints'].max() + 2)

        axs[i].set_title(locus)
        axs[i].set_xlabel('Repeat Length')

    plt.suptitle(f'Genotype Distributions of {locus} across diseases')
    plt.tight_layout()
    plt.show()





### DIFFERENCE GRAPHING FUNCTIONS ###

def graphDiff(loc, disease):
    path = os.path.join(diff_folder, f'{disease}_diff.csv')
    dat = pd.read_csv(path, usecols=[loc])
    dat[loc].hist(bins=100)
    plt.title(f'{disease} ; {loc} Distribution')
    plt.xlabel('Tumor - Normal')
    plt.ylabel('Frequency')
    # Comment out not to save the figure
    # plt.savefig(f'figs/{disease}_{loc}.png')
    plt.show()

def graphDiffLociHelper(loc, disease, ax):
    path = os.path.join(diff_folder, f'{disease}_diff.csv')
    dat = pd.read_csv(path, usecols=[loc])
    ax.hist(dat[loc], bins=100)
    ax.set_title(f'{loc} Distribution')
    ax.set_xlabel('Tumor - Normal')
    ax.set_ylabel('Frequency')
    

def graphLociDiffs(loci, disease):
    num_loci = len(loci)
    num_cols = 4
    num_rows = math.ceil(num_loci / num_cols)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(30, num_rows*5)) # Adjust the figure size as per your requirement
    axs = axs.ravel() # flatten the array of axes

    for idx, loc in enumerate(loci):
        ax = axs[idx]
        graphDiffLociHelper(loc, disease, ax)

    plt.title(f'{disease} Differences (Tumor-Normal) Distributions')
    plt.tight_layout()  # To prevent overlap of subplots
    plt.show()






def graphLoci(loc, diseases):
    # Assuming diseases is a list of your diseases
    num_diseases = len(diseases)
    num_cols = 5
    num_rows = math.ceil(num_diseases / num_cols)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(30, num_rows*5)) # Adjust the figure size as per your requirement
    axs = axs.ravel() # flatten the array of axes

    for idx, disease in enumerate(diseases):
        ax = axs[idx]
        try:
            graphLocus(loc, disease, ax)
        except:
            print(f'Could not graph {disease}')
            continue

    plt.tight_layout()  # To prevent overlap of subplots
    plt.show()

def graphLocus(loc, disease, ax):
    path = os.path.join(diff_folder, f'{disease}_diff.csv')
    dat = pd.read_csv(path, usecols=[loc])
    ax.hist(dat[loc], bins=100)
    ax.set_title(f'{disease} Distribution')
    ax.set_xlabel('Tumor - Normal')
    ax.set_ylabel('Frequency')

