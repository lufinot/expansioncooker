from datetime import datetime
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
from dbscan1d.core import DBSCAN1D
import pyranges as pr

import argparse
import os
import logging
import warnings


# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

def split_to_bed(df: pd.DataFrame, colname) -> pd.DataFrame:
    """ Splits the given column into 3 columns with titles Chromosome, Start, End. """
    split_data = df[colname].str.split(pat=':', n=1, expand=True)
    df['Chromosome'] = split_data[0]
    start_end_data = split_data[1].str.split(pat='-', n=1, expand=True)
    df['Start'] = start_end_data[0]
    df['End'] = start_end_data[1]
    return df


def get_COSMIC_regions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Get COSMIC regions from the dataframe.

    Args:
        df: Dataframe with locations of mutations.
    Returns:
        df: Dataframe with added COSMIC annotations
    """
    df_bed_format = split_to_bed(df, 'ReferenceRegion')

    # Convert pandas DataFrames to PyRanges objects
    pyranges_df = pr.PyRanges(df_bed_format)

    # Load the COSMIC data
    path = os.path.join(script_dir, 'refFiles', 'COSMIC.csv')
    cosmic_data = pd.read_csv(path)
    pyranges_cosmic = pr.PyRanges(cosmic_data)

    # Join the data with COSMIC regions using a left join
    joined_pyranges = pyranges_df.join(pyranges_cosmic, how='left')
    joined_df = joined_pyranges.df

    # Create a 'COSMIC_loc' column combining 'Start' and 'End'
    joined_df['COSMIC_loc'] = joined_df['Start_b'].astype(str) + '-' + joined_df['End_b'].astype(str)

    # Drop unnecessary columns
    columns_to_drop = ['Start', 'End', 'Start_b', 'End_b', 'Chromosome']
    joined_df.drop(columns=columns_to_drop, inplace=True)

    # Rename columns for clarity
    columns_to_rename = {
        'Tumour Types(Somatic)': 'COSMIC_TumorType', 
        'Tissue Type': 'COSMIC_TissueType',
        'Gene Symbol': 'COSMIC_GeneSymbol', 
        'Tier': 'COSMIC_Tier'
    }
    joined_df.rename(columns=columns_to_rename, inplace=True)

    return joined_df


def calculate_wilcoxon_pvals(df: pd.DataFrame) -> np.ndarray:
    """
    Calculate Wilcoxon p-values for each column in a dataframe.

    Args:
        df: Dataframe with rows as samples and cols as regions.
    Returns:
        pvals: Array of p-values for each region.
    """
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("once")
        pvals = wilcoxon(df, nan_policy='omit', zero_method='pratt', axis=0)[1]
        for warning in w:
            if str(warning.message) == 'Exact p-value calculation does not work if there are zeros. Switching to normal approximation.':
                logging.warning('Exact p-value calculation switched to normal approximation due to presence of zeros.')
            elif str(warning.message) == 'Sample size too small for normal approximation.':
                logging.warning('Sample size is too small for normal approximation.')
    _, pvals_corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    return pvals, pvals_corrected



def process_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove columns with all zeros and sample_id column if present.

    Args:
        df: Dataframe with rows as samples and cols as regions.
    Returns:
        df: Processed dataframe.
    """
    if 'sample_id' in df.columns:
        df.drop(columns=['sample_id'], axis=1, inplace=True)
    df = df.loc[:, (df.fillna(0) != 0).any(axis=0)]
    return df


def add_motif_info(df: pd.DataFrame) -> pd.DataFrame:
    """ 
    Add motif information to the dataframe, from the locus_structures.csv file.

    Args:
        df: Dataframe with rows as samples and cols as regions.
    Returns:
        df: Dataframe with motif information added.
    """
    path = os.path.join(script_dir, 'refFiles', 'locus_structures.csv')
    loc_to_motif = pd.read_csv(path)
    loc_to_motif['ReferenceRegion'] = loc_to_motif['ReferenceRegion'].str.lstrip('chr')
    df = df.merge(loc_to_motif[['ReferenceRegion', 'LocusStructure']], on='ReferenceRegion', how='left')
    return df


def cluster_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Call cluster function to get cluster features for each locus. 

    Args:
        df: Dataframe with rows as samples and cols as regions.
    Returns:
        cluster_df: With Cols ['num_clusters', 'cluster_means', 'cluster_sds', 'out3', 'out5']
    """
    result = df.apply(lambda x: cluster_and_outliers(x), axis=0).T
    result.columns = ['num_clusters', 'cluster_means', 'cluster_sds', 'out3', 'out5']
    return result

def cluster_and_outliers(x: pd.Series) -> list:
    x = x.dropna()
    n = len(x)

    # min samples is 2/3 of sqrt(num samples) or 4, whichever is larger
    n = max(int(np.sqrt(n)) * 2 / 3, 4)
    sx = x.div(x.std())
    sx = sx.values
    db = DBSCAN1D(eps=2, min_samples=n)
    labels = db.fit_predict(x.values)

    outliers = sx[labels == -1]

    # get counts of outliers above 3 and 5 sd
    out3 = len(outliers[outliers > 3])
    out5 = len(outliers[outliers > 5])

    # if no clusters, return 
    if len(np.unique(labels)) == 1:
        return [0, [], [], out3, out5]

    # get mean and sd of clusters
    means = []
    sds = []
    for label in np.unique(labels)[1:]:
        means.append(np.mean(x[labels == label]))
        sds.append(np.std(x[labels == label]))

    # remove from means and sd list the cluster with the mean closest to 0
    df = pd.DataFrame({
        'means': means,
        'sds': sds
    })
    df['abs_means'] = df['means'].abs()
    df = df.sort_values('abs_means')
    # if df.iloc[0]['abs_means'] < 3:
    #     df = df.iloc[1:, :]
 
    means = df['means'].values.tolist()  # Convert np.array to list
    sds = df['sds'].values.tolist()  # Convert np.array to list

    return [len(means), means, sds, out3, out5]  # Return a list instead of a tuple

    
def process_and_extract_features(df) -> pd.DataFrame:
    logging.debug("Processing and extracting features...")
    
    df = process_df(df)

    features_df = pd.DataFrame({'ReferenceRegion': df.columns})

    features_df['raw_pvals'], features_df['corrected_pvals'] = calculate_wilcoxon_pvals(df)
    logging.info("Calculated Wilcoxon p-values.")

    clusts = cluster_features(df)
    clusts = clusts.reset_index().rename(columns={'index': 'ReferenceRegion'})
    features_df = features_df.merge(clusts, how = 'left', on='ReferenceRegion')
    logging.info("Extracted cluster features.")

    props = df.astype(bool).sum(axis=0)/(df.shape[0]+(2*np.sqrt(df.shape[0])) + 3)
    features_df['prop_nonzero'] = props.values

    # count number of non-null values
    num = df.shape[0]
    features_df['counts'] = num - df.isnull().sum(axis=0).values
    logging.debug("Calculated proportion of non-zero values.")

    features_df['std'] = df.std().values
    logging.debug("Calculated standard deviation.")

    features_df = add_motif_info(features_df)
    logging.debug("Added motif information.")

    features_df = get_COSMIC_regions(features_df)
    logging.debug("Added COSMIC regions.")

    return features_df

def process_features(input_df: pd.DataFrame, name: str, outdir: str) -> None:

    feats_df = process_and_extract_features(input_df)

    output_path = os.path.join(outdir, f"{name}_feats.csv")
    feats_df.to_csv(output_path, index=False)
    logging.info(f"Feats saved to {output_path}")


def init_argparse():
    parser = argparse.ArgumentParser(description='Create features from ExpansionCooker output.')
    parser.add_argument('input', metavar='Diff_File', type=str, help='Location of the Expansion Cooker difference file')
    parser.add_argument('--name', '-n', help='Prefix for output file (default same as input file)')
    parser.add_argument('--outdir', '-o', default='', help='Output directory for the features. (default: script running directory)')
    return parser


def main():
    LOG_LEVEL = os.getenv('LOG_LEVEL') or 'info'
    log_dict = {'debug': logging.DEBUG, 'info': logging.INFO, 'warning': logging.WARNING, 
                'error': logging.ERROR, 'critical': logging.CRITICAL}
    log_level = log_dict.get(LOG_LEVEL.lower(), logging.INFO)

    slurm_job_id = os.getenv('SLURM_JOBID')
    if not slurm_job_id:
        slurm_job_id = datetime.now().strftime('%Y%m%d_%H%M')

    log_dir = 'ExpansionCookerLogs'
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    logging.basicConfig(filename=os.path.join(log_dir, f'{slurm_job_id}_FeatureExtractor.log'), 
                        level=log_level,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    
    parser = init_argparse()
    args = parser.parse_args()

    if not os.path.exists(args.input):
        logging.error(f"{args.input} does not exist.")
        return

    df = pd.read_csv(args.input, index_col=0)
    name = args.name or os.path.basename(args.input).split('.')[0]
    process_features(df, name, args.outdir)
 
if __name__ == '__main__':
    main()