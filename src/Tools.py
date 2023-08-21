import ExpansionFeatureExtractor as EHF
import pandas as pd
import numpy as np

def get_features(diffs):
    df = EHF.process_df(diffs)
    features_df =  pd.DataFrame({'ReferenceRegion': df.columns})

    features_df = EHF.add_motif_info(features_df)

    features_df['wilcox_pvals'] = EHF.calculate_wilcoxon_pvals(df)


    clusts = EHF.cluster_features(df)
    clusts = clusts.reset_index().rename(columns={'index': 'ReferenceRegion'})
    features_df = features_df.merge(clusts, how = 'left', on='ReferenceRegion')


    # get normalized non-zero proportion
    props = df.astype(bool).sum(axis=0)/(df.shape[0]+(2*np.sqrt(df.shape[0])))
    features_df['prop_nonzero'] = props.values

    # get standard deviation
    features_df['std'] = df.std().values
    return features_df

def getLocusWithOneCluster(features_df):
    one_clusts = features_df[features_df['num_clusters'] == 1]
    one_clusts = one_clusts[one_clusts['cluster_means'].apply(lambda x: (len(x) == 1 and 3 <= np.abs(x[0])))]
    return one_clusts


def loadGenotypes(cancer):
    ca_case = pd.read_csv(f'../data/genotypes/{cancer}_case_tidy.csv')
    ca_control = pd.read_csv(f'../data//genotypes/{cancer}_control_tidy.csv')
    return ca_case, ca_control

def loadDiffs(cancer):
    ca_diffs = pd.read_csv(f'../data/diffs/{cancer}_diffs.csv')
    return ca_diffs

def splitReferenceRegion(df):
    df['Chromosome'], rest = df['ReferenceRegion'].str.split(':', 1).str
    df['Start'], df['End'] = rest.str.split('-', 1).str