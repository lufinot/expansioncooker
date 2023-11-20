
# Outputs 
### Genotypes Information
- `{Name}_case.csv`: Case Genotypes.
- `{Name}_control.csv`: Control Genotypes.
- `{Name}_diff.csv`: Differences between case and control genotypes.
    - These are the genotypes that are considered to have a significant difference between the case and control. Rows are donors and columns are the loci.

### Tracking Information
- `{Name}_tracking.csv`: Problematic loci tracking.
    - This is a CSV containing the loci that were not included in the differences. This could be due to low coverage, low quality, or other issues.
    - Has the columns {'donor_id', 'ReferenceRegion', 'motif', 'issue'}. 

### Feature Information
- `{Name}_feats.csv`: CSV file containing extracted features from the differences.
    - (Optional) if '-f' is specified.
    - This is a CSV containing the features extracted from the differences. Rows are loci and columns are the features.
The features extracted from the differences are as follows:
- Measures of skewness
    - Pvalues from wilcoxon signed rank test
- Non-Zero Cluster information
    - Clusters of differences in genotypes found using 1d clustering
- Outlier Information
    - Outliers in the differences
- Aditional Information
    - Counts of supported differences
    - Overlap with [COSMIC](https://cancer.sanger.ac.uk/cosmic) genes