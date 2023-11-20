
# Outputs 
***

- `{Name}_case.csv`: Case Genotypes.
- `{Name}_control.csv`: Control Genotypes.
- `{Name}_diff.csv`: Differences between case and control genotypes.
    - These are the genotypes that are considered to have a significant difference between the case and control. Rows are donors and columns are the loci.

- `{Name}_tracking.csv`: Problematic loci tracking.
    - This is a CSV containing the loci that were not included in the differences. This could be due to low coverage, low quality, or other issues.
    - Has the columns {'donor_id', 'ReferenceRegion', 'motif', 'issue'}. 

- `{Name}_feats.csv`: CSV file containing extracted features from the differences.
    - (Optional) if '-f' is specified.
    - This is a CSV containing the features extracted from the differences. Rows are loci and columns are the features.