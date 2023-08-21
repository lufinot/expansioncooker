This tool is designed to process Expansion Hunter output for the analysis of paired genotype differences. It facilitates the comparison of paired genotypes through the filtering of high-quality differences and feature extraction.

## Features:

1. **Process Expansion Hunter Outputs**: The tool processes paired genotype data from Expansion Hunter's JSON outputs.
2. **Log Differences**: Logs differences in a structured format, making it easier to identify significant variations.
3. **Multiprocessing Support**: Utilizes all available CPU cores for faster processing.
4. **Feature Extraction**: Can interface with a feature extraction tool to further analyze and characterize the differences identified.

## Prerequisites:

- Python 3.x

## Usage:

	python ExpansionCooker.py <RawDir> <Manifest> -n <DiseaseName> -o <OutputDir> [-f]

### Arguments:

- `RawDir`: Directory containing Expansion Hunter output JSONs.
- `Manifest`: Manifest file with case and control object ids.
- `-n` or `--name`: Disease name for output files.
- `-o` or `--outdir`: Output directory (default is current directory).
- `-f` or `--feats`: Flag to indicate if features should be created from the output. Default is False.

### Example:

```
python ExpansionCooker.py /path/to/raw_eh /path/to/manifest.csv -n Alzheimer -o /path/to/output_dir -f
```


## Outputs:

1. `DiseaseName_case.csv`: CSV file containing the case data.
2. `DiseaseName_control.csv`: CSV file containing the control data.
3. `DiseaseName_diff.csv`: CSV file containing the differences between the case and control data.
4. `DiseaseName_tracking.csv`: CSV file containing tracking information.
5. (Optional, if `-f` is specified) `DiseaseName_feats.csv`: CSV file containing extracted features from the differences.

## Notes:

- Install packages from requirments.txt
- Check details for input specification
  
## Future Work:

- Improvement of filtering for higher resolution
- Automation of result ranking based on produced features

