This tool is designed to process Expansion Hunter output for the analysis of paired genotype differences. It facilitates the comparison of paired genotypes through the filtering of high-quality differences and feature extraction.

## Features:

1. **Process Expansion Hunter Outputs**: Filters genotypes based on paired genotype quality and produces differences.
2. **Multiprocessing Support**: Utilizes all available CPU cores for faster processing.
3. **Feature Extraction**: Feature extraction from produced differences for further analysis.
4. **Tracking & Logging**: Tracks problamatic cases and logs all events.


## Usage:

	python ExpansionCooker.py <RawDir> <Manifest> -n <DiseaseName> -o <OutputDir> [-f]

### Arguments:

- `RawDir`: Directory containing Expansion Hunter output JSONs.
- `Manifest`: Manifest file with case and control object ids.
- `-n` or `--name`: Disease name for output files.
- `-o` or `--outdir`: Output directory (default is current directory).
- `-f` or `--feats`: Flag to indicate if features should be created from the output. Default is False.


1. `Name_case.csv`: CSV file containing the case data.
2. `Name_control.csv`: CSV file containing the control data.
3. `Name_diff.csv`: CSV file containing the differences between the case and control data.
4. `Name_tracking.csv`: CSV file containing tracking information.
5. (Optional, if `-f` is specified) `Name_feats.csv`: CSV file containing extracted features from the differences.

Logging information is saved in ExpansionCookerLogs directory in script running directory.

## Notes:

- Install packages from requirments.txt
- Check usage details for input specification and examples.
  


