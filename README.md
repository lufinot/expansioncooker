# Expansion Cooker

This tool is designed to process Expansion Hunter output for the analysis of paired genotype differences. It facilitates the comparison of paired genotypes through the filtering of high-quality differences and feature extraction.

## Features:

1. **Process Expansion Hunter Outputs**: Filters genotypes based on paired genotype quality and produces differences.
2. **Multiprocessing Support**: Utilizes all available CPU cores for faster processing.
3. **Feature Extraction**: Feature extraction from produced differences for further analysis.
4. **Tracking & Logging**: Tracks problamatic cases and logs all events.


## Usage:

	1. Create Environment, install packages from requirements.txt
	2. Create config file (optional)
	3. Set environment variables for logging (optional)
	4. Run ExpansionCooker.py 
		
		python ExpansionCooker.py <RawDir> <Manifest> -n <DiseaseName> -o <OutputDir> [-f]

### Arguments:

- `RawDir`: Directory containing Expansion Hunter output JSONs.
- `Manifest`: Manifest file with case and control object ids. (Must have columns: `{donor_id,case_object_id,control_object_id}`).
- `-n` or `--name`: Disease name for output files.
- `-o` or `--outdir`: Output directory (Default: current directory).
- `-f` or `--feats`: Flag to indicate if features should be created from the output (Default: False).

### Outputs
1. `Name_case.csv`: CSV file containing the case data.
2. `Name_control.csv`: CSV file containing the control data.
3. `Name_diff.csv`: CSV file containing the differences between the case and control data.
4. `Name_tracking.csv`: CSV file containing tracking information.
5. (Optional, if `-f` is specified) `Name_feats.csv`: CSV file containing extracted features from the differences.

Logging information is saved to specified log directory or ExpansionCookerLogs. 

## Additional Parameters

### Enviormnent Variables
- 'DISEASE' : Disease name for logging filename (default datetime)
- 'LOG_LEVEL' : Logging level (default INFO)

### Config File

- 'LOG_DIR' : Directory for logging (default ExpansionCookerLogs in script running directory)
- 'MIN_COV' : Minimum coverage for a locus to be considered (default 6)
- 'HIGH_COV' : Minimum coverage for a locus to be considered well sampled (default 24)
- 'MAX_WIDTH' : Maximum confidence interval width for locus genotype (default 2)

## Notes:

- Refer to example manifest file for formatting
- Example config file with defaults is provided, but not required

  


