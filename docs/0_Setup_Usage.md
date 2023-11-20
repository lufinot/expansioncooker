# Setup 

1. Ensure you have python 3.6 or higher installed.
2. Create a python environment and install the required packages from requirements.txt.
    python3 venv -m <env_name>
    source <env_name>/bin/activate
    pip install -r requirements.txt
3. Create a config file (optional).
    Example config file with defaults is provided, but not required.
4. Set environment variables for logging (optional).
    export DISEASE=<disease_name> (default datetime)
    export LOG_LEVEL=<log_level> (default INFO)

## Usage

    python ExpansionCooker.py <RawDir> <Manifest> -n <DiseaseName> -o <OutputDir> [-f]

Using the example folder provided:

    python ExpansionCooker.py example/jsons/ example/manifest.csv -n Test -o example/out/ -f
    
Due to the nature and size of the data, the example files are blank, and outputs are not shown.

### Arguments:

- `RawDir`: Directory containing Expansion Hunter output JSONs.
- `Manifest`: Manifest file with case and control object ids.
- `-n` or `--name`: Disease name for output files.
- `-o` or `--outdir`: Output directory (Default: current directory).
- `-f` or `--feats`: Flag to indicate if features should be created from the output (Default: False).

## Outputs
Information on the outputs can be found in the [Outputs](link) file.

Logging information is saved to specified log directory or ExpansionCookerLogs. 

# Input Files Details

Check the example folder for example input files.

### Manifest File
The manifest file is a CSV file with the following columns: `{donor_id,case_object_id,control_object_id}`.
- donor_id: unique identifier for the donor. 
- case_object_id and control_object_id: prefixes of the jsons for the case and control Expansion hunter output for the donor. 

### Expansion Hunter Output
The output JSONs directly from Expansion Hunter.

### Config File Parameters

- 'LOG_DIR' : Directory for logging (default ExpansionCookerLogs in script running directory)
- 'MIN_COV' : Minimum coverage for a locus to be considered (default 6)
- 'HIGH_COV' : Minimum coverage for a locus to be considered well sampled (default 24)
- 'MAX_WIDTH' : Maximum confidence interval width for locus genotype (default 2)




