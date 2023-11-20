# Expansion Cooker

This tool is designed to process [Expansion Hunter](https://github.com/Illumina/ExpansionHunter/tree/master) output for the analysis of paired genotype differences. It facilitates the comparison of paired genotypes by identifying high-quality differences supported in the data and extracting significant features for downstream analysis.

## Features:

1. **Process Expansion Hunter Outputs**: Filters genotypes based on paired genotype quality and produces differences.
2. **Multiprocessing Support**: Utilizes all available CPU cores for faster processing.
3. **Feature Extraction**: Feature extraction from produced differences for further analysis.
4. **Tracking & Logging**: Tracks problematic cases and logs all events.


## Usage:
Setup guide, usage, and output information can be found in the [docs](docs/).

## Methods:
### Significant Difference Identification
A rough overview of the process for determing wether genotypes can be considered is as follows:
1. If the total read count is high for both case & control, trust the Expansion Hunter reported genotype.
2. If the total read count for case OR control is very low, skip.
3. If there is a big difference in read counts in case vs control, check support of reported genotypes.
	3.1 Check what genotypes in control are supported by the control reads.
	3.2 Check what genotypes in case + control are supported by case reads.
	3.3 Pair up the supported genotypes.
4. Otherwise, use Expansion Hunters confidence intervals to determine differences.

### Feature Extraction
Features of the differences which may be useful for idenitfying interesting loci are extracted. More information on the features can be found in the[outputs file](docs/1_Outputs.md).

  


