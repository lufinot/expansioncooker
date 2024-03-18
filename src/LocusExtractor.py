import os
import ijson
import pandas as pd
import argparse

def extract_locus_data(file_path, locus_name):
    """
    Extracts data for a specific locus from a JSON file without loading the entire file into memory.

    Parameters:
    file_path (str): The path to the JSON file.
    locus_name (str): The name of the locus to extract.

    Returns:
    dict: Data for the specified locus, or None if the locus is not found.
    """
    with open(file_path, 'rb') as file:  # Use 'rb' mode for binary file access
        # ijson will parse the file incrementally
        found = False
        for prefix, event, value in ijson.parse(file):
            current_variant_id = ''
            # Check for the specific path in the JSON structure
            if (prefix, event) == ('LocusResults.item.Variants.item.VariantId', 'string'):
                current_variant_id = value
            if prefix.endswith(current_variant_id + '.Genotype'):
                genotype = value
            if prefix.endswith(current_variant_id + '.CountsOfFlankingReads'):
                flanking_counts = value
            if prefix.endswith(current_variant_id + '.CountsOfInrepeatReads'):
                inrepeat_counts = value
            if prefix.endswith(current_variant_id + '.CountsOfSpanningReads'):
                spanning_counts = value
            if prefix.endswith('.ReferenceRegion') and value == locus_name:
                found = True
                break

    if found:
        return {
            "CountsOfFlankingReads": flanking_counts,
            "CountsOfInrepeatReads": inrepeat_counts,
            "CountsOfSpanningReads": spanning_counts,
            "Genotype": genotype
        }
    else:
        return None

def process_donor(donor, raw_eh_dir, locus):
    donor_id = donor['donor_id']
    # logging.info(f'Processing {donor_id}.')
    file_path_case = os.path.join(raw_eh_dir, f"{donor['case_object_id']}.json")
    file_path_control = os.path.join(raw_eh_dir, f"{donor['control_object_id']}.json")

       # Test if the files exist
    case_exists = os.path.isfile(file_path_case)
    control_exists = os.path.isfile(file_path_control)


    try:
        case_counts = extract_locus_data(file_path_case, locus)
        control_counts = extract_locus_data(file_path_control, locus)
    except Exception as e:
        # logging.error(f'Could not decode JSON for {donor_id} Error: {str(e)}')
        return 
    
    return donor_id, case_counts, control_counts, donor['case_object_id'], donor['control_object_id']


def extract_genotypes_counts(manifest_path, raw_eh_dir, output_dir, locus):
    
    # Load the manifest
    manifest = pd.read_csv(manifest_path)
    print(manifest.head())

    # Ensure output directory exists
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    results = manifest.apply(lambda x: process_donor(x, raw_eh_dir, locus), axis=1)

    donors, case_counts, control_counts, case_file, control_file = zip(*results)
    df = pd.DataFrame({'donor_id': donors, 'case_counts': case_counts, 'control_counts': control_counts, 'case_file': case_file, 'control_file': control_file})
    output_path = os.path.join(output_dir, locus + 'genotype_counts.csv')
    df.to_csv(output_path, index=False)
    df = pd.read_csv(output_path)
    # Function to extract counts from the string representation of dictionaries
    def extract_counts(counts_str):
        # Convert the string representation of dictionary into an actual dictionary
        counts_dict = eval(counts_str)
        # Extract and return the desired parts
        return counts_dict['CountsOfFlankingReads'], counts_dict['CountsOfSpanningReads'], counts_dict['CountsOfInrepeatReads'], counts_dict['Genotype']

    # Apply the function to the case_counts and control_counts columns and expand the results into new columns
    df[['case_flank', 'case_span', 'case_irr', 'case_genotype']] = df.apply(lambda x: extract_counts(x['case_counts']), axis=1, result_type='expand')
    df[['control_flank', 'control_span', 'control_irr', 'control_genotype']] = df.apply(lambda x: extract_counts(x['control_counts']), axis=1, result_type='expand')

    # Rename file name columns for clarity
    df.rename(columns={'case_file': 'case_file_name', 'control_file': 'control_file_name'}, inplace=True)

    # Dropping the original counts columns for clarity, assuming they are no longer needed
    df.drop(['case_counts', 'control_counts'], axis=1, inplace=True)

        # logging.info('Finished Saving DataFrames.')

    df.to_csv(output_path, index=False)

    return df


def init_argparse():
    parser = argparse.ArgumentParser(description='Process Expansion Hunter output for analysis of paired genotype differences.')
    parser.add_argument('raw_eh', metavar='RawDir', type=str, help='Directory with Expansion Hunter output JSONs.')
    parser.add_argument('manifest', metavar='Manifest', type=str, help='Manifest file with case and control object ids.')
    parser.add_argument('--locus', '-l', required=True, help='Locus of Interest.')
    parser.add_argument('--outdir', '-o', required=True, help='Output directory (default .).')
    return parser


def main():
    parser = init_argparse()
    args = parser.parse_args()
    extract_genotypes_counts(args.manifest, args.raw_eh, args.outdir, args.locus)



if __name__ == "__main__":
    main()
    
  