import argparse
import pandas as pd
import orjson
import json5
import os
import numpy as np
import logging
from datetime import datetime
import logging.handlers
from ExpansionFeatureExtractor import process_features
import re
from collections import Counter
import random

from multiprocessing import Pool, cpu_count
from functools import partial

slurm_cpus = os.getenv('SLURM_CPUS_PER_TASK')
if slurm_cpus:
    cpu_count = int(slurm_cpus)
else:
    cpu_count = cpu_count()


# Load configuration values from the default location if the file exists
config_values = {}
if os.path.exists('./config.json5'):
    with open('./config.json5', 'r') as file:
        config_values = json5.load(file)

LOG_DIR = config_values.get("LOG_DIR", 'ExpansionCookerLogs')
MIN_READS = config_values.get("MIN_READS", 6)
HIGH_COV = config_values.get("HIGH_COV", 24)
MAX_WIDTH = config_values.get("MAX_WIDTH", 2)

# Logging parameters, ugly setup since weird to pass values from script call with the multiprocessing
#       Needs enviormental variables to be set for proper logging (Might be better soln)
DISEASE = os.getenv('DISEASE') or datetime.now().strftime('%Y%m%d_%H%M')
LOG_LEVEL = os.getenv('LOG_LEVEL') or 'info'

log_dict = {'debug': logging.DEBUG, 'info': logging.INFO, 'warning': logging.WARNING, 
            'error': logging.ERROR, 'critical': logging.CRITICAL}
log_level = log_dict.get(LOG_LEVEL.lower(), logging.INFO)

if not os.path.exists(LOG_DIR):
    os.makedirs(LOG_DIR)

logging.basicConfig(filename=os.path.join(LOG_DIR, f'{DISEASE}_ExpansionCooker.log'), 
                    level=log_level,
                    format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')


def is_wide(ci):
    # Determines if a confidence interval is wide based on the MAX_WIDTH value.
    ci = list(map(int, ci.split('-')))
    return ci[1] - ci[0] > MAX_WIDTH

def decide_genotype_order(case, control):
    # Swaps the order of the genotypes to pair matching genotypes.
    if case[1] == control[0]:
        temp = case[1]
        case[1] = case[0]
        case[0] = temp


    return case, control

def getPairs(case_ci, control_ci, case_genotypes, control_genotypes):
    # Pairs together non-wide intervals, assuming there one or two wide intervals present.
    if is_wide(case_ci[0]):
        if is_wide(control_ci[1]):
            return [case_genotypes[1], control_genotypes[0]], [case_ci[0], control_ci[1]]
        else:
            return [case_genotypes[1], control_genotypes[1]], [case_ci[0], control_ci[0]]
    else:
        if is_wide(control_ci[0]):
            return [case_genotypes[0], control_genotypes[1]], [case_ci[1], control_ci[0]]
        else:
            return [case_genotypes[0], control_genotypes[0]], [case_ci[1], control_ci[1]]


class GenotypeChecker:
    """ 
    Class to check if a genotype is supported by the reads.

    :param genotypes: List of genotypes to check.
    :param spanning_reads: String containing the counts of spanning reads.
    :param flanking_reads: String containing the counts of flanking reads.
    """
    def __init__(self, genotypes, spanning_reads, flanking_reads):
        self.genotypes = genotypes
        self.spanning_reads_dict = self._parse_reads_to_dict(spanning_reads)
        self.flanking_reads_list = sorted(self.parse_counts(flanking_reads))
        self.supported_genotypes = []
        self.tot_spanning = sum(self.spanning_reads_dict.values())
        self.tot_flanking = len(self.flanking_reads_list)

    def get_lowest_genotype(self):
        return min(self.supported_genotypes)
    
    def _parse_reads_to_dict(self, reads):
        reads = self.parse_counts(reads)
        return Counter(reads)

    def add_genotypes(self, genotypes):
        # append to front of list
        self.genotypes = genotypes + self.genotypes
        
    def num_reads(self):
        return self.tot_spanning + (self.tot_flanking / 4)

    def _count_flanking_below_threshold(self, threshold):
        count = 0
        for read in self.flanking_reads_list:
            if read <= threshold:
                count += 1
            else:
                break  # since the list is sorted, we can break once we pass the threshold
        return count
    
    def parse_counts(self, s):
        """
        Parse the counts of reads.

        :param s: String containing the counts of reads.
        """
        
        # Adjusted the regex to ensure proper matching
        pairs = re.findall(r'\((\d+),\s*(\d+)\)', s)
        
        # Convert strings to integers and expand based on counts
        result = []
        for length, count in pairs:
            result.extend([int(length)] * int(count))
            
        return result

    def _get_spanning_reads_near_genotype(self, genotype):
        """
        Helper method to get spanning reads near the given genotype.
        
        :param genotype: The genotype to check.
        :return: List of spanning reads near the genotype.
        """

        return self.spanning_reads_dict[genotype - 1] + self.spanning_reads_dict[genotype] + self.spanning_reads_dict[genotype + 1]

    def _get_flanking_reads_below_genotype(self, genotype):
        """
        Helper method to get flanking reads below the given genotype.
        Returns proportion of reads under given genotype based on next greatest genotype.
        
        :param genotype: The genotype to check.
        :return: List of flanking reads below the genotype.
        """

        return [read for read in self.flanking_reads_list if read <= genotype]

    def _check_support(self, genotype):
        """
        Private method to check if a genotype is supported by the reads.
        
        :param genotype: The genotype to check.
        :return: Boolean indicating if the genotype is supported.
        """
        SPANNING_UPPER_THRESHOLD = max(4, self.tot_spanning / 3)
        SPANNING_LOWER_THRESHOLD = max(2, self.tot_spanning / 6)
        FLANKING_THRESHOLD = max(6, self.tot_flanking / 4)
        
        close_spanning_reads = self._get_spanning_reads_near_genotype(genotype)
        
        # if genotype shows up > 4 times in spanning reads
        if close_spanning_reads > SPANNING_UPPER_THRESHOLD:
            self._remove_supporting_reads(genotype)
            return True
        
        # elif genotype shows up 2<4 times in spanning reads
        elif SPANNING_LOWER_THRESHOLD < close_spanning_reads <= SPANNING_UPPER_THRESHOLD:
            below_genotype_flanking_reads = self._get_flanking_reads_below_genotype(genotype)
            
            # if more than THRESHOLD reads below genotype
            if len(below_genotype_flanking_reads) > FLANKING_THRESHOLD:
                self._remove_supporting_spanning_reads(genotype)
                self._remove_supporting_flanking_reads(genotype)
                return True
        
        return False
    
    def _remove_supporting_reads(self, genotype):
        """
        Private method to remove reads that support the given genotype.
        """
        self._remove_supporting_spanning_reads(genotype)
        self._remove_supporting_flanking_reads(genotype)
    
    def _remove_supporting_flanking_reads(self, genotype):
        """
        Private method to remove reads that support the given genotype from flanking reads.
        """
        # make new array with reads above genotype
        above = [read for read in self.flanking_reads_list if read > genotype]
        below = [read for read in self.flanking_reads_list if read <= genotype]

        # Calcualte proportion of reads below genotype that belong to gneotype being removed
        b = max(above) if len(above) > 0 else genotype
        b += 0.01
        r = genotype / b
        prop = 1 - (r / (r+1))

        # remove prop that belonged to genotype
        below = random.sample(below, int(len(below) * prop))
        
        self.flanking_reads_list = above + below
        

    def _remove_supporting_spanning_reads(self, genotype):
        """
        Private method to remove reads that support the given genotype from spanning reads.
        """

        self.spanning_reads_dict[genotype - 1] /= 2
        self.spanning_reads_dict[genotype] /= 2
        self.spanning_reads_dict[genotype + 1] /= 2

    def check_genotype(self, genotype):
        """
        Checks if the given genotype is supported by the reads. Removes the supporting reads if so
        
        :param genotype: The genotype to check.
        :return: Boolean indicating if the genotype is supported.
        """
        return self._check_support(genotype)


    def identify_supported_genotypes(self):
        """
        Identify and return the genotypes that are supported by the reads.
        
        :return: List of supported genotypes.
        """
        for genotype in self.genotypes:
            if self._check_support(genotype):
                self.supported_genotypes.append(genotype)
                self._remove_supporting_reads(genotype)
        return self.supported_genotypes or []
    
def append_genotype_data(case_genotypes, control_genotypes, donor_id, ReferenceRegion,
                        data_frames):
    if len(case_genotypes) == 2:
        case_genotypes, control_genotypes = decide_genotype_order(case_genotypes, control_genotypes)

    for i in range(len(case_genotypes)):
        data_frames['case_df'].append({'donor_id': donor_id + f'_{i}', 'ReferenceRegion': ReferenceRegion, 'value': case_genotypes[i]})
        data_frames['control_df'].append({'donor_id': donor_id + f'_{i}', 'ReferenceRegion': ReferenceRegion, 'value': control_genotypes[i]})
        data_frames['diff_df'].append({'donor_id': donor_id + f'_{i}', 'ReferenceRegion': ReferenceRegion, 'value': case_genotypes[i] - control_genotypes[i]})

def track_issue(tracking_df, donor_id, ReferenceRegion, motif, issue, count=2):
    for i in range(count):
        tracking_df.append({'donor_id': donor_id, 
                            'ReferenceRegion': ReferenceRegion,
                            'motif': motif,
                            'issue': issue})

            
def process_locus(donor_id, data_case, data_control, data_frames):
    """ 
    Function to process a single locus from the Expansion Hunter output. For each variant in the locus, it will check if the reported genotype is supported by the reads. If not, it will try to find a supported genotype. If it cannot find a supported genotype, it will log the locus as problematic. 

    Logic for checking if a genotype is supported by the reads:
    1. If the total read count is very high, trust the reported genotype.
    2. If the total read count for case OR control is very low, skip.
    3. If there is a big difference in read counts, use the support checking method.
        3.1 Check what genotypes in control are supported in control reads
        3.2 Check what genotypes in case + control are supported by case reads
        3.3 Pair up the supported genotypes
    4. Otherwise, use the confidence interval approach.

    Keyword arguments:
    donor_id -- The donor id.
    data_case -- The data for the case.
    data_control -- The data for the control.
    data_frames -- The data frames to append the genotypes to.


    """

    allele_count = data_case['AlleleCount']
    high_cov = HIGH_COV
    min_reads = MIN_READS

    # Check if allele count is 1, if so, half the thresholds to maintain consistency in required reads
    if allele_count == 1:
        high_cov = HIGH_COV / 2
        min_reads = MIN_READS / 2

    # For each variant in the locus look at the genotypes, number of reads, and confidence intervals to try and detemrine wether the reported genotype is properly supported by the reads
    for variant in set(data_case['Variants']):
        case = data_case['Variants'][variant]
        control = data_control['Variants'][variant]
        motif = case.get('RepeatUnit')
        ReferenceRegion = data_case['Variants'][variant]['ReferenceRegion']

        try: 
            control_genotypes = list(map(int, control.get('Genotype').split('/')))
            case_genotypes = list(map(int, case.get('Genotype').split('/')))
        except AttributeError as a:
            continue
        
        # make genotype checker objects for case and cotrol
        case_genotype_checker = GenotypeChecker(case_genotypes, case['CountsOfSpanningReads'], case['CountsOfFlankingReads'])
        control_genotype_checker = GenotypeChecker(control_genotypes, control['CountsOfSpanningReads'], control['CountsOfFlankingReads'])

        case_num_reads = case_genotype_checker.num_reads()
        control_num_reads = control_genotype_checker.num_reads()

        # if very high read count, trust Expansion Hunter's Original genotypes
        if case_num_reads > high_cov and control_num_reads > high_cov:
            append_genotype_data(case_genotypes, control_genotypes, donor_id, ReferenceRegion, 
                                 data_frames)
            continue

        # if very low read count, skip (REVISIT)
        if case_num_reads < min_reads:
            track_issue(data_frames['tracking_df'], donor_id, ReferenceRegion, motif, 'case_low_reads', count=allele_count)
            continue

        if control_num_reads < min_reads:
            track_issue(data_frames['tracking_df'], donor_id, ReferenceRegion, motif, 'control_low_reads', count=allele_count)
            continue

        
        # if there is a big percentage difference in read counts, use support checking method
        diff = abs(case_num_reads - control_num_reads)/min(case_num_reads, control_num_reads)
        if diff > 0.30:
            # (REVISIT) Potential bias towards expansions, ideally this function doesnt start with contorl and then check case, 
            #   should find way to make it start based on amount of support or something else IDRK

            # Find the supported genotypes for control, method in GenotypeChecker class
            checked_control_genotypes = control_genotype_checker.identify_supported_genotypes()

            # If no supported control genotypes, log to tracking and move to next variant
            if len(checked_control_genotypes) == 0:
                track_issue(data_frames['tracking_df'], donor_id, ReferenceRegion, motif, 'control_low_reads', count=allele_count)
                continue

            # Add the supported genotypes for control to the case genotype checker, 
            #      so that the control gentoypes are checked for in control.
            #      This is done to ensure that is there is a zero difference it is detected, minimizing false positives
            case_genotype_checker.add_genotypes(checked_control_genotypes)
            checked_case_genotypes = case_genotype_checker.identify_supported_genotypes()


            if len(checked_case_genotypes) == 0:
                track_issue(data_frames['tracking_df'], donor_id, ReferenceRegion, motif, 'case_low_reads', count=allele_count)
                continue

            # Supported control assumed to be final here.
            final_case_genotypes = []

            # Look for matching supported genotypes, add to final_case_genotypes if found
            for genotype in checked_control_genotypes:
                if genotype in checked_case_genotypes:
                    checked_case_genotypes.remove(genotype)
                    final_case_genotypes.append(genotype)
                
            # If not matching, add supported case genotype until we have the same number of genotypes or 
            #   exhaust the list of supported case genotypes
            i = 0
            checked_case_genotypes = sorted(checked_case_genotypes)
            while len(final_case_genotypes) < len(checked_control_genotypes) and i < len(checked_case_genotypes):
                final_case_genotypes.append(checked_case_genotypes[i])
                i += 1
                    
            
            # make sure we have the same number of genotypes
            #   (REVISE) currently just takes the first genotype, should be more robust
            if len(checked_control_genotypes) < len(final_case_genotypes):
                if checked_control_genotypes[0] in final_case_genotypes:
                    checked_control_genotypes = [checked_control_genotypes[0]]
                else:
                    checked_control_genotypes = [control_genotypes[0]]
            elif len(checked_control_genotypes) > len(final_case_genotypes):
                if final_case_genotypes[0] in checked_control_genotypes:
                    checked_control_genotypes = [final_case_genotypes[0]]
                else:
                    checked_control_genotypes = [control_genotypes[0]]

            if len(final_case_genotypes) < allele_count:
                track_issue(data_frames['tracking_df'], donor_id, ReferenceRegion, motif, 'low_reads', count=1)

            append_genotype_data(final_case_genotypes, checked_control_genotypes, donor_id, ReferenceRegion, 
                                data_frames)
            continue

        
        # otherwise, use CI approach
        ci_approach(allele_count, donor_id, case, control, data_frames)
        

def ci_approach(allele_count, donor_id, case, control, data_frames):
    """
    Perform confidence interval approach for allele count analysis. Checks for wide confidence intervals and tries to pair up genotypes with intervals within specified range.

    Args:
        allele_count (int): The allele count.
        donor_id (str): The donor ID.
        case (dict): The case data.
        control (dict): The control data.
        data_frames (dict): The data frames to append the genotypes to.

    Returns:
        None
    """
    ReferenceRegion = case['ReferenceRegion']
    # get values
    case_ci = case.get('GenotypeConfidenceInterval')
    control_ci = control.get('GenotypeConfidenceInterval')
    motif = case.get('RepeatUnit')

    if allele_count == 1:
        if is_wide(case_ci) or is_wide(control_ci):
            track_issue(data_frames['tracking_df'], donor_id, ReferenceRegion, motif, 'wide_interval', count=1)
            return

        append_genotype_data([int(case.get('Genotype'))], [int(control.get('Genotype'))], donor_id, ReferenceRegion, 
                            data_frames)
        return
    if allele_count == 2:

        case_ci = case_ci.split('/')
        control_ci = control_ci.split('/')
        case_genotypes = case.get('Genotype')
        control_genotypes = control.get('Genotype')
        case_genotypes = list(map(int, case_genotypes.split('/')))
        control_genotypes = list(map(int, control_genotypes.split('/')))


        # make any genotypes with a wide confidence interval nan
        if is_wide(case_ci[0]):
            case_genotypes[0] = np.nan
        if is_wide(case_ci[1]):
            case_genotypes[1] = np.nan
        if is_wide(control_ci[0]):
            control_genotypes[0] = np.nan
        if is_wide(control_ci[1]):
            control_genotypes[1] = np.nan

        tot_wide = np.isnan(case_genotypes).sum() + np.isnan(control_genotypes).sum()

        if tot_wide == 0:
            append_genotype_data(case_genotypes, control_genotypes, donor_id, ReferenceRegion,
                                 data_frames)
            return

        # if 3 out of 4 values are nan, skip
        if tot_wide >= 3 or (is_wide(case_ci[0]) and is_wide(case_ci[1])) or (is_wide(control_ci[0]) and is_wide(control_ci[1])):
            track_issue(data_frames['tracking_df'], donor_id, ReferenceRegion, motif, 'wide_interval', count=2)
            return
        
        goodPair, badCi = getPairs(case_ci, control_ci, case_genotypes, control_genotypes)
        append_genotype_data([goodPair[0]], [goodPair[1]], donor_id, ReferenceRegion, 
                            data_frames)
        track_issue(data_frames['tracking_df'], donor_id, ReferenceRegion, motif, 'wide_interval', count=1)
        

def process_donor(donor, raw_eh_dir):
    """ Process all loci for a single donor.

    :param donor: donor information from the manifest, including the case and control object ids and the donor id
    :param raw_eh_dir: directory with Expansion Hunter output JSONs

    :return: DataFrames with the genotypes for the case and control, the difference between the genotypes, and a DataFrame with tracking information and the donor id.
    """
    donor_id = donor['donor_id']
    logging.info(f'Processing {donor_id}.')
    file_path_case = os.path.join(raw_eh_dir, f"{donor['case_object_id']}.json")
    file_path_control = os.path.join(raw_eh_dir, f"{donor['control_object_id']}.json")

    data_frames = {
        'case_df': [],
        'control_df': [],
        'diff_df': [],
        'tracking_df': []
    }

    # Test if the files exist
    case_exists = os.path.isfile(file_path_case)
    control_exists = os.path.isfile(file_path_control)

    if not case_exists and not control_exists:
        logging.error(f'Missing both files for{donor_id}, both files do not exist: {file_path_case}, {file_path_control}')
    elif not case_exists:
        logging.error(f'Missing case for {donor_id}: {file_path_case}')
    elif not control_exists:
        logging.error(f'Missing control for {donor_id}: {file_path_control}')
        
    if not case_exists or not control_exists:
        return data_frames['case_df'], data_frames['control_df'], data_frames['diff_df'], data_frames['tracking_df'], donor_id

    # Load the JSONs
    with open(file_path_case, 'r') as file_case, open(file_path_control, 'r') as file_control:
        try:
            data_case = orjson.loads(file_case.read())
            data_control = orjson.loads(file_control.read())
        except Exception as e:
            logging.error(f'Could not decode JSON for {donor_id} Error: {str(e)}')
            return data_frames['case_df'], data_frames['control_df'], data_frames['diff_df'], data_frames['tracking_df'], donor_id

        # Process each locus getting the genotypes for the case and control, adding to the DataFrames
        for locus in set(data_case['LocusResults']):
            process_locus(donor_id, data_case['LocusResults'][locus], data_control['LocusResults'][locus], data_frames)

    logging.info(f'Finished {donor_id} succesfully.')
    return data_frames['case_df'], data_frames['control_df'], data_frames['diff_df'], data_frames['tracking_df'], donor_id

def extract_genotypes_diffs(manifest_path, disease_name, raw_eh_dir, output_dir):
    
    # Load the manifest
    manifest = pd.read_csv(manifest_path)


    # Ensure output directory exists
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    with Pool(processes=cpu_count) as pool:
        func = partial(process_donor, raw_eh_dir=raw_eh_dir)
        results = pool.map(func, manifest.to_dict('records'))


    logging.info('Finished processing files, combining results.')
    # Aggregate results
    case_df_list, control_df_list, diff_df_list, df_tracking_list, donor_ids = zip(*results)
    case_df = [item for sublist in case_df_list for item in sublist]
    control_df = [item for sublist in control_df_list for item in sublist]
    diff_df = [item for sublist in diff_df_list for item in sublist]
    df_tracking = [item for sublist in df_tracking_list for item in sublist]

    # Convert lists to DataFrames and pivot
    case_df = pd.DataFrame(case_df).pivot(index='donor_id', columns='ReferenceRegion', values='value')
    control_df = pd.DataFrame(control_df).pivot(index='donor_id', columns='ReferenceRegion', values='value')
    diff_df = pd.DataFrame(diff_df).pivot(index='donor_id', columns='ReferenceRegion', values='value')
    df_tracking = pd.DataFrame(df_tracking)
    
    # log proportion of problematic loci
    logging.info(f'Proportion of problematic loci: {len(df_tracking)/(diff_df.shape[1] * diff_df.shape[0])}')

    logging.info(f'Saving DataFrames.')

    # Save the DataFrames
    case_df.to_csv(os.path.join(output_dir, f'{disease_name}_case.csv'))
    control_df.to_csv(os.path.join(output_dir, f'{disease_name}_control.csv'))
    diff_df.to_csv(os.path.join(output_dir, f'{disease_name}_diff.csv'))
    df_tracking.to_csv(os.path.join(output_dir, f'{disease_name}_tracking.csv'))

    logging.info('Finished Saving DataFrames.')

    return diff_df

def init_argparse():
    parser = argparse.ArgumentParser(description='Process Expansion Hunter output for analysis of paired genotype differences.')
    parser.add_argument('raw_eh', metavar='RawDir', type=str, help='Directory with Expansion Hunter output JSONs.')
    parser.add_argument('manifest', metavar='Manifest', type=str, help='Manifest file with case and control object ids.')
    parser.add_argument('--name', '-n', required=True, help='Disease name for output files.')
    parser.add_argument('--outdir', '-o', required=True, help='Output directory (default .).')
    parser.add_argument('--feats', '-f', default=False, action='store_true', help='Create features from the output? (Default: False)')
    return parser


def main():
    parser = init_argparse()
    args = parser.parse_args()
    diffs = extract_genotypes_diffs(args.manifest, args.name, args.raw_eh, args.outdir)

    if args.feats:
        logging.info('Creating features from the output.')
        process_features(diffs, args.name, args.outdir)  # Pass the diffs dataframe to the refactored function

    logging.info('Finished.')


if __name__ == "__main__":
    main()
    
  