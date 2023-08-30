import argparse
import pandas as pd
import orjson
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

MIN_READS = 6
HIGH_COV = 24
MAX_WIDTH = 4


LOG_LEVEL = os.getenv('LOG_LEVEL') or 'info'
log_dict = {'debug': logging.DEBUG, 'info': logging.INFO, 'warning': logging.WARNING, 
            'error': logging.ERROR, 'critical': logging.CRITICAL}
log_level = log_dict.get(LOG_LEVEL.lower(), logging.INFO)

dis_id = os.getenv('DISEASE')
if not dis_id:
    dis_id = datetime.now().strftime('%Y%m%d_%H%M')


log_dir = os.getenv('LOG_DIR') or 'ExpansionCookerLogs'
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

logging.basicConfig(filename=os.path.join(log_dir, f'{dis_id}_ExpansionCooker.log'), 
                    level=log_level,
                    format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')


def is_wide(ci):
    ci = list(map(int, ci.split('-')))
    return ci[1] - ci[0] > MAX_WIDTH

def decide_genotype_order(case, control):

    if case[1] == control[0]:
        temp = case[1]
        case[1] = case[0]
        case[0] = temp


    return case, control

def getPairs(case_ci, control_ci, case_genotypes, control_genotypes):
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
                         local_case_df, local_control_df, local_diff_df):
    if len(case_genotypes) == 2:
        case_genotypes, control_genotypes = decide_genotype_order(case_genotypes, control_genotypes)

    for i in range(len(case_genotypes)):
        local_case_df.append({'donor_id': donor_id + f'_{i}', 'ReferenceRegion': ReferenceRegion, 'value': case_genotypes[i]})
        local_control_df.append({'donor_id': donor_id + f'_{i}', 'ReferenceRegion': ReferenceRegion, 'value': control_genotypes[i]})
        local_diff_df.append({'donor_id': donor_id + f'_{i}', 'ReferenceRegion': ReferenceRegion, 'value': case_genotypes[i] - control_genotypes[i]})


def process_locus(donor_id, data_case, data_control, local_case_df, local_control_df, local_diff_df, local_df_tracking):
    allele_count = data_case['AlleleCount']
    high_cov = HIGH_COV
    min_reads = MIN_READS
    if allele_count == 1:
        high_cov = HIGH_COV / 2
        min_reads = MIN_READS / 2
    for variant in set(data_case['Variants']):
        case = data_case['Variants'][variant]
        control = data_control['Variants'][variant]
        ReferenceRegion = data_case['Variants'][variant]['ReferenceRegion']

        try: 
            control_genotypes = list(map(int, control.get('Genotype').split('/')))
            case_genotypes = list(map(int, case.get('Genotype').split('/')))
        except AttributeError as a:
            continue
        
        # make genotype checker objects for case and cotrol
        case_genotype_checker = GenotypeChecker(case_genotypes, case.get('CountsOfSpanningReads'), case.get('CountsOfFlankingReads'))
        control_genotype_checker = GenotypeChecker(control_genotypes, control.get('CountsOfSpanningReads'), control.get('CountsOfFlankingReads'))

        case_num = case_genotype_checker.num_reads()
        control_num = control_genotype_checker.num_reads()

        # if very high read count, trust Egor's genotypes
        if case_num > high_cov and control_num > high_cov:
            append_genotype_data(case_genotypes, control_genotypes, donor_id, ReferenceRegion, 
                                 local_case_df, local_control_df, local_diff_df)
            continue

        # if very low read count, skip (REVISIT)
        if case_num < min_reads:
            # log to tracking and continue
            local_df_tracking.append({'donor_id': donor_id, 
                            'ReferenceRegion': ReferenceRegion,
                            'motif': case.get('RepeatUnit'),
                            'issue': 'case_low_reads'})
            local_df_tracking.append({'donor_id': donor_id, 
                            'ReferenceRegion': ReferenceRegion,
                            'motif': case.get('RepeatUnit'),
                            'issue': 'case_low_reads'})
            continue

        if control_num < min_reads:
            local_df_tracking.append({'donor_id': donor_id, 
                            'ReferenceRegion': ReferenceRegion,
                            'motif': case.get('RepeatUnit'),
                            'issue': 'control_low_reads'})
            local_df_tracking.append({'donor_id': donor_id, 
                            'ReferenceRegion': ReferenceRegion,
                            'motif': case.get('RepeatUnit'),
                            'issue': 'case_low_reads'})
            continue

        
        # if there is a big difference in read counts, use support checking method (REVISIT, potential bias towards expansions )
        diff = abs(case_num - control_num)/min(case_num, control_num)
        if diff > 0.40:
            checked_control_genotypes = control_genotype_checker.identify_supported_genotypes()

            if len(checked_control_genotypes) == 0:
                local_df_tracking.append({'donor_id': donor_id, 
                                'ReferenceRegion': ReferenceRegion,
                                'motif': case.get('RepeatUnit'),
                                'issue': 'control_low_reads'})
                local_df_tracking.append({'donor_id': donor_id, 
                                'ReferenceRegion': ReferenceRegion,
                                'motif': case.get('RepeatUnit'),
                                'issue': 'control_low_reads'})
                continue
            control_is_first = checked_control_genotypes[0] == control_genotypes[0]
            case_genotype_checker.add_genotypes(checked_control_genotypes)
            checked_case_genotypes = case_genotype_checker.identify_supported_genotypes()

            if len(checked_case_genotypes) == 0:
                local_df_tracking.append({'donor_id': donor_id, 
                                'ReferenceRegion': ReferenceRegion,
                                'motif': case.get('RepeatUnit'),
                                'issue': 'case_low_reads'})
                local_df_tracking.append({'donor_id': donor_id, 
                                'ReferenceRegion': ReferenceRegion,
                                'motif': case.get('RepeatUnit'),
                                'issue': 'case_low_reads'})
                continue
            final_case_genotypes = []
            
            i = 0
            while i < len(checked_control_genotypes):
                if checked_control_genotypes[i] in checked_case_genotypes:
                    checked_case_genotypes.remove(checked_control_genotypes[i])
                    final_case_genotypes.append(checked_control_genotypes[i])
                i += 1
            i = 0
            checked_case_genotypes = sorted(checked_case_genotypes)
            while len(final_case_genotypes) < len(checked_control_genotypes) and i < len(checked_case_genotypes):
                if control_is_first:
                    final_case_genotypes.append(checked_case_genotypes[i])

                i += 1
                    
            # make sure we have the same number of genotypes
            if len(checked_control_genotypes) < len(final_case_genotypes):
                final_case_genotypes = final_case_genotypes[:len(checked_control_genotypes)]
            elif len(checked_control_genotypes) > len(final_case_genotypes):
                checked_control_genotypes = checked_control_genotypes[:len(final_case_genotypes)]
            
            append_genotype_data(final_case_genotypes, checked_control_genotypes, donor_id, ReferenceRegion, 
                                 local_case_df, local_control_df, local_diff_df)
            continue

        
        # otherwise, use CI approach
        ci_approach(allele_count, donor_id, case, control, local_case_df, local_control_df, local_diff_df, local_df_tracking)
        

def ci_approach(allele_count, donor_id, case, control, local_case_df, local_control_df, local_diff_df, local_df_tracking):
    ReferenceRegion = case['ReferenceRegion']
    # get values
    case_ci = case.get('GenotypeConfidenceInterval')
    control_ci = control.get('GenotypeConfidenceInterval')


    if allele_count == 1:
        if is_wide(case_ci) or is_wide(control_ci):
            local_df_tracking.append({'donor_id': donor_id, 
                                'ReferenceRegion': ReferenceRegion,
                                'motif': case.get('RepeatUnit'),
                                'control_ci': control_ci, 
                                'case_ci': case_ci})
            return

        local_case_df.append({'donor_id': donor_id + '_0', 'ReferenceRegion': ReferenceRegion, 'value': case.get('Genotype')})
        local_control_df.append({'donor_id': donor_id + '_0', 'ReferenceRegion': ReferenceRegion, 'value': control.get('Genotype')})
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
                                 local_case_df, local_control_df, local_diff_df)
            return

        # if 3 out of 4 values are nan, skip
        if tot_wide >= 3 or (is_wide(case_ci[0]) and is_wide(case_ci[1])) or (is_wide(control_ci[0]) and is_wide(control_ci[1])):
            local_df_tracking.append({'donor_id': donor_id, 
                                'ReferenceRegion': ReferenceRegion, 
                                'motif': case.get('RepeatUnit'),
                                'control_ci': control_ci[0], 
                                'case_ci': case_ci[0]})
            local_df_tracking.append({'donor_id': donor_id, 
                                'ReferenceRegion': ReferenceRegion, 
                                'motif': case.get('RepeatUnit'),
                                'control_ci': control_ci[1], 
                                'case_ci': case_ci[1]})
            return
        
        goodPair, badCi = getPairs(case_ci, control_ci, case_genotypes, control_genotypes)
        local_case_df.append({'donor_id': donor_id + '_0', 'ReferenceRegion': ReferenceRegion, 'value': goodPair[0]})
        local_control_df.append({'donor_id': donor_id + '_0', 'ReferenceRegion': ReferenceRegion, 'value': goodPair[1]})
        local_diff_df.append({'donor_id': donor_id + '_0', 'ReferenceRegion': ReferenceRegion, 'value': goodPair[0] - goodPair[1]})
        local_df_tracking.append({'donor_id': donor_id, 
                            'ReferenceRegion': ReferenceRegion,
                            'motif': case.get('RepeatUnit'),
                            'control_ci': badCi[1], 
                            'case_ci': badCi[0]})
def process_donor(donor, raw_eh_dir):
    donor_id = donor['donor_id']
    logging.info(f'Processing {donor_id}.')
    file_path_case = os.path.join(raw_eh_dir, f"{donor['case_object_id']}.json")
    file_path_control = os.path.join(raw_eh_dir, f"{donor['control_object_id']}.json")

    local_case_df = []
    local_control_df = []
    local_diff_df = []
    local_df_tracking = []

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
        return local_case_df, local_control_df, local_diff_df, local_df_tracking, donor_id

    with open(file_path_case, 'r') as file_case, open(file_path_control, 'r') as file_control:
        try:
            data_case = orjson.loads(file_case.read())
            data_control = orjson.loads(file_control.read())
        except Exception as e:
            logging.error(f'Could not decode JSON for {donor_id} Error: {str(e)}')
            return local_case_df, local_control_df, local_diff_df, local_df_tracking, donor_id

        for locus in set(data_case['LocusResults']):
            process_locus(donor_id, data_case['LocusResults'][locus], data_control['LocusResults'][locus], local_case_df, local_control_df, local_diff_df, local_df_tracking)

    logging.info(f'Finished {donor_id} succesfully.')
    return local_case_df, local_control_df, local_diff_df, local_df_tracking, donor_id


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
    
  