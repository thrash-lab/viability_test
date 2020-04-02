import numpy as np
import pandas as pd
from utility import get_ci
import argparse
import multiprocessing
from multiprocessing import Pool
from timeit import default_timer as timer

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def setup_args():
    # Create argument parser
    parser = argparse.ArgumentParser('viability_test')

    parser.add_argument("-b", "--n_bootstraps", help="The number of bootstraps to run per experiment", type=int, default=9999)
    parser.add_argument("-p", "--threads", help="The number of threads to use for bootstrapping (will use all by default)", type=int, default=-1)


    subparsers = parser.add_subparsers(title='Please select a function:',
                                       dest='subcommand',
                                       description='valid subcommands')

    parser_estimate_viability = subparsers.add_parser('estimate_viability', help='Estimate viability of a culture from a DTE experiment')

    # Positional mandatory arguments
    parser_estimate_viability.add_argument("-w", "--wells", help="The number of wells to simulate", type=int, required=True)
    parser_estimate_viability.add_argument("-i", "--inoculum", help="The number of cells added to each well", type=float, required=True)
    parser_estimate_viability.add_argument("-o", "--num_observed", help="The number of observed positive wells", type=int, required=True)
    parser_estimate_viability.add_argument("-r", "--rel_abund", help="The relative abundance of a particular taxon in your inoculum", type=float, default=1)

    parser_predict_wells = subparsers.add_parser('predict_wells', help='Predict the number of wells you are likely to observe for a given taxon in a DTE experiment')
    parser_predict_wells.add_argument("-w", "--wells", help="The number of wells to simulate", type=int, required=True)
    parser_predict_wells.add_argument("-i", "--inoculum",
                                       help="The number of cells added to each well", type=float,
                                       required=True)
    parser_predict_wells.add_argument('-r', '--rel_abund',
                                      help='The relative abundance of a particular taxon in your inoculum', type=float, default=1)
    parser_predict_wells.add_argument('-v', '--viability',
                                      help='The estimated viability of a particular taxon in your inoculum', type=float,
                                      default=1)


    parser_bulk_estimate = subparsers.add_parser('bulk', help='Perform a bulk estimate of viabilities from a TSV input')
    parser_bulk_estimate.add_argument('-i', '--input', help='A TSV file containing rows with 4 columns defined as: wells, inoculum, rel_abund, viability')
    parser_bulk_estimate.add_argument('-o', '--output', help='The output file')
    return parser.parse_args()


def calculate_dte(inoculum_size, proportion_of_taxa, viability=1, num_wells=96):
    """
    Main method to calculate DTE probabilities
    :param inoculum_size: The number of cells added to each well
    :param proportion_of_taxa: The relative abundance of the target taxa in the inoculum
    :param viability: The viability of the target taxa
    :param num_wells: The number of wells being inoculated
    :return: the number of positive wells,
            the number of pure wells,
            the number of positive wells containing the target taxon,
            the number of pure wells containing the target taxon
    """
    inocula = np.random.poisson(inoculum_size, num_wells)
    taxon_wells = np.random.binomial(inocula, proportion_of_taxa)
    positive_wells = np.sum(inocula >= 1)
    single_wells = np.sum(inocula == 1)
    # This captures all the wells with greater than 0 cells where the number
    # of cells in it are all from the selected taxon
    positive_wells_with_taxon = np.sum(taxon_wells >= 1)
    taxon_pure_wells = np.sum((taxon_wells == inocula) & (inocula > 0))
    if viability < 1:
        pure_wells_to_test = taxon_wells[(taxon_wells == inocula) & (inocula > 0)]
        viable_pure_wells = np.random.binomial(pure_wells_to_test, viability)
        taxon_pure_wells = np.sum(viable_pure_wells >= 1)
        unviable_taxon_pure = np.sum(viable_pure_wells == 0)
        positive_wells = positive_wells - unviable_taxon_pure
    return positive_wells, single_wells, positive_wells_with_taxon, taxon_pure_wells


def bootstrap_dte(inoculum, rel_abund, viability, num_wells):
    positive_wells = np.zeros(number_of_experiments)
    single_wells = np.zeros(number_of_experiments)
    taxon_positive_wells = np.zeros(number_of_experiments)
    taxon_pure_wells = np.zeros(number_of_experiments)

    for i in range(number_of_experiments):
        positive_wells[i], single_wells[i], taxon_positive_wells[i], taxon_pure_wells[i] = calculate_dte(inoculum,
                                                                                                         rel_abund,
                                                                                                         viability=viability,
                                                                                                         num_wells=num_wells)

    taxon_pure_med, taxon_pure_low, taxon_pure_high = get_ci(taxon_pure_wells, as_string=False)
    taxon_positive_med, taxon_positive_low, taxon_positive_high = get_ci(taxon_positive_wells, as_string=False)
    positive_med, positive_low, positive_high = get_ci(positive_wells, as_string=False)
    single_med, single_low, single_high = get_ci(single_wells, as_string=False)

    return pd.Series([inoculum, rel_abund, viability,
                      num_wells, number_of_experiments,
                      positive_med, positive_low, positive_high,
                      single_med, single_low, single_high,
                      taxon_positive_med, taxon_positive_low, taxon_positive_high,
                      taxon_pure_med, taxon_pure_low, taxon_pure_high],
                     index=['cells_per_well', 'rel_abund', 'viability',
                            'num_wells', 'number_of_experiments',
                            'positive_well_med', 'positive_well_95pc_low', 'positive_well_95pc_high',
                            'single_med', 'single_low', 'single_high',
                            'taxon_positive_med', 'taxon_positive_95pc_low', 'taxon_positive_95pc_high',
                            'pure_well_med', 'pure_well_95pc_low', 'pure_well_95pc_high'])


def process_viability(items):
    """
    Utility method for multiprocessing
    :param items: The item to be processed
    :return: the value of viability if it can explain the observed wells, or 0 if not.
    """
    s = bootstrap_dte(items[0], items[1], items[2], items[3])
    taxon_pure_low = s.pure_well_95pc_low
    taxon_pure_high = s.pure_well_95pc_high
    if items[4] >= taxon_pure_low and items[4] <= taxon_pure_high:
        return items[2]
    else:
        return 0

def predict_wells(inoculum, rel_abund, viability, num_wells):
    """
    Function for predicting outcome of a DTE experiment
    :param inoculum: The size of the inoculum
    :param rel_abund: The relative abundance of your taxon of interest
    :param viability: The estimated viability of your taxon of interest
    :param num_wells: The number of simulated wells.
    :return: Nothing.
    """
    print(f'Predicting outcome of a DTE experiment using {num_threads} processors and {number_of_experiments} bootstraps')
    result = bootstrap_dte(inoculum, rel_abund, viability, num_wells)
    print(f'''{bcolors.OKGREEN}
    You simulated a DTE where you inoculated {num_wells} wells with an average of {inoculum} cells per well.
    Your estimated relative abundance for your taxon of interest was:       {rel_abund:1.3%}
    Your estimated viability of your taxon of interest was:                 {viability:1.3%}

    After {number_of_experiments} simulations:
    The median number of positive (not taxon-specific) wells was:           {result.positive_well_med:.0f} ({result.positive_well_95pc_low:.0f}-{result.positive_well_95pc_high:.0f}, 95%CI)
    The median number of wells inoculated with one cell was:                {result.single_med:.0f} ({result.single_low:.0f}-{result.single_high:.0f}, 95%CI)
    The median number of wells containing your taxon of interest was:       {result.taxon_positive_med:.0f} ({result.taxon_positive_95pc_low:.0f}-{result.taxon_positive_95pc_high:.0f}, 95%CI)
    The median number of pure wells for your taxon of interest was:         {result.pure_well_med:.0f} ({result.pure_well_95pc_low:.0f}-{result.pure_well_95pc_high:.0f}, 95%CI)
    {bcolors.ENDC}''')



def estimate_viability(inoculum, num_wells, observed_pure, rel_abund=1, test_min=0, test_max=1, initial_test_step=0.01):
    """
    Estimates the viability of a culture
    :param inoculum: The mean number of cells added to a well
    :param num_wells: The number of wells inoculated in the experiment
    :param observed_pure: The number of observed positive wells following DTE
    :param rel_abund: The relative abundance of the taxon in the culture
    :param test_min: The minimum viability to test
    :param test_max: The maximum viability to test
    :param initial_test_step: The size of the steps taken to test viability
    :return: min,max viability
    """

    viability_range = np.arange(start=test_min, stop=test_max, step=initial_test_step)
    item_list = [(inoculum, rel_abund, v, num_wells, observed_pure) for v in viability_range]
    with Pool(num_threads) as p:
        viability_values = list(p.imap(process_viability, item_list))

    cleaned_viability_values = [x for x in viability_values if x != 0]
    try:
        min = np.min(cleaned_viability_values)
        max = np.max(cleaned_viability_values)

        new_min = np.max([0, min - initial_test_step])
        new_max = np.min([1, max + initial_test_step])
        new_step = initial_test_step / 10
        #print(
        #    f'''Minimum and maximum values are {min} and {max}, respectively - refining:\nTesting values between {new_min:.3f} and {new_max:.3f} in increments of {new_step}''')

        viability_range = np.arange(start=new_min, stop=new_max, step=new_step)
        item_list = [(inoculum, rel_abund, v, num_wells, observed_pure) for v in viability_range]

        with Pool(num_threads) as p:
            viability_values = list(p.imap(process_viability, item_list))
        cleaned_viability_values = [x for x in viability_values if x != 0]

        min = np.min(cleaned_viability_values)
        max = np.max(cleaned_viability_values)
    except ValueError:
        # To get to here, the list of viability values is full of zeros, which means
        # the range hasn't been found. This is most likely due to extremely low viability
        # So, we set new_min to 0 and new_max to 0.1 and try at steps of 0.0001
        new_step = initial_test_step
        while len(cleaned_viability_values) == 0 and new_step > 0.00001:
            new_min = 0
            new_max = new_step
            new_step = new_step / 10
            #print(f'Trying values between {new_min} and {new_max} with increments of {new_step}')
            viability_range = np.arange(start=new_min, stop=new_max, step=new_step)
            item_list = [(inoculum, rel_abund, v, num_wells, observed_pure) for v in viability_range]

            with Pool(num_threads) as p:
                viability_values = list(p.imap(process_viability, item_list))
            cleaned_viability_values = [x for x in viability_values if x != 0]
        if len(cleaned_viability_values) > 0:
            min = np.min(cleaned_viability_values)
            max = np.max(cleaned_viability_values)
        else:
            #print(f'Gave up - estimated viability is really rather low.')
            min = np.NaN
            max = np.NaN

    return min, max

def estimate_viability_wrapper(inoculum, wells, num_observed, rel_abund, bulk=False):
    start = timer()

    #print(f'Estimating viability using {num_threads} processors and {number_of_experiments} bootstraps')
    result = bootstrap_dte(inoculum, rel_abund, 1, wells)
    if result.pure_well_95pc_high *1.1 < num_observed:
        print(f'''
            {bcolors.FAIL}Your number of observed wells ({num_observed}) is at least 10% greater than the 95%CI maximum if viability were 100% ({result.pure_well_95pc_high:.0f}).
            Consequently, it makes little sense to try and estimate viability >100%, so I will set min,max values to Inf.
            I advise caution in interpreting such extreme values.{bcolors.ENDC}''')
        min=np.Inf
        max=np.Inf
    else:
        min, max = estimate_viability(inoculum, wells, num_observed, rel_abund)

    end = timer()
    clock_time = end - start

    print(f'''
            You simulated a DTE where you inoculated {wells} wells with an average of {inoculum} cells per well.
            A range of decreasing viability was tested using {number_of_experiments} bootstraps per experiment.
            {bcolors.OKGREEN}Your estimated relative abundance for your taxon of interest was:       {rel_abund:1.3%}
            You observed:                                                           {num_observed} wells of interest

            If viability of your taxon had been 100%, you would have expected:      {result.pure_well_med:.0f} ({result.pure_well_95pc_low:.0f}-{result.pure_well_95pc_high:.0f}, 95%CI)

            The viability estimates that explain your observed counts is between:   {min:1.3%} - {max:1.3%}
            {bcolors.ENDC}(the whole process took {clock_time:1.3f} seconds)
            ''')

    return pd.Series(data=[result.pure_well_med,  result.pure_well_95pc_low, result.pure_well_95pc_high, min, max],
                    index=['pure_well_med', 'pure_well_95pc_low', 'pure_well_95pc_high', 'viability_95pc_low', 'viabiilty_95pc_high'])



def bulk_estimate(args):
    df = pd.read_csv(args.input, sep='\t')
    results = []
    for index, row in df.iterrows():
        results.append(estimate_viability_wrapper(row['inoculum'],
                                                    row['wells'],
                                                    row['num_observed'],
                                                    row['rel_abund'], bulk=True))
    results_df = pd.DataFrame(results)
    out_df = df.join(results_df)
    out_df['within_range']  = (out_df.num_observed >= out_df.pure_well_95pc_low) & (out_df.num_observed<=out_df.pure_well_95pc_high)
    out_df['deviance'] = None
    out_df.loc[~out_df.within_range, 'deviance'] = out_df.num_observed - out_df.pure_well_med
    out_df.to_csv(args.output, sep='\t', index=False)

def main():

    args = setup_args()
    global number_of_experiments
    number_of_experiments = args.n_bootstraps

    global number_of_wells_to_simulate

    global num_threads
    num_threads = args.threads

    if num_threads == -1:
        num_threads = multiprocessing.cpu_count()

    #This tests the viability
    if args.subcommand == 'estimate_viability':
        number_of_wells_to_simulate = args.wells
        estimate_viability_wrapper(args.inoculum, args.wells, args.num_observed, args.rel_abund)
    elif args.subcommand =='predict_wells':
        predict_wells(args.inoculum, args.rel_abund, args.viability, args.wells)
    elif args.subcommand =='bulk':
        bulk_estimate(args)


if __name__=="__main__":
    main()
