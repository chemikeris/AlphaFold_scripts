#! /usr/bin/env python3
"""Remove files and data for models that were not top ranked.
"""

import sys
import os
import json
import glob
import argparse
import logging


def main(arguments):
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('results_directory')
    arguments_parser.add_argument(
        '--debug', action='store_true',
        help='use debug mode with more logging'
        )
    arguments_parser.add_argument(
        '--only-pkl-file', action='store_true',
        help='remove only pkl file for lower rank models'
        )
    keep = 1
    
    args = arguments_parser.parse_args()
    if args.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    logging.basicConfig(format='%(levelname)s: %(message)s', level=log_level)

    ranking_json = os.path.join(args.results_directory, 'ranking_debug.json')
    try:
        with open(ranking_json) as f:
            ranking_data = json.load(f)
    except FileNotFoundError as err:
        logging.error('Problem reading ranking_debug.json: %s', err)
        return 1
    order = ranking_data['order']
    for good_ranked in order[:keep]:
        logging.info('Keeping top ranked model files: %s', good_ranked)
    for i, badly_ranked in enumerate(order[keep:]):
        logging.info('Removing files related to %s', badly_ranked)
        wildcard = os.path.join(args.results_directory, '*%s*' % badly_ranked)
        if args.only_pkl_file:
            wildcard += '.pkl'
        remove_files = glob.glob(wildcard)
        for fname in remove_files:
            logging.info('Removing file %s', fname)
            os.remove(fname)
        if not args.only_pkl_file:
            ranked_fname = 'ranked_%s.pdb' % str(keep+i)
            logging.info('Removing file %s', ranked_fname)
            try:
                os.remove(os.path.join(args.results_directory, ranked_fname))
            except FileNotFoundError:
                logging.warning('File %s not found!', ranked_fname)


if __name__ == '__main__':
    sys.exit(main(sys.argv))

