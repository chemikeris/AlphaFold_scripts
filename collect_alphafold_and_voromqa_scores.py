#! /usr/bin/env python3
"""Collect AlphaFold modeling and VoroMQA scores, assuming that they are there.
"""

import sys
import os
import logging
import csv
import argparse


def get_alphafold_scores_for_model(model_no, results_directory, af_model):
    "Get summary of calculated scores"
    logging.info('Getting AlphaFold scores for model %s.', model_no)
    filename = os.path.join(
        results_directory,
        f'result_model_{model_no}_{af_model}_pred_0.af_scores'
        )
    plddt = None
    ptm = None
    iptm = None
    with open(filename) as f:
        for line in f:
            if line.startswith('pLDDT'):
                plddt = line.split()[1]
            if line.startswith('pTM'):
                ptm = line.split()[1]
            if line.startswith('ipTM'):
                iptm = line.split()[1]
    return plddt, ptm, iptm


def read_voromqa_scores_for_model(model_no, results_directory, af_model):
    "Read VoroMQA scores"
    logging.info('Reading VoroMQA scores for model %s', model_no)
    filename = os.path.join(
        results_directory,
        f'relaxed_model_{model_no}_{af_model}_pred_0.voromqa'
        )
    with open(filename) as f:
        voromqa_line_parts = f.read().splitlines()[-1].split()
        voromqa = voromqa_line_parts[1]
        try:
            voromqa_energy = voromqa_line_parts[9]
        except IndexError:
            voromqa_energy = None
    return voromqa, voromqa_energy


def collect_scores(directory, protein, multimer):
    "Collect scores for all AlphaFold models"
    results_directory = os.path.join(directory, protein)
    models = range(1, 6)
    if multimer:
        af_model = 'multimer_v3'
    else:
        af_model = 'ptm'
    logging.info('Getting AlphaFold scores for protein %s', protein)
    all_scores = []
    for m in models:
        scores = {}
        scores['Protein'] = protein
        scores['Model'] = m
        plddt, ptm, iptm = get_alphafold_scores_for_model(
            m, results_directory, af_model
            )
        scores['pLDDT'] = plddt
        scores['pTM'] = ptm
        if iptm:
            scores['ipTM'] = iptm
        voromqa, voromqa_energy = read_voromqa_scores_for_model(
               m, results_directory, af_model
               )
        scores['voromqa'] = voromqa
        if voromqa_energy:
            scores['voromqa_energy'] = voromqa_energy
        all_scores.append(scores)
    return all_scores


def print_scores(scores):
    "Print collected scores"
    w = csv.DictWriter(
        sys.stdout, fieldnames=scores[0].keys(), delimiter='\t',
        lineterminator='\n')
    w.writeheader()
    for s in scores:
        w.writerow(s)


def main(arguments):
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('directory')
    arguments_parser.add_argument(
        'protein',
        help='Protein subdirectory name '\
            '(if it is "all", all subdirectories are analyzed)'
        )
    arguments_parser.add_argument(
        '--multimer', action='store_true', help='parse multimer results'
        )
    args = arguments_parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.protein == 'all':
        logging.info('Collecting scores for all proteins in given directory.')
        proteins = os.listdir(args.directory)
    else:
        proteins = [args.protein]
    scores = []
    for p in proteins:
        scores += collect_scores(args.directory, p, args.multimer)
    print_scores(scores)


if __name__ == '__main__':
    sys.exit(main(sys.argv))

