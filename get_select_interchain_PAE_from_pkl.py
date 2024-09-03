#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A script to select custom inter-chain PAE from AlphaFold 2 pkl file.

Originally written by Lukas Valančauskas.
"""

import argparse
import json
import string
import csv
import sys
import logging
import os

from Bio import SeqIO
import numpy as np

from analyze_alphafold_pickle import ModelData, NumPyEncoder


def fasta_reader(file):
    "Read fasta file and return header, sequence tuples"
    ret = []
    with open(file, "r", encoding="UTF-8") as f_handle:
        for sr_tuple in SeqIO.FastaIO.SimpleFastaParser(f_handle):
            ret.append(sr_tuple)
    return ret


def get_chain_descriptions(fasta_data, chains):
    "Get chain descriptions that correspond to given letters"
    alphabet = string.ascii_uppercase + string.ascii_lowercase
    descriptions = []
    for chain_letter in chains:
        i = alphabet.index(chain_letter)
        descriptions.append(fasta_data[i][0])
    return descriptions
    

def get_chain_v_chain_PAE_res_pkl(
        pkl_model_data, seqs, query_chains, target_chains, gather_all=False):
    "Get chain vs chain PAE from AlphaFold pkl."
    def get_ich_pae(start0, stop0, start1, stop1, pae):
        "chain v chain interchain pae"
        ich_byrows = pae[start0:stop0, start1:stop1]
        ich_byrows = ich_byrows.reshape(1, ich_byrows.size)
        ich_bycols = pae[start1:stop1, start0:stop0]
        ich_bycols = ich_bycols.reshape(1, ich_bycols.size)
        ich_pae = np.concatenate((ich_bycols, ich_byrows), axis=1)
        return ich_pae

    def get_ch_pae(start, stop, pae):
        ch_pae = pae[start:stop, start:stop]
        ch_pae = ch_pae.reshape(1, ch_pae.size)
        return ch_pae

    def process_pair_pae(name0, name1, seqs, pae):
        idx0 = [n for n, i in enumerate(seqs) if i[0] == name0]
        if not idx0:
            msg = f"Couldnt find target seq {name0} among seqs: "
            msg += ",".join([i[0] for i in seqs])
            logging.error(msg)
            raise RuntimeError
        idx0 = idx0[0]
        start0 = [len(i[1]) for n, i in enumerate(seqs) if n < idx0]
        start0 = sum(start0)
        stop0 = start0 + len(seqs[idx0][1])
        idx1 = [n for n, i in enumerate(seqs) if i[0] == name1]
        if not idx1:
            msg = f"Couldnt find target seq {name1} among seqs: "
            msg += ",".join([i[0] for i in seqs])
            print(msg)
            raise RuntimeError
        idx1 = idx1[0]
        start1 = [len(i[1]) for n, i in enumerate(seqs) if n < idx1]
        start1 = sum(start1)
        stop1 = start1 + len(seqs[idx1][1])
        ch_pae = get_ch_pae(start0, stop0, pae)
        ich_pae = get_ich_pae(start0, stop0, start1, stop1, pae)
        return ch_pae, ich_pae

    pae = pkl_model_data.predicted_aligned_error
    pairs = []
    for query in query_chains:
        for target in target_chains:
            pairs.append((query, target))
    pair_res = [process_pair_pae(i[0], i[1], seqs, pae) for i in pairs]
    ch_pae = np.concatenate([i[0] for i in pair_res], axis=1)
    ich_pae = np.concatenate([i[1] for i in pair_res], axis=1)
    # chain pae
    if ch_pae.size == 0:
        logging.error(f"Couldn't select PAE result for file {pkl_f}!")
        return {}
    if ich_pae.size == 0:
        logging.error(f"Failed to select inter-chain PAE for file {pkl_f}")
        return {}
    res = {
        "pkl_file": os.path.basename(pkl_model_data.file_name),
        "ch_PAE_median": np.median(ch_pae),
        "ch_PAE_mean": np.mean(ch_pae),
        "ich_PAE_median": np.median(ich_pae),
        "ich_PAE_mean": np.mean(ich_pae),
    }
    if gather_all:
        res["ptm"] = pkl_model_data.ptm
        res["iptm"] = pkl_model_data.iptm
        res["plddt_median"] = np.median(pkl_model_data.plddt)
    res["PAE_query_chains"] = query_chains
    res["PAE_target_chains"] = target_chains
    return res


def print_csv(results_dict, skip_header=False):
    "Output results to CSV"
    writer = csv.DictWriter(sys.stdout, fieldnames=results_dict.keys())
    for item_to_join in ('PAE_query_chains', 'PAE_target_chains'):
        results_dict[item_to_join] = ':'.join(results_dict[item_to_join])
    if not skip_header:
        writer.writeheader()
    writer.writerow(results_dict)


def print_json(results):
    "Output results to JSON"
    print(json.dumps(results, cls=NumPyEncoder))


def main():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('pkl_file')
    arguments_parser.add_argument(
        '--fasta', required=True, help='fasta file which was used for modeling'
        )
    arguments_parser.add_argument(
        '--chains1', required=True,
        help='comma separated chain names of first subunit'
        )
    arguments_parser.add_argument(
        '--chains2', required=True,
        help='comma separated chain names of second subunit'
        )
    arguments_parser.add_argument(
        '--output-csv', help='output CSV instead of JSON',
        action='store_true'
        )
    arguments_parser.add_argument(
        '--no-header', action='store_true', help='skip CSV header row'
        )
    args = arguments_parser.parse_args()

    # Reading input data.
    model_data = ModelData(args.pkl_file, multimer=True)
    sequences = fasta_reader(args.fasta)
    query_chains = get_chain_descriptions(sequences, args.chains1.split(','))
    target_chains = get_chain_descriptions(sequences, args.chains2.split(','))
    gather_all = True  # also extract ptm, iptm and plddt values, rather than only PAE
    # Processing.
    res = get_chain_v_chain_PAE_res_pkl(
        pkl_model_data=model_data,
        seqs=sequences,
        query_chains=query_chains,
        target_chains=target_chains,
        gather_all=gather_all,
    )
    # Output.
    if args.output_csv:
        print_csv(res, args.no_header)
    else:
        print_json(res)


if __name__ == "__main__":
    main()
