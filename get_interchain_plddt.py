#! /usr/bin/env python3
"""Get pLDDT values of residues forming inter-chain interfaces.

Depends on Voronota-JS package.
"""

import sys
import os
import argparse
import subprocess
import logging
import errno
import csv
import string

from io import StringIO

from analyze_alphafold_pickle import ModelData
from get_select_interchain_PAE_from_pkl import fasta_reader

VORONOTA_CONTACTS_SCRIPT = 'voronota-js-fast-iface-contacts'


def voronota_js_script_present():
    "Check if necessary script voronota-js-fast-iface-contacts is available"
    try:
        subprocess.run(VORONOTA_CONTACTS_SCRIPT,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        return True
    except FileNotFoundError as err:
        if err.errno == errno.ENOENT:
            logging.error('Voronota script %s not found.\n',
                          VORONOTA_CONTACTS_SCRIPT)
        else:
            raise
        return False


def run_voronota_inter_chain_contacts_script(pdb_file):
    "Calculate inter-chain contacts using Voronota-JS fast contacts script"
    voronota_command = [
        VORONOTA_CONTACTS_SCRIPT,
        '--input', pdb_file,
        '--expand-ids',
        '--coarse-grained',
        ]
    try:
        s = subprocess.run(
            voronota_command, capture_output=True, encoding='utf-8'
            )
    except subprocess.CalledProcessError as err:
        logging.error('Could not calculate inter-chain contacts: %s', err)
        return None
    if s.stderr:
        logging.error(
            'Voronota-JS contacts script proceeded with errors: %s', s.stderr
            )
        return None
    return s.stdout


def get_interface_residues_from_contacts(voronota_js_contacts_str):
    "Get interface residue numbers from Voronota-JS contacts"
    header, contacts_data = voronota_js_contacts_str.split('\n', 1)
    fieldnames = header.split()
    contacts = csv.DictReader(
        StringIO(contacts_data), delimiter='\t', fieldnames=fieldnames
        )
    interface_residues = {}
    interface_residue_areas = {}
    for contact in contacts:
        chain_1 = contact['ID1_chainID']
        residue_1 = contact['ID1_resSeq']
        chain_2 = contact['ID2_chainID']
        residue_2 = contact['ID2_resSeq']
        area = float(contact['area'])
        for c, r in ((chain_1, residue_1), (chain_2, residue_2)):
            residue_no = int(r)
            residue_descriptor = (c, residue_no,)
            try:
                interface_residues[c].add(residue_no)
            except KeyError:
                interface_residues[c] = set([residue_no])
            try:
                interface_residue_areas[residue_descriptor] += area
            except KeyError:
                interface_residue_areas[residue_descriptor] = area
    return interface_residues, interface_residue_areas


def all_residues_in_format_of_interface_residues(sequences_data):
    "Fake interface residues from all sequences for testing"
    fake_interface_residues = {}
    fake_interface_areas = {}
    alphabet = string.ascii_uppercase + string.ascii_lowercase
    for description, sequence in sequences_data:
        letter, alphabet = alphabet[0], alphabet[1:]
        fake_interface_residues[letter] = []
        for i in range(len(sequence)):
            fake_interface_residues[letter].append(i+1)
            fake_interface_areas[(letter, i+1,)] = 1
    return fake_interface_residues, fake_interface_areas


def get_plddt_values_for_residues(
        model_data, sequences_data, interface_residues, interface_areas):
    "Get values in pkl's pLDDT array for the interface residues"
    plddt_values = []
    weighted_plddt_values = []
    alphabet = string.ascii_uppercase + string.ascii_lowercase
    current_sequence_start_position = 0
    for description, sequence in sequences_data:
        letter, alphabet = alphabet[0], alphabet[1:]
        try:
            interesting_positions = interface_residues[letter]
        except KeyError:
            logging.debug('Chain "%s" is not in the interface.', letter)
            current_sequence_start_position += len(sequence)
            continue
        for residue_no in interesting_positions:
            residue_descriptor = (letter, residue_no,)
            plddt = model_data.plddt[
                current_sequence_start_position + residue_no - 1
                ]
            weighted_plddt = plddt * interface_areas[residue_descriptor]
            plddt_values.append(plddt)
            weighted_plddt_values.append(weighted_plddt)
            logging.debug(
                'Residue %s in chain %s: pLDDT %s, interface area %s',
                residue_no, letter, plddt, interface_areas[residue_descriptor]
                )
        current_sequence_start_position += len(sequence)
    return plddt_values, weighted_plddt_values


def main(arguments):
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('pkl_file')
    arguments_parser.add_argument('pdb_file')
    arguments_parser.add_argument(
        '--fasta', required=True, help='fasta file which was used for modeling'
        )
    arguments_parser.add_argument(
        '--test', action='store_true', help='use test mode'
        )
    arguments_parser.add_argument(
        '--debug', action='store_true', help='use debug logging'
        )
    args = arguments_parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    if not voronota_js_script_present():
        print(__doc__)
        return 1
    # Reading input data.
    model_data = ModelData(args.pkl_file, multimer=True)
    contacts = run_voronota_inter_chain_contacts_script(args.pdb_file)
    sequences_data = fasta_reader(args.fasta)
    # Processing
    if args.test:
        interface_residues, interface_areas = \
            all_residues_in_format_of_interface_residues(sequences_data)
    else:
        interface_residues, interface_areas = \
            get_interface_residues_from_contacts(contacts)
    plddt_values, weighted_plddt_values = get_plddt_values_for_residues(
        model_data, sequences_data, interface_residues, interface_areas
        )
    # Output:
    #   model file
    #   interface pLDDT
    #   interface pLDDT weighted by areas
    #   global pLDDT
    print(args.pdb_file,
          sum(plddt_values)/len(plddt_values),
          sum(weighted_plddt_values)/sum(interface_areas.values()),
          model_data.global_plddt
          )


if __name__ == '__main__':
    sys.exit(main(sys.argv))

