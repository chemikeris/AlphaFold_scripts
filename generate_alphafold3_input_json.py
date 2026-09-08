#!/usr/bin/env python3
"""
Create AlphaFold3 input JSON (version 3) from separate FASTA files
for proteins, DNA, and RNA, with stoichiometry implemented via
unique chain IDs (A–Z).
"""

import os
import argparse
import json
import logging
import sys
from pathlib import Path
from typing import List, Tuple, Dict, Iterator
import string
import random

# -----------------------------------------------------------------------------
# FASTA parsing
# -----------------------------------------------------------------------------
def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    records = []
    current_id = None
    buf = []

    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(buf)))
                current_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)

    if current_id is not None:
        records.append((current_id, "".join(buf)))

    return records

# -----------------------------------------------------------------------------
# Stoichiometry parsing (numbers only)
# -----------------------------------------------------------------------------
def parse_numeric_stoichiometry(stoich: str, n_items: int) -> List[int]:
    if not stoich:
        return [1] * n_items

    parts = stoich.split(":")
    if len(parts) != n_items:
        raise ValueError(
            f"Stoichiometry length ({len(parts)}) does not match "
            f"number of sequences ({n_items})."
        )

    try:
        counts = [int(p) for p in parts]
    except ValueError:
        raise ValueError("Stoichiometry must contain integers only.")

    if any(c < 1 for c in counts):
        raise ValueError("Stoichiometry values must be >= 1.")

    return counts

# -----------------------------------------------------------------------------
# Chain ID generator (A–Z)
# -----------------------------------------------------------------------------
def chain_id_generator() -> Iterator[str]:
    for c in string.ascii_uppercase:
        yield c

    for first in string.ascii_uppercase:
        for second in string.ascii_uppercase:
            yield first + second

    raise RuntimeError("Exceeded maximum number of chains (702).")

def custom_chain_id_generator(prefix) -> Iterator[str]:
    chain_ids = chain_id_generator()
    for c in chain_ids:
        yield prefix+c

# -----------------------------------------------------------------------------
# Generate random seeds
# -----------------------------------------------------------------------------
def generate_random_seeds(n: int ) -> List[int]:
    return list(set(random.sample(range(4294967295), n)))

# -----------------------------------------------------------------------------
# Entity builders
# -----------------------------------------------------------------------------
def protein_entity(
        chain_ids: List[str],
        seq: str,
        no_templates: bool,
        msa_path: str
        ) -> Dict:
    templates_info = [] if no_templates else None
    return_value = {
        "protein": {
            "id": chain_ids,
            "sequence": seq,
            'templates': templates_info
            }
        }
    if msa_path:
        return_value["protein"]["pairedMsaPath"] = \
            os.path.join(msa_path, "pairedMsa.a3m")
        return_value["protein"]["unpairedMsaPath"] = \
            os.path.join(msa_path, "unpairedMsa.a3m")
    return return_value

def dna_entity(chain_ids: List[str], seq: str) -> Dict:
    return {"dna": {"id": chain_ids, "sequence": seq}}

def rna_entity(chain_ids: List[str], seq: str, msa_path: str) -> Dict:
    return_value = {"rna": {"id": chain_ids, "sequence": seq}}
    if msa_path:
        return_value['rna']['unpairedMsaPath'] = msa_path
    return return_value

def ligand_entity_from_ccd(chain_ids: str, ccd: str) -> Dict:
    return {"ligand": {"id": chain_ids, "ccdCodes": ccd}}

def ligand_entity_from_smiles(chain_ids: List[str], smiles: str) -> Dict:
    return {"ligand": {"id": chain_ids, "smiles": smiles}}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Create AlphaFold3 JSON input (version 3)."
    )

    parser.add_argument("--proteins", type=Path, help="Protein FASTA file")
    parser.add_argument("--protein-stoich", type=str, default="",
                        help="Protein stoichiometry (e.g. 2:1)")
    parser.add_argument('--no-protein-templates', action='store_true',
                        help='Do not use templates for proteins')
    parser.add_argument('--protein-msas-dir', type=str, default="", nargs='*',
                        help="Directories with precomputed paired and "\
                             "unpaired protein MSAs"
                        )

    parser.add_argument("--dna", type=Path, help="DNA FASTA file")
    parser.add_argument("--dna-stoich", type=str, default="",
                        help="DNA stoichiometry (e.g. 1:2)")

    parser.add_argument("--rna", type=Path, help="RNA FASTA file")
    parser.add_argument("--rna-stoich", type=str, default="",
                        help="RNA stoichiometry (e.g. 1:1)")
    parser.add_argument("--rna-msa-path", type=str, nargs="*", default="",
                        help="RNA MSA A3M file(s)")

    parser.add_argument("--ligand-ccd", type=str, default=None,
                        help="Ligand PDB CCD codes (comma separated)")
    parser.add_argument("--ligand-ccd-stoich", type=str, default="",
                        help="Ligand from CCD stoichiometry")
    parser.add_argument("--ligand-ccd-prefix", type=str, default="",
                        help="Prefixes for CCD ligand chain names")

    parser.add_argument("--user-ccd-path", type=str, default="",
                        help="Path to custom CCD file")

    parser.add_argument("--ligand-smiles", type=str, default=None,
                        help="Ligand smiles (comma separated)")
    parser.add_argument("--ligand-smiles-stoich", type=str, default="",
                        help="Ligand from SMILES stoichiometry")

    parser.add_argument('--number-of-seeds', type=int, default=1,
                        help='Number of random seeds (default: 1)')

    parser.add_argument("--name", type=str, default="af3_job",
                        help="Job name")
    parser.add_argument("--debug", action="store_true",
                        help="Enable debug logging")

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # Logging
    # -------------------------------------------------------------------------
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s"
    )
    logger = logging.getLogger("af3-json")

    sequences = []
    chain_ids = chain_id_generator()

    # -------------------------------------------------------------------------
    # Proteins
    # -------------------------------------------------------------------------
    if args.proteins:
        logger.info(f"Reading protein FASTA: {args.proteins}")
        proteins = parse_fasta(args.proteins)
        counts = parse_numeric_stoichiometry(
            args.protein_stoich, len(proteins)
        )
        if args.protein_msas_dir:
            logger.info("Using pre-calculated protein MSAs")
            if len(args.protein_msas_dir) != len(proteins):
                logger.error(
                    "Number of protein MSA directories does not match number "\
                    "of protein sequences!"
                    )
                sys.exit(1)
            protein_msa_paths = args.protein_msas_dir
        else:
            protein_msa_paths = [None] * len(proteins)

        for (_, seq), count, p_msa in zip(proteins, counts, protein_msa_paths):
            ids = [next(chain_ids) for _ in range(count)]
            sequences.append(
                protein_entity(ids, seq, args.no_protein_templates, p_msa)
                )

    # -------------------------------------------------------------------------
    # DNA
    # -------------------------------------------------------------------------
    if args.dna:
        logger.info(f"Reading DNA FASTA: {args.dna}")
        dnas = parse_fasta(args.dna)
        counts = parse_numeric_stoichiometry(
            args.dna_stoich, len(dnas)
        )

        for (_, seq), count in zip(dnas, counts):
            ids = [next(chain_ids) for _ in range(count)]
            sequences.append(dna_entity(ids, seq))

    # -------------------------------------------------------------------------
    # RNA
    # -------------------------------------------------------------------------
    if args.rna:
        logger.info(f"Reading RNA FASTA: {args.rna}")
        rnas = parse_fasta(args.rna)
        if args.rna_msa_path:
            logger.info("Using pre-calculated RNA MSAs")
            if len(args.rna_msa_path) != len(rnas):
                logger.error(
                    "Number of RNA MSA files does not match number of RNA "\
                    "sequences!"
                    )
                sys.exit(1)
            rna_msa_paths = args.rna_msa_path
        else:
            rna_msa_paths = [None] * len(rnas)

        counts = parse_numeric_stoichiometry(
            args.rna_stoich, len(rnas)
        )

        for (_, seq), count, msa_path in zip(rnas, counts, rna_msa_paths):
            ids = [next(chain_ids) for _ in range(count)]
            sequences.append(rna_entity(ids, seq, msa_path))

    # -------------------------------------------------------------------------
    # Ligand CCD
    # -------------------------------------------------------------------------
    if args.ligand_ccd:
        logger.info("Reading ligands from given CCD codes")
        ligand_ccd_codes = args.ligand_ccd.split(',')
        if args.ligand_ccd_prefix:
            ligand_ccd_prefixes = args.ligand_ccd_prefix.split(',')
            if len(ligand_ccd_codes) != len(ligand_ccd_prefixes):
                logging.error(
                    "Number of CCD codes does not match number of given "\
                    "prefixes!"
                )
                sys.exit(1)
            else:
                use_custom_ligand_chain_names = True
        else:
            use_custom_ligand_chain_names = False
        counts = parse_numeric_stoichiometry(
            args.ligand_ccd_stoich, len(ligand_ccd_codes)
        )
        for i, ligand_ccd in enumerate(ligand_ccd_codes):
            count = counts[i]
            if use_custom_ligand_chain_names:
                names = custom_chain_id_generator(ligand_ccd_prefixes[i])
                ids = [next(names) for _ in range(count)]
            else:
                ids = [next(chain_ids) for _ in range(count)]
            sequences.append(ligand_entity_from_ccd(ids, [ligand_ccd]))

    # -------------------------------------------------------------------------
    # Ligand SMILES
    # -------------------------------------------------------------------------
    if args.ligand_smiles:
        logger.info("Reading ligand SMILES")
        smiles = args.ligand_smiles.split(',')
        counts = parse_numeric_stoichiometry(
            args.ligand_smiles_stoich, len(smiles)
        )
        for sm, count in zip(smiles, counts):
            chain = [next(chain_ids) for _ in range(count)]
            sequences.append(ligand_entity_from_smiles(chain, sm))

    if not sequences:
        logger.error("No sequences provided.")
        sys.exit(1)

    af3_json = {
        "name": args.name,
        "sequences": sequences,
        "dialect": "alphafold3",
        "version": 3,
        "modelSeeds": generate_random_seeds(args.number_of_seeds),
    }
    if args.user_ccd_path:
        af3_json['userCCDPath'] = args.user_ccd_path

    logger.info("Writing AlphaFold3 JSON to stdout")
    json.dump(af3_json, sys.stdout, indent=2)
    sys.stdout.write("\n")

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()

