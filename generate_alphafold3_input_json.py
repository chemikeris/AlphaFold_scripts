#!/usr/bin/env python3
"""
Create AlphaFold3 input JSON (version 3) from separate FASTA files
for proteins, DNA, and RNA, with stoichiometry implemented via
unique chain IDs (A–Z).
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import List, Tuple, Dict, Iterator
import string

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
    raise RuntimeError("Exceeded maximum number of chains (26).")

# -----------------------------------------------------------------------------
# Entity builders
# -----------------------------------------------------------------------------
def protein_entity(chain_ids: List[str], seq: str) -> Dict:
    return {"protein": {"id": chain_ids, "sequence": seq}}

def dna_entity(chain_ids: List[str], seq: str) -> Dict:
    return {"dna": {"id": chain_ids, "sequence": seq}}

def rna_entity(chain_ids: List[str], seq: str) -> Dict:
    return {"rna": {"id": chain_ids, "sequence": seq}}

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

    parser.add_argument("--dna", type=Path, help="DNA FASTA file")
    parser.add_argument("--dna-stoich", type=str, default="",
                        help="DNA stoichiometry (e.g. 1:2)")

    parser.add_argument("--rna", type=Path, help="RNA FASTA file")
    parser.add_argument("--rna-stoich", type=str, default="",
                        help="RNA stoichiometry (e.g. 1:1)")

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
        prots = parse_fasta(args.proteins)
        counts = parse_numeric_stoichiometry(
            args.protein_stoich, len(prots)
        )

        for (_, seq), count in zip(prots, counts):
            ids = [next(chain_ids) for _ in range(count)]
            sequences.append(protein_entity(ids, seq))

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
        counts = parse_numeric_stoichiometry(
            args.rna_stoich, len(rnas)
        )

        for (_, seq), count in zip(rnas, counts):
            ids = [next(chain_ids) for _ in range(count)]
            sequences.append(rna_entity(ids, seq))

    if not sequences:
        logger.error("No sequences provided.")
        sys.exit(1)

    af3_json = {
        "name": args.name,
        "sequences": sequences,
        "dialect": "alphafold3",
        "version": 3,
        "modelSeeds": [0],
    }

    logger.info("Writing AlphaFold3 JSON to stdout")
    json.dump(af3_json, sys.stdout, indent=2)
    sys.stdout.write("\n")

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()

