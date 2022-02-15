#! /bin/bash

set -o nounset
set -o errexit
set -o pipefail

conda_directory="$HOME/miniconda3"
alphafold_installation_directory="$HOME/software/alphafold"
databases_directory="/tmp/alphafold"
input_fasta=""
output_directory=""
db_preset=""
model_preset=""
max_template_date="2200-02-02"
dry_run=false

usage="A script to run AlphaFold2 from DeepMind.

Command line arguments:
-C: Miniconda directory where AlphaFold2 is installed (current: $conda_directory);
-I: AlphaFold installation directory (current: $alphafold_installation_directory);
-D: directory path to AlphaFold databases (current: $databases_directory);
-i: input Fasta file;
-o: output directory (default: directory of input file);
-d: AlphaFold databases preset (reduced_dbs or full_dbs);
-p: AlphaFold model preset (monomer_ptm, multimer);
-T: maximum template date (current: $max_template_date);
-t: dry run: do not run AlphaFold, only test configuration;
-h: show this message.
"

if [ $# -eq 0 ]; then
    echo "$usage"
    exit 0
fi
while getopts 'I:D:i:o:d:p:T:th' opt; do
    case $opt in
        C) conda_directory="$OPTARG";;
        I) alphafold_installation_directory="$OPTARG";;
        D) databases_directory="$OPTARG";;
        i) input_fasta="$OPTARG";;
        o) output_directory="$OPTARG";;
        d) db_preset="$OPTARG";;
        p) model_preset="$OPTARG";;
        T) max_template_date="$OPTARG";;
        t) dry_run=true;;
        h) echo "$usage"; exit 0;;
    esac
done

if [ "$output_directory" == "" ]; then
    output_directory=`dirname $input_fasta`/${model_preset}_${db_preset}
fi
if [ "$db_preset" == "full_dbs" ]; then
    bfd_argument="--bfd_database_path=$databases_directory/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    small_bfd_argument=""
    uniclust30_argument="--uniclust30_database_path=$databases_directory/uniclust30/uniclust30_2018_08/uniclust30_2018_08"
elif [ "$db_preset" == "reduced_dbs" ]; then
    small_bfd_argument="--small_bfd_database_path=$databases_directory/small_bfd/bfd-first_non_consensus_sequences.fasta"
    bfd_argument=""
    uniclust30_argument=""
else
    echo "ERROR: unknown AlphaFold databases preset setting!"
    echo
    echo "$usage"
    exit 1
fi
if [ "$model_preset" == "monomer_ptm" ]; then
    uniprot_argument=""
    pdb_seqres_argument=""
    pdb70_argument="--pdb70_database_path=$databases_directory/pdb70/pdb70"
elif [ "$model_preset" == "multimer" ]; then
   uniprot_argument="--uniprot_database_path=$databases_directory/uniprot/uniprot.fasta"
   pdb_seqres_argument="--pdb_seqres_database_path=$databases_directory/pdb_seqres/pdb_seqres.txt"
   pdb70_argument=""
else
    echo "ERROR: unknown AlphaFold model preset!"
    echo
    echo "$usage"
    exit 1
fi

echo "Running AlphaFold for $input_fasta."
echo "Using databases preset '$db_preset', AlphaFold model preset '$model_preset'."
echo "Output directory: $output_directory."
echo "Using Miniconda from $conda_directory, AlphaFold from $alphafold_installation_directory, and databases from $databases_directory."
if $dry_run; then
    exit 0
fi

. $conda_directory/bin/activate alphafold2
python $alphafold_installation_directory/run_alphafold.py \
    --fasta_paths=$input_fasta \
    --output_dir=$output_directory \
    --data_dir=$databases_directory \
    --db_preset=$db_preset \
    --model_preset=$model_preset \
    $bfd_argument \
    $small_bfd_argument \
    $uniprot_argument \
    $pdb_seqres_argument \
    $pdb70_argument \
    $uniclust30_argument \
    --uniref90_database_path=$databases_directory/uniref90/uniref90.fasta \
    --mgnify_database_path=$databases_directory/mgnify/mgy_clusters_2018_12.fa \
    --max_template_date=$max_template_date \
    --template_mmcif_dir=$databases_directory/pdb_mmcif/mmcif_files/ \
    --obsolete_pdbs_path=$databases_directory/pdb_mmcif/obsolete.dat

