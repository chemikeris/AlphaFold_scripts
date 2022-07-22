#! /bin/bash

use_pdb70=true
alphafold_databases_directory="/data/alphafold_dbs"
output_directory=""
templates=""
conda_directory="/data/miniconda3"

usage="Prepare a library of custom templates for structure modeling.

Usage:
-i, --alphafold-dbs-directory: directory with AlphaFold databases (default: $alphafold_databases_directory)
-o, --output-directory: directory where custom templates DB should be written
    --multimer: use full PDB database (for AlphaFold-multimer, default is PDB70)
-t, --templates: comma separated list of templates
-C, --conda-path: specify custom MiniConda directory path (default: $conda_directory)
-h, --help: show this message
"

# Parsing command line options.
if [ $# -eq 0 ]; then
    echo "$usage"
    exit 0
fi
while [[ $# > 0 ]]; do
    opt="$1"
    arg="$2"
    shift
    case $opt in
        -h|--help)
            echo "$usage"
            exit 0
            ;;
        --multimer)
            use_pdb70=false
            ;;
        -i|--alphafold-dbs-directory)
            alphafold_databases_directory="$arg"
            shift
            ;;
        -o|--output-directory)
            output_directory="$arg"
            shift
            ;;
        -t|--templates)
            templates="$arg"
            shift
            ;;
        -C|--conda-path)
            conda_directory="$arg"
            shift
            ;;
    esac
done

set -o nounset
set -o errexit
set -o pipefail

if  [[ "$output_directory" == "" ]] || [[ "$templates" == "" ]]; then
    echo "$usage"
    echo "Please specify the output directory and list of templates!"
    exit 1
fi

echo "Preparing custom templates from $alphafold_databases_directory."
echo "Copying templates information for $templates."
templates=${templates/,/" "}

if [ ! -d $output_directory ]; then
    mkdir -p $output_directory
fi
cd $output_directory

if $use_pdb70; then
    echo "Copying profiles from PDB70 database, which is default for monomeric AlphaFold."
    source_directory="$alphafold_databases_directory/pdb70"
    mkdir -p pdb70
    . $conda_directory/bin/activate alphafold
    tmp_dir=`mktemp -d`
    mkdir -p $tmp_dir/a3m $tmp_dir/cs219 $tmp_dir/hhm
    for t in $templates; do
        echo "Getting data for $t"
        for suffix in a3m cs219 hhm; do
            ffindex_get $source_directory/pdb70_$suffix.ffdata $source_directory/pdb70_$suffix.ffindex $t > $tmp_dir/$suffix/$t
        done
    done
    echo "Joining profiles to a ffindex database."
    for suffix in a3m cs219 hhm; do
        ffindex_build -as pdb70/pdb70_$suffix.ffdata pdb70/pdb70_$suffix.ffindex $tmp_dir/$suffix
    done
    rm -rf $tmp_dir
else
    echo "Collecting sequences from the PDB database, which is default for multimeric AlphaFold."
    mkdir -p pdb_seqres
    pdb_seqres_file="pdb_seqres/pdb_seqres.txt"
    rm -f $pdb_seqres_file
    for t in $templates; do
        echo "Getting data for $t"
        grep -A1 $t $alphafold_databases_directory/$pdb_seqres_file >> $pdb_seqres_file
    done

fi
echo "Collecting structure files for templates."
mmcifs_directory="pdb_mmcif/mmcif_files"
mmcifs_source="$alphafold_databases_directory/$mmcifs_directory"
mkdir -p $mmcifs_directory
for t in $templates; do
    pdb_id=`echo "$t" | awk -F'_' '{print $1}' | tr '[:upper:]' '[:lower:]'`
    rsync -i --ignore-existing "$mmcifs_source/$pdb_id.cif" $mmcifs_directory
done

