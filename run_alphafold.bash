#! /bin/bash

set -o nounset
set -o errexit
set -o pipefail

conda_directory="/data/miniconda3"
alphafold_installation_directory="/data/alphafold"
databases_directory="/data/alphafold_dbs"
input_fasta=""
output_directory=""
db_preset=""
model_preset=""
max_template_date="2200-02-02"
gpu_present=true
precomputed_msas=false
analyze_results=false
calculate_voromqa=false
dry_run=false

usage="A script to run AlphaFold2 from DeepMind.

Command line arguments:
-C: Miniconda directory where AlphaFold2 is installed (current: $conda_directory);
-I: AlphaFold installation directory (current: $alphafold_installation_directory);
-D: directory path to AlphaFold databases (current: $databases_directory);
-S: analyze AlphaFold prediction results using scripts from given scripts directory;
-i: input Fasta file;
-o: output directory (default: directory of input file);
-d: AlphaFold databases preset (reduced_dbs, full_dbs, full_dbs_newer_uniprot);
-p: AlphaFold model preset (monomer_ptm, multimer);
-T: maximum template date (current: $max_template_date);
-m: use precomputed MSAs;
-G: indicate if GPU is not possible;
-V: calculate VoroMQA scores for relaxed models;
-t: dry run: do not run AlphaFold, only test configuration;
-h: show this message.

To run using SLURM on GPU server, use this sbatch command:
sbatch --ntasks=1 --gres=gpu:1 --cpus-per-task=24 ...
"

if [ $# -eq 0 ]; then
    echo "$usage"
    exit 0
fi
while getopts 'C:I:D:S:i:o:d:p:T:thGmV' opt; do
    case $opt in
        C) conda_directory="$OPTARG";;
        I) alphafold_installation_directory="$OPTARG";;
        D) databases_directory="$OPTARG";;
        S) analyze_results=true; scripts_directory="$OPTARG";;
        i) input_fasta="$OPTARG";;
        o) output_directory="$OPTARG";;
        d) db_preset="$OPTARG";;
        p) model_preset="$OPTARG";;
        T) max_template_date="$OPTARG";;
        m) precomputed_msas=true;;
        G) gpu_present=false;;
        t) dry_run=true;;
        V) calculate_voromqa=true;;
        h) echo "$usage"; exit 0;;
    esac
done

if [ "$output_directory" == "" ]; then
    output_directory=`dirname $input_fasta`/${model_preset}_${db_preset}
fi
if [ "$db_preset" == "full_dbs" ]; then
    bfd_argument="--bfd_database_path=$databases_directory/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    small_bfd_argument=""
    uniref30_argument="--uniref30_database_path=$databases_directory/uniref30/UniRef30_2021_03"
elif [ "$db_preset" == "full_dbs_newer_uniprot" ]; then
    bfd_argument="--bfd_database_path=$databases_directory/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    small_bfd_argument=""
    uniref30_argument="--uniref30_database_path=$databases_directory/uniref30/UniRef30_2023_02"
    db_preset='full_dbs'
elif [ "$db_preset" == "reduced_dbs" ]; then
    small_bfd_argument="--small_bfd_database_path=$databases_directory/small_bfd/bfd-first_non_consensus_sequences.fasta"
    bfd_argument=""
    uniref30_argument=""
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
    num_multimer_models_argument=""
    pkl_analysis_multimer_setting=""
    voromqa_inter_chain_argument=""
elif [ "$model_preset" == "multimer" ]; then
    uniprot_argument="--uniprot_database_path=$databases_directory/uniprot/uniprot.fasta"
    pdb_seqres_argument="--pdb_seqres_database_path=$databases_directory/pdb_seqres/pdb_seqres.txt"
    pdb70_argument=""
    num_multimer_models_argument="--num_multimer_predictions_per_model 1"
    pkl_analysis_multimer_setting="--multimer"
    voromqa_inter_chain_argument="--score-inter-chain"
else
    echo "ERROR: unknown AlphaFold model preset!"
    echo
    echo "$usage"
    exit 1
fi

echo "Running AlphaFold for $input_fasta."
echo "Using databases preset '$db_preset', AlphaFold model preset '$model_preset'."
if $analyze_results; then
    echo "Running AlphaFold prediction analysis scripts from $scripts_directory."
fi
if $calculate_voromqa; then
    if which voronota-voromqa > /dev/null; then
        echo "Calculating VoroMQA for predictions."
    else
        echo "VoroMQA executable 'voronota-voromqa' not found in $PATH, skipping VoroMQA calculations."
        calculate_voromqa=false
    fi
fi
echo "Output directory: $output_directory."
input_basename=`basename $input_fasta`
output_subdirectory=$output_directory/${input_basename%.*}
echo "Results directory: $output_subdirectory."
echo "Using Miniconda from $conda_directory, AlphaFold from $alphafold_installation_directory, and databases from $databases_directory."
if $gpu_present; then
    echo "Relaxing using GPU."
    gpu_relax_argument="--use_gpu_relax"
    if [ ! -z ${CUDA_VISIBLE_DEVICES:-} ]; then
        export NVIDIA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES
    fi
    export TF_FORCE_UNIFIED_MEMORY=1
    export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
else
    echo "Relaxing using CPU."
    gpu_relax_argument="--nouse_gpu_relax"
fi
if $precomputed_msas; then
    echo "Using precomputed MSAs."
    precomputed_msas_argument="--use_precomputed_msas"
else
    echo "Running sequence searches for MSAs."
    precomputed_msas_argument=""
fi
if $dry_run; then
    exit 0
fi

echo "Running AlphaFold predictions."
. $conda_directory/bin/activate alphafold
python $alphafold_installation_directory/run_alphafold.py \
    $gpu_relax_argument $precomputed_msas_argument \
    --fasta_paths=$input_fasta \
    --output_dir=$output_directory \
    --data_dir=$databases_directory \
    --db_preset=$db_preset \
    --model_preset=$model_preset \
    $num_multimer_models_argument \
    $bfd_argument \
    $small_bfd_argument \
    $uniprot_argument \
    $pdb_seqres_argument \
    $pdb70_argument \
    $uniref30_argument \
    --uniref90_database_path=$databases_directory/uniref90/uniref90.fasta \
    --mgnify_database_path=$databases_directory/mgnify/mgy_clusters_2022_05.fa \
    --max_template_date=$max_template_date \
    --template_mmcif_dir=$databases_directory/pdb_mmcif/mmcif_files/ \
    --obsolete_pdbs_path=$databases_directory/pdb_mmcif/obsolete.dat
if $analyze_results; then
    echo "Analyzing prediction results."
    find $output_subdirectory -type f -name result*.pkl | while read f; do
        out=${f/.pkl/.af_scores}
        python $scripts_directory/analyze_alphafold_pickle.py $f --save-plots $pkl_analysis_multimer_setting | tee $out
    done
fi
if $calculate_voromqa; then
    echo "Calculating VoroMQA scores for relaxed models."
    find $output_subdirectory -type f -name relaxed*.pdb | while read f; do
        out=${f/.pdb/.voromqa}
        voronota-voromqa -i $f $voromqa_inter_chain_argument --print-header | tee $out
    done
fi
