#! /bin/bash

usage="Collect model PDB files from a given directory based on given pattern.

Usage:
$0 input_directory output_directory (--pattern comma_separated_filename_patterns)
"

pattern=''
positional_arguments=()

# Parsing command line options.
while [[ $# > 0 ]]; do
    opt="$1"
    arg="$2"
    shift
    case $opt in
        -h|--help)
            echo "$usage"
            exit 0
            ;;
        --pattern)
            pattern=$arg
            shift
            ;;
        *)
            positional_arguments+=("$opt")
            ;;
    esac
done
if [ ${#positional_arguments[@]} -lt 2 ]; then
    echo "$usage"
    exit
fi

input_directory=${positional_arguments[0]}
output_directory=${positional_arguments[1]}

set -o nounset
set -o errexit

if [ "$input_directory" == '' ] || [ "$output_directory" == '' ]; then
    echo "$usage"
    exit
fi

# Input.
echo "Collecting necessary files from $input_directory"
pdb_files=$(find $input_directory -type f -name *.pdb)
if [[ "$pattern" != '' ]]; then
    pattern_file=$(mktemp)
    trap "rm $pattern_file" EXIT
    echo "$pattern" | sed 's|,|\n|g' > $pattern_file
    pdb_files=$(echo "$pdb_files" | grep -f $pattern_file)
fi
num_files=$(echo "$pdb_files" | wc -l)
echo "Found $num_files files that match given pattern."
# Output.
if [ ! -d $output_directory ]; then
    mkdir -p $output_directory
fi
echo "$pdb_files" | while read fname; do
    new_fname=$(
        echo "$fname" \
        | sed "s|^"$input_directory"||" \
        | sed 's|/|__|g'
    )
    echo "Copying $fname to $new_fname"
    cp $fname $output_directory/$new_fname
done

