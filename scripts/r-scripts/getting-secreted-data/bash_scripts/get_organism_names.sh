#!/bin/bash

# Function
function get_organism_name() {
	# echo "Reading file" $1"..."
	filename="${1%.*}"
	cat $1 | grep -E "^>" | awk -v orgname=$filename '{ print orgname, $1 }' | sed 's/.pep.all >/,/g'
}


# --------------------------------------------------------------------------
# READ ORGANISM NAMES
# --------------------------------------------------------------------------

# Set working directory
cd $1

# Output subdirectory

# List of FASTA files
fa_files=(*".fa")

# Process all fles
for ((i = 0; i < ${#fa_files[@]}; i++)); do
	get_organism_name ${fa_files[$i]} >> ../$1_organism_names.csv
done
