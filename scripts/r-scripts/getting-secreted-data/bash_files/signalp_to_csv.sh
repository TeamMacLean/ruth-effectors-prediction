#!/bin/bash

# $1 euk|gram-|gram+
# $2 directory

# SRCs:
# https://stackoverflow.com/questions/19075671/how-do-i-use-shell-variables-in-an-awk-script
# https://stackoverflow.com/questions/6116994/padding-with-sprintf
# https://stackoverflow.com/questions/21668471/bash-script-create-array-of-all-files-in-a-directory

# Function
function split_file () {
	echo "Splitting file" $1"..."
	filename="${1%.*}"
	extension="${1##*.}"
	awk -v pattern=$filename -v outputdir=$2 -v ext=$extension '
		BEGIN {
			n_seq=0;
		} /^>/ {
			if (n_seq%500==0) {
				file = sprintf(outputdir"/"pattern"_%06d."ext, n_seq);
			}
			print >> file;
			n_seq++; next;
		} {
			print >> file;
		}
	' < $1

	echo "✓ Split finished."
	echo
}

function predict_signalp () {
	echo "Predicting signalp on file" $2"..."
	filename="${2%.*}"
	extension="${2##*.}"
	signalp -t $1 $2 | grep -B 1 -A 3 "Prediction:" | sed ':a;N;s/\n/,/g;ba' | sed 's/--,>/\n>/g' | sed 's/,,*$//g' > $3/$filename"_results.csv"

	echo "✓ Prediction finished."
	echo
}

# --------------------------------------------------------------------------
# SPLIT FASTA FILES
# --------------------------------------------------------------------------

echo "-------------------------------------------------"
echo "Splitting files"
echo "-------------------------------------------------"

# Set working directory
cd $2

# Output subdirectory
splitdir="split"
rm -rf $splitdir
if [ ! -d $splitdir ]; then
	mkdir $splitdir
	echo "Created" $splitdir "subdirectory."
fi

# List of FASTA files
fa_files=(*".fa")

# Process all fles
for ((i = 0; i < ${#fa_files[@]}; i++)); do
	split_file ${fa_files[$i]} $splitdir
done

# --------------------------------------------------------------------------
# PREDICT SIGNALP
# --------------------------------------------------------------------------

echo "-------------------------------------------------"
echo "Predicting signalp"
echo "-------------------------------------------------"

# Set working directory
cd $splitdir

# Output subdirectory
signalpdir="signalp_results"
rm -rf $signalpdir
if [ ! -d $signalpdir ]; then
	mkdir $signalpdir
	echo "Created" $signalpdir "subdirectory."
fi

# List of FASTA files
fa_files=(*".fa")

# Process all fles
for ((i = 0; i < ${#fa_files[@]}; i++)); do
	predict_signalp $1 ${fa_files[$i]} $signalpdir
done


# --------------------------------------------------------------------------
# CLEANUP
# --------------------------------------------------------------------------

cd ..
cat $splitdir/$signalpdir/*.csv > ../$2_signalp_results.csv
rm -rf $splitdir



