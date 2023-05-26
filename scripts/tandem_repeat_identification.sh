#!/bin/bash

# Input parameters
assembly_file=$1

# Extract the base name of the assembly file without extension for further usage
assembly_base_name=$(basename "$assembly_file" .fa)
search_id="windowed100K_${assembly_base_name}"

# If tandem repeats are prevalent and multimegabase, it's more convenient to segment the chromosomes into shorter sequences
# 'makewindows' function creates windows along the genome, 'getfasta' retrieves the sequence for these windows
# The output is passed to 'seqkit' for further processing, the result is saved in a file with 'windowed100K_' prefix
bedtools makewindows -g ${assembly_file}.fai -w 100000 | bedtools getfasta -fi $assembly_file -bed - | seqkit seq > "windowed100K_${assembly_base_name}.fa"

# Path to singularity binary. Singularity is used here to run the TRF (Tandem Repeats Finder) program within a virtual machine
singularity_cmd="/software/singularity-v3.6.4/bin/singularity run -B /lustre:/lustre /lustre/scratch123/tol/teams/blaxter/projects/tol-nemotodes/sw/.singularity/tetools_1.3.1.sif trf"


# Execute TRF command with parameters
$singularity_cmd ${search_id}.fa 2 7 7 80 10 100 600 -f -h -d -m

# Rename output files for clarity
mv ${search_id}.fa.2.7.7.80.10.100.600.dat ${search_id}.trf.dat
mv ${search_id}.fa.2.7.7.80.10.100.600.mask ${search_id}.trf.mask

# Extract and format the output, only selecting rows that give the coordinates of the repeats
# Add sequence coordinates to the resulting file. The result is saved in a tab-separated values file (.tsv)
grep -P '^[S0-9]' ${search_id}.trf.dat | cut -d ' ' -f 1,2,3,4,6,7,8,13,14 | awk '{if($1 ~ /Sequence/){chr=$2} else {print chr, $0}}' | tr [:blank:] '\t' > ${search_id}.trf.tsv

# View the largest repeats, those that have a score (in this case, the result of 2 for match -7 for mismatch and -7 for indel in the coordinate) greater than 80000 
# bioawk -t '$8 > 80000' ${search_id}.trf.tsv | less