threads=1                           # Number of threads for parallel processing
ident_cutoff=0.9                    # Identity cutoff for clustering
cluster_fasta=arhod_reps_clustrd.fasta  # Output file name for clustered FASTA sequences
concat_reps=arhod_reps_for_vsearch.fasta  # Input file name for concatenated FASTA sequences

seqkit tab2fx arhod_reps_for_vsearch.tsv > $concat_reps  # Convert TSV file to FASTA format using seqkit

/software/team301/vsearch-2.21.1/bin/vsearch --cluster_size $concat_reps \  # Execute VSEARCH tool for clustering
    --threads $threads \          # Specify the number of threads for parallel processing
    --id $ident_cutoff \          # Set the identity cutoff for clustering
    --strand both \               # Consider both strands for clustering
    --sizein \                    # Input file contains sequence size information
    --fasta_width 0 \             # Disable line wrapping in the output FASTA file
    --relabel cluster \           # Relabel the sequences with cluster names
    --uc all.preclustered.uc \    # Output file for cluster information in UC format
    --centroids $cluster_fasta    # Output file for clustered centroid sequences in FASTA format


# Execute the RepeatClassifier tool within a Singularity container
# Mount the /lustre directory from the host into the Singularity container
# Use the tetools_1.3.1.sif Singularity image
# Provide the clustered FASTA file ($cluster_fasta) as input to RepeatClassifier
/software/singularity-v3.6.4/bin/singularity run -B /lustre:/lustre /lustre/scratch123/tol/teams/blaxter/projects/tol-nemotodes/sw/.singularity/tetools_1.3.1.sif RepeatClassifier -consensi $cluster_fasta


# Perform a blastn search
# Use the $assem variable as the subject (database) for the search
# Use the $cluster_fasta file as the query for the search
# Specify the output format as tabular (format '6')
# Save the output to arhod_reps_clustrd.tsv file
# Set the e-value threshold to 0.01
blastn -subject $assem -query $cluster_fasta -outfmt '6 qaccver saccver pident length qlen slen qstart qend sstart send evalue bitscore qcovs' -out arhod_reps_clustrd.tsv -evalue 0.01