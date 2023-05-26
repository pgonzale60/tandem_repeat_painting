# Genomic Tandem Repeat Identification and Analysis

This repository contains scripts for identifying multimegabase tandem repeats in genome sequences, reducing redundancy, locating the occurrences of these monomers genome-wide, and creating visualizations for better understanding of the genome-wide distribution of these repeats.

## Features

- Multimegabase Tandem Repeat Identification: Identify and extract tandem repeat sequences from genomic data.
- Redundancy Reduction: Use vsearch to minimize redundancy in the identified repeats.
- Genome-wide Monomer Location: Find the occurrences of these monomers throughout the entire genome.
- Visualization: Create plots with chromosomes as panels and color-coded repeat monomers for easy identification of identical repeats in different parts of the genome.

## Dependencies

- R
- vsearch
- blastn
- bioawk

## Usage

1. Clone this repository to your local machine using `https://github.com/<YourUserName>/Genomic-Tandem-Repeat-Identification-and-Analysis.git`
2. Navigate to the directory of the cloned repository. (`cd Genomic-Tandem-Repeat-Identification-and-Analysis`)
3. Run the scripts with your genome sequence file as an argument.
4. Use the generated outputs for further analysis or visualization.

## Scripts

### tandem_repeat_identification.sh

Identify multimegabase tandem repeats in genome sequences.

Usage: `bash tandem_repeat_identification.sh genome_sequence.fasta`

### redundancy_reduction.R

Reduce redundancy in the identified repeats using vsearch.

Usage: `Rscript redundancy_reduction.R repeats.windowed.trf.tsv`

### monomer_location.sh

Find the locations of these monomers throughout the entire genome.

Usage: `bash monomer_location.sh genome_sequence.fasta reps_for_vsearch.tsv`

### visualization.R

Create plots with chromosomes as panels and color-coded repeat monomers.

Usage: `Rscript visualization.R blasted_clustrd_reps.tsv.gz`

## Contributing

Contributions are welcome! Please read the [contribution guidelines](CONTRIBUTING.md) first.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact

For help and feedback, please feel free to contact the repository owner.

This README was last updated on 26th May 2023.
