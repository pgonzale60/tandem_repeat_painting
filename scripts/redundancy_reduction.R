library(GenomicRanges)   # Load the GenomicRanges library for genomic range operations
library(tidyverse)       # Load the tidyverse library for data manipulation

# Find faidx index files ".primary.fa.gz.fai" to get sequence sizes
fai_files <- list.files("analyses/genome_features/sequence_sizes/",
                        full.names = TRUE, pattern = ".primary.fa.gz.fai")

# Set names for the files based on the extracted names from the file paths, correspond to assembly ID
names(fai_files) <- make.names(sub(".+//(.+).primary.fa.gz.fai", "\\1", fai_files))

# Read the sequence sizes from the files and combine them into a single data frame
seq_sizes <- map_df(fai_files, read_tsv,
                    col_names = c("Sequence", "size", "L1", "L2", "L3"),
                    col_types = c("ciiii"),
                    .id = "assembly") %>%
  select(assembly, Sequence, size) %>%
  filter(!grepl("MT", Sequence)) %>%
  mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

# Create a GRanges object from the seq_sizes data frame, specifying the columns to use for chromosome, start, and end positions
gnm <- mutate(seq_sizes, chrom = multispecies_sequence, start = 1, end = size) %>%
  select(chrom, start, end) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = TRUE)

# Read a TSV file containing data on AuaRhod1.1 TRF (tandem repeat finder) results
aua_reps <- read_tsv("~/Downloads/windowed100K_nxAuaRhod1.1.trf.tsv",
                     col_names = c("chrom", "start", "end", "unit_size",
                                   "num_copies", "perc_ident", "perc_indels",
                                   "score", "entropy", "sequence")) %>%
  mutate(prestart = as.integer(sub(".+:(.+)-.+", "\\1", chrom)),
         start = start + prestart,
         end = end + prestart,
         chrom = sub(":.+", "", chrom))

# Filter AuaRhod1.1 TRF data for rows where the score is greater than 80000
# Group the filtered data by chromosome and select the row with the minimum unit_size for each group
# Ungroup the data and create a new column seq_for_clust by concatenating the sequence with itself
# Select the chrom and seq_for_clust columns and write them to a TSV file
filter(aua_reps, score > 80000) %>%
  group_by(chrom) %>%
  slice_min(unit_size) %>%
  ungroup() %>%
  mutate(seq_for_clust = paste0(sequence, sequence)) %>%
  select(chrom, seq_for_clust) %>%
  write_tsv("~/Downloads/arhod_reps_for_vsearch.tsv", col_names = FALSE)