library(GenomicRanges)   # Load the GenomicRanges library for genomic range operations
library(tidyverse)       # Load the tidyverse library for data manipulation

#### Genomic ranges

# Create Seqinfo object with sequence names and lengths
gnm_gr <- Seqinfo(seqnames = seq_sizes$multispecies_sequence, seqlengths = seq_sizes$size,
                  isCircular = rep(F, nrow(seq_sizes)), genome = "Auanema")

# Create GRanges object for TRF data
trf_gr <- mutate(rep_for_plot, chrom = paste0("nxAuaRhod1.1_", saccver),
                 start = ifelse(sstart < send, sstart, send),
                 end = ifelse(sstart > send, sstart, send),
                 rep_id = qaccver) %>%
  select(chrom, start, end, rep_id) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = T)

# Set the window size for binning
window_size <- 100000

# Create tiled genome of 1 kb bins
aua_1kb <- tileGenome(gnm_gr, tilewidth = window_size, cut.last.tile.in.chrom = T)

### Get coverage signal as Rle object

### Get average coverage in each bin

clust_abundance <- tibble(seqnames = character(), start = integer(), end = integer(),
                          width = integer(), strand = factor(), binned_cov = double(),
                          clust_id = character())

# Calculate average coverage in each bin for each cluster
for (clust_id in unique(trf_gr$rep_id)) {
  gr1_cov <- GenomicRanges::coverage(GenomicRanges::reduce(trf_gr[trf_gr$rep_id == clust_id]))
  clust_abundance <- GenomicRanges::binnedAverage(aua_1kb, gr1_cov, "binned_cov") %>%
    as_tibble() %>%
    mutate(clust_id = clust_id) %>%
    rbind(clust_abundance)
}

# Combine break_sites and clust_abundance data for plotting
trf_plt <- bind_rows(dplyr::select(break_sites, seqnames = multispecies_sequence, diminution_pos),
                     clust_abundance) %>%
  mutate(chr = sub("[^_]+[_\\.]1_SUPER_(.+)", "chr\\1", seqnames),
         clust_id = sub("cluster", "", clust_id),
         clust_id = factor(clust_id, levels = order(clust_id))) %>%
  filter(!grepl("MT|unloc|scaff", seqnames)) %>%
  ggplot(aes(x = start, y = binned_cov, fill = clust_id, diminution_pos)) +
  facet_grid(chr ~ .) +
  geom_bar(position = "stack", stat = "identity", width = window_size, alpha = 0.8) +
  geom_vline(aes(xintercept = diminution_pos), linetype = "dotted",
             color = "black", size = 1) +
  scale_x_continuous(labels = label_number_si()) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(0, 1)) +
  ylab("fraction per 100 Kb window") +
  xlab("Position along chromosome") +
  theme_bw()