# Required Libraries
library(dplyr)
library(ggtree)
library(ggplot2)
library(treeio)

# Set input directory (location of all input files) - Should not be changed
setwd("/projects/intro2gds/I2GDS2025/G4_Viruses/mitchell/individual_project/inputs")

# Set output directory variable for saving
indir <- "PATH/TO/DESIRED/OUTPUT"

# Taxonomic lineage input file
tax <- read.csv("tax_map.csv", header = T, sep = ",")

# Sample metadata input file
metadata <- read.delim("sequence_info.tsv", sep = "\t", stringsAsFactors = FALSE)

# Tree file input
contree_24557 <- read.tree("24-557_aln.fna.contree")

# Tree input variable for each sample along with sample name
trees <- list(
   "24-557 (Host: Mum)" = contree_24557)

# Correct tip label naming on a sample-by-sample process
extras_list <- list(
  "24-557 (Host: Mum)" = data.frame(
    accession_tips = c("Chrysanthemum", "Carlavirus"),
    accession_species = c(
      "24-557 (Host: Mum) - Chrysanthemum virus B contig",
      "24-557 (Host: Mum) - Chrysanthemum virus R contig"
    ),
    stringsAsFactors = FALSE
  )
#  Add second sample if necessary
)


# Storage for plots
plots <- list()

# Set sample_id variable
sample_ids <- c("24-557") # Add more samples as necessary

# Build sample_id bolding pattern (Sample1|Sample2|Sample3|etc.)
bold_pattern <- paste(sample_ids, collapse = "|")

# Modify tip labels/plot trees for each sample in trees list + extras_list
for (sample_id in names(trees)) {

  contree <- trees[[sample_id]]

  # Extract tree tip labels
  tip_labels <- contree$tip.label
  tip_labels <- sub("^([A-Z]{1,4}_?[0-9]+\\.[0-9]+).*", "\\1", tip_labels) # Strips tip label down to only accession number
  accession_tips <- grep("^[A-Z]{1,4}_?[0-9]+\\.[0-9]+$", tip_labels, value = TRUE) # Sets only accession number tips (reference genome related tips) to a separate variable for easier editing

  # Filter tip labels and add species information from taxonomic lineage file
  filtered_tax <- tax %>%
    filter(Accession %in% accession_tips) %>%
    transmute(accession_species = paste(Accession, Species.name, sep = " ")) # Adds species information to accession numbers that are present in the taxonomic lineage file "tax"

  # Sort updated tip labels by species name
  filtered_tax_sorted <- filtered_tax %>% arrange(accession_species)

  # Prepare tip labels to combine with the new sample tip labels (i.e., Chrysanthemum > 24-557 (Host: Mum) - Chrysanthemum virus B contig)
  accession_tips_df <- as.data.frame(sort(accession_tips))
  tip_df <- bind_cols(accession_tips_df, filtered_tax_sorted) %>%
    rename(accession_tips = 1)

  # Select the correct sample tip labels for each sample
  extras_df <- extras_list[[sample_id]]

  # Add Geographic location to each accession tip label + species resulting in: AccessionNumber Species (Geographic location)
  other_df <- tip_df %>% rename(Accession = accession_tips)
  combined_df <- left_join(other_df, metadata, by = "Accession") %>%
    mutate(
      formatted_location = gsub(":(\\S)", ": \\1", `Geographic.Location`),
      location_parens = ifelse(is.na(formatted_location) | formatted_location == "",
                               "",
                               paste0("(", formatted_location, ")")),
      Metadata_Combined = paste(`accession_species`, location_parens, sep = " ")
    )
  
  # Prepare tip_df to be used to rename the original tip labels with the newly formatted tip labels
  tip_df <- bind_cols(tip_df$accession_tips, combined_df$Metadata_Combined) %>%
    rename(accession_tips = ...1,
           accession_species = ...2)
  # Adds the sample tip labels which are specified prior to the function (i.e., Chrysanthemum > 24-557 (Host: Mum) - Chrysanthemum virus B contig)
  full_df <- bind_rows(tip_df, extras_df)

  # Rename tree tip labels with the newly formatted tip labels
  contree <- rename_taxa(contree, full_df, key = "accession_tips", value = "accession_species")

  # Plot the tree for tree data
  contree_plot <- ggtree(contree, branch.length = "none")
  
  # Add tip label bolding if the sample has a specified sample ID within the tip labels (sample related tips will be easier to see for comparison to reference genomes)
  contree_data <- contree_plot$data %>%
    mutate(label_bold = ifelse(grepl(bold_pattern, label), "bold", "plain"))

  # Plot consensus tree data with bolding and bootstrap values by circles of varying sizes
  p <- ggtree(contree_data) +
    geom_tiplab(aes(label = label, fontface = label_bold), align = TRUE, size = 5) +
    xlim(0, 30) +
    labs(title = paste("Phylogenetic Tree: Sample", sample_id)) +
    geom_point(
      data = subset(contree_plot$data, !is.na(label) & isTip == FALSE),
      aes(x = x, y = y, size = as.numeric(label)),
      color = "steelblue",
      alpha = 0.6
    ) +
    scale_size_continuous(name = "Bootstrap", range = c(1, 6)) +
    theme(legend.position = c(0.80, 0.5))

  # Store plot for current sample_id
  plots[[sample_id]] <- p
  
  # Save current plot as PDF (the file usually has incorrect formatting due to how ggtree displays/saves plots; it is best to print the plot, view it in full screen for the image to resize, then take a screenshot to save it - I was unable to find a fix for this issue)
  ggsave(
    filename = file.path(indir, paste0("phylo_tree-", sample_id, ".pdf")),
    plot = p,
    width = 10, height = 10, dpi = 300,
    bg = "white"
   )
}

# Display first sample plot
print(plots[[1]])
