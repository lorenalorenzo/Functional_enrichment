# run_enrichment.R
# Generalized GO enrichment script using topGO

library(topGO)
library(tidyverse)
library(optparse)

# ---------------------
# Command-line options
# ---------------------
option_list <- list(
  make_option(c("-g", "--genes"), type="character", default=NULL,
              help="Candidate gene list file (one gene per line)", metavar="character"),

  make_option(c("--annotation_source"), type = "character", default = "custom",
            help = "Source of annotation: 'custom' or 'ensembl'", metavar = "character"),

  make_option(c("--ensembl_dataset"), type = "character", default = "fcatus_gene_ensembl",
            help = "Dataset name in Ensembl (e.g., fcatus_gene_ensembl, hsapiens_gene_ensembl)", metavar = "character"),

  make_option(c("--annotation"), type="character", default=NULL,
              help="Gene-to-GO annotation file (only used if annotation_source is 'custom')", metavar="character"),

  make_option(c("-o", "--out"), type="character", default="go_enrichment", 
              help="Output file prefix [default %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$genes)) {
  stop("You must provide --genes file.")
}
if (opt$annotation_source == "custom" && is.null(opt$annotation)) {
  stop("You must provide --annotation file when using annotation_source 'custom'.")
}

# ---------------------
# Load input files
# ---------------------
candidateGenes <- readLines(opt$genes)

if (opt$annotation_source == "ensembl") {
  library(biomaRt)

  message("🔄 Downloading annotation from Ensembl...")
  ensembl <- useMart("ensembl", dataset = opt$ensembl_dataset)
  attributes <- c("ensembl_gene_id", "external_gene_name", "go_id", "gene_biotype")

  ensembl_to_go <- getBM(attributes = attributes, mart = ensembl) %>%
    filter(gene_biotype == "protein_coding")

  gene2GO <- split(ensembl_to_go$go_id, ensembl_to_go$ensembl_gene_id)
  allGenesList <- unique(ensembl_to_go$ensembl_gene_id)

  # Map for gene name lookup
  gene_id_to_name <- ensembl_to_go %>%
    select(ensembl_gene_id, external_gene_name) %>%
    distinct() %>%
    deframe()

} else {
  gene2GO_df <- read.delim(opt$annotation, header = TRUE, sep = "\t")
  gene2GO <- setNames(strsplit(gene2GO_df$go_terms, ";"), gene2GO_df$gene_id)
  allGenesList <- unique(gene2GO_df$gene_id)

  # Use gene_id as name if no external names are given
  gene_id_to_name <- setNames(allGenesList, allGenesList)
}

geneList <- setNames(as.numeric(allGenesList %in% candidateGenes), allGenesList)

# ---------------------
# Create topGO object
# ---------------------
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSel = function(p) p == 1,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)

# ---------------------
# Overrepresentation test
# ---------------------
over_test <- runTest(GOdata, statistic = "fisher")

# ---------------------
# Output results
# ---------------------
p_values <- score(over_test)
result_table <- GenTable(GOdata, Fisher=over_test, topNodes=length(p_values), numChar=1000) %>%
  as_tibble()

# Get candidate genes per GO term
term_genes <- genesInTerm(GOdata)
gene_lists <- lapply(result_table$GO.ID, function(go_id) {
  intersect(term_genes[[go_id]], candidateGenes)
})

# Convert to names if available
gene_lists_named <- lapply(gene_lists, function(ids) {
  sapply(ids, function(id) ifelse(!is.na(gene_id_to_name[id]), gene_id_to_name[id], id))
})

result_table$Genes <- sapply(gene_lists_named, function(g) paste(g, collapse = ";"))

# Filter significant
significant_table <- result_table %>%  
  filter(Fisher < 0.05 & Significant > 1)

# Save output
write.table(result_table, paste0(opt$out, "_over.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(significant_table, paste0(opt$out, "_significant.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nEnrichment analysis completed. Results saved to:\n")
cat(paste0(" - ", opt$out, "_over.tsv\n"))
cat(paste0(" - ", opt$out, "_significant.tsv\n"))
