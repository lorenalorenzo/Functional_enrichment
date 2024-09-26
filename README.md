Enrichment in introgressed regions
================
Lorena Lorenzo Fernández
26-Sep-2024

## Testing for over and under representation on introgressed regions between Iberian and Eurasian lynxes

**DISCLAIMER: This is part of Enrico’s proyect in detecting and
characterising introgression in lynxes. Go to ebazzicalupo GitHub page
for more information.**

Here, we are analyzing the excess or lack of functions in GO terms in a
list of genes obtained from the introgressed regions. The data consist
of:

1.  List of genes introgressed from Eurasian to Iberian
    (“ab.geneids.txt”)
2.  List of genes introgressed from Iberian to Eurasian
    (“ba.geneids.txt”)
3.  List of genes introgressed in both directions (“bi.geneids.txt”)
4.  Table with the annotation of the genes and the corresponding GO
    terms (“LYRU_2A.FA.genego_table.tsv”)

``` r
# Load necessary libraries
library(topGO)
library(tidyverse)
```

``` r
# Custom method to perform Fisher's exact test for underrepresentation (less frequent GO terms)
if (!isGeneric("GOFisherUnder")) {
  setGeneric("GOFisherUnder", function(object) standardGeneric("GOFisherUnder"))
}

setMethod("GOFisherUnder", "classicCount", function(object) {
  contMat <- contTable(object)
  if (all(contMat == 0)) {
    p.value <- 1
  } else {
    p.value <- fisher.test(contMat, alternative = "less")$p.value
  }
  return(p.value)
})
```

``` r
# Define the gene sets to analyze
gene_sets <- list(
  ab = "ab.geneids.txt",
  ba = "ba.geneids.txt",
  bi = "bi.geneids.txt"
)
# Read the gene-to-GO mappings from a .tsv file
gene2GO_df <- read.delim("LYRU2_2A.FA.genego_table.tsv", header = TRUE, sep = "\t")

# Convert the gene2GO data frame into a named list where each gene ID points to a vector of GO terms
gene2GO <- setNames(strsplit(gene2GO_df$go_terms, ","), gene2GO_df$gene_id)

# Get the list of all unique genes from the gene2GO mapping
allGenesList <- unique(gene2GO_df$gene_id)

# Loop through each gene set
for (set_name in names(gene_sets)) {
  
  # Read the list of candidate genes for the current gene set
  candidateGenes <- readLines(gene_sets[[set_name]])
  
  # Create a named vector (geneList) for all genes:
  # 1 for candidate genes, 0 for non-candidate genes
  geneList <- setNames(as.numeric(allGenesList %in% candidateGenes), allGenesList)
  
  # Create the topGOdata object for both overrepresentation and underrepresentation tests
  GOdata <- new(
    "topGOdata",
    ontology = "BP",  # Choose "BP", "MF", or "CC" as appropriate
    allGenes = geneList,  # The named vector of all genes with 1 for candidates, 0 for non-candidates
    geneSel = function(p) p == 1,  # Function to select candidate genes (where p == 1)
    annot = annFUN.gene2GO,  # Annotation function for gene2GO mappings
    gene2GO = gene2GO  # Your gene-to-GO mappings
  )
  
  # Perform Overrepresentation Test using the "classic" algorithm and default Fisher's exact test
  resultFisherOver <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  # Set up a new test for underrepresentation using the custom GOFisherUnder function
  test.stat <- new("classicCount", testStatistic = GOFisherUnder, name = "Fisher's exact test for underrepresentation")
  
  # Perform Underrepresentation Test using the "classic" algorithm and custom GOFisherUnder function
  resultFisherUnder <- getSigGroups(GOdata, test.stat)
  
  # Get all GO terms used in the analysis
  allGO <- usedGO(GOdata)
  
#  # Generate a table of results for overrepresentation
#  resultsOver <- GenTable(GOdata, classicFisher = resultFisherOver, orderBy = "classicFisher", #topNodes = length(allGO))
#  
#  # Generate a table of results for underrepresentation
#  resultsUnder <- GenTable(GOdata, classicFisher = resultFisherUnder, orderBy = "classicFisher"#, topNodes = length(allGO))
  
  
 result_over <- GenTable(GOdata, Fisher=resultFisherOver, topNodes=resultFisherOver@geneData[2], numChar=1000) %>% 
        as_tibble() %>% 
        mutate(p.adj = round(p.adjust(as.numeric(gsub("<", "", Fisher)), method="BH"), 15)) %>%
        filter(p.adj<0.05) 
 
  result_under <- GenTable(GOdata, Fisher=resultFisherUnder, topNodes=resultFisherUnder@geneData[2], numChar=1000) %>% 
        as_tibble() %>% 
        mutate(p.adj = round(p.adjust(as.numeric(gsub("<", "", Fisher)), method="BH"), 15)) %>%
        filter(p.adj<0.05) 

  # Save the results to files
  write.table(result_over, paste0("significant_overrepresented_GO_terms_", set_name, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(result_under, paste0("significant_underrepresented_GO_terms_", set_name, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
```
