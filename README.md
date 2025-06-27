# GO Enrichment Analysis Pipeline

This repository contains a generalized and reusable pipeline for performing Gene Ontology (GO) enrichment analysis using the R package **topGO**.

Originally designed for [Enrico's lynx introgression project](https://github.com/Enricobazzi/Lynxtrogression_v2), this pipeline has been modularized to be applicable to any gene set.

---

## 📦 Requirements

- R (≥ 4.0 recommended)
- R packages: `topGO`, `tidyverse`, `optparse` and `biomaRt` if used with ensembl annotation

Automatically installed and load when running the run_enrichment.R script

---

## 📥 Input Files

### 1. `candidate_genes.txt`
A plain text file listing one gene ID per line:

```
GeneA
GeneB
GeneC
...
```

### 2. `gene_to_go.tsv` *(only if using `--annotation_source custom`, see below)*
A tab-delimited file with two columns:
- `gene_id`: the gene name
- `go_terms`: semicolon-separated GO terms

You can directly use the [gene_to_go_gff3.py](https://github.com/lorenalorenzo/Functional_enrichment/blob/main/scripts/gene_to_go_gff3.py) Python script to extract the `gene_to_go.tsv` file. The script inspects the gff3 looking for transcript rows and takes the Parent ID (gene_id) and the GO ontology terms associated (go_terms).

Usage:
```bash
python gene_to_go_gff3.py <input.gff3> gene_to_go.tsv
```

This file is ignored if `--annotation_source ensembl` is used, in which case you must specify `--ensembl_dataset` to select the species.

Example:
```
gene_id    go_terms
GeneA      GO:0008150,GO:0003674
GeneB      GO:0009987
```

---

## 🚀 Running the Pipeline

From the project root, run the script via Rscript:

### Option A: Use custom/local annotation

```bash
Rscript scripts/run_enrichment.R \
  --genes data/candidate_genes.txt \
  --annotation data/gene_to_GO.tsv \
  --annotation_source custom \
  --out results/my_project
```

### Option B: Use Ensembl annotation 

```bash
Rscript scripts/run_enrichment.R \
  --genes data/candidate_genes.txt \
  --annotation_source ensembl \
  --ensembl_dataset fcatus_gene_ensembl \
  --out results/my_project_ensembl
```

Use the `--ensembl_dataset` argument to specify the species dataset from Ensembl. Example datasets:
- `fcatus_gene_ensembl` — domestic cat (Felis catus)
- `hsapiens_gene_ensembl` — human
- `mmusculus_gene_ensembl` — mouse

You can find the full list using the listDatasets(ensembl) function of biomaRt package.

⚠️ When using `--annotation_source ensembl`, you do **not** need to specify a `--annotation` file. The script will download the gene-to-GO mappings directly from Ensembl using the `biomaRt` package.

In this case, the script will fetch GO terms for **Felis catus** from Ensembl using the `biomaRt` package.


---

## 📤 Output Files

- `results/my_project_over.tsv`: Full enrichment results for all GO terms
- `results/my_project_significant.tsv`: Filtered results with Fisher < 0.05 and >1 significant gene

Each row includes:
- `GO.ID`: GO term identifier
- `Term`: Description of the GO term
- `Significant`: Number of candidate genes associated with the term
- `Fisher`: Raw p-value from Fisher's exact test
- `Genes`: Semicolon-separated list of candidate genes contributing to the term, using gene names when available (otherwise Ensembl IDs)
---

## 🧪 Notes

- This script uses the "Biological Process" (BP) ontology by default.
  - You can change it to "MF" (Molecular Function) or "CC" (Cellular Component) by modifying the line `ontology = "BP"` inside the script.
- GO terms are considered significantly enriched if:
  - Fisher's exact test p-value < 0.05, **and**
  - More than 1 candidate gene is annotated to the term.
- The output includes a `Genes` column with contributing gene names (or Ensembl IDs if names are unavailable).
- When using Ensembl annotations, only **protein-coding** genes are considered.