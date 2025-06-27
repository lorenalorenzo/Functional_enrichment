import sys

def extract_gene_to_go(gff3_file, output_file):
    gene_to_go = {}

    with open(gff3_file, 'r') as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            if parts[2] != "transcript":
                continue

            attr_field = parts[8]

            # Get the gene ID from Parent=
            parent_id = None
            for field in attr_field.split(';'):
                if field.startswith("Parent="):
                    parent_id = field.replace("Parent=", "").strip()
                    break

            if not parent_id or "Ontology_term=" not in attr_field:
                continue

            # Extract everything after Ontology_term=
            go_part = attr_field.split("Ontology_term=")[-1]
            go_terms = [go.strip() for go in go_part.split(';') if go.startswith("GO:")]

            if go_terms:
                if parent_id in gene_to_go:
                    gene_to_go[parent_id].update(go_terms)
                else:
                    gene_to_go[parent_id] = set(go_terms)

    with open(output_file, 'w') as out:
        out.write("gene_id\tgo_terms\n")  # header
        for gene_id, go_list in gene_to_go.items():
            out.write(f"{gene_id}\t{','.join(sorted(go_list))}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_gene_to_go_from_gff3.py input.gff3 output.tsv")
        sys.exit(1)

    extract_gene_to_go(sys.argv[1], sys.argv[2])