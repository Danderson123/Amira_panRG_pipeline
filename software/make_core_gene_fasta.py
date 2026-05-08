import glob
from tqdm import tqdm
import os

allele_dir = "/hps/nobackup/iqbal/dander/Amira_panRG_pipeline_test/Enterococcus_faecium_panRG/qced_unaligned_gene_sequences"
core_gene_file = "/hps/nobackup/iqbal/dander/Amira_panRG_pipeline_test/Enterococcus_faecium_panRG/core_genes.txt"
outfile = "/hps/nobackup/iqbal/dander/Amira_panRG_pipeline_test/Enterococcus_faecium_panRG/core_genes.fasta"

with open(core_gene_file) as i:
    core_genes = i.read().split("\n")
core_alleles = []
for c in tqdm(core_genes):
    with open(os.path.join(allele_dir, f"{c}.fasta")) as i:
        alleles = i.read().split(">")[1:]
    for a in alleles:
        seq = "".join(a.split("\n")[1:])
        core_alleles.append(f">{c}\n{seq}")
        # take only the first allele
        break
with open(outfile, "w") as o:
    o.write("\n".join(core_alleles))