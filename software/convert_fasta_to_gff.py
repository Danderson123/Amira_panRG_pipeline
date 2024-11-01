import argparse
from Bio import SeqIO
import os
from tqdm import tqdm

# Initialize the ArgumentParser
parser = argparse.ArgumentParser(description="Convert a FASTA file to a GFF3 file.")

# Add the argument for the FASTA file
parser.add_argument("fasta_file", help="The path to the FASTA file.")
parser.add_argument("output_file", help="The path to the output gff file.")

# Parse the arguments
args = parser.parse_args()

# GFF3 format headers
header = '##gff-version 3\n'

junk_sequences = []
# Open the GFF3 file
with open(args.output_file, 'w') as f:
    f.write(header)

    for record in tqdm(SeqIO.parse(args.fasta_file, "fasta")):
        # skip if any chars are not "A", "C", "T" or "G"
        if any(base not in ["A", "C", "T", "G"] for base in str(record.seq)):
            junk_sequences.append(record.id)
            continue
        # Each line is a separate CDS feature.
        # Format: seqid source type start end score strand phase attributes
        line = '{}\t.\tCDS\t1\t{}\t.\t+\t0\tID={};Name={};gene={}\n'.format(record.id, len(record.seq), record.id, record.id, record.id)
        f.write(line)

    # Include the sequences after a ##FASTA line
    f.write("##FASTA\n")
    for record in SeqIO.parse(args.fasta_file, "fasta"):
        f.write(">{}\n{}\n".format(record.id, record.seq))

with open(os.path.join(os.path.dirname(args.output_file), "junk_alleles.txt"), "w") as o:
    o.write("\n".join(junk_sequences))
print("GFF3 file has been written successfully.")

