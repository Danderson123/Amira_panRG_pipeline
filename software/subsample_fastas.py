import os
import random
import shutil
import glob
from tqdm import tqdm

random.seed(2025)

def subsample_fastas(input_dir, output_dir, sample_size=500):
    # Ensure the input directory exists
    if not os.path.isdir(input_dir):
        raise ValueError(f"Input directory does not exist: {input_dir}")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Get all FASTA files in the input directory
    fasta_files = glob.glob(os.path.join(input_dir, "*", "*.fna"))

    if not fasta_files:
        raise ValueError("No FASTA files found in the input directory.")

    # Determine the number of files to sample
    sample_size = min(sample_size, len(fasta_files))

    # Randomly select files
    sampled_files = random.sample(fasta_files, sample_size)

    # Copy the selected files to the output directory
    for file in tqdm(sampled_files):
        src = os.path.join(input_dir, file)
        dst = os.path.join(output_dir, os.path.basename(file).replace(".fna", ".fa"))
        shutil.copy(src, dst)

    print(f"Sampled {sample_size} FASTA files to {output_dir}")

if __name__ == "__main__":
    # Example usage:
    input_directory = "/hps/nobackup/iqbal/dander/Amira_panRG_pipeline_test/data/reference_assemblies/ncbi_dataset/data"
    output_directory = "/hps/nobackup/iqbal/dander/Amira_panRG_pipeline_test/data/reference_assemblies/metagenome"
    subsample_fastas(input_directory, output_directory, sample_size=100)
