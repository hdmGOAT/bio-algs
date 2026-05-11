from pathlib import Path

from bioalgs.genome_assembly import genome_path

path = Path(__file__).parent / "datasets" / "dataset_30182_3.txt"

kmers = path.read_text().split()

with open("output.txt", "w") as f:
    f.write(genome_path(kmers))
