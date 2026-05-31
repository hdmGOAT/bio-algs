def compute_nx(lengths, x):
    """
    Compute Nx (e.g. N50, N75, NGA50) from a list of contig or aligned block lengths.
    """
    lengths = sorted(lengths)
    total = sum(lengths)
    threshold = total * (x / 100.0)

    running = 0
    for i in range(len(lengths) - 1, -1, -1):
        running += lengths[i]
        if running >= threshold:
            return lengths[i]


def main():
    contig_lens = [20, 20, 30, 30, 60, 60, 80, 100, 200]

    n50 = compute_nx(contig_lens, 50)
    n75 = compute_nx(contig_lens, 75)

    print("Original assembly:")
    print(f"  N50  = {n50}")
    print(f"  N75  = {n75}")
    print()

    aligned_blocks = [20, 20, 30, 30, 50, 50, 60, 60, 80, 200]

    nga50 = compute_nx(aligned_blocks, 50)
    nga75 = compute_nx(aligned_blocks, 75)

    print("After misassembly (aligned blocks):")
    print(f"  NGA50 = {nga50}")
    print(f"  NGA75 = {nga75}")


if __name__ == "__main__": main()
