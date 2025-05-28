
import pysam

import multiprocessing
import argparse

def process_sequence(bam_file, sequence_name):
    """
    Process all reads from a given sequence name and compute the RE tag ratio for each CB tag.
    Args:
        bam_file (str): Path to the BAM file.
        sequence_name (str): Sequence name to fetch reads for.
    Returns:
        dict: A dictionary with CB tags as keys and the ratio of N over (N + E) as values.
    """
    re_counts = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(sequence_name):
            if read.has_tag('CB') and read.has_tag('RE'):
                cb_tag = read.get_tag('CB')
                re_tag = read.get_tag('RE')
                if cb_tag not in re_counts:
                    re_counts[cb_tag] = {"N": 0, "E": 0}
                if re_tag == "N":
                    re_counts[cb_tag]["N"] += 1
                elif re_tag == "E":
                    re_counts[cb_tag]["E"] += 1

    # Compute ratios
    ratios = {}
    for cb_tag, counts in re_counts.items():
        n = counts["N"]
        e = counts["E"]
        ratios[cb_tag] = n / (n + e) if (n + e) > 0 else 0
    return ratios

def worker(bam_file, sequence_name):
    """
    Worker function for multiprocessing to handle one sequence.
    Args:
        bam_file (str): Path to the BAM file.
        sequence_name (str): Sequence name to process.
    Returns:
        dict: A dictionary summarizing ratios for the given sequence.
    """
    return process_sequence(bam_file, sequence_name)

def parallel_process_bam(bam_file, num_threads):
    """
    Parallel process a BAM file to compute the ratio of N over (N + E) for each CB tag using sequences.
    Args:
        bam_file (str): Path to the BAM file.
        num_threads (int): Number of parallel threads to use.
    Returns:
        dict: A dictionary with CB tags as keys and the computed ratio as values.
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Get the list of sequence names (references) in the BAM file
        sequence_names = bam.references

    # Use multiprocessing to process each sequence
    with multiprocessing.Pool(num_threads) as pool:
        results = pool.starmap(worker, [(bam_file, seq) for seq in sequence_names])
    
    # Merge results
    merged_ratios = {}
    for result in results:
        merged_ratios.update(result)

    return merged_ratios

def main():
    parser = argparse.ArgumentParser(description="Summarize RE tag ratios for each CB tag in a BAM file.")
    parser.add_argument("bam_file", help="Path to the input BAM file.")
    parser.add_argument("output_file", help="Path to the output TSV file.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use (default: 4).")
    
    args = parser.parse_args()
    bam_file = args.bam_file
    output_file = args.output_file
    num_threads = args.threads
    
    # Compute ratios
    ratios = parallel_process_bam(bam_file, num_threads)
    
    # Write results to output file
    with open(output_file, "w") as out:
        out.write("barcode\tnuclear_fraction\n")  # Header
        for cb_tag, ratio in ratios.items():
            out.write(f"{cb_tag}\t{ratio:.4f}\n")
    
    print(f"Results written to {output_file}")

if __name__ == "__main__":
    main()
