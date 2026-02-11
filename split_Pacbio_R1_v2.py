import gzip
import os, sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import csv
from tqdm import tqdm

def reverse_complement(seq):
    """Generate inverse complementary sequences"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))

def hamming_distance(seq1, seq2):
    """Calculate the Hamming distance, allowing a maximum of mismatches of `mismatch_distance`."""
    if len(seq1) != len(seq2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def process_fastq_single_end(r1_path, mapping_file, output_dir, mismatch_distance):
    """
    Split the R1 FASTQ file based on the inverse complementary sequences at the beginning and end of the TSV file.
    Specify the number of mismatches.
    """
    os.makedirs(output_dir, exist_ok=True)
    os.system(f"cp {mapping_file} {output_dir}")

    # === Read TSV mapping file ===
    barcode_patterns = {}
    with open(mapping_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) >= 3:
                filename = row[0]
                start_seq = row[1].upper()
                end_seq = row[2].upper()
                end_seq_rc = reverse_complement(end_seq)
                start_seq_rc = reverse_complement(start_seq)
                barcode_patterns[filename] = {
                    'start': start_seq,
                    'end_rc': end_seq_rc,
                    'start_len': len(start_seq),
                    'end_len': len(end_seq_rc),
                    'start_rc': start_seq_rc,
                    'end': end_seq
                }

    records_in_memory = defaultdict(list)
    batch_size = 5000000
    processed_count = 0
    matched_count = 0
    unmatched_count = 0

    r1_open = gzip.open if r1_path.endswith('.gz') else open

    with r1_open(r1_path, 'rt') as r1_in:
        r1_iterator = SeqIO.parse(r1_in, 'fastq')

        for record in tqdm(r1_iterator, desc="Processing reads"):
            seq = str(record.seq).upper()
            seq_len = len(seq)
            matched = False

            for filename, pattern in barcode_patterns.items():
                start_seq = pattern['start']
                end_seq_rc = pattern['end_rc']
                start_len = pattern['start_len']
                end_len = pattern['end_len']
                start_seq_rc = pattern['start_rc']
                end_seq = pattern['end']

                if seq_len < start_len + end_len:
                    continue

                start_match = hamming_distance(seq[:start_len], start_seq) <= mismatch_distance
                end_match = hamming_distance(seq[-end_len:], end_seq_rc) <= mismatch_distance

                start_match_2 = hamming_distance(seq[:end_len], end_seq) <= mismatch_distance
                end_match_2 = hamming_distance(seq[-start_len:], start_seq_rc) <= mismatch_distance

                if (start_match and end_match) or (start_match_2 and end_match_2):
                    records_in_memory[filename].append(record)
                    matched_count += 1
                    matched = True
                    break

            if not matched:
                records_in_memory['unmatched'].append(record)
                unmatched_count += 1

            processed_count += 1

            if processed_count % batch_size == 0:
                write_batches(records_in_memory, output_dir)
                for key in list(records_in_memory.keys()):
                    if key != 'unmatched':
                        records_in_memory[key].clear()

        write_batches(records_in_memory, output_dir, final_write=True)

    print(f"Processing complete! A total of {processed_count} reads were processed.")
    print(f"Matched reads: {matched_count}")
    print(f"Unmatched reads: {unmatched_count}")
    print(f"Match rate: {matched_count/processed_count*100:.2f}%")

def write_batches(records, output_dir, final_write=False):
    """Write records in memory in batches to a file"""
    for filename, record_list in records.items():
        if not record_list and not final_write:
            continue

        # === If it is unmatched, then output separately as follows: unmatched_R1.fastq.gz ===
        if filename == 'unmatched':
            output_path = os.path.join(output_dir, "unmatched_R1.fastq.gz")
        else:
            output_path = os.path.join(output_dir, f"{filename}.fastq.gz")

        with gzip.open(output_path, 'at') as out_file:
            SeqIO.write(record_list, out_file, 'fastq')

        if filename != 'unmatched' or final_write:
            record_list.clear()

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python split_Pacbio_R1_v2.py <R1 file> <TSV mapping file> <output directory> <distance>")
        print("TSV file format: First column = output file name, second column = first , third column = last ")
        sys.exit(1)

    r1_file = sys.argv[1]
    tsv_file = sys.argv[2]
    output_directory = sys.argv[3]
    distance = int(sys.argv[4])
    process_fastq_single_end(r1_file, tsv_file, output_directory, distance)
