#!/usr/bin/env python3
"""
Trim sequences to region between positions 1233 and 1682 (inclusive) using Biopython.
"""

from Bio import SeqIO

def trim_sequences(input_fasta, output_fasta, start_pos, end_pos):
    """
    Trim sequences to specified region.
    Positions are 1-based (as typically used in biology).
    Biopython uses 0-based indexing, so we convert.
    """
    
    # Convert to 0-based indexing for Python
    start_idx = start_pos - 1  # 1233 -> 1232
    end_idx = end_pos          # 1682 -> 1682 (inclusive end in slice)
    
    sequences_processed = 0
    sequences_trimmed = 0
    sequences_too_short = 0
    
    trimmed_sequences = []
    
    print(f"Trimming sequences to positions {start_pos}-{end_pos} (1-based, inclusive)")
    print(f"Python slice indices: [{start_idx}:{end_idx}]")
    print(f"Expected trimmed length: {end_pos - start_pos + 1} bp\n")
    
    with open(input_fasta, 'r') as fin:
        for record in SeqIO.parse(fin, "fasta"):
            sequences_processed += 1
            
            seq_length = len(record.seq)
            
            # Check if sequence is long enough
            if seq_length < end_pos:
                sequences_too_short += 1
                # For short sequences, take what we can
                if seq_length > start_idx:
                    trimmed_seq = record.seq[start_idx:]
                    record.seq = trimmed_seq
                    trimmed_sequences.append(record)
                    sequences_trimmed += 1
                continue
            
            # Trim the sequence
            trimmed_seq = record.seq[start_idx:end_idx]
            record.seq = trimmed_seq
            trimmed_sequences.append(record)
            sequences_trimmed += 1
            
            # Progress indicator
            if sequences_processed % 1000 == 0:
                print(f"Processed {sequences_processed} sequences...")
    
    # Write trimmed sequences to output file
    with open(output_fasta, 'w') as fout:
        SeqIO.write(trimmed_sequences, fout, "fasta")
    
    print(f"\n{'='*60}")
    print(f"TRIMMING SUMMARY")
    print(f"{'='*60}")
    print(f"Total sequences processed: {sequences_processed}")
    print(f"Sequences trimmed and kept: {sequences_trimmed}")
    print(f"Sequences too short (< {end_pos} bp): {sequences_too_short}")
    
    if trimmed_sequences:
        # Show length distribution
        lengths = [len(rec.seq) for rec in trimmed_sequences]
        expected_length = end_pos - start_pos + 1
        full_length_count = sum(1 for l in lengths if l == expected_length)
        
        print(f"\nTrimmed sequence lengths:")
        print(f"  Full length ({expected_length} bp): {full_length_count}")
        print(f"  Shorter (original seq too short): {len(lengths) - full_length_count}")
        print(f"  Min length: {min(lengths)} bp")
        print(f"  Max length: {max(lengths)} bp")

def main():
    input_fasta = "/Users/sarajordan/Desktop/Gobal_Measles_DTA_analysis/Filtered_sequences_with_metadata_no_N.fasta"
    output_fasta = "/Users/sarajordan/Desktop/Gobal_Measles_DTA_analysis/Trimmed_sequences_1233-1682.fasta"
    
    # Positions 1233 to 1682 (1-based, inclusive)
    start_position = 1233
    end_position = 1682
    
    trim_sequences(input_fasta, output_fasta, start_position, end_position)
    
    print(f"\nOutput written to: {output_fasta}")

if __name__ == "__main__":
    main()
