#!/usr/bin/env python3
"""
Script to process measles sequences:
1. Add metadata (sampling date, country, sequence length, genotype) to sequence headers
2. Remove sequences without complete sampling dates (but include YYYY-MM dates with day set to 15)
3. Remove sequences with blank country field
"""

import sys
from collections import defaultdict

def parse_metadata(metadata_file):
    """Parse the metadata TSV file and create a dictionary keyed by accession."""
    metadata = {}
    
    with open(metadata_file, 'r') as f:
        header = f.readline().strip().split('\t')
        
        # Find column indices
        acc_idx = header.index('accessionVersion')
        genotype_idx = header.index('genotype')
        date_idx = header.index('sampleCollectionDate')
        country_idx = header.index('geoLocCountry')
        length_idx = header.index('length')
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) > max(acc_idx, genotype_idx, date_idx, country_idx, length_idx):
                accession = fields[acc_idx]
                genotype = fields[genotype_idx] if fields[genotype_idx] else 'unknown'
                date = fields[date_idx] if fields[date_idx] else ''
                country = fields[country_idx] if fields[country_idx] else ''
                length = fields[length_idx] if fields[length_idx] else ''
                
                metadata[accession] = {
                    'genotype': genotype,
                    'date': date,
                    'country': country,
                    'length': length
                }
    
    return metadata

def process_date(date_str):
    """
    Process the date string:
    - Return None if empty or incomplete (not YYYY-MM-DD or YYYY-MM format)
    - If YYYY-MM format, add '-15' to make it YYYY-MM-15
    - If YYYY-MM-DD format, return as is
    """
    if not date_str:
        return None
    
    parts = date_str.split('-')
    
    # Complete date YYYY-MM-DD
    if len(parts) == 3:
        return date_str
    
    # Year-month only YYYY-MM - add day 15
    elif len(parts) == 2:
        return f"{date_str}-15"
    
    # Incomplete date (only year or other format)
    else:
        return None

def create_new_header(accession, metadata_dict):
    """
    Create a new header with format: >accession/date/country/length/genotype
    All fields separated by '/' with no whitespaces.
    Replace any spaces in fields with underscores.
    """
    if accession not in metadata_dict:
        return None
    
    meta = metadata_dict[accession]
    
    # Replace spaces with underscores in all fields
    date = meta['date'].replace(' ', '_') if meta['date'] else ''
    country = meta['country'].replace(' ', '_') if meta['country'] else ''
    length = meta['length'].replace(' ', '_') if meta['length'] else ''
    genotype = meta['genotype'].replace(' ', '_') if meta['genotype'] else ''
    
    # Create header
    new_header = f">{accession}/{date}/{country}/{length}/{genotype}"
    return new_header

def process_fasta(input_fasta, output_fasta, metadata_dict):
    """
    Process FASTA file:
    - Add metadata to headers
    - Filter based on date and country requirements
    """
    sequences_processed = 0
    sequences_kept = 0
    sequences_removed_no_date = 0
    sequences_removed_no_country = 0
    
    with open(input_fasta, 'r') as fin, open(output_fasta, 'w') as fout:
        current_accession = None
        current_sequence = []
        skip_current = False
        
        for line in fin:
            line = line.strip()
            
            if line.startswith('>'):
                # Process previous sequence if exists
                if current_accession and not skip_current and current_sequence:
                    fout.write(f"{new_header}\n")
                    fout.write('\n'.join(current_sequence) + '\n')
                    sequences_kept += 1
                
                # Parse new header
                current_accession = line[1:].split()[0]  # Get accession (first part after >)
                current_sequence = []
                skip_current = False
                sequences_processed += 1
                
                # Check if accession exists in metadata
                if current_accession not in metadata_dict:
                    skip_current = True
                    sequences_removed_no_country += 1
                    continue
                
                meta = metadata_dict[current_accession]
                
                # Check country (must not be blank)
                if not meta['country']:
                    skip_current = True
                    sequences_removed_no_country += 1
                    continue
                
                # Process and check date
                processed_date = process_date(meta['date'])
                if processed_date is None:
                    skip_current = True
                    sequences_removed_no_date += 1
                    continue
                
                # Update the date in metadata with processed date
                meta['date'] = processed_date
                
                # Create new header
                new_header = create_new_header(current_accession, metadata_dict)
                
            else:
                # Accumulate sequence lines
                if not skip_current and line:
                    current_sequence.append(line)
        
        # Process last sequence
        if current_accession and not skip_current and current_sequence:
            fout.write(f"{new_header}\n")
            fout.write('\n'.join(current_sequence) + '\n')
            sequences_kept += 1
    
    # Print summary
    print(f"\n=== Processing Summary ===")
    print(f"Total sequences processed: {sequences_processed}")
    print(f"Sequences kept: {sequences_kept}")
    print(f"Sequences removed (no/incomplete date): {sequences_removed_no_date}")
    print(f"Sequences removed (no country): {sequences_removed_no_country}")
    print(f"Total removed: {sequences_removed_no_date + sequences_removed_no_country}")

def main():
    # File paths
    metadata_file = "/Users/sarajordan/Desktop/Gobal_Measles_DTA_analysis/metadata_global_measles_data.tsv"
    input_fasta = "/Users/sarajordan/Desktop/Gobal_Measles_DTA_analysis/Unaligned_unfiltered_global_sequences.fasta"
    output_fasta = "/Users/sarajordan/Desktop/Gobal_Measles_DTA_analysis/Filtered_sequences_with_metadata.fasta"
    
    print("Step 1: Parsing metadata...")
    metadata = parse_metadata(metadata_file)
    print(f"Loaded metadata for {len(metadata)} sequences")
    
    print("\nStep 2: Processing FASTA file...")
    process_fasta(input_fasta, output_fasta, metadata)
    
    print(f"\nOutput written to: {output_fasta}")

if __name__ == "__main__":
    main()
