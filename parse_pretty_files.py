#!/usr/bin/env python3

import os
import glob
import re
import argparse

def parse_pretty_files(data_folder="data", output_file="parsed_results.tsv"):
    """
    Parse all *.pretty files in the data folder and extract blast hit information.
    
    Args:
        data_folder (str): Path to folder containing .pretty files
        output_file (str): Path to output TSV file
    """
    
    all_hits = []
    
    pretty_files = glob.glob(os.path.join(data_folder, "*.pretty"))
    
    for file_path in pretty_files:
        print(f"Processing {file_path}...")
        hits = parse_single_pretty_file(file_path)
        all_hits.extend(hits)
    
    write_output(all_hits, output_file)
    print(f"Parsed {len(all_hits)} hits from {len(pretty_files)} files")
    print(f"Results written to {output_file}")

def parse_single_pretty_file(file_path):
    """
    Parse a single .pretty file and extract all blast hits.
    
    Args:
        file_path (str): Path to the .pretty file
        
    Returns:
        list: List of dictionaries containing hit information
    """
    
    hits = []
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Split content by "newSP" markers
    sections = re.split(r'_\|\*\|.*?newSP.*?\|\*\|_', content)
    
    # Skip the first section (before first newSP)
    for section in sections[1:]:
        hit_info = parse_hit_section(section)
        if hit_info:
            hits.append(hit_info)
    
    return hits

def parse_hit_section(section):
    """
    Parse a single hit section to extract required information.
    
    Args:
        section (str): Text section for one blast hit
        
    Returns:
        dict: Dictionary containing extracted information
    """
    
    hit_info = {
        'Blastx evalue': '',
        'Query ID': '',
        'Protein Description': '',
        'Count': '',
        'Taxonomy': '',
        'Tax ID': '',
        'Rep ID': '',
        'Target': '',
        'Chromosome': '',
        'UGA-SECIS': '',
        'Free Energy': ''
    }
    
    lines = section.split('\n')
    
    for line in lines:
        line = line.strip()
        
        # Extract Blastx evalue
        if line.startswith('Blastx evalue:'):
            hit_info['Blastx evalue'] = line.split(':', 1)[1].strip()
        
        # Extract Query name and split into components
        elif line.startswith('Query name:'):
            query_name = line.split(':', 1)[1].strip()
            parse_query_name(query_name, hit_info)
        
        # Extract Target
        elif line.startswith('Target:'):
            target_path = line.split(':', 1)[1].strip()
            # Extract just the filename
            hit_info['Target'] = os.path.basename(target_path)
        
        # Extract Chromosome
        elif line.startswith('Chromosome:'):
            chrom_info = line.split(':', 1)[1].strip()
            hit_info['Chromosome'] = chrom_info
        
        # Extract UGA-SECIS distance
        elif 'UGA-SECIS:' in line:
            match = re.search(r'UGA-SECIS:\s*(\d+)', line)
            if match:
                hit_info['UGA-SECIS'] = match.group(1)
        
        # Extract Free Energy
        elif 'Free Energy =' in line:
            match = re.search(r'Free Energy\s*=\s*([-\d.]+)', line)
            if match:
                hit_info['Free Energy'] = match.group(1)
    
    # Only return hit if we found at least the basic info
    if hit_info['Blastx evalue'] and hit_info['Query ID']:
        return hit_info
    
    return None

def parse_query_name(query_name, hit_info):
    """
    Parse the query name string into separate components.
    
    Example: "lcl|UniRef50_Q50KB1 Protein disulfide-isomerase-like protein EhSep2 n=18 Tax=Eukaryota TaxID=2759 RepID=SEP2_EMIHU"
    
    Args:
        query_name (str): The full query name string
        hit_info (dict): Dictionary to store parsed components
    """
    
    # Split by spaces but handle special cases
    parts = query_name.split()
    
    if len(parts) == 0:
        return
    
    # First part is always the Query ID
    hit_info['Query ID'] = parts[0]
    
    # Find positions of special markers
    n_index = -1
    tax_index = -1
    taxid_index = -1
    repid_index = -1
    
    for i, part in enumerate(parts):
        if part.startswith('n='):
            n_index = i
        elif part.startswith('Tax='):
            tax_index = i
        elif part.startswith('TaxID='):
            taxid_index = i
        elif part.startswith('RepID='):
            repid_index = i
    
    # Extract each component
    if n_index > 0:
        # Protein description is everything between Query ID and n=
        hit_info['Protein Description'] = ' '.join(parts[1:n_index])
        hit_info['Count'] = parts[n_index]
    elif len(parts) > 1:
        # If no n= found, put everything after Query ID as protein description
        hit_info['Protein Description'] = ' '.join(parts[1:])
    
    if tax_index >= 0:
        hit_info['Taxonomy'] = parts[tax_index]
    
    if taxid_index >= 0:
        hit_info['Tax ID'] = parts[taxid_index]
    
    if repid_index >= 0:
        hit_info['Rep ID'] = parts[repid_index]

def write_output(hits, output_file):
    """
    Write parsed hits to output TSV file.
    
    Args:
        hits (list): List of hit dictionaries
        output_file (str): Path to output file
    """
    
    headers = ['Blastx evalue', 'Query ID', 'Protein Description', 'Count', 'Taxonomy', 'Tax ID', 'Rep ID', 'Target', 'Chromosome', 'UGA-SECIS', 'Free Energy']
    
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(headers) + '\n')
        
        # Write data
        for hit in hits:
            row = [hit.get(header, '') for header in headers]
            f.write('\t'.join(row) + '\n')

def main():
    """
    Main function to handle command line arguments and run the parser.
    """
    parser = argparse.ArgumentParser(
        description='Parse *.pretty files and extract blast hit information',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
examples:
  python parse_pretty_files.py
  python parse_pretty_files.py -d /path/to/data -o results.tsv
  python parse_pretty_files.py --data-folder data --output-file my_results.tsv
        '''
    )
    
    parser.add_argument(
        '-d', '--data-folder',
        default='data',
        help='Path to folder containing *.pretty files (default: data)'
    )
    
    parser.add_argument(
        '-o', '--output-file',
        default='parsed_results.tsv',
        help='Path to output TSV file (default: parsed_results.tsv)'
    )
    
    args = parser.parse_args()
    
    # Check if data folder exists
    if not os.path.exists(args.data_folder):
        print(f"Error: Data folder '{args.data_folder}' does not exist.")
        return 1
    
    # Check if data folder contains any .pretty files
    pretty_files = glob.glob(os.path.join(args.data_folder, "*.pretty"))
    if not pretty_files:
        print(f"Error: No *.pretty files found in '{args.data_folder}'.")
        return 1
    
    print(f"Found {len(pretty_files)} *.pretty files in '{args.data_folder}'")
    
    try:
        parse_pretty_files(args.data_folder, args.output_file)
        return 0
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    exit(main())