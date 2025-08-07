#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import glob

def get_all_targets_from_data_folder(data_folder="data"):
    """
    Get all target filenames from .pretty files in the data folder.
    
    Args:
        data_folder (str): Path to folder containing .pretty files
        
    Returns:
        set: Set of all target filenames
    """
    all_targets = set()
    
    # Find all .pretty files
    pretty_files = glob.glob(os.path.join(data_folder, "*.pretty"))
    
    if not pretty_files:
        print(f"Warning: No .pretty files found in '{data_folder}'")
        return all_targets
    
    print(f"Scanning {len(pretty_files)} .pretty files to find all targets...")
    
    for pretty_file in pretty_files:
        # Extract target name from the pretty filename
        base_name = os.path.basename(pretty_file)
        
        # Handle different file naming patterns
        if base_name == 'Emilania_huxleyi.output.pretty':
            target_name = 'Emihu1_scaffold.fa'
        elif base_name.endswith('.mainGenome.fasta.output.output.pretty'):
            # Special case: Prymnesium_platychrysis_sp_CCMP1217.mainGenome.fasta.output.output.pretty
            target_name = base_name.replace('.output.output.pretty', '')
        elif base_name.endswith('.mainGenome.output.output.pretty'):
            # Special case: Chrysochromulina_leadbeateri_var_UIO393.mainGenome.output.output.pretty  
            target_name = base_name.replace('.mainGenome.output.output.pretty', '.mainGenome.fasta')
        elif base_name.endswith('.output.output.pretty'):
            target_name = base_name.replace('.output.output.pretty', '.mainGenome.fasta')
        elif base_name.endswith('.output.pretty'):
            if base_name == 'Calcidiscus_leptoporus_var_RCC1130.output.pretty':
                target_name = 'Calcidiscus_leptoporus_var_RCC1130.mainGenome.fasta'
            else:
                target_name = base_name.replace('.output.pretty', '.mainGenome.fasta')
        else:
            target_name = base_name + '.mainGenome.fasta'
        
        all_targets.add(target_name)
    
    return all_targets

def generate_hit_matrix(input_file, output_file="hit_matrix.tsv", output_format="tsv", data_folder="data"):
    """
    Generate a matrix counting hits of each Query ID against each Target.
    
    Args:
        input_file (str): Path to parsed results TSV file
        output_file (str): Path to output matrix file
        output_format (str): Output format ('tsv' or 'csv')
        data_folder (str): Path to data folder containing .pretty files
    """
    
    try:
        # Read the input TSV file
        print(f"Reading data from {input_file}...")
        df = pd.read_csv(input_file, sep='\t')
        
        # Check if required columns exist
        required_columns = ['Query ID', 'Target']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")
        
        print(f"Found {len(df)} total hits")
        print(f"Unique Query IDs: {df['Query ID'].nunique()}")
        print(f"Unique Targets in hits: {df['Target'].nunique()}")
        
        # Get all possible targets from data folder
        all_targets = get_all_targets_from_data_folder(data_folder)
        print(f"Total targets in data folder: {len(all_targets)}")
        
        # Identify targets with zero hits
        targets_with_hits = set(df['Target'].unique())
        targets_with_zero_hits = all_targets - targets_with_hits
        if targets_with_zero_hits:
            print(f"Targets with zero hits: {len(targets_with_zero_hits)}")
            for target in sorted(targets_with_zero_hits):
                print(f"  - {target}")
        
        # Create a pivot table (matrix) counting hits
        print("Generating hit matrix...")
        hit_matrix = df.pivot_table(
            index='Query ID', 
            columns='Target', 
            aggfunc='size',  # Count occurrences
            fill_value=0     # Fill missing combinations with 0
        )
        
        # Add columns for targets with zero hits
        all_targets_sorted = sorted(all_targets)
        hit_matrix = hit_matrix.reindex(columns=all_targets_sorted, fill_value=0)
        
        # Clean up target names by removing common suffixes
        def clean_target_name(target_name):
            """Remove common suffixes from target names for cleaner display."""
            if target_name.endswith('.mainGenome.fasta'):
                return target_name.replace('.mainGenome.fasta', '')
            elif target_name.endswith('.fa'):
                return target_name.replace('.fa', '')
            return target_name
        
        # Rename columns to cleaner names
        clean_column_names = {col: clean_target_name(col) for col in hit_matrix.columns}
        hit_matrix = hit_matrix.rename(columns=clean_column_names)
        
        # Sort the matrix by Query ID
        hit_matrix = hit_matrix.sort_index()
        
        # Write the matrix to file
        if output_format.lower() == 'csv':
            hit_matrix.to_csv(output_file)
        else:
            hit_matrix.to_csv(output_file, sep='\t')
        
        print(f"Matrix saved to {output_file}")
        print(f"Matrix dimensions: {hit_matrix.shape[0]} Query IDs × {hit_matrix.shape[1]} Targets")
        print(f"Targets with hits: {len(targets_with_hits)}")
        print(f"Targets with zero hits: {len(targets_with_zero_hits)}")
        
        # Print summary statistics
        total_hits = hit_matrix.sum().sum()
        non_zero_entries = (hit_matrix > 0).sum().sum()
        print(f"Total hits in matrix: {total_hits}")
        print(f"Non-zero entries: {non_zero_entries}")
        print(f"Matrix density: {non_zero_entries / (hit_matrix.shape[0] * hit_matrix.shape[1]):.2%}")
        
        return hit_matrix
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return None
    except pd.errors.EmptyDataError:
        print(f"Error: Input file '{input_file}' is empty.")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None

def print_matrix_summary(matrix, top_n=10):
    """
    Print a summary of the matrix including top Query IDs and Targets.
    
    Args:
        matrix (pd.DataFrame): The hit matrix
        top_n (int): Number of top items to display
    """
    
    if matrix is None:
        return
    
    print(f"\n=== Matrix Summary ===")
    
    # Top Query IDs by total hits
    query_totals = matrix.sum(axis=1).sort_values(ascending=False)
    print(f"\nTop {top_n} Query IDs by total hits:")
    for i, (query_id, total) in enumerate(query_totals.head(top_n).items(), 1):
        print(f"  {i:2d}. {query_id}: {total} hits")
    
    # Top Targets by total hits
    target_totals = matrix.sum(axis=0).sort_values(ascending=False)
    print(f"\nTop {top_n} Targets by total hits:")
    for i, (target, total) in enumerate(target_totals.head(top_n).items(), 1):
        print(f"  {i:2d}. {target}: {total} hits")

def main():
    """
    Main function to handle command line arguments and run the matrix generation.
    """
    parser = argparse.ArgumentParser(
        description='Generate a hit matrix from parsed pretty file results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
examples:
  python generate_hit_matrix.py parsed_results.tsv
  python generate_hit_matrix.py -i parsed_results.tsv -o my_matrix.tsv
  python generate_hit_matrix.py -i results.tsv -o matrix.csv -f csv
  python generate_hit_matrix.py -i results.tsv --summary 5
        '''
    )
    
    parser.add_argument(
        'input_file',
        nargs='?',
        default='parsed_results.tsv',
        help='Path to input TSV file (default: parsed_results.tsv)'
    )
    
    parser.add_argument(
        '-i', '--input',
        dest='input_file_alt',
        help='Path to input TSV file (alternative way to specify input)'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='hit_matrix.tsv',
        help='Path to output matrix file (default: hit_matrix.tsv)'
    )
    
    parser.add_argument(
        '-f', '--format',
        choices=['tsv', 'csv'],
        default='tsv',
        help='Output format: tsv or csv (default: tsv)'
    )
    
    parser.add_argument(
        '-s', '--summary',
        type=int,
        default=10,
        help='Number of top items to show in summary (default: 10)'
    )
    
    parser.add_argument(
        '--no-summary',
        action='store_true',
        help='Skip printing matrix summary'
    )
    
    parser.add_argument(
        '-d', '--data-folder',
        default='data',
        help='Path to folder containing .pretty files (default: data)'
    )
    
    args = parser.parse_args()
    
    # Determine input file (prioritize alternative flag if provided)
    input_file = args.input_file_alt if args.input_file_alt else args.input_file
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        return 1
    
    try:
        # Generate the matrix
        matrix = generate_hit_matrix(input_file, args.output, args.format, args.data_folder)
        
        if matrix is not None:
            # Print summary unless disabled
            if not args.no_summary:
                print_matrix_summary(matrix, args.summary)
            
            return 0
        else:
            return 1
            
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    exit(main())