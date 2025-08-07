#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def simplify_query_name(query_name):
    """
    Simplify query names by extracting protein information.
    
    Args:
        query_name (str): Full query identifier
        
    Returns:
        str: Simplified query name
    """
    if query_name.startswith('lcl|UniRef50_'):
        return query_name.replace('lcl|UniRef50_', '')
    return query_name

def create_bubble_plot(input_file, output_file=None, figsize=None, 
                      max_bubbles=1000, simplify_names=True, min_hits=1):
    """
    Create a bubble plot showing hits between queries and genomes.
    
    Args:
        input_file (str): Path to hit matrix TSV file
        output_file (str): Path to save the plot (optional)
        figsize (tuple): Figure size as (width, height)
        max_bubbles (int): Maximum number of bubbles to show (for performance)
        simplify_names (bool): Simplify query names
        min_hits (int): Minimum number of hits to show a bubble
    """
    
    try:
        # Read the hit matrix
        print(f"Reading hit matrix from {input_file}...")
        df = pd.read_csv(input_file, sep='\t', index_col=0)
        
        print(f"Matrix dimensions: {df.shape[0]} queries × {df.shape[1]} genomes")
        
        # Filter out zero hits and apply minimum threshold
        nonzero_hits = df[df >= min_hits].stack()
        nonzero_hits = nonzero_hits.dropna()
        
        print(f"Found {len(nonzero_hits)} query-genome pairs with ≥{min_hits} hits")
        
        # If too many bubbles, sample or filter to top hits
        if len(nonzero_hits) > max_bubbles:
            print(f"Too many bubbles ({len(nonzero_hits)}), showing top {max_bubbles} by hit count")
            nonzero_hits = nonzero_hits.nlargest(max_bubbles)
        
        # Create coordinate arrays
        queries = []
        genomes = []
        hit_counts = []
        
        for (query, genome), hits in nonzero_hits.items():
            queries.append(query)
            genomes.append(genome)
            hit_counts.append(hits)
        
        # Get unique queries and genomes for positioning
        unique_queries = df.index.tolist()
        unique_genomes = df.columns.tolist()
        
        # Create position mappings
        query_positions = {query: i for i, query in enumerate(unique_queries)}
        genome_positions = {genome: i for i, genome in enumerate(unique_genomes)}
        
        # Convert to coordinates
        x_coords = [genome_positions[genome] for genome in genomes]
        y_coords = [query_positions[query] for query in queries]
        
        # Normalize bubble sizes (radius proportional to sqrt of hits for better visual scaling)
        max_hits = max(hit_counts)
        min_hits_val = min(hit_counts)
        bubble_sizes = [(np.sqrt(hits) / np.sqrt(max_hits)) * 200 + 10 for hits in hit_counts]
        
        # Create colors for different queries
        n_queries = len(unique_queries)
        colors = plt.cm.Set3(np.linspace(0, 1, min(n_queries, 12)))  # Use Set3 colormap
        if n_queries > 12:
            # For many queries, cycle through colors
            colors = plt.cm.tab20(np.linspace(0, 1, 20))
        
        query_colors = [colors[query_positions[query] % len(colors)] for query in queries]
        
        # Calculate dynamic figure size based on number of queries and genomes
        if figsize is None:
            # Base size with scaling for number of queries/genomes
            width = max(20, len(unique_genomes) * 0.5)
            height = max(15, len(unique_queries) * 0.3)
            # Cap maximum size to prevent extremely large plots
            width = min(width, 50)
            height = min(height, 40)
            figsize = (width, height)
        
        # Create the plot
        plt.figure(figsize=figsize)
        
        # Create scatter plot
        plt.scatter(x_coords, y_coords, s=bubble_sizes, c=query_colors, 
                   alpha=0.6, edgecolors='black', linewidth=0.5)
        
        # Set labels and title
        plt.xlabel('Genomes', fontsize=12)
        plt.ylabel('Queries', fontsize=12)
        plt.title(f'Query-Genome Hit Matrix Bubble Plot\n({len(nonzero_hits)} query-genome pairs)', fontsize=14)
        
        # Simplify query names if requested
        if simplify_names:
            display_queries = [simplify_query_name(q) for q in unique_queries]
        else:
            display_queries = unique_queries
        
        # Set ticks and labels - always show all query names
        plt.yticks(range(len(unique_queries)), display_queries, fontsize=max(6, min(12, 200 // len(unique_queries))))
        
        # For genomes, show all if reasonable, otherwise subset
        if len(unique_genomes) <= 100:
            plt.xticks(range(len(unique_genomes)), unique_genomes, 
                      rotation=45, ha='right', fontsize=max(6, min(10, 150 // len(unique_genomes))))
        else:
            genome_step = max(1, len(unique_genomes) // 50)
            plt.xticks(range(0, len(unique_genomes), genome_step), 
                      [unique_genomes[i] for i in range(0, len(unique_genomes), genome_step)], 
                      rotation=45, ha='right', fontsize=8)
        
        # Add grid
        plt.grid(True, alpha=0.3)
        
        # Create comprehensive legend for bubble sizes
        if max_hits > min_hits_val:
            # Create 5 representative sizes
            legend_sizes = [
                min_hits_val,
                min_hits_val + (max_hits - min_hits_val) // 4,
                min_hits_val + (max_hits - min_hits_val) // 2,
                min_hits_val + 3 * (max_hits - min_hits_val) // 4,
                max_hits
            ]
        else:
            legend_sizes = [max_hits]
        
        legend_bubbles = []
        legend_labels = []
        
        for size in legend_sizes:
            bubble_size = (np.sqrt(size) / np.sqrt(max_hits)) * 200 + 10
            legend_bubbles.append(plt.scatter([], [], s=bubble_size, 
                                            c='gray', alpha=0.6, edgecolors='black', linewidth=0.5))
            legend_labels.append(f'{size} hit{"s" if size != 1 else ""}')
        
        # Add size legend with better positioning
        size_legend = plt.legend(legend_bubbles, legend_labels, 
                               title="Circle Size = Hit Count", loc='upper left', 
                               bbox_to_anchor=(1.02, 1), fontsize=10, title_fontsize=11)
        plt.gca().add_artist(size_legend)
        
        # Add color legend (for top queries only to avoid clutter)
        if len(unique_queries) <= 20:
            color_handles = []
            color_labels = []
            for i, query in enumerate(unique_queries[:min(10, len(unique_queries))]):
                color_handles.append(plt.scatter([], [], c=[colors[i % len(colors)]], 
                                               s=50, alpha=0.6, edgecolors='black', linewidth=0.5))
                display_name = simplify_query_name(query) if simplify_names else query
                color_labels.append(display_name)
            
            plt.legend(color_handles, color_labels, title="Queries (sample)", 
                      loc='upper left', bbox_to_anchor=(1.02, 0.7), fontsize=8)
        
        plt.tight_layout()
        
        # Save or show the plot
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_file}")
        else:
            plt.show()
        
        # Print summary statistics
        print(f"\n=== Bubble Plot Summary ===")
        print(f"Total bubbles shown: {len(nonzero_hits)}")
        print(f"Hit count range: {min_hits_val} - {max_hits}")
        print(f"Matrix coverage: {len(nonzero_hits)}/{df.shape[0] * df.shape[1]} ({100*len(nonzero_hits)/(df.shape[0] * df.shape[1]):.2f}%) non-zero entries")
        
        # Top query-genome pairs
        print(f"\nTop 10 query-genome pairs by hits:")
        top_pairs = nonzero_hits.nlargest(10)
        for i, ((query, genome), hits) in enumerate(top_pairs.items(), 1):
            query_display = simplify_query_name(query) if simplify_names else query
            print(f"  {i:2d}. {query_display} × {genome}: {hits} hits")
        
        return nonzero_hits
        
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None

def main():
    """
    Main function to handle command line arguments and create the bubble plot.
    """
    parser = argparse.ArgumentParser(
        description='Create a bubble plot showing hits between queries and genomes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
examples:
  python bubble_plot_hits.py hit_matrix.tsv
  python bubble_plot_hits.py -i matrix.tsv -o bubble_plot.png --max-bubbles 500
  python bubble_plot_hits.py --size 25 20 --min-hits 5
  python bubble_plot_hits.py --no-simplify -o large_plot.pdf
        '''
    )
    
    parser.add_argument(
        'input_file',
        nargs='?',
        default='hit_matrix.tsv',
        help='Path to hit matrix TSV file (default: hit_matrix.tsv)'
    )
    
    parser.add_argument(
        '-i', '--input',
        dest='input_file_alt',
        help='Path to hit matrix TSV file (alternative way to specify)'
    )
    
    parser.add_argument(
        '-o', '--output',
        help='Path to save the plot (PNG, PDF, SVG, etc.)'
    )
    
    parser.add_argument(
        '--size',
        type=float,
        nargs=2,
        default=None,
        help='Figure size as width height (default: auto-calculated based on data size)'
    )
    
    parser.add_argument(
        '--max-bubbles',
        type=int,
        default=1000,
        help='Maximum number of bubbles to show (default: 1000)'
    )
    
    parser.add_argument(
        '--min-hits',
        type=int,
        default=1,
        help='Minimum number of hits to show a bubble (default: 1)'
    )
    
    parser.add_argument(
        '--no-simplify',
        action='store_true',
        help='Do not simplify query names'
    )
    
    args = parser.parse_args()
    
    # Determine input file
    input_file = args.input_file_alt if args.input_file_alt else args.input_file
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        return 1
    
    try:
        # Create the bubble plot
        result = create_bubble_plot(
            input_file=input_file,
            output_file=args.output,
            figsize=tuple(args.size) if args.size else None,
            max_bubbles=args.max_bubbles,
            simplify_names=not args.no_simplify,
            min_hits=args.min_hits
        )
        
        return 0 if result is not None else 1
        
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    exit(main())