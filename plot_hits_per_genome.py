#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import sys

def simplify_genome_name(genome_name):
    """
    Simplify genome names by extracting species/strain information.
    
    Args:
        genome_name (str): Full genome filename
        
    Returns:
        str: Simplified genome name
    """
    # Remove common suffixes
    name = genome_name.replace('.mainGenome.fasta', '').replace('.fa', '')
    
    # Split by underscores and take first few meaningful parts
    parts = name.split('_')
    
    if len(parts) >= 2:
        # Take genus and species, plus strain/var if present
        if 'var' in parts:
            var_idx = parts.index('var')
            if var_idx + 1 < len(parts):
                return f"{parts[0]} {parts[1]} {parts[var_idx + 1]}"
        elif 'sp' in parts:
            sp_idx = parts.index('sp')
            if sp_idx + 1 < len(parts):
                return f"{parts[0]} {parts[1]} {parts[sp_idx + 1]}"
        else:
            return f"{parts[0]} {parts[1]}"
    
    return name

def plot_hits_per_genome(input_file, output_file=None, plot_type="bar", figsize=(12, 8), 
                        top_n=None, rotate_labels=True, simplify_names=True, plot_metric="hits"):
    """
    Plot the total number of hits or unique queries per genome from the hit matrix.
    
    Args:
        input_file (str): Path to hit matrix TSV file
        output_file (str): Path to save the plot (optional)
        plot_type (str): Type of plot ('bar', 'horizontal', 'pie')
        figsize (tuple): Figure size as (width, height)
        top_n (int): Show only top N genomes by hits/queries
        rotate_labels (bool): Rotate x-axis labels
        simplify_names (bool): Simplify genome names
        plot_metric (str): What to plot - 'hits' or 'queries'
    """
    
    try:
        # Read the hit matrix
        print(f"Reading hit matrix from {input_file}...")
        df = pd.read_csv(input_file, sep='\t', index_col=0)
        
        if plot_metric == "queries":
            # Calculate unique queries per genome (count non-zero entries in each column)
            metric_per_genome = (df > 0).sum(axis=0).sort_values(ascending=False)
            metric_name = "unique queries"
            y_label = "Number of Unique Queries"
            title_suffix = "Unique Queries"
        else:
            # Calculate total hits per genome (sum of each column)
            metric_per_genome = df.sum(axis=0).sort_values(ascending=False)
            metric_name = "hits"
            y_label = "Total Number of Hits"
            title_suffix = "Hits"
        
        print(f"Found {len(metric_per_genome)} genomes with total of {metric_per_genome.sum()} {metric_name}")
        
        # Filter to top N if specified
        if top_n is not None:
            metric_per_genome = metric_per_genome.head(top_n)
            print(f"Showing top {len(metric_per_genome)} genomes")
        
        # Simplify genome names if requested
        if simplify_names:
            genome_names = [simplify_genome_name(name) for name in metric_per_genome.index]
        else:
            genome_names = metric_per_genome.index
        
        # Create the plot
        plt.figure(figsize=figsize)
        
        if plot_type == "bar":
            bars = plt.bar(range(len(metric_per_genome)), metric_per_genome.values, 
                          color=plt.cm.viridis(range(len(metric_per_genome))))
            plt.xlabel('Genome')
            plt.ylabel(y_label)
            plt.title(f'Total Number of {title_suffix} per Genome')
            plt.xticks(range(len(metric_per_genome)), genome_names, 
                      rotation=45 if rotate_labels else 0, ha='right')
            
            # Add value labels on bars
            for i, (bar, value) in enumerate(zip(bars, metric_per_genome.values)):
                plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                        str(value), ha='center', va='bottom', fontsize=8)
        
        elif plot_type == "horizontal":
            bars = plt.barh(range(len(metric_per_genome)), metric_per_genome.values,
                           color=plt.cm.viridis(range(len(metric_per_genome))))
            plt.ylabel('Genome')
            plt.xlabel(y_label)
            plt.title(f'Total Number of {title_suffix} per Genome')
            plt.yticks(range(len(metric_per_genome)), genome_names)
            
            # Add value labels on bars
            for i, (bar, value) in enumerate(zip(bars, metric_per_genome.values)):
                plt.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                        str(value), ha='left', va='center', fontsize=8)
        
        elif plot_type == "pie":
            # For pie chart, show only non-zero values
            non_zero_metric = metric_per_genome[metric_per_genome > 0]
            if top_n is None or len(non_zero_metric) <= 20:
                plt.pie(non_zero_metric.values, labels=genome_names[:len(non_zero_metric)], 
                       autopct='%1.1f%%', startangle=90)
                plt.title(f'Distribution of {title_suffix} per Genome')
            else:
                print("Too many genomes for pie chart. Using bar chart instead.")
                plot_type = "bar"
                return plot_hits_per_genome(input_file, output_file, "bar", figsize, 
                                          top_n, rotate_labels, simplify_names, plot_metric)
        
        plt.tight_layout()
        
        # Save or show the plot
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_file}")
        else:
            plt.show()
        
        # Print summary statistics
        print(f"\n=== Summary Statistics ===")
        print(f"Total {metric_name} across all genomes: {metric_per_genome.sum()}")
        print(f"Number of genomes with {metric_name}: {(metric_per_genome > 0).sum()}")
        print(f"Average {metric_name} per genome: {metric_per_genome.mean():.1f}")
        print(f"Median {metric_name} per genome: {metric_per_genome.median():.1f}")
        print(f"Max {metric_name} in a single genome: {metric_per_genome.max()}")
        print(f"Min {metric_name} in a single genome: {metric_per_genome.min()}")
        
        # Top genomes
        print(f"\nTop 10 genomes by {metric_name}:")
        for i, (genome, count) in enumerate(metric_per_genome.head(10).items(), 1):
            simple_name = simplify_genome_name(genome) if simplify_names else genome
            print(f"  {i:2d}. {simple_name}: {count} {metric_name}")
        
        return metric_per_genome
        
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None

def main():
    """
    Main function to handle command line arguments and run the plotting.
    """
    parser = argparse.ArgumentParser(
        description='Plot total number of hits or unique queries per genome from hit matrix',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
examples:
  python plot_hits_per_genome.py hit_matrix.tsv
  python plot_hits_per_genome.py -i matrix.tsv -o plot.png --top 15
  python plot_hits_per_genome.py -t horizontal --size 10 6
  python plot_hits_per_genome.py -t pie --top 10 --no-simplify
  python plot_hits_per_genome.py --metric queries -o unique_queries.png
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
        '-t', '--type',
        choices=['bar', 'horizontal', 'pie'],
        default='bar',
        help='Type of plot: bar, horizontal, or pie (default: bar)'
    )
    
    parser.add_argument(
        '--size',
        type=float,
        nargs=2,
        default=[12, 8],
        help='Figure size as width height (default: 12 8)'
    )
    
    parser.add_argument(
        '--top',
        type=int,
        help='Show only top N genomes by hit/query count'
    )
    
    parser.add_argument(
        '--no-rotate',
        action='store_true',
        help='Do not rotate x-axis labels (for bar plots)'
    )
    
    parser.add_argument(
        '--no-simplify',
        action='store_true',
        help='Do not simplify genome names'
    )
    
    parser.add_argument(
        '--metric',
        choices=['hits', 'queries'],
        default='hits',
        help='What to plot: hits (total hit count) or queries (unique query count) (default: hits)'
    )
    
    args = parser.parse_args()
    
    # Determine input file
    input_file = args.input_file_alt if args.input_file_alt else args.input_file
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        return 1
    
    try:
        # Create the plot
        result = plot_hits_per_genome(
            input_file=input_file,
            output_file=args.output,
            plot_type=args.type,
            figsize=tuple(args.size),
            top_n=args.top,
            rotate_labels=not args.no_rotate,
            simplify_names=not args.no_simplify,
            plot_metric=args.metric
        )
        
        return 0 if result is not None else 1
        
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    exit(main())