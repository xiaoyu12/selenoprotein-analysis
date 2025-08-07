#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import numpy as np

def simplify_query_name(query_name):
    """
    Simplify query names by extracting protein information.
    
    Args:
        query_name (str): Full query identifier
        
    Returns:
        str: Simplified query name
    """
    # Remove the lcl|UniRef50_ prefix for cleaner display
    if query_name.startswith('lcl|UniRef50_'):
        return query_name.replace('lcl|UniRef50_', '')
    return query_name

def plot_hits_per_query(input_file, output_file=None, plot_type="bar", figsize=(15, 8), 
                       top_n=None, rotate_labels=True, simplify_names=True, plot_metric="hits"):
    """
    Plot the total number of hits or unique genomes per query from the hit matrix.
    
    Args:
        input_file (str): Path to hit matrix TSV file
        output_file (str): Path to save the plot (optional)
        plot_type (str): Type of plot ('bar', 'horizontal', 'pie', 'combined')
        figsize (tuple): Figure size as (width, height)
        top_n (int): Show only top N queries by hits/genomes
        rotate_labels (bool): Rotate x-axis labels
        simplify_names (bool): Simplify query names
        plot_metric (str): What to plot - 'hits', 'genomes', or 'combined'
    """
    
    try:
        # Read the hit matrix
        print(f"Reading hit matrix from {input_file}...")
        df = pd.read_csv(input_file, sep='\t', index_col=0)
        
        # Calculate both metrics
        hits_per_query = df.sum(axis=1)
        genomes_per_query = (df > 0).sum(axis=1)
        
        if plot_metric == "genomes":
            # Use unique genomes per query
            metric_per_query = genomes_per_query.sort_values(ascending=False)
            metric_name = "unique genomes"
            y_label = "Number of Unique Genomes Hit"
            title_suffix = "Unique Genomes"
        elif plot_metric == "combined":
            # For combined plots, sort by hits
            metric_per_query = hits_per_query.sort_values(ascending=False)
            metric_name = "hits and genomes"
            y_label = "Count"
            title_suffix = "Hits and Unique Genomes"
        else:
            # Use total hits per query
            metric_per_query = hits_per_query.sort_values(ascending=False)
            metric_name = "hits"
            y_label = "Total Number of Hits"
            title_suffix = "Hits"
        
        print(f"Found {len(metric_per_query)} queries with total of {metric_per_query.sum()} {metric_name}")
        
        # Filter to top N if specified
        if top_n is not None:
            top_queries = metric_per_query.head(top_n).index
            metric_per_query = metric_per_query.loc[top_queries]
            hits_per_query = hits_per_query.loc[top_queries]
            genomes_per_query = genomes_per_query.loc[top_queries]
            print(f"Showing top {len(metric_per_query)} queries")
        
        # Simplify query names if requested
        if simplify_names:
            query_names = [simplify_query_name(name) for name in metric_per_query.index]
        else:
            query_names = metric_per_query.index
        
        # Create the plot
        if plot_metric == "combined" or plot_type == "combined":
            # Create combined subplot for both metrics
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(figsize[0], figsize[1]*1.2))
            
            # Plot 1: Total hits
            hits_data = hits_per_query.loc[metric_per_query.index]
            bars1 = ax1.bar(range(len(hits_data)), hits_data.values, 
                           color=plt.cm.viridis(np.linspace(0, 1, len(hits_data))))
            ax1.set_ylabel('Total Number of Hits')
            ax1.set_title('Total Hits per Query')
            ax1.set_xticks(range(len(hits_data)))
            ax1.set_xticklabels(query_names, rotation=45 if rotate_labels else 0, ha='right')
            
            # Add value labels on bars
            for i, (bar, value) in enumerate(zip(bars1, hits_data.values)):
                ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                        str(value), ha='center', va='bottom', fontsize=8)
            
            # Plot 2: Unique genomes
            genomes_data = genomes_per_query.loc[metric_per_query.index]
            bars2 = ax2.bar(range(len(genomes_data)), genomes_data.values, 
                           color=plt.cm.plasma(np.linspace(0, 1, len(genomes_data))))
            ax2.set_xlabel('Query ID')
            ax2.set_ylabel('Number of Unique Genomes Hit')
            ax2.set_title('Unique Genomes Hit per Query')
            ax2.set_xticks(range(len(genomes_data)))
            ax2.set_xticklabels(query_names, rotation=45 if rotate_labels else 0, ha='right')
            
            # Add value labels on bars
            for i, (bar, value) in enumerate(zip(bars2, genomes_data.values)):
                ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                        str(value), ha='center', va='bottom', fontsize=8)
            
        else:
            plt.figure(figsize=figsize)
            
            if plot_type == "bar":
                bars = plt.bar(range(len(metric_per_query)), metric_per_query.values, 
                              color=plt.cm.plasma(range(len(metric_per_query))))
                plt.xlabel('Query ID')
                plt.ylabel(y_label)
                plt.title(f'Total Number of {title_suffix} per Query')
                plt.xticks(range(len(metric_per_query)), query_names, 
                          rotation=45 if rotate_labels else 0, ha='right')
                
                # Add value labels on bars
                for i, (bar, value) in enumerate(zip(bars, metric_per_query.values)):
                    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                            str(value), ha='center', va='bottom', fontsize=8)
            
            elif plot_type == "horizontal":
                bars = plt.barh(range(len(metric_per_query)), metric_per_query.values,
                               color=plt.cm.plasma(range(len(metric_per_query))))
                plt.ylabel('Query ID')
                plt.xlabel(y_label)
                plt.title(f'Total Number of {title_suffix} per Query')
                plt.yticks(range(len(metric_per_query)), query_names)
                
                # Add value labels on bars
                for i, (bar, value) in enumerate(zip(bars, metric_per_query.values)):
                    plt.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                            str(value), ha='left', va='center', fontsize=8)
            
            elif plot_type == "pie":
                # For pie chart, show only non-zero values and limit to reasonable number
                non_zero_metric = metric_per_query[metric_per_query > 0]
                if len(non_zero_metric) > 20:
                    print(f"Too many queries ({len(non_zero_metric)}) for pie chart. Using bar chart instead.")
                    plot_type = "bar"
                    return plot_hits_per_query(input_file, output_file, "bar", figsize, 
                                              min(20, top_n) if top_n else 20, rotate_labels, simplify_names, plot_metric)
                else:
                    plt.pie(non_zero_metric.values, labels=query_names[:len(non_zero_metric)], 
                           autopct='%1.1f%%', startangle=90)
                    plt.title(f'Distribution of {title_suffix} per Query')
        
        plt.tight_layout()
        
        # Save or show the plot
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_file}")
        else:
            plt.show()
        
        # Print summary statistics
        if plot_metric == "combined":
            print(f"\n=== Combined Summary Statistics ===")
            print(f"Total hits across all queries: {hits_per_query.sum()}")
            print(f"Total unique query-genome pairs: {genomes_per_query.sum()}")
            print(f"Number of queries: {len(metric_per_query)}")
            print(f"Average hits per query: {hits_per_query.mean():.1f}")
            print(f"Average unique genomes per query: {genomes_per_query.mean():.1f}")
            
            print(f"\nTop 10 queries by total hits:")
            for i, (query, count) in enumerate(hits_per_query.head(10).items(), 1):
                simple_name = simplify_query_name(query) if simplify_names else query
                genomes_count = genomes_per_query.loc[query]
                print(f"  {i:2d}. {simple_name}: {count} hits, {genomes_count} genomes")
        else:
            print(f"\n=== Summary Statistics ===")
            print(f"Total {metric_name} across all queries: {metric_per_query.sum()}")
            print(f"Number of queries with {metric_name}: {(metric_per_query > 0).sum()}")
            print(f"Average {metric_name} per query: {metric_per_query.mean():.1f}")
            print(f"Median {metric_name} per query: {metric_per_query.median():.1f}")
            print(f"Max {metric_name} for a single query: {metric_per_query.max()}")
            print(f"Min {metric_name} for a single query: {metric_per_query.min()}")
            
            # Top queries
            print(f"\nTop 10 queries by {metric_name}:")
            for i, (query, count) in enumerate(metric_per_query.head(10).items(), 1):
                simple_name = simplify_query_name(query) if simplify_names else query
                print(f"  {i:2d}. {simple_name}: {count} {metric_name}")
        
        return metric_per_query
        
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
        description='Plot total number of hits or unique genomes per query from hit matrix',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
examples:
  python plot_hits_per_query.py hit_matrix.tsv
  python plot_hits_per_query.py -i matrix.tsv -o plot.png --top 20
  python plot_hits_per_query.py -t horizontal --size 12 10
  python plot_hits_per_query.py --metric genomes -o unique_genomes.png
  python plot_hits_per_query.py --metric combined -o combined_plot.png --top 15
  python plot_hits_per_query.py --top 15 --no-simplify
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
        choices=['bar', 'horizontal', 'pie', 'combined'],
        default='bar',
        help='Type of plot: bar, horizontal, pie, or combined (default: bar)'
    )
    
    parser.add_argument(
        '--size',
        type=float,
        nargs=2,
        default=[15, 8],
        help='Figure size as width height (default: 15 8)'
    )
    
    parser.add_argument(
        '--top',
        type=int,
        help='Show only top N queries by hit/genome count'
    )
    
    parser.add_argument(
        '--no-rotate',
        action='store_true',
        help='Do not rotate x-axis labels (for bar plots)'
    )
    
    parser.add_argument(
        '--no-simplify',
        action='store_true',
        help='Do not simplify query names'
    )
    
    parser.add_argument(
        '--metric',
        choices=['hits', 'genomes', 'combined'],
        default='combined',
        help='What to plot: hits (total hit count), genomes (unique genome count), or combined (both) (default: combined)'
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
        result = plot_hits_per_query(
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