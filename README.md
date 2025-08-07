# Selenoprotein Analysis Tools

A comprehensive toolkit for analyzing and visualizing selenoprotein BLAST hit data from genomic searches. This project provides scripts to parse BLAST output files, generate hit matrices, and create various visualizations to understand the distribution of selenoproteins across different genomes and queries.

## Features

- **BLAST Output Parsing**: Parse `.pretty` formatted BLAST output files
- **Hit Matrix Generation**: Create comprehensive hit matrices from parsed data
- **Multiple Visualization Types**:
  - Bar plots for hits per query/genome
  - Combined plots showing both hits and unique genome counts
  - Interactive bubble plots showing query-genome relationships
- **Flexible Output**: Support for various image formats (PNG, PDF, SVG)
- **Customizable Analysis**: Configurable filtering, sorting, and display options

## Installation

### From Source

```bash
git clone <repository-url>
cd selenoprotein-analysis
pip install -e .
```

### Dependencies

- Python >= 3.8
- pandas >= 1.3.0
- matplotlib >= 3.5.0
- seaborn >= 0.11.0
- numpy >= 1.20.0

## Quick Start

### 1. Parse BLAST Output Files

```bash
# Parse all .pretty files in the data/ directory
parse-blast-hits data/ -o parsed_results.tsv

# Or using the script directly
python parse_pretty_files.py data/ -o parsed_results.tsv
```

### 2. Generate Hit Matrix

```bash
# Create a hit matrix from parsed results
generate-hit-matrix parsed_results.tsv -o hit_matrix.tsv

# Or using the script directly
python generate_hit_matrix.py parsed_results.tsv -o hit_matrix.tsv
```

### 3. Create Visualizations

```bash
# Plot hits per query (default: combined view)
plot-hits-per-query hit_matrix.tsv -o query_plot.png

# Plot hits per genome
plot-hits-per-genome hit_matrix.tsv -o genome_plot.png

# Create bubble plot
bubble-plot-hits hit_matrix.tsv -o bubble_plot.png
```

## Scripts Overview

### `parse_pretty_files.py`
Parses BLAST `.pretty` output files and extracts hit information into a structured TSV format.

**Usage:**
```bash
python parse_pretty_files.py [data_folder] -o output.tsv
```

**Features:**
- Processes all `.pretty` files in a directory
- Extracts query, subject, and hit count information
- Handles various BLAST output formats

### `generate_hit_matrix.py`
Converts parsed BLAST results into a query-genome hit matrix.

**Usage:**
```bash
python generate_hit_matrix.py input.tsv -o matrix.tsv [options]
```

**Options:**
- `--min-hits`: Filter entries with minimum hit count
- `--top-queries`: Keep only top N queries by hit count
- `--top-genomes`: Keep only top N genomes by hit count

### `plot_hits_per_query.py`
Creates visualizations showing hit distribution per query.

**Usage:**
```bash
python plot_hits_per_query.py matrix.tsv [options]
```

**Options:**
- `--metric {hits,genomes,combined}`: What to plot (default: combined)
- `--type {bar,horizontal,pie}`: Plot type (default: bar)
- `--top N`: Show only top N queries
- `--size W H`: Figure dimensions

**Plot Types:**
- `hits`: Total hit count per query
- `genomes`: Number of unique genomes hit per query
- `combined`: Both metrics in subplot layout

### `plot_hits_per_genome.py`
Creates visualizations showing hit distribution per genome.

**Usage:**
```bash
python plot_hits_per_genome.py matrix.tsv [options]
```

**Features:**
- Similar options to query plotting
- Focuses on genome-centric analysis

### `bubble_plot_hits.py`
Creates bubble plots showing query-genome relationships.

**Usage:**
```bash
python bubble_plot_hits.py matrix.tsv [options]
```

**Features:**
- **Grid Layout**: Genomes on X-axis, queries on Y-axis
- **Bubble Size**: Proportional to hit count
- **Color Coding**: Different colors for each query
- **Auto-scaling**: Figure size adapts to data dimensions
- **Comprehensive Legend**: Maps bubble size to hit counts

**Options:**
- `--max-bubbles N`: Limit bubbles for performance (default: 1000)
- `--min-hits N`: Minimum hits to show bubble (default: 1)
- `--size W H`: Custom figure dimensions (default: auto-calculated)

## Example Workflow

```bash
# 1. Parse BLAST output files
python parse_pretty_files.py data/ -o results.tsv

# 2. Generate hit matrix with filtering
python generate_hit_matrix.py results.tsv -o matrix.tsv --min-hits 5

# 3. Create comprehensive visualizations
python plot_hits_per_query.py matrix.tsv --metric combined -o query_analysis.png --top 20
python bubble_plot_hits.py matrix.tsv -o interaction_map.png --min-hits 10

# 4. Analyze genome-specific patterns
python plot_hits_per_genome.py matrix.tsv -o genome_analysis.png --top 15
```

## Data Format

### Input: BLAST Pretty Files
The tools expect BLAST output in "pretty" format (human-readable text format).

### Intermediate: Parsed Results TSV
```
Query	Subject	Hits
lcl|UniRef50_Q8WZ42	Emiliana_huxleyi	15
lcl|UniRef50_Q8WZ42	Calcidiscus_leptoporus	8
...
```

### Output: Hit Matrix TSV
```
Query/Genome	Emiliana_huxleyi	Calcidiscus_leptoporus	...
lcl|UniRef50_Q8WZ42	15	8	...
lcl|UniRef50_P12345	0	12	...
...
```

## Configuration

All scripts support:
- **Flexible input/output**: Various file formats and naming schemes
- **Query name simplification**: Removes `lcl|UniRef50_` prefixes for cleaner display
- **Performance optimization**: Handles large datasets efficiently
- **Customizable visualization**: Colors, sizes, labels, and layouts

## Development

### Setup Development Environment

```bash
pip install -e ".[dev]"
```

### Running Tests

```bash
pytest tests/
```

### Code Formatting

```bash
black .
isort .
```

### Type Checking

```bash
mypy .
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Run the test suite and linting tools
6. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this tool in your research, please cite:

```
[Your Citation Here]
```