# BPE-GloVe DNA Sequence Analysis Project

This project applies Byte Pair Encoding (BPE) and GloVe-style co-occurrence analysis to DNA sequences, enabling the discovery of meaningful patterns and relationships in genomic data.

## Overview

The project consists of two main components:
1. **BPE Processor**: Applies Byte Pair Encoding to DNA sequences to create meaningful subword tokens
2. **GloVe Co-occurrence Matrix Generator**: Computes token co-occurrence statistics similar to the Stanford GloVe algorithm

### Key Results (from `cbe.fa` analysis)
- **Vocabulary Size**: 99 unique BPE tokens discovered
- **Window Size**: 5 (default sliding window for co-occurrence)
- **Co-occurrence Relationships**: 857 unique token pair relationships identified

## Project Structure

### Main Scripts

- **`bpe_processor.py`**: Implements Byte Pair Encoding for DNA sequences
  - Reads FASTA format DNA sequences
  - Learns BPE codes from the input sequences
  - Encodes sequences using learned BPE tokens
  - Outputs compression statistics and encoded sequences

- **`glove_cooccurrence.py`**: Generates co-occurrence matrices for BPE tokens
  - Builds vocabulary from BPE-encoded sequences
  - Computes token co-occurrence statistics using a sliding window
  - Generates sparse and dense matrix representations
  - Calculates marginal frequencies for statistical analysis

- **`session_logger.py`**: Utility for logging processing sessions
  - Tracks command history and outputs
  - Creates timestamped session logs

### Input Data (`iData/`)

Contains various biological sequence and annotation files:
- **`cbe.fa`**: FASTA format DNA sequences
- **`Ebox.txt`**: E-box motif sequences
- **`NEEandREPseq.txt`**: Regulatory element sequences
- **`expre12.tab`**, **`factorexpdts2.tab`**: Expression data tables
- **`factorinfodts.txt`**: Transcription factor information
- **`topbot2Int.txt`**: Interaction data
- Various other sequence and annotation files

### Output Files

#### `output/` Directory
- **`results.csv`**: Primary BPE encoding results
  - Contains: sequence_id, original_sequence, bpe_encoded, original_length, encoded_tokens, compression_ratio
  - Shows how each DNA sequence is tokenized using BPE

- **`results2.csv`**: Additional BPE analysis results

- **`results2_cooccurrence.csv`**: Token co-occurrence statistics
  - Format: token1, token2, cooccurrence_count
  - Shows how often pairs of BPE tokens appear together within the window size

- **`results2_dense_matrix.csv`**: Dense co-occurrence matrix
  - Full matrix representation of token relationships
  - Rows and columns represent tokens, values are co-occurrence counts

- **`results2_marginal_frequencies.csv`**: Marginal frequency statistics
  - Contains frequency information for each token
  - Used for statistical analysis and normalization

- **`results2_vocabulary.json`**: Vocabulary mapping
  - Maps between tokens and their numeric IDs
  - Essential for interpreting the matrix outputs

#### `output2/` Directory
Contains the most recent run results with the same file types as above.

## Usage

### BPE Processing
```bash
python bpe_processor.py -i iData/cbe.fa -o output/results.csv --num-merges 1000
```

### Co-occurrence Analysis
```bash
python glove_cooccurrence.py -i output/results.csv -o output/results2 --window-size 5
```

## Output Interpretation

- **Compression Ratio**: Higher values indicate better tokenization efficiency
- **Co-occurrence Counts**: Higher counts between tokens suggest functional relationships
- **Vocabulary Size**: Indicates the diversity of discovered sequence patterns

## Experiments

### BPE Analysis on cbe.fa

**Input**: `iData/cbe.fa` (314 DNA sequences)

**Command**:
```bash
python bpe_processor.py iData/cbe.fa --merges 1000 --show-tokens --output output/results.csv
```

**Results Summary**:
- **Total unique tokens discovered**: 173
- **Average compression ratio**: 15.31x
- **Output file**: `output/results.csv`

**Top 10 Most Frequent Tokens**:
1. `ATGGCTTTTATATT</w>`: 11 occurrences
2. `GCATTATTTTCCGC</w>`: 8 occurrences  
3. `AAAGAATAACCCAA</w>`: 8 occurrences
4. `GCGGAATTTCCTG</w>`: 7 occurrences
5. `GCGTACTTTCCGC</w>`: 7 occurrences
6. `GGCCGGAAATTCCCC</w>`: 6 occurrences
7. `AGGGAAACCCCAAT</w>`: 5 occurrences
8. `ATGGGAAAACCGGAA</w>`: 5 occurrences
9. `ACGCGGAAATTCCACGGC</w>`: 5 occurrences
10. `ACGCGGAAATTCCGCAGC</w>`: 5 occurrences

## Applications

This analysis can be used for:
- Discovering regulatory motifs in DNA sequences
- Identifying conserved sequence patterns
- Preparing sequence data for machine learning models
- Analyzing relationships between sequence elements