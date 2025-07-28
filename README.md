# BPE-GloVe DNA Sequence Analysis Project

This project applies Byte Pair Encoding (BPE) and GloVe-style co-occurrence analysis to DNA sequences, enabling the discovery of meaningful patterns and relationships in genomic data.

## Overview

The project consists of two main components:
1. **BPE Processor**: Applies Byte Pair Encoding to DNA sequences to create meaningful subword tokens
2. **GloVe Co-occurrence Matrix Generator**: Computes token co-occurrence statistics similar to the Stanford GloVe algorithm

### Key Results (from `cbe.fa` analysis using `standard_bpe.py`)
- **Input Sequences**: 314 DNA sequences
- **Vocabulary Size**: 173 unique BPE tokens discovered
- **Learned Merges**: 661 BPE merge operations
- **Compression Ratio**: 22.31x (7005 characters → 314 tokens)
- **Window Size**: 5 (default sliding window for co-occurrence)
- **Co-occurrence Relationships**: 857 unique token pair relationships identified

## Project Structure

### Main Scripts

- **`bpe_processor.py`**: Original BPE implementation for DNA sequences
  - Reads FASTA format DNA sequences
  - Learns BPE codes from the input sequences
  - Encodes sequences using learned BPE tokens
  - Outputs compression statistics and encoded sequences

- **`standard_bpe.py`**: Standard BPE implementation with enhanced features
  - Implements the classic BPE algorithm for DNA sequences
  - Character-level tokenization with proper corpus processing
  - Multiple output formats including alternating case visualization
  - Comprehensive command-line interface with flexible I/O options

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

#### `output/` Directory (Standard BPE Results)
- **`corpus.txt`**: Processed corpus with start/end tags for each sequence
  - DNA sequences formatted as: `<s>SEQUENCE</s> <s>SEQUENCE2</s> ...`
  - Used as input for BPE learning process

- **`encoded_results.txt`**: BPE-encoded sequences (space-separated tokens)
  - Shows the tokenized representation of all sequences
  - Each sequence encoded as learned BPE tokens

- **`alternating_case.txt`**: Visual token boundary inspection file
  - Concatenated tokens with alternating UPPER/lower case
  - Helps identify token segmentation patterns and context effects

- **`token_frequencies.csv`**: Complete token frequency analysis
  - Format: token, frequency, count  
  - All 173 unique tokens with their occurrence statistics

#### Legacy Output Files
- **`results.csv`**: Legacy BPE encoding results (from `bpe_processor.py`)
- **`results2_cooccurrence.csv`**: Token co-occurrence statistics  
- **`results2_dense_matrix.csv`**: Dense co-occurrence matrix
- **`results2_marginal_frequencies.csv`**: Marginal frequency statistics
- **`results2_vocabulary.json`**: Vocabulary mapping

#### `output2/` Directory
Contains the most recent run results with the same file types as above.

## Usage

### Standard BPE Processing
```bash
# Basic usage with FASTA input
python standard_bpe.py -i iData/NEEandREPseq.txt --merges 1000

# Auto-optimize merge count for best efficiency
python standard_bpe.py -i iData/cbe.fa --auto-optimize --output-tokens output/optimized_tokens.csv

# Full feature usage with all outputs
python standard_bpe.py -i iData/NEEandREPseq.txt \
  --merges 1000 \
  --output-tokens output/token_frequencies.csv \
  --output-corpus output/corpus.txt \
  --output-encoded output/encoded_results.txt \
  --output-alternating output/alternating_case.txt
```

### Legacy BPE Processing
```bash
python bpe_processor.py -i iData/cbe.fa -o output/results.csv --num-merges 1000
```

### Co-occurrence Analysis
```bash
python glove_cooccurrence.py -i output/results.csv -o output/results2 --window-size 5
```

## Auto-Optimization Feature

The `standard_bpe.py` script includes an **auto-optimization** feature (`--auto-optimize`) that automatically finds the optimal number of BPE merges by testing different values and maximizing the efficiency score (compression ratio ÷ vocabulary size).

### How Auto-Optimization Works
1. **Tests multiple merge counts**: 0, 5, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 500
2. **Calculates efficiency score**: compression_ratio ÷ vocabulary_size for each
3. **Selects optimal**: The merge count with the highest efficiency score
4. **Balances trade-offs**: Finds the sweet spot between compression and vocabulary size

This is crucial because more merges don't always mean better results - too many merges can create an oversized vocabulary that's expensive to store and process.

## Alternating Case Token Visualization

The `standard_bpe.py` script includes a unique **alternating case visualization** feature (`--output-alternating`) that helps identify how BPE tokenization segments sequences in different contexts.

### How It Works
- Tokens are concatenated without spaces
- Even-indexed tokens appear in **UPPERCASE**
- Odd-indexed tokens appear in **lowercase**
- This creates a visual pattern that makes token boundaries immediately apparent

### Why This Is Useful
BPE tokenization can vary depending on context, meaning the same substring might be segmented differently in different sequences. The alternating case visualization helps identify these variations.

#### Example Problem
Consider the sequence `GGGAATTCCC`:
- In the token frequency file, it appears **5 times** as a standalone token
- But grepping the original FASTA file shows it occurs **17 times** total
- The discrepancy occurs because:
  - **7 instances** got merged into longer tokens like `GGGAATTCCCTTCCGC`
  - **5 instances** remained as standalone `GGGAATTCCC` tokens  
  - **5 instances** got segmented differently based on their flanking sequences

#### Using Alternating Case to Debug
```bash
# Generate alternating case output
python standard_bpe.py -i iData/NEEandREPseq.txt --output-alternating debug_tokens.txt

# Search case-insensitively to find all instances
grep -i "gggaattccc" debug_tokens.txt
```

The alternating case makes it easy to see:
- Complete tokens: `GGGAATTCCC` or `gggaattccc`
- Split tokens: `GGGaa` followed by `ttccc` 
- Embedded tokens: `ATGGGAATTCCCttccgc`

This visualization is essential for understanding how BPE's context-dependent tokenization affects sequence analysis.

## Output Interpretation

- **Compression Ratio**: Higher values indicate better tokenization efficiency
- **Co-occurrence Counts**: Higher counts between tokens suggest functional relationships
- **Vocabulary Size**: Indicates the diversity of discovered sequence patterns

## Experiments

### Standard BPE Analysis on cbe.fa

**Input**: `iData/cbe.fa` (314 DNA sequences)

**Command**:
```bash
python standard_bpe.py -i iData/cbe.fa --merges 1000 \
  --output-tokens output/token_frequencies.csv \
  --output-corpus output/corpus.txt \
  --output-encoded output/encoded_results.txt \
  --output-alternating output/alternating_case.txt
```

**Results Summary**:
- **Input**: 314 DNA sequences, 7005 total characters
- **BPE Learning**: 661 merge operations learned (stopped early)  
- **Final vocabulary**: 173 unique tokens
- **Compression**: 22.31x (7005 characters → 314 tokens)
- **Tokenization**: Each sequence becomes one complete token, but sequences vary in length from short (13 bases) to longer (19+ bases)

**Top 10 Most Frequent Tokens**:
1. `<s>ATGGCTTTTATATT</s></w>`: 11 occurrences (3.50%)
2. `<s>GCATTATTTTCCGC</s></w>`: 8 occurrences (2.55%)  
3. `<s>AAAGAATAACCCAA</s></w>`: 8 occurrences (2.55%)
4. `<s>GCGGAATTTCCTG</s></w>`: 7 occurrences (2.23%)
5. `<s>GCGTACTTTCCGC</s></w>`: 7 occurrences (2.23%)
6. `<s>GGCCGGAAATTCCCC</s></w>`: 6 occurrences (1.91%)
7. `<s>AGGGAAACCCCAAT</s></w>`: 5 occurrences (1.59%)
8. `<s>ATGGGAAAACCGGAA</s></w>`: 5 occurrences (1.59%)
9. `<s>ACGCGGAAATTCCACGGC</s></w>`: 5 occurrences (1.59%)
10. `<s>ACGCGGAAATTCCGCAGC</s></w>`: 5 occurrences (1.59%)

### Standard BPE Analysis on NEEandREPseq.txt

**Input**: `iData/NEEandREPseq.txt` (120 longer DNA sequences)

**Results Summary** (More realistic BPE behavior):
- **Input**: 120 sequences, 86,705 total characters  
- **BPE Learning**: 1000 merge operations completed
- **Final vocabulary**: 929 unique tokens
- **Compression**: 4.74x (86,705 characters → 18,303 tokens)
- **Tokenization**: Meaningful subword patterns like `TGC`, `ATT`, `AGC`

This demonstrates how BPE performs better with longer, more diverse sequences.

### Auto-Optimized BPE Analysis on cbe.fa

**Input**: `iData/cbe.fa` (314 DNA sequences)

**Command**:
```bash
python standard_bpe.py -i iData/cbe.fa --auto-optimize --output-tokens output/optimized_tokens.csv
```

**Auto-Optimization Process**:
The algorithm tested merge counts from 0 to 500 and found the optimal efficiency score at **10 merges**.

**Optimized Results Summary**:
- **Optimal merges (k)**: 10 
- **Learned vocabulary size**: 12 tokens
- **Compression ratio**: 1.79x (7005 characters → 3917 tokens)
- **Efficiency score**: 0.1490 (compression ÷ vocabulary size)
- **Generated files**: `output/optimized_tokens.csv`, `output/encoded_results.txt`

**Top 5 Most Frequent Tokens**:
1. `G`: 510 occurrences (13.02%)
2. `A`: 487 occurrences (12.43%)  
3. `C`: 484 occurrences (12.36%)
4. `AA`: 395 occurrences (10.08%)
5. `T`: 390 occurrences (9.96%)

**Key Insight**: The optimal solution uses only 12 tokens (individual DNA bases + common pairs like AA, CC, GG, TT + markup tokens) rather than 173 unique sequence tokens. This provides a much more efficient vocabulary while maintaining reasonable compression.

## Applications

This analysis can be used for:
- Discovering regulatory motifs in DNA sequences
- Identifying conserved sequence patterns
- Preparing sequence data for machine learning models
- Analyzing relationships between sequence elements