#!/usr/bin/env python3
"""
GloVe-style Co-occurrence Matrix Generator for BPE Tokens
Based on the Stanford GloVe algorithm for computing token co-occurrence statistics.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict, Counter
import json


class CooccurrenceMatrix:
    def __init__(self, window_size=5, min_count=1):
        self.window_size = window_size
        self.min_count = min_count
        self.token_to_id = {}
        self.id_to_token = {}
        self.token_counts = Counter()
        self.cooccurrence_counts = defaultdict(float)
        self.vocab_size = 0
    
    def build_vocabulary(self, token_sequences):
        """Build vocabulary from token sequences."""
        print(f"Building vocabulary...")
        
        # Count all tokens
        for sequence in token_sequences:
            for token in sequence:
                self.token_counts[token] += 1
        
        # Filter by minimum count
        filtered_tokens = {token: count for token, count in self.token_counts.items() 
                          if count >= self.min_count}
        
        # Create token-to-id mapping
        self.vocab_size = len(filtered_tokens)
        for i, token in enumerate(sorted(filtered_tokens.keys())):
            self.token_to_id[token] = i
            self.id_to_token[i] = token
        
        print(f"Vocabulary size: {self.vocab_size} tokens (min_count >= {self.min_count})")
        return self.vocab_size
    
    def compute_cooccurrence(self, token_sequences):
        """Compute co-occurrence statistics using sliding window."""
        print(f"Computing co-occurrence matrix with window size {self.window_size}...")
        
        total_pairs = 0
        
        for seq_idx, sequence in enumerate(token_sequences):
            if (seq_idx + 1) % 50 == 0:
                print(f"Processed {seq_idx + 1}/{len(token_sequences)} sequences...")
                
            # Convert tokens to IDs, skip unknown tokens
            token_ids = []
            for token in sequence:
                if token in self.token_to_id:
                    token_ids.append(self.token_to_id[token])
            
            # Sliding window co-occurrence
            for i, center_id in enumerate(token_ids):
                # Define window boundaries
                start = max(0, i - self.window_size)
                end = min(len(token_ids), i + self.window_size + 1)
                
                # Count co-occurrences within window
                for j in range(start, end):
                    if i != j:  # Don't count self-occurrence
                        context_id = token_ids[j]
                        distance = abs(i - j)
                        
                        # Weight by inverse distance (closer tokens have higher weight)
                        weight = 1.0 / distance
                        
                        # Store symmetric co-occurrence
                        pair = (min(center_id, context_id), max(center_id, context_id))
                        self.cooccurrence_counts[pair] += weight
                        total_pairs += 1
        
        print(f"Total co-occurrence pairs: {len(self.cooccurrence_counts)}")
        print(f"Total weighted counts: {total_pairs}")
    
    def get_matrix_stats(self):
        """Get statistics about the co-occurrence matrix."""
        print("Computing matrix statistics...")
        
        total_entries = self.vocab_size * self.vocab_size
        non_zero_entries = len(self.cooccurrence_counts) * 2  # symmetric
        sparsity = (1 - non_zero_entries / total_entries) * 100
        
        print(f"Matrix shape: {self.vocab_size}x{self.vocab_size}")
        print(f"Non-zero entries: {non_zero_entries}")
        print(f"Sparsity: {sparsity:.2f}%")
        
        return {
            'shape': (self.vocab_size, self.vocab_size),
            'non_zero': non_zero_entries,
            'sparsity': sparsity
        }
    
    def save_vocabulary(self, output_file):
        """Save vocabulary mapping to JSON file."""
        vocab_data = {
            'token_to_id': self.token_to_id,
            'id_to_token': self.id_to_token,
            'token_counts': dict(self.token_counts),
            'vocab_size': self.vocab_size,
            'window_size': self.window_size,
            'min_count': self.min_count
        }
        
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(vocab_data, f, indent=2, ensure_ascii=False)
        
        print(f"Vocabulary saved to: {output_file}")
    
    def save_cooccurrence_csv(self, output_file, top_n=None):
        """Save co-occurrence matrix as CSV (sparse format)."""
        print(f"Saving co-occurrence matrix to CSV...")
        
        # Sort by count (descending)
        sorted_pairs = sorted(self.cooccurrence_counts.items(), 
                             key=lambda x: x[1], reverse=True)
        
        if top_n:
            sorted_pairs = sorted_pairs[:top_n]
            print(f"Saving top {top_n} co-occurrence pairs...")
        
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['token1', 'token2', 'token1_id', 'token2_id', 'cooccurrence_count']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            
            for (i, j), count in sorted_pairs:
                writer.writerow({
                    'token1': self.id_to_token[i],
                    'token2': self.id_to_token[j],
                    'token1_id': i,
                    'token2_id': j,
                    'cooccurrence_count': round(count, 4)
                })
        
        print(f"Co-occurrence matrix saved to: {output_file}")
        return len(sorted_pairs)
    
    def save_dense_matrix_csv(self, output_file, max_tokens=50):
        """Save dense co-occurrence matrix for top tokens (for visualization)."""
        print(f"Saving dense matrix for top {max_tokens} tokens...")
        
        # Get top tokens by frequency
        top_tokens = sorted(self.token_counts.items(), key=lambda x: x[1], reverse=True)[:max_tokens]
        top_token_ids = [self.token_to_id[token] for token, _ in top_tokens if token in self.token_to_id]
        
        # Create dense submatrix (list of lists instead of numpy)
        matrix = [[0.0 for _ in range(len(top_token_ids))] for _ in range(len(top_token_ids))]
        
        for i, id1 in enumerate(top_token_ids):
            for j, id2 in enumerate(top_token_ids):
                pair = (min(id1, id2), max(id1, id2))
                if pair in self.cooccurrence_counts:
                    matrix[i][j] = self.cooccurrence_counts[pair]
        
        # Save as CSV
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            
            # Header with token names
            header = ['token'] + [self.id_to_token[id_] for id_ in top_token_ids]
            writer.writerow(header)
            
            # Matrix rows
            for i, id1 in enumerate(top_token_ids):
                row = [self.id_to_token[id1]] + matrix[i]
                writer.writerow(row)
        
        print(f"Dense matrix saved to: {output_file}")
        return len(top_token_ids)


def load_bpe_results(csv_file):
    """Load BPE results from CSV file."""
    print(f"Loading BPE results from: {csv_file}")
    
    token_sequences = []
    
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            bpe_encoded = row['bpe_encoded']
            tokens = bpe_encoded.split()
            token_sequences.append(tokens)
    
    print(f"Loaded {len(token_sequences)} sequences")
    return token_sequences


def main():
    parser = argparse.ArgumentParser(
        description='Generate GloVe-style co-occurrence matrix from BPE tokens'
    )
    parser.add_argument('bpe_file', help='Input BPE results CSV file')
    parser.add_argument('--window-size', '-w', type=int, default=5,
                       help='Context window size (default: 5)')
    parser.add_argument('--min-count', '-c', type=int, default=1,
                       help='Minimum token count threshold (default: 1)')
    parser.add_argument('--output-dir', '-o', default='output',
                       help='Output directory (default: output)')
    parser.add_argument('--top-pairs', '-t', type=int, default=10000,
                       help='Save top N co-occurrence pairs (default: 10000)')
    parser.add_argument('--dense-size', '-d', type=int, default=50,
                       help='Size of dense matrix for top tokens (default: 50)')
    
    args = parser.parse_args()
    
    # Load BPE token sequences
    token_sequences = load_bpe_results(args.bpe_file)
    
    if not token_sequences:
        print("No token sequences found!")
        sys.exit(1)
    
    # Initialize co-occurrence matrix
    cooc_matrix = CooccurrenceMatrix(
        window_size=args.window_size,
        min_count=args.min_count
    )
    
    # Build vocabulary
    vocab_size = cooc_matrix.build_vocabulary(token_sequences)
    
    if vocab_size == 0:
        print("No tokens meet minimum count threshold!")
        sys.exit(1)
    
    # Compute co-occurrence statistics
    cooc_matrix.compute_cooccurrence(token_sequences)
    
    # Get matrix statistics
    stats = cooc_matrix.get_matrix_stats()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Save results
    base_name = os.path.splitext(os.path.basename(args.bpe_file))[0]
    
    # Save vocabulary
    vocab_file = os.path.join(args.output_dir, f"{base_name}_vocabulary.json")
    cooc_matrix.save_vocabulary(vocab_file)
    
    # Save sparse co-occurrence matrix
    cooc_file = os.path.join(args.output_dir, f"{base_name}_cooccurrence.csv")
    num_pairs = cooc_matrix.save_cooccurrence_csv(cooc_file, top_n=args.top_pairs)
    
    # Save dense matrix for visualization
    dense_file = os.path.join(args.output_dir, f"{base_name}_dense_matrix.csv")
    dense_size = cooc_matrix.save_dense_matrix_csv(dense_file, max_tokens=args.dense_size)
    
    # Print summary
    print(f"\nCo-occurrence Analysis Complete!")
    print("=" * 50)
    print(f"Input sequences: {len(token_sequences)}")
    print(f"Vocabulary size: {vocab_size}")
    print(f"Window size: {args.window_size}")
    print(f"Matrix shape: {stats['shape'][0]}x{stats['shape'][1]}")
    print(f"Matrix sparsity: {stats['sparsity']:.2f}%")
    print(f"Co-occurrence pairs saved: {num_pairs}")
    print(f"Dense matrix size: {dense_size}x{dense_size}")
    print(f"\nOutput files:")
    print(f"  - Vocabulary: {vocab_file}")
    print(f"  - Co-occurrence: {cooc_file}")
    print(f"  - Dense matrix: {dense_file}")


if __name__ == "__main__":
    main()