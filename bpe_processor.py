#!/usr/bin/env python3
"""
BPE (Byte Pair Encoding) processor for DNA sequence files.
Takes a FASTA file as input and applies BPE to the sequences.
"""

import argparse
import sys
import os
import csv
from collections import defaultdict, Counter
import re


class BPEProcessor:
    def __init__(self, num_merges=1000):
        self.num_merges = num_merges
        self.bpe_codes = []
        self.vocab = Counter()
    
    def get_stats(self, vocab):
        """Count frequency of consecutive pairs in vocabulary."""
        pairs = defaultdict(int)
        for word, freq in vocab.items():
            symbols = word.split()
            for i in range(len(symbols) - 1):
                pairs[(symbols[i], symbols[i + 1])] += freq
        return pairs
    
    def merge_vocab(self, pair, v_in):
        """Merge the most frequent pair in the vocabulary."""
        v_out = {}
        bigram = re.escape(' '.join(pair))
        p = re.compile(r'(?<!\S)' + bigram + r'(?!\S)')
        for word in v_in:
            new_word = p.sub(''.join(pair), word)
            v_out[new_word] = v_in[word]
        return v_out
    
    def prepare_sequences(self, sequences):
        """Prepare sequences for BPE by treating each sequence as a word in corpus."""
        vocab = Counter()
        for seq in sequences:
            # Each sequence is a "word" with begin/end markers and spaces between chars
            spaced_seq = '<s> ' + ' '.join(seq) + ' </s>'
            vocab[spaced_seq] += 1
        return vocab
    
    def learn_bpe(self, sequences):
        """Learn BPE codes from the input sequences."""
        print(f"Learning BPE with {self.num_merges} merges...")
        
        # Prepare vocabulary
        self.vocab = self.prepare_sequences(sequences)
        print(f"Initial vocabulary size: {len(self.vocab)}")
        
        # Learn BPE merges
        for i in range(self.num_merges):
            pairs = self.get_stats(self.vocab)
            if not pairs:
                break
            
            best_pair = max(pairs, key=pairs.get)
            self.vocab = self.merge_vocab(best_pair, self.vocab)
            self.bpe_codes.append(best_pair)
            
            if (i + 1) % 100 == 0:
                print(f"Completed {i + 1} merges...")
        
        print(f"Learned {len(self.bpe_codes)} BPE codes")
    
    def apply_bpe(self, sequence):
        """Apply learned BPE codes to a sequence."""
        # Start with spaced sequence with begin/end markers
        word = '<s> ' + ' '.join(sequence) + ' </s>'
        
        # Apply each BPE merge in order
        for pair in self.bpe_codes:
            bigram = re.escape(' '.join(pair))
            p = re.compile(r'(?<!\S)' + bigram + r'(?!\S)')
            word = p.sub(''.join(pair), word)
        
        return word.split()
    
    def get_unique_tokens(self):
        """Get all unique tokens from the final vocabulary."""
        tokens = set()
        for word in self.vocab:
            tokens.update(word.split())
        return sorted(tokens)


def parse_fasta(file_path):
    """Parse FASTA file and return sequences with their headers."""
    sequences = []
    headers = []
    current_seq = ""
    current_header = ""
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(current_seq)
                        headers.append(current_header)
                    current_header = line[1:]  # Remove '>'
                    current_seq = ""
                else:
                    current_seq += line
            
            # Add last sequence
            if current_seq:
                sequences.append(current_seq)
                headers.append(current_header)
    
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
    
    return headers, sequences


def save_results_to_csv(headers, sequences, bpe_results, output_file):
    """Save BPE results to CSV file."""
    print(f"\nSaving results to: {output_file}")
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['sequence_id', 'original_sequence', 'bpe_encoded', 'original_length', 'encoded_tokens', 'compression_ratio']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        
        for header, seq, bpe_tokens in zip(headers, sequences, bpe_results):
            compression_ratio = len(seq) / len(bpe_tokens) if len(bpe_tokens) > 0 else 0
            
            writer.writerow({
                'sequence_id': header,
                'original_sequence': seq,
                'bpe_encoded': ' '.join(bpe_tokens),
                'original_length': len(seq),
                'encoded_tokens': len(bpe_tokens),
                'compression_ratio': round(compression_ratio, 2)
            })
    
    print(f"Results saved successfully!")


def save_marginal_frequencies(bpe_results, output_file):
    """Save token frequency analysis (marginal frequencies) to CSV file."""
    print(f"\nGenerating marginal frequencies...")
    
    # Count all token frequencies
    token_counts = Counter()
    total_tokens = 0
    
    for bpe_tokens in bpe_results:
        for token in bpe_tokens:
            token_counts[token] += 1
            total_tokens += 1
    
    # Calculate frequencies and probabilities
    token_stats = []
    for token, count in token_counts.most_common():
        frequency = count / total_tokens if total_tokens > 0 else 0
        token_stats.append({
            'token': token,
            'count': count,
            'frequency': round(frequency, 6),
            'percentage': round(frequency * 100, 2)
        })
    
    # Save to CSV
    print(f"Saving marginal frequencies to: {output_file}")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['token', 'count', 'frequency', 'percentage']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        writer.writerows(token_stats)
    
    print(f"Marginal frequencies saved! Found {len(token_stats)} unique tokens.")
    print(f"Top 5 most frequent tokens:")
    for i, stats in enumerate(token_stats[:5], 1):
        print(f"  {i}. '{stats['token']}': {stats['count']} occurrences ({stats['percentage']}%)")


def main():
    parser = argparse.ArgumentParser(description='Apply BPE to DNA sequences in FASTA format')
    parser.add_argument('input_file', help='Input FASTA file path')
    parser.add_argument('--merges', '-m', type=int, default=1000, 
                       help='Number of BPE merges to perform (default: 1000)')
    parser.add_argument('--show-tokens', '-t', action='store_true',
                       help='Show learned tokens/vocabulary')
    parser.add_argument('--show-codes', '-c', action='store_true',
                       help='Show learned BPE merge codes')
    parser.add_argument('--output', '-o', default='output/results.csv',
                       help='Output CSV file path (default: output/results.csv)')
    parser.add_argument('--marginal-freq', '-f', action='store_true',
                       help='Generate marginal frequencies file')
    
    args = parser.parse_args()
    
    # Parse input file
    print(f"Reading sequences from: {args.input_file}")
    headers, sequences = parse_fasta(args.input_file)
    print(f"Found {len(sequences)} sequences")
    
    if not sequences:
        print("No sequences found in file.")
        sys.exit(1)
    
    # Initialize BPE processor
    bpe = BPEProcessor(num_merges=args.merges)
    
    # Learn BPE from sequences
    bpe.learn_bpe(sequences)
    
    # Show learned codes if requested
    if args.show_codes:
        print("\nLearned BPE merge codes:")
        for i, (a, b) in enumerate(bpe.bpe_codes, 1):
            print(f"{i:4d}: '{a}' + '{b}' -> '{a + b}'")
    
    # Show vocabulary if requested
    if args.show_tokens:
        tokens = bpe.get_unique_tokens()
        print(f"\nLearned vocabulary ({len(tokens)} tokens):")
        for token in tokens:
            print(f"  '{token}'")
    
    # Apply BPE to all sequences
    print(f"\nApplying BPE to sequences...")
    bpe_results = []
    
    for header, seq in zip(headers, sequences):
        bpe_tokens = bpe.apply_bpe(seq)
        bpe_results.append(bpe_tokens)
    
    # Save results to CSV
    save_results_to_csv(headers, sequences, bpe_results, args.output)
    
    # Generate marginal frequencies if requested
    if args.marginal_freq:
        freq_output = args.output.replace('.csv', '_marginal_frequencies.csv')
        save_marginal_frequencies(bpe_results, freq_output)
    
    # Show summary results
    print(f"\nBPE Results Summary:")
    print("=" * 60)
    
    total_original_chars = sum(len(seq) for seq in sequences)
    total_encoded_tokens = sum(len(tokens) for tokens in bpe_results)
    avg_compression = total_original_chars / total_encoded_tokens if total_encoded_tokens > 0 else 0
    
    print(f"Total sequences processed: {len(sequences)}")
    print(f"Total original characters: {total_original_chars}")
    print(f"Total encoded tokens: {total_encoded_tokens}")
    print(f"Average compression ratio: {avg_compression:.2f}x")
    print(f"Results saved to: {args.output}")
    
    # Show first few examples if requested
    if not args.show_tokens and not args.show_codes:
        print(f"\nFirst 5 sequences (use -t or -c for more details):")
        print("-" * 60)
        for i, (header, seq, tokens) in enumerate(zip(headers[:5], sequences[:5], bpe_results[:5])):
            compression = len(seq) / len(tokens) if len(tokens) > 0 else 0
            print(f">{header}")
            print(f"Original: {seq}")
            print(f"BPE:      {' '.join(tokens)}")
            print(f"Tokens:   {len(tokens)} (compression: {compression:.2f}x)")
            print("-" * 40)


if __name__ == "__main__":
    main()