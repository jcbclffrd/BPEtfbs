#!/usr/bin/env python3
"""
Standard BPE (Byte Pair Encoding) processor for DNA sequence files.
Implements the classic BPE algorithm as described in the original paper.
"""

import argparse
import sys
import os
from collections import defaultdict, Counter
import re


class StandardBPE:
    def __init__(self, num_merges=1000):
        self.num_merges = num_merges
        self.bpe_codes = []
        self.word_freqs = Counter()
        self.vocab = set()
    
    def get_word_tokens(self, word):
        """Split word into characters with end-of-word marker."""
        # For DNA sequences, split into individual characters
        return list(word) + ['</w>']
    
    def get_stats(self, vocab):
        """Count frequency of consecutive symbol pairs in vocabulary."""
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
    
    def learn_bpe(self, corpus_text):
        """Learn BPE codes from corpus text."""
        print(f"Learning BPE with {self.num_merges} merges...")
        
        # Split corpus into words (sequences)
        words = corpus_text.strip().split()
        
        # Initialize vocabulary with character-level tokens
        vocab = {}
        for word in words:
            word_tokens = ' '.join(self.get_word_tokens(word))
            vocab[word_tokens] = vocab.get(word_tokens, 0) + 1
        
        print(f"Initial vocabulary size: {len(vocab)}")
        print(f"Total words in corpus: {len(words)}")
        
        # Learn BPE merges
        for i in range(self.num_merges):
            pairs = self.get_stats(vocab)
            if not pairs:
                print(f"No more pairs to merge after {i} merges")
                break
            
            best_pair = max(pairs, key=pairs.get)
            vocab = self.merge_vocab(best_pair, vocab)
            self.bpe_codes.append(best_pair)
            
            if (i + 1) % 100 == 0:
                print(f"Completed {i + 1} merges...")
        
        print(f"Learned {len(self.bpe_codes)} BPE codes")
        
        # Store final vocabulary
        self.vocab = set()
        for word in vocab:
            self.vocab.update(word.split())
    
    def apply_bpe(self, word):
        """Apply learned BPE codes to a single word."""
        word_tokens = ' '.join(self.get_word_tokens(word))
        
        # Apply each BPE merge in order
        for pair in self.bpe_codes:
            bigram = re.escape(' '.join(pair))
            p = re.compile(r'(?<!\S)' + bigram + r'(?!\S)')
            word_tokens = p.sub(''.join(pair), word_tokens)
        
        return word_tokens.split()
    
    def encode_corpus(self, corpus_text):
        """Encode entire corpus using learned BPE."""
        words = corpus_text.strip().split()
        encoded_words = []
        
        for word in words:
            encoded_tokens = self.apply_bpe(word)
            encoded_words.extend(encoded_tokens)
        
        return encoded_words
    
    def get_token_frequencies(self, encoded_tokens):
        """Get frequency counts of all tokens."""
        return Counter(encoded_tokens)


def parse_fasta_to_corpus(fasta_file):
    """Parse FASTA file and create corpus with start/end tags."""
    sequences = []
    current_seq = ""
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(current_seq)
                    current_seq = ""
                else:
                    current_seq += line
            
            # Add last sequence
            if current_seq:
                sequences.append(current_seq)
    
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
    
    # Create corpus with start/end tags
    corpus_words = []
    for seq in sequences:
        # Add start and end tags to each sequence
        tagged_seq = f"<s>{seq}</s>"
        corpus_words.append(tagged_seq)
    
    return ' '.join(corpus_words)


def save_corpus(corpus_text, output_file):
    """Save corpus text to file."""
    print(f"Saving corpus to: {output_file}")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(corpus_text)
    
    print(f"Corpus saved successfully!")


def load_corpus(corpus_file):
    """Load corpus text from file."""
    try:
        with open(corpus_file, 'r', encoding='utf-8') as f:
            return f.read()
    except FileNotFoundError:
        print(f"Error: Corpus file '{corpus_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading corpus file: {e}")
        sys.exit(1)


def save_encoded_results(encoded_tokens, output_file):
    """Save encoded results to file."""
    print(f"Saving encoded results to: {output_file}")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(' '.join(encoded_tokens))
    
    print(f"Encoded results saved successfully!")


def save_alternating_case_tokens(encoded_tokens, output_file):
    """Save encoded tokens with alternating case (caps/lower) for visual inspection."""
    print(f"Saving alternating case tokens to: {output_file}")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    alternating_tokens = []
    for i, token in enumerate(encoded_tokens):
        if i % 2 == 0:  # Even indices: uppercase
            alternating_tokens.append(token.upper())
        else:  # Odd indices: lowercase
            alternating_tokens.append(token.lower())
    
    # Concatenate all tokens (no spaces between them)
    concatenated_text = ''.join(alternating_tokens)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(concatenated_text)
    
    print(f"Alternating case tokens saved successfully!")


def save_token_frequencies(token_freqs, output_file):
    """Save token frequencies to file."""
    print(f"Saving token frequencies to: {output_file}")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("token,frequency,count\n")
        total_tokens = sum(token_freqs.values())
        
        for token, count in token_freqs.most_common():
            frequency = count / total_tokens if total_tokens > 0 else 0
            f.write(f'"{token}",{frequency:.6f},{count}\n')
    
    print(f"Token frequencies saved! Found {len(token_freqs)} unique tokens.")


def main():
    parser = argparse.ArgumentParser(description='Standard BPE processor for DNA sequences')
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input-fasta', help='Input FASTA file path')
    input_group.add_argument('--input-corpus', help='Input corpus text file path')
    
    # Output options
    parser.add_argument('--output-corpus', help='Output corpus file path')
    parser.add_argument('--output-encoded', default='output/encoded_results.txt',
                       help='Output encoded results file (default: output/encoded_results.txt)')
    parser.add_argument('--output-tokens', help='Output token frequencies file')
    parser.add_argument('--output-alternating', help='Output tokens with alternating case for visual inspection')
    
    # BPE parameters
    parser.add_argument('--merges', '-m', type=int, default=1000,
                       help='Number of BPE merges to perform (default: 1000)')
    
    # Display options
    parser.add_argument('--show-vocab', action='store_true',
                       help='Show learned vocabulary')
    parser.add_argument('--show-codes', action='store_true',
                       help='Show learned BPE merge codes')
    
    args = parser.parse_args()
    
    # Load or create corpus
    if args.input_fasta:
        print(f"Reading FASTA file: {args.input_fasta}")
        corpus_text = parse_fasta_to_corpus(args.input_fasta)
        
        if args.output_corpus:
            save_corpus(corpus_text, args.output_corpus)
    else:
        print(f"Loading corpus from: {args.input_corpus}")
        corpus_text = load_corpus(args.input_corpus)
    
    print(f"Corpus length: {len(corpus_text)} characters")
    print(f"Corpus words: {len(corpus_text.split())} words")
    
    # Initialize and train BPE
    bpe = StandardBPE(num_merges=args.merges)
    bpe.learn_bpe(corpus_text)
    
    # Show learned codes if requested
    if args.show_codes:
        print("\nLearned BPE merge codes:")
        for i, (a, b) in enumerate(bpe.bpe_codes, 1):
            print(f"{i:4d}: '{a}' + '{b}' -> '{a + b}'")
    
    # Show vocabulary if requested
    if args.show_vocab:
        print(f"\nLearned vocabulary ({len(bpe.vocab)} tokens):")
        for token in sorted(bpe.vocab):
            print(f"  '{token}'")
    
    # Encode corpus
    print("\nEncoding corpus with learned BPE...")
    encoded_tokens = bpe.encode_corpus(corpus_text)
    
    # Save encoded results
    save_encoded_results(encoded_tokens, args.output_encoded)
    
    # Save alternating case tokens if requested
    if args.output_alternating:
        save_alternating_case_tokens(encoded_tokens, args.output_alternating)
    
    # Calculate and show statistics
    original_chars = len(corpus_text.replace(' ', ''))
    encoded_tokens_count = len(encoded_tokens)
    compression_ratio = original_chars / encoded_tokens_count if encoded_tokens_count > 0 else 0
    
    print(f"\nBPE Results Summary:")
    print("=" * 60)
    print(f"Original characters: {original_chars}")
    print(f"Encoded tokens: {encoded_tokens_count}")
    print(f"Compression ratio: {compression_ratio:.2f}x")
    print(f"Vocabulary size: {len(bpe.vocab)}")
    print(f"Learned merges: {len(bpe.bpe_codes)}")
    
    # Save token frequencies if requested
    if args.output_tokens:
        token_freqs = bpe.get_token_frequencies(encoded_tokens)
        save_token_frequencies(token_freqs, args.output_tokens)
        
        print(f"\nTop 10 most frequent tokens:")
        for i, (token, count) in enumerate(token_freqs.most_common(10), 1):
            freq = count / len(encoded_tokens)
            print(f"  {i:2d}. '{token}': {count} occurrences ({freq:.4f})")


if __name__ == "__main__":
    main()