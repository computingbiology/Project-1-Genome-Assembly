import numpy as np
import random

nucleotides_dict = {0:"A", 1:"T", 2:"G", 3:"C"}


def generate_random_sequence(N, rng):
    """
        Create a random sequence of given length
    """
    synthetic_sequence = []
    for i in range(N):
        randnum = rng.choice(np.arange(len(nucleotides_dict)))
        nucleotide = nucleotides_dict[randnum]
        synthetic_sequence.append(nucleotide)
    return synthetic_sequence


def write_to_file(k_mers, k, filename, metadata):
    """
        Write given k-mer reads to file
    """
     
    # Shuffle k-mer list before giving to students
    random.shuffle(k_mers)
    
    # Generate random phred scores
    phred_scores = []
    for score in np.random.choice(np.arange(20, 50), size=len(k_mers)*k):
        phred_scores.append(chr(score + 33))

    # Write generated k-mer reads to a file in FASTQ format
    i = 0

    with open(filename, "w") as f:
        for k_mer in k_mers:
            k = len(k_mer)
            # FASTQ format:
            # Metadata
            f.write(metadata + "\n")
            # Read
            f.write(k_mer + "\n")
            # '+' on line
            f.write("+\n")
            # Phred33 scores
            f.write("".join(phred_scores[i:i+k]) + "\n")
            i += k


def generate_sythetic_data(sequence_len, k, seed, filename="TeleTubby.fastq"):
    """
        Generate synthetic data and write to file
    """
    rng = np.random.default_rng(seed)

    # Generate random sequence
    sequence = generate_random_sequence(sequence_len, rng)

    # Get k-mers for given sequence
    k_mers = get_k_mers(sequence, k)
    
    # Write to file (with random scores)
    metadata = "@TeleTubby Genome: Project 1"
    write_to_file(k_mers, k, filename, metadata)
    
    # Return sequence (for reference)
    return "".join(sequence)


def read_fastq(filepath):
    """
        Read and parse fastq file
    """
    sequences = []
    qualities = []
    with open(filepath) as f:
        while True:
            meta_data = f.readline()
            seq = f.readline().rstrip()
            ph = f.readline()
            qual = f.readline().rstrip()
            
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


# Generate k-mers for a given sequence
def get_k_mers(sequence, k_):
    """
        Extract k-mers of given length from sequence
    """
    k_mers = []
    for i in range(len(sequence)-k_+1):
        k_mer = sequence[i:k_+i]
        k_mers.append(''.join(k_mer))
    return k_mers


def viz_graph(n, e):
    """ visualize graph"""
    dot_str = 'digraph "Assembly graph" {\n'
    for src, dst in e:
        dot_str += '{} -> {} ;\n'.format(src, dst)
    return dot_str + '}\n'


def match_score(seq1, seq2):
    """
        Print out match-score between two sequences
    """
    if len(seq1) != len(seq2):
        return 0
    score = 0
    for s1, s2 in zip(seq1, seq2):
        score += (s1 == s2)
    return score / len(seq1)