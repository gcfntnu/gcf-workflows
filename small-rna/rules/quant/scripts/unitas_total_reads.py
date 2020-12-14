import sys
import glob

def parse_fasta(iterator):
    if not isinstance(iterator, str):
        iterator = '\n'.join(iterator)
    sequence_blocks = [x.strip() for x in iterator.split('>') if x.strip()]
    name = None
    sequences = {}
    for sequence_block in sequence_blocks:
        sequence_block = sequence_block.splitlines()
        name = sequence_block[0].strip()
        sequence = [x.strip() for x in sequence_block[1:] if x.strip()]
        sequence = ''.join(sequence)
        sequences[name] = sequence
    return sequences

def load_fasta(pathname):
    sequences = {}
    for filename in glob.glob(pathname):
        with open(filename, 'rb') as f:
            sequences.update(parse_fasta(f))
    return sequences

def count_sequences(sequences):
    c = collections.defaultdict(int)
    for num, seq in sequences.items():
        c[seq] += int(num)
    return counts

