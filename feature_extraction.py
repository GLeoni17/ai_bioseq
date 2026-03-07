from Bio import SeqIO
import numpy as np
import os
from pathlib import Path

def save_kmers_organized(results: list[np.ndarray], sequences: list[str]):
    Path("data").mkdir(exist_ok=True)

    all_data = np.concatenate([r.reshape(-1, 1) for r in results])
    np.savetxt("data/todos_kmers.csv", all_data, delimiter=",", fmt="%.2f")
    
    for i, (kmers, seq) in enumerate(zip(results, sequences)):
        nome_arq = f"data/seq_{i+1:04d}_len{len(seq)}_kmers{len(kmers)}.csv"
        np.savetxt(nome_arq, kmers.reshape(-1, 1), delimiter=",", fmt="%.2f")
    
    with open("data/metadados.txt", "w") as f:        
        for i, (kmers, seq) in enumerate(zip(results, sequences)):
            f.write(f"Seq {i+1:04d}: {len(seq)} bases → {len(kmers)} kmers")


def dna_to_numbers(sequences: list[str]) -> list[np.ndarray]:
    base_map = {
        'A': 0, 'C': 1, 'T': 2, 'G': 3,
        'N': 1.5, 'R': 1.5, 'Y': 1.5, 'S': 1.0, 'W': 2.0,
        'K': 2.5, 'M': 0.5, 'B': 2.0, 'D': 1.67, 'H': 1.0, 'V': 1.33
    }
    
    results = []
    for seq in sequences:
        seq_array = np.array([base_map.get(b.upper(), 0) for b in seq])
        kmers = []
        for i in range(len(seq) - 3):
            kmer = seq_array[i:i+4]
            num = kmer[0]*64 + kmer[1]*16 + kmer[2]*4 + kmer[3]
            kmers.append(num)
        results.append(np.array(kmers))

    save_kmers_organized(results, sequences)

def read_sequences():
    ref_seq = list(SeqIO.parse("data/cv_ref.fasta", "fasta"))
    lst_seqs = list(SeqIO.parse("data/cv_sequences.fasta", "fasta"))

    dna_to_numbers(ref_seq + lst_seqs)

if __name__ == "__main__":
    read_sequences()