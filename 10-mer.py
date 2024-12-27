#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 14:35:56 2024

@author: batuhancifer
"""

from Bio import SeqIO

def extract_unique_nmers(sequence, n=10):
    """
    Verilen bir sekans için benzersiz n'li gruplar (n-mers) oluşturur.

    Args:
        sequence (str): Protein dizisi.
        n (int): N-mer uzunluğu (default 10).

    Returns:
        set: Benzersiz n-mer gruplarının kümesi.
    """
    nmers = {sequence[i:i + n] for i in range(len(sequence) - n + 1)}
    return nmers

def process_fasta_file(input_file, n=10):
    """
    Bir FASTA dosyasını okuyarak benzersiz n'li grupları çıkarır.

    Args:
        input_file (str): FASTA dosyasının yolu.
        n (int): N-mer uzunluğu (default 10).

    Returns:
        set: Benzersiz n-mer gruplarının kümesi.
    """
    unique_nmers = set()
    for record in SeqIO.parse(input_file, "fasta"):
        unique_nmers.update(extract_unique_nmers(str(record.seq), n))
    return unique_nmers

def write_fasta(unique_motifs, output_file):
    """
    Benzersiz motifleri FASTA formatında sadece motif olarak yazar.

    Args:
        unique_motifs (set): Benzersiz motiflerin kümesi.
        output_file (str): Çıkış için FASTA dosyası yolu.
    """
    with open(output_file, "w") as fasta_file:
        for motif in unique_motifs:
            fasta_file.write(f"{motif}\n")

def main():
    # Girdi dosyaları
    positive_file = "/Users/batuhancifer/Desktop/positive_proteins.fasta"
    negative_file = "/Users/batuhancifer/Desktop/negative_proteins.fasta"

    # Çıkış dosyaları (isimlere '10mer' ekleniyor)
    positive_output_file = "/Users/batuhancifer/Desktop/unique_positive_10mer_motifs.fasta"
    negative_output_file = "/Users/batuhancifer/Desktop/unique_negative_10mer_motifs.fasta"

    # Pozitif protein dosyasını işleme
    print("Pozitif protein dosyası işleniyor...")
    positive_unique_motifs = process_fasta_file(positive_file)
    write_fasta(positive_unique_motifs, positive_output_file)

    # Negatif protein dosyasını işleme
    print("Negatif protein dosyası işleniyor...")
    negative_unique_motifs = process_fasta_file(negative_file)
    write_fasta(negative_unique_motifs, negative_output_file)

    print("Motifler başarıyla işlendi ve dosyalar oluşturuldu.")

if __name__ == "__main__":
    main() 