from Bio import SeqIO
import os

# Dosya yolları
positive_file = "/Users/batuhancifer/Desktop/positive_proteins.fasta"
negative_file = "/Users/batuhancifer/Desktop/negative_proteins.fasta"
full_sequences_file = "/Users/batuhancifer/Downloads/output_proteins.fasta"
output_positive_file = "/Users/batuhancifer/Desktop/matched_positive_sequences.fasta"
output_negative_file = "/Users/batuhancifer/Desktop/matched_negative_sequences.fasta"

# Tüm sekansları yükle
def load_sequences(file_path):
    return list(SeqIO.parse(file_path, "fasta"))

# Kısmı sekansları tam sekanslarla kıyasla ve eşleşenleri bul
def find_matches(partial_seqs, full_seqs):
    matches = []
    for partial in partial_seqs:
        for full in full_seqs:
            if str(partial.seq).replace("\n", "") in str(full.seq):
                matches.append((partial, full))
    return matches

# Sekansları yükle
positive_seqs = load_sequences(positive_file)
negative_seqs = load_sequences(negative_file)
full_seqs = load_sequences(full_sequences_file)

# Pozitif ve negatif eşleşmeleri bul
positive_matches = find_matches(positive_seqs, full_seqs)
negative_matches = find_matches(negative_seqs, full_seqs)

# Pozitif eşleşmeleri yeni bir dosyaya yaz
with open(output_positive_file, "w") as output:
    for partial, full in positive_matches:
        header = f">{partial.id}_matches_{full.id} | {full.description}"
        output.write(f"{header}\n{partial.seq}\n")

# Negatif eşleşmeleri yeni bir dosyaya yaz
with open(output_negative_file, "w") as output:
    for partial, full in negative_matches:
        header = f">{partial.id}_matches_{full.id} | {full.description}"
        output.write(f"{header}\n{partial.seq}\n")

print(f"Pozitif eşleşmeler {output_positive_file} dosyasına kaydedildi.")
print(f"Negatif eşleşmeler {output_negative_file} dosyasına kaydedildi.")