from Bio import SeqIO
from collections import Counter
import csv

def extract_nmers(sequence, n=9):
    """
    Verilen bir sekans için 9'lu gruplar (n-mers) oluşturur.

    Args:
        sequence (str): Protein dizisi.
        n (int): N-mer uzunluğu (default 9).

    Returns:
        list: Sekanstan çıkarılan n-mer gruplarının listesi.
    """
    nmers = [sequence[i:i + n] for i in range(len(sequence) - n + 1)]
    return nmers

def process_fasta_file(input_file, n=9):
    """
    Bir FASTA dosyasını okuyarak 9'lu grupları çıkarır ve tekrar sayılarını hesaplar.

    Args:
        input_file (str): FASTA dosyasının yolu.
        n (int): N-mer uzunluğu (default 9).

    Returns:
        Counter: N-mer gruplarının ve tekrar sayılarının Counter nesnesi.
    """
    nmers = []
    for record in SeqIO.parse(input_file, "fasta"):
        nmers.extend(extract_nmers(str(record.seq), n))
    return Counter(nmers)

def write_output(sorted_motifs, output_fasta, output_csv):
    """
    Çıkarılan motifleri FASTA ve CSV formatında yazar.

    Args:
        sorted_motifs (list): Motiflerin tekrar sayılarına göre sıralı listesi.
        output_fasta (str): Çıkış için FASTA dosyası yolu.
        output_csv (str): Çıkış için CSV dosyası yolu.
    """
    # FASTA formatına yazma
    with open(output_fasta, "w") as fasta_file:
        for idx, (motif, count) in enumerate(sorted_motifs):
            fasta_file.write(f">{motif}_count_{count}\n{motif}\n")

    # CSV formatına yazma
    with open(output_csv, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Motif", "Count"])  # Başlık satırı
        for motif, count in sorted_motifs:
            writer.writerow([motif, count])

def main():
    # Girdi dosyaları
    positive_file = "/Users/batuhancifer/Desktop/positive_proteins.fasta"
    negative_file = "/Users/batuhancifer/Desktop/negative_proteins.fasta"

    # Çıkış dosyaları
    positive_output_fasta = "/Users/batuhancifer/Desktop/positive_motifs.fasta"
    negative_output_fasta = "/Users/batuhancifer/Desktop/negative_motifs.fasta"
    positive_output_csv = "/Users/batuhancifer/Desktop/positive_motifs.csv"
    negative_output_csv = "/Users/batuhancifer/Desktop/negative_motifs.csv"

    # Pozitif protein dosyasını işleme
    print("Pozitif protein dosyası işleniyor...")
    positive_motifs = process_fasta_file(positive_file)
    sorted_positive_motifs = sorted(positive_motifs.items(), key=lambda x: x[1], reverse=True)

    # Negatif protein dosyasını işleme
    print("Negatif protein dosyası işleniyor...")
    negative_motifs = process_fasta_file(negative_file)
    sorted_negative_motifs = sorted(negative_motifs.items(), key=lambda x: x[1], reverse=True)

    # Çıktıları yazma
    write_output(sorted_positive_motifs, positive_output_fasta, positive_output_csv)
    write_output(sorted_negative_motifs, negative_output_fasta, negative_output_csv)

    print("Motifler başarıyla işlendi ve dosyalar oluşturuldu.")

if __name__ == "__main__":
    main()