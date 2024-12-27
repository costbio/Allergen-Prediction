from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# FASTA formatındaki dosyayı okuma
def read_fasta_file(file_path):
    """
    FASTA formatındaki protein başlıklarını ve dizilimlerini okuyarak bir sözlük döndürür.
    """
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        header = record.id  # FASTA başlığı
        sequence = str(record.seq)  # Sekans
        sequences[header] = sequence
    return sequences

# Proteinleri gruplandırma
def group_proteins(sequences):
    """
    Protein başlıklarına göre gruplandırır. Her gruptaki ilk protein pozitif, geri kalanlar negatif olur.
    """
    grouped_data = {}
    for header, sequence in sequences.items():
        # Başlık "group_index" formatında değilse hatayı önlemek için kontrol eklenir
        if "_" not in header or len(header.split("_")) != 2:
            print(f"Warning: Invalid header format '{header}', skipping...")
            continue
        
        group_id, index = header.split("_")  # Başlık örneği: "001_1"
        if group_id not in grouped_data:
            grouped_data[group_id] = {'positive': sequence, 'negative': []}
        else:
            grouped_data[group_id]['negative'].append(sequence)
    return grouped_data

# Grupları FASTA formatında yazma
def write_to_fasta(grouped_proteins, positive_file, negative_file):
    """
    Gruplandırılmış protein verilerini pozitif ve negatif olarak iki ayrı FASTA dosyasına yazar.
    """
    # Pozitif proteinleri yaz
    positive_records = []
    for group, data in grouped_proteins.items():
        positive_records.append(SeqRecord(seq=data['positive'], id=f"{group}_positive", description=""))
    
    SeqIO.write(positive_records, positive_file, "fasta")

    # Negatif proteinleri yaz
    negative_records = []
    for group, data in grouped_proteins.items():
        for idx, seq in enumerate(data['negative']):
            negative_records.append(SeqRecord(seq=seq, id=f"{group}_negative_{idx+1}", description=""))
    
    SeqIO.write(negative_records, negative_file, "fasta")

# Ana fonksiyon
def main():
    # FASTA dosya yolu
    file_path = "/Users/batuhancifer/Downloads/trainingdataset.fa"
    
    # Proteinleri oku
    protein_sequences = read_fasta_file(file_path)
    
    # Proteinleri gruplandır
    grouped_proteins = group_proteins(protein_sequences)

    # Gruplandırılmış veriyi ekrana yazdır
    for group, data in grouped_proteins.items():
        print(f"Group: {group}")
        print(f"  Positive: {data['positive'][:30]}...")  # Pozitif proteinin ilk 30 karakteri
        print(f"  Negative: {[seq[:30] for seq in data['negative']]}")  # Negatif proteinlerin ilk 30 karakteri

    # FASTA dosyası yolları
    positive_output_file = "/Users/batuhancifer/Desktop/positive_proteins.fasta"
    negative_output_file = "/Users/batuhancifer/Desktop/negative_proteins.fasta"

    # Veriyi FASTA formatında yaz
    write_to_fasta(grouped_proteins, positive_output_file, negative_output_file)

    print(f"Pozitif protein veritabanı oluşturuldu: {positive_output_file}")
    print(f"Negatif protein veritabanı oluşturuldu: {negative_output_file}")

# Kodun çalıştırılması
if __name__ == "__main__":
    main()