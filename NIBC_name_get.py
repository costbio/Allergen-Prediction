#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 15:34:29 2024

@author: batuhancifer
"""
import os
import time
import pandas as pd
from Bio import Entrez, SeqIO

# Entrez ayarları
Entrez.email = "alibatuhancifer@hotmail.com"

# Girdi ve çıktı dosyaları
input_file = "/Users/batuhancifer/Downloads/Browse the Database.csv"
output_fasta = "/Users/batuhancifer/Downloads/output_proteins.fasta"

# Veriyi oku
df = pd.read_csv(input_file)

def fetch_protein_data(accession):
    """
    NCBI'den protein adı, organizması ve sekansını çeker.
    """
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        time.sleep(1)  # Rate limit için bekleme süresi
        
        protein_name = record.description.split(",")[0]
        organism = record.annotations.get("organism", "Unknown Organism")
        sequence = str(record.seq)
        
        return protein_name, organism, sequence
    
    except Exception as e:
        print(f"ERROR: Could not fetch data for {accession}. Reason: {e}")
        return None, None, None

# FASTA formatında dosyaya yaz
with open(output_fasta, "w") as fasta_file:
    for idx, accession in enumerate(df['Accession']):
        print(f"Fetching data for: {accession}")
        protein_name, organism, sequence = fetch_protein_data(accession)
        
        if protein_name and organism and sequence:
            header = f">{accession}_gi|{idx+1}|_{protein_name}_{organism}_&_id={idx+1:.3f}"
            fasta_file.write(f"{header}\n{sequence}\n")

print(f"Protein verileri '{output_fasta}' dosyasına kaydedildi.")
