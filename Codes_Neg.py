import time
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
###############################################################################


# print(mhc_list)

#/home/sevval/desktop/netMHCpan-4.1

# netMHCpan proje dizini
project_pwd = "/home/sevval/desktop/netMHCpan-4.1"
NETMHCpan_pwd="/home/sevval/desktop/netMHCpan-4.1/Linux_x86_64"
binary_path = os.path.join(project_pwd, "Linux_x86_64", "bin", "netMHCpan")

env = os.environ.copy()
env["PROJECTPWD"] = project_pwd
env["NETMHCpan"] = NETMHCpan_pwd
env["TMPDIR"] = os.path.join(project_pwd, "tmp")
# tmp klasörü varsa oluştur
os.makedirs(env["TMPDIR"], exist_ok=True)
###############################################################################
# Dosya yolunu belirt
file_path = "Data/MHCs.txt"

# Dosyadaki her satırı liste olarak oku ve boşlukları temizle
with open(file_path, "r") as f:
    mhc_list = [line.strip() for line in f if line.strip()]

# Hedef klasörün yolu
folder_path = "Data/PeptitNeg/"

# Klasördeki dosyaları listele (sadece dosyalar)
files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
# print(files)
###############################################################################


###############################################################################
lengths = "8,10"

###############################################################################

###############################################################################
# Tek işi çalıştıran fonksiyon
def run_job(fasta, allele):
    fasta_file = f"./Data/PeptitNeg/{fasta}"
    output_file = f"./Output/PeptitNeg/{fasta}@{allele}.xls"
    cmd = [
        binary_path,
        "-f", fasta_file,
        "-a", allele,
        "-l", lengths,
        "-xls",
        "-xlsfile", output_file
    ]
    try:
        subprocess.run(cmd, check=True, env=env, capture_output=True, text=True)
        return f"[✓] {fasta} + {allele}"
    except subprocess.CalledProcessError as e:
        return f"[✗] {fasta} + {allele}:\n{e.stderr or e.stdout}"

# Parametreleri açacak yardımcı fonksiyon
def run_job_wrapper(args):
    return run_job(*args)



start_time = time.time()


tasks = [(fasta, allele) for fasta in files for allele in mhc_list]

with ProcessPoolExecutor() as executor:
    results = executor.map(run_job_wrapper, tasks)

for r in results:
    print(r)



end_time = time.time()
elapsed = end_time - start_time
print(f"Geçen süre: {elapsed:.2f} saniye")


#Geçen süre: 271.72 saniye





###############################################################################









