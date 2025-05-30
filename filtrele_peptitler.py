import pandas as pd
import os
import sys

# --- Ayarlar ---
klasor_yolu = 'PeptitPos' # İçinde Excel/TSV dosyalarınızın olduğu klasör
rank_sutunu_ana_htarif = 'El_Rank' 
aranacak_uzantilar = ('.xls', '.xlsx')
# Yeni Ayar: Tüm sonuçların toplanacağı tek Excel dosyasının adı
birlesik_cikti_dosyasi_adi = 'Tum_Filtrelenmis_Peptidler.xlsx' 
# Bu dosya, script'in çalıştığı dizine kaydedilecektir.

# --- İşleme Fonksiyonu ---
def siniflandir_peptid(rank_value):
    try:
        if isinstance(rank_value, str):
            rank_value = rank_value.replace(',', '.')
        rank_f = float(rank_value)
        
        if rank_f < 0.500:
            return "Strong binding peptides"
        elif 0.500 <= rank_f <= 2.000:
            return "Weak binding peptides"
        else:
            return None
    except (ValueError, TypeError):
        return None

# --- Ana İşlem ---
print(f"'{klasor_yolu}' klasöründeki dosyalar işleniyor...")

if not os.path.isdir(klasor_yolu):
    print(f"HATA: '{klasor_yolu}' adında bir klasör bulunamadı!")
    sys.exit(1)

try:
    dosya_adlari = [f for f in os.listdir(klasor_yolu) if f.endswith(aranacak_uzantilar) and not f.startswith('~') and not f.endswith('_filtered.xlsx')]
except Exception as e:
    print(f"HATA: '{klasor_yolu}' klasörü okunurken bir sorun oluştu: {e}")
    sys.exit(1)

if not dosya_adlari:
    print(f"HATA: '{klasor_yolu}' klasöründe '{aranacak_uzantilar}' uzantılı dosya bulunamadı.")
    sys.exit(1)

toplam_dosya = len(dosya_adlari)
print(f"Toplam {toplam_dosya} dosya bulundu. İşlem başlatılıyor...")
print("-" * 40)

dosya_sayaci = 0
basarili_dosya_sayaci = 0 # Hatasız işlenen dosya sayısı
hatali_dosya_sayaci = 0   # İşlenirken hata alınan dosya sayısı
bulunan_peptid_dosya_sayaci = 0 # İçinde uygun peptid bulunan dosya sayısı

tum_filtrelenmis_peptid_listesi = [] # Tüm dosyalardan gelen filtrelenmis verileri toplamak için

for dosya_adi in dosya_adlari:
    dosya_sayaci += 1
    tam_dosya_yolu = os.path.join(klasor_yolu, dosya_adi)
    # Bireysel çıktı dosyası oluşturma kısmı kaldırıldı.
    # cikti_dosya_adi = os.path.splitext(dosya_adi)[0] + '_filtered.xlsx'
    # cikti_dosya_yolu = os.path.join(klasor_yolu, cikti_dosya_adi)

    print(f"({dosya_sayaci}/{toplam_dosya}) '{dosya_adi}' işleniyor...", end="")
    df = None 

    try:
        file_extension = os.path.splitext(tam_dosya_yolu)[1].lower()
        print(f" (Uzantı: {file_extension})", end="")

        if file_extension == '.xls':
            try:
                df = pd.read_csv(tam_dosya_yolu, sep='\s+', skiprows=[0, 2], engine='python', dtype=str, na_filter=False)
                print(" [TSV (skiprows=[0,2], sep=\\s+) olarak okundu]", end="")
            except Exception as e_tsv:
                print(f" [TSV HATA: {type(e_tsv).__name__}: {e_tsv}]", end="")
        elif file_extension == '.xlsx':
            try:
                df = pd.read_excel(tam_dosya_yolu, engine='openpyxl')
                print(" [openpyxl ile okundu]", end="")
            except Exception as e_xlsx:
                 print(f" [openpyxl HATA: {type(e_xlsx).__name__}]", end="")
        else:
            try:
                df = pd.read_excel(tam_dosya_yolu, engine=None)
                print(" [engine=None bilinmeyen uzantı]", end="")
            except Exception as e_unknown:
                print(f" [engine=None bilinmeyen uzantı HATA: {type(e_unknown).__name__}]", end="")
        
        if df is None: 
            raise ValueError("Dosya okuma başarısız (df None kaldı).")
        
        if df.empty:
            print(" [Uyarı: Dosya okundu ancak boş geldi.]")
            basarili_dosya_sayaci += 1 # Boş da olsa hatasız işlendi.
            continue # Bu dosyadan peptid eklenmeyecek, sonraki dosyaya geç.

        # --- Sütun Adı Bulma Mantığı ---
        gecerli_rank_sutunu = None
        df_columns_str = [str(col).strip() for col in df.columns] 
        df.columns = df_columns_str 

        priority_rank_names_lower = [rank_sutunu_ana_htarif.lower(), 'തിരഞ്ഞെടുക്കൽ_rank', '%rank', 'rank_el']
        for col_str in df_columns_str:
            if col_str.lower() in priority_rank_names_lower:
                gecerli_rank_sutunu = col_str
                print(f" [Sütun: '{gecerli_rank_sutunu}' (öncelikli)]", end="")
                break
        
        if gecerli_rank_sutunu is None:
            secondary_keywords = ['rank', 'score']
            excluded_general_terms = ['core', 'icore']
            for col_str in df_columns_str:
                col_lower = col_str.lower()
                is_excluded = any(term in col_lower for term in excluded_general_terms)
                if not is_excluded:
                    if any(keyword in col_lower for keyword in secondary_keywords):
                        gecerli_rank_sutunu = col_str
                        print(f" [Sütun: '{gecerli_rank_sutunu}' (ikincil 'rank'/'score')]", end="")
                        break
        
        if gecerli_rank_sutunu is None:
            general_alternatives_lower = ['core', 'icore', 'prediction score', 
                                          'binding affinity (nm)', 'affinity(nm)', 'nm', 'affinity',
                                          'presentazione rank']
            for col_str in df_columns_str:
                if col_str.lower() in general_alternatives_lower:
                    gecerli_rank_sutunu = col_str
                    print(f" [Sütun: '{gecerli_rank_sutunu}' (genel alt.)]", end="")
                    break
            
            if gecerli_rank_sutunu is None: 
                 for col_str in df_columns_str:
                    if col_str.startswith('HLA-') and ':' in col_str:
                        gecerli_rank_sutunu = col_str
                        print(f" [Sütun: '{gecerli_rank_sutunu}' (HLA paterni)]", end="")
                        break
                            
        if gecerli_rank_sutunu is None:
             raise ValueError(f"Rank sütunu bulunamadı. Aranan: '{rank_sutunu_ana_htarif}'. Dosyadaki sütunlar: {df_columns_str}")
        
        # --- Sınıflandırma ---
        try:
            if gecerli_rank_sutunu not in df.columns:
                 raise ValueError(f"Belirlenen rank sütunu '{gecerli_rank_sutunu}' dataframe'de bulunamadı.")
            df['Binding_Strength'] = df[gecerli_rank_sutunu].apply(siniflandir_peptid)
        except Exception as e_apply:
            problematic_values = []
            if gecerli_rank_sutunu in df: 
                 problematic_values = df[gecerli_rank_sutunu].unique()[:5]
            raise ValueError(f"'{gecerli_rank_sutunu}' sütununa sınıflandırma uygulanırken hata: {e_apply}. Sütundaki değerler: {problematic_values}")

        filtrelenmis_df = df.dropna(subset=['Binding_Strength']).copy() # .copy() ile SettingWithCopyWarning önlemi

        if not filtrelenmis_df.empty:
            filtrelenmis_df.loc[:, 'KaynakDosya'] = dosya_adi # Kaynak dosya adını yeni bir sütun olarak ekle
            tum_filtrelenmis_peptid_listesi.append(filtrelenmis_df)
            print(f" [OK - {len(filtrelenmis_df)} peptid bulundu ve listeye eklendi]")
            bulunan_peptid_dosya_sayaci += 1
        else:
            print(" [OK - Uygun peptid yok]")
        
        basarili_dosya_sayaci += 1 # Dosya hatasız işlendi

    except Exception as e: 
        error_type = type(e).__name__
        print(f" [DOSYA İŞLEME HATA ({error_type}): {e}]")
        hatali_dosya_sayaci += 1

# --- Döngü sonrası birleştirme ve kaydetme ---
print("-" * 40)
if tum_filtrelenmis_peptid_listesi:
    print(f"{bulunan_peptid_dosya_sayaci} adet dosyadan uygun peptid(ler) bulundu.")
    birlesik_df = pd.concat(tum_filtrelenmis_peptid_listesi, ignore_index=True)
    toplam_bulunan_peptid = len(birlesik_df)
    print(f"Toplam {toplam_bulunan_peptid} uygun peptid tüm dosyalardan toplandı.")

    try:
        # Birleşik dosyayı script'in çalıştığı ana dizine kaydet
        birlesik_df.to_excel(birlesik_cikti_dosyasi_adi, index=False, engine='openpyxl')
        print(f"Tüm filtrelenmiş peptidler '{birlesik_cikti_dosyasi_adi}' dosyasına başarıyla kaydedildi.")
        
        summary_stats = birlesik_df['Binding_Strength'].value_counts()
        print("\nBirleşik Çıktıdaki Peptid Sınıflandırma Özeti:")
        print(summary_stats)
        if 'KaynakDosya' in birlesik_df:
            print(f"\nBu peptidler {birlesik_df['KaynakDosya'].nunique()} farklı kaynak dosyadan gelmektedir.")

    except Exception as e_save:
        print(f"HATA: Birleşik Excel dosyası '{birlesik_cikti_dosyasi_adi}' kaydedilemedi: {e_save}")
else:
    print("Hiçbir dosyada uygun peptid bulunamadığı için birleşik çıktı dosyası oluşturulmadı.")

print("-" * 40)
print("Genel İşlem Özeti:")
print(f"Toplam İncelenen Dosya: {dosya_sayaci}")
print(f"Hatasız İşlenen Dosya Sayısı: {basarili_dosya_sayaci}")
print(f"İçinde Uygun Peptid Bulunan Dosya Sayısı: {bulunan_peptid_dosya_sayaci}")
print(f"İşlenirken Hata Alınan Dosya Sayısı: {hatali_dosya_sayaci}")