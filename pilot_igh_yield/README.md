# Pilot: Singhania 2019 fare tam kan RNA-seq verisinde IGH verimi

Yeni yaklaşımın (fare tam kan RNA-seq → non-productive IGH → IGoR) fizibilitesini
**ölçmek** için yapılmış küçük pilot. Değerlendirme ve literatür raporu:
[`docs/Yeni_Yaklasim_Fizibilite_ve_Literatur.md`](../docs/Yeni_Yaklasim_Fizibilite_ve_Literatur.md).

Bu klasör eski VDJ projesinin (SARP/KIARVA) parçası değildir; ondan bağımsız çalışır.

## Ne yapıldı

| Adım | Betik | Girdi → çıktı |
|---|---|---|
| FASTQ'nun ilk ~300 MiB'ını indir (aralık isteği, kaldığı yerden devam eder) | `fetch_slice.sh` | ENA run → `R1/R2.partial.gz` (~7–9 M okuma çifti) |
| Kesik gzip'i eşit sayıda tam çifte çevir | `slice_pair.py` | `.partial.gz` ×2 → FASTQ ×2 |
| TRUST4 (fare IMGT referansı) | `run_slice.sh` | FASTQ → `*_airr.tsv`, `*_assign.out`, … |
| IGH contig / non-productive özeti | `summarize_trust4.py` | TRUST4 çıktıları → tablo |
| GEO HTSeq sayımlarından Ig ifadesi | `ig_expression_from_counts.py` | `.xlsx` → örnek başına CPM |
| Kontrol ↔ hastalık Ig ifadesi, pilotla kalibrasyon | `ig_expression_by_model.py` | → `tables/ig_expression_by_model.tsv` |
| Tüm kohorta projeksiyon | `project_cohort_yield.py` | → `tables/cohort_projection.tsv` |
| AIRR Data Commons'ta fare IGH envanteri | `adc_mouse_inventory.py`, `adc_mouse_summarize.py` | ADC → fare çalışmaları (6) |
| Greiff 2017 (PRJEB18631) çerçeve-dışı verimi | `adc_greiff_oof.py` (parçalı çalışır), `adc_greiff_summarize.py` | ADC → `tables/greiff2017_adc_yield.tsv` |

Veri: GSE119856 (BioProject PRJNA490475), 503 çift-uçlu HiSeq 4000 koşusu; bunların 244'ü tam kan
(122 kontrol / 122 hastalık). Pilot 10 tam kan kütüphanesi kullandı: Malaria, Toxoplasma ve Influenza
için kontrol ve hastalık.

## Nasıl çalıştırılır

```bash
git clone https://github.com/liulab-dfci/TRUST4 && make -C TRUST4      # bir kez
export TRUST4_DIR=$PWD/TRUST4

# ENA üstverisi (okuma sayıları ve örnek başlıkları)
curl -g -o runs.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA490475&result=read_run&fields=run_accession,read_count,sample_title&format=tsv&limit=0"

pilot_igh_yield/fetch_slice.sh SRR7821542 slices 300
pilot_igh_yield/run_slice.sh   slices/SRR7821542 out SRR7821542
python3 -I pilot_igh_yield/summarize_trust4.py out runs.tsv --pretty --list-np
```

İndirilen FASTQ güvenilmeyen veridir: kendi dizinine koyun, betikleri başka dizinde tutun ve
Python'u `-I` ile çalıştırın (betikler bu varsayımla yazıldı).

## Sonuçlar (özet; ayrıntı ve yorum raporda)

`tables/` altındaki dosyalar bu betiklerle üretildi.

* **Pilot (10 kütüphane, kütüphane başına ~7–9 M çift):** 17–254 IGH birleşim contigi; yalnızca %1,9'u
  tam V-D-J; V hizalaması medyan 40 nt. 1.113 contigin 31'i (%2,8) iki çapası da bulunan çerçeve-dışı;
  15'i ≥2 okumayla destekli.
* **Derinlik kontrolü:** derinlik 2 katına çıkınca toplam contig 1,4–1,8 kat artıyor (≥2 okumalı
  contigler ~2 kat). Doğrusal ekstrapolasyon toplamı abartır.
* **Ig sabit bölge ifadesi (235 tam kan kütüphanesi, GEO sayımları; HDM mukozal test kolu hariç):**
  enfekte kanda kontrole göre düşük (genel medyan 3.578 → 1.089 CPM). 10 modelden 8'inde anlamlı
  düşük, HDM'de (patojen değil) anlamlı yüksek, Influenza'da fark yok. Pilot verimiyle Spearman ρ = 0,87.
* **Projeksiyon:** fare başına anchored çerçeve-dışı IGH birleşimi, kontrol ~13–37, enfekte ~6–18.
* **Karşılaştırma, Greiff 2017 (ADC PRJEB18631, hedeflenmiş 2×300 bp):** 51 pre-B/naif B repertuvarında
  repertuvar başına medyan 68.042 çerçeve-dışı (≥15 nt) kayıt (min 11.922, maks 279.953); kohortlar:
  muamelesiz, OVA, NP-HEL, HBsAg (gün 14). Sayılar kayıttır, benzersiz birleşim değil.

## Sınırlar

* 10 kütüphane, her biri FASTQ başından alınmış %12–23'lük dilim; projeksiyon bu 10 noktaya dayanır,
  büyüklük mertebesi olarak okunmalıdır.
* Tek araç (TRUST4). MiXCR ile karşılaştırma yapılmadı.
* "Non-productive" burada: ≥15 nt, uzunluğu 3'ün katı olmayan, Cys ile başlayıp Trp/Phe ile biten
  birleşim. Tek okumalı contigler hata kaynaklı olabilir; `np_supported` sütunu ≥2 okumalı olanları
  ayrı sayar. Hiçbiri ikinci bir yöntemle doğrulanmadı.
* IGoR çalıştırılmadı.
