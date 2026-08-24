# VDJ Rekombinasyon Coğrafi Bias Analizi

Hoolehan ve ark. 2022 (NAR) makalesindeki SARP-seq RSS aktivite skorlarını
kullanarak, farklı kıtasal popülasyonlarda D segmentlerini flanke eden RSS
dizilerinin RAG1/2 rekombinasyon verimliliği açısından anlamlı bir farklılık
gösterip göstermediğini test eder.

## Çalıştırma

```bash
pip install -r requirements.txt
python run_analysis.py --n-per-group 410
```

Sonuçlar `results/` klasörüne CSV olarak, insan-okunur özet `results/report.txt`
dosyasına yazılır. İlk çalıştırma internetten gerçek veri indirdiği için birkaç
dakika sürebilir; sonraki çalıştırmalar `.cache/` sayesinde hızlıdır.

## Veri kaynakları ve neden bu tasarım

Üç gerçek, halka açık kaynak kullanılıyor:

1. **SARP skorları** — Hoolehan et al. 2022, NAR 50(20):11696-11711,
   Supplementary Dataset S1 (Europe PMC üzerinden canlı indiriliyor). 4096
   olası 9-mer'in (heptamer'in son 4 bazı + spacer'ın ilk 2 bazı) 1879'u için
   gerçek, üç tekrarlı ölçüm ortalaması.
2. **Popülasyon allel frekansları** — KIARVA (kiarva.scilifelab.se),
   Corcoran et al. 2026, Immunity, 1000 Genomes'tan 2486 kişi, 25 popülasyon.
   Gerçek API'sinden (`/api/data/frequencies/superpopulations`) her D geni
   aleli için süper-popülasyon başına gerçek frekans + örneklem büyüklüğü (n).
3. **Gerçek D-RSS dizileri** — VDJbase genomik API'si
   (Rodriguez et al. 2023, Nat Commun, IGenotyper uzun-okuma verisi, ~102
   birey). Her D aleli için gerçek 5' ve 3' RSS (heptamer+spacer+nonamer).

### Neden "410 gerçek kişi" yerine simülasyon?

Bu üç kaynağı ayrı ayrı denedik:

- **KIARVA** 2486 kişilik gerçek genotipleme verisine sahip ama kişi-bazlı
  ham genotip tablosunu (makalenin Table S2'si) herkese açık API üzerinden
  vermiyor — sadece allel + popülasyon frekansı sunuyor. Ayrıca RSS flanking
  dizisini de vermiyor, sadece D geninin kodlayan kısmını veriyor.
- **VDJbase**'de gerçek kişi-bazlı RSS dizisi var ama toplam ulaşılabilir
  birey sayısı ~102 (Avrupa ~46, Afrika kökenli ~18, Güney+Doğu Asya ~18) —
  410/grup hedefinin çok altında; bu N ile grup karşılaştırması istatistiksel
  güç açısından zayıf olurdu.

Bu yüzden pipeline şunu yapıyor: her kıtasal grup için 410 "pseudo-birey",
KIARVA'nın **gerçek** süper-popülasyon allel frekanslarından Hardy-Weinberg
varsayımıyla (iki alel bağımsız çekilir) simüle ediliyor, sonra her simüle
bireyin alellerine VDJbase'den gelen **gerçek** RSS dizisi atanıyor. Yani
örneklem büyüklüğü hedefine (410) ulaşılıyor, ama bunu sağlayan olasılık
dağılımının kendisi gerçek, ölçülmüş popülasyon genetiği verisi — uydurma
sayı yok. Bu, popülasyon genetiğinde standart bir tekniktir (allel
frekansından genotip simülasyonu), ama "410 isimli, gerçekten sekanslanmış
kişi" ile karıştırılmamalı; `results/` çıktısı ve konsol mesajları bunu
her zaman "pseudo-birey" olarak etiketler.

Grup eşlemesi: Africa=AFR, Europe=EUR, Asia=EAS+SAS (n-ağırlıklı ortalama).
AMR (Amerika kökenli, karışık) analiz dışı bırakıldı çünkü kıtasal
"yerli" popülasyon karşılaştırmasını bulanıklaştırıyor.

### Bilimsel sınırlama

Yüksek SARP skoru = o RSS'nin RAG1/2 tarafından daha verimli rekombine
edilmesi, yani o D segmentinin repertuvara daha sık girmesi demektir. Bu
pipeline sadece "popülasyonlar arasında bu verimlilikte istatistiksel
anlamlı fark var mı" sorusunu test eder; "bu fark filanca patojene karşı
avantaj sağlıyor" iddiası ayrı bir biyolojik hipotezdir ve bu analizle
kanıtlanmaz, sadece motive edilebilir.
