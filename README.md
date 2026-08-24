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
dosyasına yazılır. İlk çalıştırma ~200MB'lık gerçek genotip dosyasını indirdiği
için birkaç dakika sürebilir; sonraki çalıştırmalar `.cache/` sayesinde hızlıdır.

## Veri kaynakları

Üç gerçek, halka açık kaynak kullanılıyor — hiçbiri uydurma/örnek veri değil:

1. **SARP skorları** — Hoolehan et al. 2022, NAR 50(20):11696-11711,
   Supplementary Dataset S1 (Europe PMC üzerinden canlı indiriliyor). 4096
   olası 9-mer'in (heptamer'in son 4 bazı + spacer'ın ilk 2 bazı) 1879'u için
   gerçek, üç tekrarlı ölçüm ortalaması.
2. **Gerçek kişi bazlı IGHD genotipleri** — KIARVA'nın resmi, herkese açık
   backend deposu: [`github.com/ScilifelabDataCentre/kiarva-backend`](https://github.com/ScilifelabDataCentre/kiarva-backend),
   `data/compressed/tsv_files-prod.zip` içindeki `1KGP_long_genotypes.tsv`.
   Bu, KIARVA'nın canlı API'sinin sunduğu **aynı** veri (repo README'sinde
   "the currently public data that is exposed in our production environment"
   diye tanımlanıyor), toplu dosya olarak. Corcoran et al. 2026 (Immunity),
   1000 Genomes'tan 2486 kişi, 25 popülasyon. Gerçek, doğrulanmış grup
   büyüklükleri: AFR=708, SAS=543, EAS=540, AMR=281, EUR=414.
   Lisans: veri CC BY-NC 4.0 (atıfla, ticari olmayan kullanım serbest), kod
   MIT — atıf: KIARVA (RRID: SCR_026682) + Corcoran et al., 10.1016/j.immuni.2026.01.026.
3. **Gerçek D-RSS dizileri** — birincil kaynak: yukarıdaki genotip
   dosyasının kendisi. Satıların yaklaşık yarısında (`_F1` gibi eklerle
   işaretli "flank-uzatılmış" okumalar) her alelin D-REGION çekirdek dizisi,
   gerçek çevreleyen genomik bazlarla (heptamer+spacer dahil) birlikte
   veriliyor — 27/28 D geninde, binlerce bireyde doğrulandı. Bu satırlarda
   çekirdek diziyi bulup çevresindeki gerçek bazları çıkararak RSS 9-mer'i
   üretiyoruz. Bir alel bu şekilde hiç yakalanmamışsa, yedek olarak
   VDJbase'in genomik API'si (Rodriguez et al. 2023, Nat Commun, IGenotyper
   uzun-okuma verisi, ~102 birey) kullanılıyor.

## 410 kişi/grup nasıl seçiliyor

Her kıtasal grup (Africa=AFR, Europe=EUR, Asia=EAS+SAS) için, o gruptaki
**gerçek** 2472 bireyden 410'u, her alt popülasyonun (örn. Avrupa için
FIN/GBR/IBS/TSI) gerçek oranı korunarak (en büyük kalan yöntemiyle
yuvarlanmış kota + rastgele örnekleme) seçiliyor. Bunlar simüle değil,
gerçekten var olan, gerçekten genotiplenmiş kişiler — kod `--seed` ile
yeniden üretilebilir şekilde seçiyor.

AMR (Amerika kökenli, karışık) analiz dışı bırakıldı çünkü kıtasal "yerli"
popülasyon karşılaştırmasını bulanıklaştırıyor.

## Bilimsel sınırlama

Yüksek SARP skoru = o RSS'nin RAG1/2 tarafından daha verimli rekombine
edilmesi, yani o D segmentinin repertuvara daha sık girmesi demektir. Bu
pipeline sadece "popülasyonlar arasında bu verimlilikte istatistiksel
anlamlı fark var mı" sorusunu test eder; "bu fark filanca patojene karşı
avantaj sağlıyor" iddiası ayrı bir biyolojik hipotezdir ve bu analizle
kanıtlanmaz, sadece motive edilebilir.

Ayrıca: bir kişinin genotipinde bir D geni için hiç alel çağrısı yoksa
(gerçek bir yapısal delesyon olabilir), o kişi o gen için skorlamaya dahil
edilmiyor (NaN olarak atlanıyor, sıfır veya ortalama değer atanmıyor).
