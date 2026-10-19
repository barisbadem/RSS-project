# Yeni yaklaşım: fare tam kan RNA-seq → non-productive IGH → IGoR

**Fizibilite değerlendirmesi ve literatür taraması** · 2026-10-19
Veri: Singhania ve ark. 2019, *Nat Commun* 10:2887 (GEO GSE119856)

Bu rapor yalnızca bir görüş değil: veriyi indirip işleyerek ölçtüm (pilot, bölüm 3). Ölçümlerin kodu ve
çıktıları [`pilot_igh_yield/`](../pilot_igh_yield/) altında.

---

## 0. Kısa cevap

**Veriden IGH çıkarma kısmı kurulabiliyor ve çalışıyor (IGoR adımını çalıştırmadım). Ama bu veriyle
istediğiniz soru, yani enfeksiyon öncesi/sonrası V(D)J rekombinasyon istatistiğinin IGoR ile
karşılaştırılması, yanıtlanamaz.** Eksik olan kod değil; dizi sayısı ve örnekleme tasarımı. Sırasıyla:

1. **Dizi sayısı yetersiz.** Tam derinliğe yansıtılınca fare başına yaklaşık **13–37 (kontrol) / 6–18
   (enfekte)** adet "çerçeve-dışı ve iki çapası da bulunan" IGH birleşimi çıkıyor. Marcou ve ark. (IGoR)
   TRB için ≥5.000 benzersiz çerçeve-dışı dizi gerektiğini gösteriyor (IGH'de parametre sayısı daha
   fazla). Bir model içinde 10 farelik kol: ~130–370 (kontrol), ~60–175 (enfekte), yani eşiğin 13–80 kat
   altı. 118 kontrol faresinin **hepsini** havuzlasanız bile 1.500–4.400 (eşiğin altında); enfekte kol
   700–2.000, üstelik havuzlamak patojen-özgü soruyu siler.
2. **"Sonra" kolunda bilgi, "önce" kolundan çok daha az.** Ig sabit-bölge ifadesi, 10 modelin 8'inde enfekte
   kanda kontrole göre anlamlı düşük (genel medyan 3.578 → 1.089 CPM; T. gondii 0,20×, Listeria 0,26×,
   B. pseudomallei akut 0,17×). Kollar eşit bilgi içermiyor.
3. **Zamanlama kandaki rekombinasyonu göstermeye uygun değil.** Örnekleme gün 1–7 (tek istisna kronik
   *B. pseudomallei*, gün 27–90). Kandaki olgun B hücrelerinin büyük kısmı enfeksiyondan önce üretilmiş;
   bu pencerede görülecek fark yeni V(D)J olayları değil, hücre bileşimi, seçilim ve trafik olur.
4. **Non-productive, tam anlamıyla "seçilimsiz" değil.** NMD ile kısmen siliniyor, B hücresi durumuna göre
   değişiyor; somatik hipermutasyon (SHM) alır; klonal genişlemeyle birlikte taşınır; plazma hücresi
   transkripsiyonu ağırlığı belirler.
5. **IGoR'da fare IGH için hazır model yok** (hazır modeller: insan TRA/TRB/IGH/IGK/IGL ve fare TRB).
   Özel model kurulabiliyor (Slabodkin 2021 kurmuş) ama ayrı bir iş.

**Literatür:** Birebir bu tasarımı (fare enfeksiyon kontrolü + bulk RNA-seq + non-productive IGH + IGoR)
bulamadım. En yakını Slabodkin ve ark. 2021: IGoR modelleri **inbred C57BL/6 farelerde bile** bireyler
arasında ayırt edilebiliyor, yani rekombinasyon istatistiğinde genetik olmayan bir bileşen var; ve modelin
**zamanla/çevreyle değişip değişmediğini** açıkça "gelecek çalışma" olarak bırakıyorlar. Soru açık ve
değerli; sorun bu veri setinin onu cevaplayacak güçte olmaması.

**Öneri (ayrıntı bölüm 7):** (a) mevcut 244 kütüphaneyle *tanımlayıcı* bir kan IGH repertuvar çalışması
yapılabilir, ama başlığı "rekombinasyon" değil "bileşim/seçilim" olmalı; (b) **asıl soru için hazır ve çok
daha uygun bir veri var: Greiff 2017 (AIRR Data Commons, PRJEB18631).** C57BL/6J pre-B (kemik iliği) ve naif B
(dalak) hücreleri, muamelesiz vs üç antijenle bağışıklanmış kohort (gün 14); repertuvar başına medyan
**68 bin çerçeve-dışı kayıt** (RNA-seq kanda fare başına 6–37). İlk test oradan yapılmalı. Enfeksiyon için
ise ya uygun halka açık veri taranmalı ya da küçük, doğru tasarlanmış bir deney yapılmalı (sıralanmış
naif/geçiş B hücreleri, gDNA, ≥10⁴–10⁵ benzersiz dizi/kol).

> "nonsense" ifadesini **non-productive (çerçeve-dışı veya stop kodonlu, NMD hedefi) IGH transkriptleri**
> olarak yorumladım. Başka bir şey kastettiyseniz söyleyin, bölüm 3–4 değişir.

---

## 1. Soru ve varsayımlar

* Hedef: enfeksiyona maruz kalmış ve kontrol farelerde, seçilimden geçmemiş (pre-selection) V, D, J gen
  kullanımını ve delesyon/insersiyon istatistiğini karşılaştırmak; "VDJ rekombinasyonu çevreye göre
  dinamik olarak değişiyor mu?" sorusuna cevap aramak.
* Yöntem fikri: tam kan RNA-seq'ten IGH non-productive transkriptleri ayıkla → IGoR ile model çıkar
  → enfeksiyon vs kontrol karşılaştır.
* Mantık sağlam: non-productive dizi protein yapmaz, dolayısıyla doğrudan seçilmez (Murugan 2012;
  Marcou 2018; sizin eski projenizin "selection-free" sayımı da aynı mantığa dayanıyor).

## 2. Veri seti: doğrulanmış gerçekler

Kaynak: makalenin Methods bölümü (tam metin), GEO üstverisi, ENA run tablosu (503 run).

* GEO **GSE119856** (SuperSeries; alt seriler GSE119848–GSE119855), BioProject PRJNA490475.
  503 RNA-seq koşusu, hepsi **çift-uçlu Illumina HiSeq 4000, 2×75 bp** (FACS ile ayrılmış CD4 T hücreleri 2×101).
  Toplam FASTQ ≈ **2,2 TB**. 244'ü tam kan (122 kontrol / 122 hastalık); gerisi akciğer, dalak, karaciğer, BAL ve sıralanmış T hücreleri.
* Fareler: **C57BL/6(J), dişi, SPF, 6–18 haftalık**; T. gondii kolunda ayrıca Ifnar⁻/⁻, Ifngr⁻/⁻ ve çift KO.
* Kan RNA'sı: Tempus, PerfectPure RNA Blood Kit, **Mouse GLOBINclear ile globin azaltma**, sonra
  **TruSeq Stranded mRNA (poly-A seçilimi)**, 200 ng; makale "en az 25 M çift okuma" diyor, gerçekte
  kütüphane başına medyan 19–87 M çift.
* Bu seride **Mtb fare örneği yok**; fare tarafında Mtb yerine aşağıdaki modeller var.

| Model (tam kan) | Örnekleme günü | n kontrol/hastalık | Medyan çift (M) kontrol/hastalık |
|---|---|---|---|
| *P. chabaudi* AS (malarya) | 6 | 12 / 13 | 44 / 65 |
| MCMV | 2 | 10 / 10 | 67 / 62 |
| *L. monocytogenes* | 3 | 10 / 10 | 65 / 59 |
| *B. pseudomallei* akut | 3 | 8 / 8 | 59 / 56 |
| *B. pseudomallei* kronik | 27, 39, 49, 65, 90 | 16 / 16 | 62 / 62 |
| *T. gondii* ME-49 | 6–7 | 10 / 7 | 47 / 40 |
| *T. gondii* (Crick: WT, Ifnar, Ifngr, çift KO) | 6–7 | 3/5, 3/5, 3/5, 3/4 | 20–31 |
| Influenza A/X-31 | 6 | 10 / 8 | 87 / 53 |
| RSV A2 | 2 | 10 / 10 | 48 / 68 |
| *C. albicans* | 1 | 10 / 6 | 65 / 49 |
| HDM alerjisi (patojen değil) | 25 | 10 / 10 (+ mukozal test 4/5) | 52 / 19 |

Tasarımla ilgili üç not:

* **Aynı hayvanda önce/sonra yok.** Kontrol ve hastalık farklı farelerdir (gruplar arası karşılaştırma).
* Veri **IFN imzası** için tasarlandı, BCR repertuvarı için değil; örnekleme günleri erken akut faz.
* Modeller farklı tesislerde/tedarikçilerde üretildi (NIMR, Crick, LSHTM, Imperial, Cardiff, Harlan);
  modeller arası havuzlama tedarikçi/yaş/tesis farkıyla karışır, karşılaştırma **model içinde** yapılmalı.

## 3. Pilot: bu veride ne kadar IGH var? (ölçüm)

**Yöntem.** 10 tam kan kütüphanesinin (Malaria, Toxoplasma, Influenza; ikişer kontrol/hastalık, Influenza
birer) FASTQ başından ~300 MiB (**7,1–9,2 M çift**, kütüphanenin %12–23'ü) indirildi. TRUST4 v1.1.11,
fare IMGT referansı (TRUST4 deposundaki `mouse_IMGT+C.fa`). Çerçeve ve stop çağrılarını junction
nükleotidinden kendim yaptım, çünkü TRUST4'ün `productive=F` bayrağı kısmi CDR3'leri de kapsıyor.
"Non-productive" = ≥15 nt, uzunluğu 3'ün katı olmayan, Cys ile başlayıp Trp/Phe ile biten birleşim.

| Kütüphane | Grup | Çift (dilim/toplam, M) | IGH contig | ≥2 okumalı | Çerçeve-dışı | ≥2 okumalı çerçeve-dışı | Tam derinlikte çerçeve-dışı (doğrusal üst sınır) |
|---|---|---|---|---|---|---|---|
| SRR7821541 | Malaria kontrol | 8,8 / 41,3 | 118 | 58 | 3 | 1 | 14 |
| SRR7821542 | Malaria kontrol | 9,2 / 40,2 | 37 | 15 | 1 | 0 | 4 |
| SRR7821545 | Malaria hastalık | 7,2 / 50,6 | 32 | 11 | 1 | 0 | 7 |
| SRR7821555 | Malaria hastalık | 7,8 / 65,0 | 50 | 14 | 2 | 0 | 17 |
| SRR7821637 | Toxoplasma kontrol | 7,3 / 40,6 | 203 | 119 | 3 | 3 | 17 |
| SRR7821644 | Toxoplasma kontrol | 7,3 / 44,1 | 231 | 100 | 7 | 5 | 42 |
| SRR7821646 | Toxoplasma hastalık | 7,1 / 50,0 | 32 | 13 | 0 | 0 | 0 |
| SRR7821647 | Toxoplasma hastalık | 7,5 / 39,5 | 17 | 3 | 0 | 0 | 0 |
| SRR7821659 | Influenza kontrol | 8,0 / 52,2 | 139 | 86 | 6 | 2 | 39 |
| SRR7821668 | Influenza hastalık | 7,5 / 40,1 | 254 | 109 | 8 | 4 | 43 |

Bulgular:

* **Dilim başına 17–254 IGH birleşim contigi** (toplam 1.113). Okuma çifti olarak IGH'ye atananlar
  dilimin %0,0006–0,057'si. Örnekler arası fark 10 katı aşıyor.
* **1.113 contigin 31'i (%2,8) çerçeve-dışı**, 15'i ≥2 okumayla destekli, 8'i çerçeve-içi stop. Bu oran,
  NMD beklentisiyle ve eski projenizin Briney pilot repertuvarındaki büyüklükle (commit mesajlarına göre
  ≈1,34 M productive, 29,6 bin gerçek çerçeve-dışı kayıt, yani ≈%2) aynı mertebede; iki sayının aynı
  repertuvara ait olduğu mesajlardan kesin değil.
* **Contigler kısmi:** yalnızca %1,9'u tam V-D-J. V'ye hizalanan kısım **medyan 40 nt** (IQR 27–94;
  %65'i ≤60 nt), contig uzunluğu medyan 156 nt. Yani V geni neredeyse yalnızca 3' ucundan tanınıyor;
  V geni düzeyinde kullanım ancak aile düzeyinde güvenilir olur. IGoR 130 bp IGH ve 60 bp TRB ile
  test edilmiş olduğundan uzunluk IGoR için sorun değil, V ayırt edilebilirliği sorun.
* **D çağrısı:** %67,8 tekil, %27,9 belirsiz (ör. `IGHD2-4*01|IGHD2-9*02`), %4,2 yok. Eski projedeki
  insan sonucuyla (D'nin %44'ü çağrılamıyor) aynı yönde: D genleri kısa ve kırpılıyor.
* **Tek bir klon desteği domine edebiliyor:** SRR7821644'te çerçeve-dışı bir contig 31 okumayla
  destekli (V çerçevesi 0, J çerçevesi 2; J %100 özdeş, gerçek bir birleşim). Okuma ağırlıklı analiz bu
  klonun etkisinde kalır; benzersiz birleşim bazında sayılmalı (eski projenizin dersi).
* **Derinlik doğrusal ölçeklenmiyor:** aynı kütüphanelerin ilk yarısıyla koşuldu; derinlik 2 katına
  çıkınca toplam contig 1,4–1,8 kat (≥2 okumalı contig ~2 kat) arttı. Tablodaki "tam derinlik"
  sütunu doğrusal ekstrapolasyon olduğundan **üst sınırdır**.

**Tüm kohorta yansıtma.** 235 tam kan kütüphanesinin (HDM mukozal test kolu hariç) GEO HTSeq sayımlarından Ig sabit-bölge ifadesi
(IGH-C CPM) hesaplandı; pilotta bu ifade TRUST4 verimiyle güçlü ilişkili (Spearman ρ = 0,87, p = 0,001,
n = 10). Kalibrasyon doğrusu ve çerçeve-dışı oranı (%2,8) ile her kütüphane için verim yansıtıldı
(alt sınır: derinlik üssü 0,5; üst sınır: doğrusal).

| | Fare başına çerçeve-dışı IGH birleşimi | 10 fare | Tüm fareler havuzlanırsa |
|---|---|---|---|
| Kontrol (118) | 13–37 | 130–370 | 1.500–4.400 |
| Enfekte (117) | 6–18 | 60–175 | 700–2.000 |

Model bazında tek kol: kontrol 16–844, enfekte 2–348 (n = 3–16; en büyük kol kronik *B. pseudomallei*,
16 fare: kontrol 293–844, enfekte 120–348). **5.000 benzersiz diziye** ulaşmak için kontrolde
134–391, enfektede 286–806 fare gerekirdi (mevcut: model başına 3–16). Ayrıntı:
[`tables/cohort_projection.tsv`](../pilot_igh_yield/tables/cohort_projection.tsv).

**Ig ifadesi, kontrol ↔ hastalık (medyan CPM, tam kan).**

| Model | n | Kontrol | Hastalık | Oran | Mann-Whitney p |
|---|---|---|---|---|---|
| *B. pseudomallei* akut | 8/8 | 4.262 | 725 | 0,17 | 0,0002 |
| *B. pseudomallei* kronik | 16/16 | 4.402 | 1.770 | 0,40 | <0,0001 |
| *C. albicans* | 10/6 | 1.471 | 215 | 0,15 | 0,003 |
| *L. monocytogenes* | 10/10 | 3.423 | 894 | 0,26 | 0,0002 |
| MCMV | 10/10 | 3.409 | 1.360 | 0,40 | 0,0002 |
| *P. chabaudi* | 12/13 | 1.814 | 854 | 0,47 | 0,0007 |
| RSV | 10/10 | 3.988 | 3.016 | 0,76 | 0,001 |
| *T. gondii* | 10/7 | 3.807 | 743 | 0,20 | 0,0001 |
| Influenza | 10/8 | 3.354 | 3.650 | 1,09 | 0,20 (fark yok) |
| HDM (patojen değil) | 10/10 | 3.920 | 5.399 | 1,38 | 0,0004 (artış) |
| **Tümü** | 118/117 | **3.578** | **1.089** | **0,30** | |

Cd79a ve Igkc de aynı yönde (ör. T. gondii 0,13× ve 0,23×), yani düşüş yalnızca ağır zincir değil, kandaki
B hücresi RNA payı. Olası nedenler (lenfopeni/yeniden dağılım, nötrofil-monosit akını ile seyrelme) bu
veriyle ayrıştırılamaz; sonuç için nedeni bilmek gerekmiyor: enfekte kolda **daha az dizi** var.

## 4. Teknik riskler

| # | Risk | Kanıt | Ağırlık |
|---|---|---|---|
| 4.1 | **Yetersiz dizi sayısı** | Bölüm 3. IGoR: ~5.000 benzersiz çerçeve-dışı dizi TRB'yi öğrenmeye yetiyor, 10⁵'te çıkarım çok doğru (Marcou 2018, tam metin). Slabodkin 2021 tüm (productive+non-productive) diziyle, 1.000 okumada bile bireyler arası model farkını gürültüden ayırabiliyor; fakat bu, hedeflenmiş 2×300 bp MiSeq repertuvar verisiydi (RNA-seq değil) | **Engelleyici** |
| 4.2 | **Non-productive tanımı** | `productive=F` kısmi CDR3'leri de içeriyor; eski proje 15 nt tabanı koymak zorunda kalmıştı (başarısız anotasyonlar). Tek okumalı contigler hata kaynaklı olabilir; pilotta çerçeve-dışı 31 contigin yalnızca 15'i ≥2 okumalı. Hiçbiri ikinci bir araçla doğrulanmadı | Yüksek |
| 4.3 | **NMD ve hücre durumu** | Pro-B hattında non-productive μ mRNA'sı ~40× düşük, ön-mRNA benzer (Eberle 2009, özet düzeyi). Gözetim B hücresi durumuna göre değişiyor; naif B'de işlenmemiş transkript birikirken aktivasyonla transkripsiyon, splicing ve NMD artıyor (Tinguely 2012); plazma hücrelerinde ek kalite kontrol (Lambert 2019 derlemesi). Srour 2016: Igκ'da V ekzonu stop kodonu ekson atlamasına yol açıp kesik hafif zincir üretiyor ve bu **bialelik rearranjmanlı B hücre klonları plazma hücresi evresinde ayıklanıyor**; yani non-productive allel plazma hücresi evresinde klon için nötr olmayabilir (ağır zincirde test edilip edilmediğini doğrulamadım). **Gruplar arası hücre bileşimi değişirse non-productive havuzun kendisi değişir** | Yüksek |
| 4.4 | **Klonal tekrar, plazma hücresi baskınlığı** | RNA-seq okuma sayısı ifade ağırlıklıdır. Non-productive allel, sağlam diğer allelle birlikte genişleyen klonla taşınır (seçilim diğer allelde). Pilotta tek contig 31 okuma. Eski proje: junction'ların yalnızca %68,8'i benzersiz, biri 34 kez. Benzersiz birleşim bazında sayılmalı | Orta |
| 4.5 | **SHM non-productive allelde de olur** | Pasif allel sistemleriyle gösterilmiş (Chaudhary & Wesemann 2018 derlemesi). IGoR yalnızca naif çerçeve-dışı dizilerle rekombinasyon istatistiği öğrenmiş; bellek çerçeve-dışı diziler SHM modeli için kullanılmış (Marcou 2018). Enfekte farede GC/plazma B hücreleri payı artar. Eski projede V özdeşliği %90–99,6 idi. Pilotta medyan V özdeşliği %100 (V'nin yalnızca 3' ucu hizalanıyor) | Orta |
| 4.6 | **Referans ve suş** | TRUST4'ün fare dosyasında 406 IGHV aleli var; C57BL/6'da ~99 fonksiyonel IGHV bildirilmiş, BALB/c'de 164, yalnızca 5'i ortak (Collins 2015). Dosya büyük olasılıkla birden çok suştan alel içeriyor; C57BL/6 referansına daraltmak gerekir | Orta |
| 4.7 | **IGoR fare IGH** | Hazır modeller: insan TRA/TRB/IGH/IGK/IGL, fare yalnızca TRB (IGoR dokümantasyonu). Özel V/D/J FASTA + CDR3 çapaları + model topolojisi gerekir; Slabodkin 2021 fare IGH için kurmuş (238 V aleli, 25 D, 4 J). IGoR'un en olası senaryosu 130 bp IGH dizilerinin %72'sinde yanlış: belirsizlik kaçınılmaz, bu yüzden olasılıksal çıkarım kullanılıyor | Orta |
| 4.8 | **Kısa okuma, V ve D belirsizliği** | V hizalaması medyan 40 nt; D çağrısının %32'si belirsiz/yok. "Hangi D genleri kullanılıyor" sorusu bu veriyle zayıf, V için aile düzeyi gerçekçi | Orta |
| 4.9 | **İstatistik tasarım** | Eşleşmemiş hayvanlar; kol başına 3–16 fare; fare başına birkaç–onlarca dizi, dolayısıyla havuzlama kaçınılmaz ve diziler bağımsız değil (aynı fareden). İnferans birimi fare olmalı (fare düzeyinde bootstrap/permütasyon). İnbred C57BL/6 farelerde bile model parametreleri bireyler arasında farklı (Slabodkin 2021), yani sıfır hipotezi "fark yok" değil "bireyler arası değişkenlik" | Yüksek |
| 4.10 | **Maliyet** | 2,2 TB; bu sandbox'ta ~8 MB/s ile günler sürer. TRUST4 ~100 sn / 8 M çift (4 iş parçacığı) → 244 kütüphane için yaklaşık 40 CPU-saat. Küme/iş istasyonu gerekir; "indir → işle → sil" akışı disk için şart | Düşük |

## 5. Kavramsal risk: "rekombinasyon dinamik mi?" ile "repertuvar değişti mi?" aynı soru değil

1. **Rekombinasyon kemik iliğinde olur.** Olgun periferik B hücresinde yeni V(D)J olayı beklenmez:
   RAG ifadesi geç evrelerde sürüyor ama immünizasyon sonrası belirgin yeniden indüksiyon görülmemiş
   (Yu 1999; GC'de RAG bildiren Han 1996 tartışmalı).
2. **Kan B hücre havuzu yavaş yenilenir.** Yetişkin farelerde dalak B hücrelerinin üçte ikisinden fazlası
   haftalar-aylar yaşıyor ve havuz yavaş yenileniyor (Förster & Rajewsky 1990); periferik B hücre
   homeostazı kemik iliği akışından kısmen bağımsız (Hao & Rajewsky 2001, başlık düzeyinde). Gün 1–7
   örneklerinde kandaki B hücreleri büyük oranda enfeksiyondan önce üretilmiştir.
3. **Enfeksiyon kemik iliğini etkiler**, ama bu rekombinasyon *kurallarının* değiştiği anlamına gelmez:
   tip I IFN erken B/T gelişimini IL-7'yi baskılayarak geri dönüşümlü engeller (Lin 1998); inflamasyon
   kemik iliğinde lenfopoezi azaltır, granülopoezi artırır (Ueda 2005); *P. chabaudi* kemik iliği B
   lenfopoezini bozar ve dalakta transitional/MZ B hücrelerinde apoptoz yapar (Bockstal 2011). Bunlar
   yeni üretilen hücrelerin **sayısını/zamanlamasını** değiştirir; V(D)J istatistiği değişirse bile bu
   etki kanda ancak haftalar sonra ve seyreltilmiş görünür.
4. **Gözlenen kan farkı büyük olasılıkla bileşim, seçilim ve trafiktir.** Ig ifadesi ve B hücresi
   belirteçleri (Cd79a) enfekte kanda ~3 kat düşüyor (bölüm 3); pilotta kütüphaneler arası IgM/IgD payı
   %42–71 arasında değişiyor ve gruplarla tutarlı bir yön yok (n = 2). Periferik V-aile temsili
   çevresel antijen ve serum Ig'ye göre seçiliyor (Freitas 1991, germ-free fareler).
5. **"Dinamik" olduğu bilinen şey gelişim, enfeksiyon değil:** fare TCRβ'de N insersiyonu embriyodan
   6. haftaya hızla artıyor (Sethna 2017; yalnızca çerçeve-dışı diziyle), fetal/neonatal fare B
   hücrelerinde N bölgesi neredeyse yok (Feeney 1990). Naif V/D kullanımı ikizlerde kalıtsal ve çevre,
   hastalık, kronik lenfosit ablasyonuyla değişmiyor (Glanville 2011); fare naif B hücrelerinde V
   kullanımının >%90'ı genetik olarak öngörülebilir (Greiff 2017).

Sonuç: sıfır hipotezi güçlü (rekombinasyon çevreyle büyük oranda değişmiyor). Pozitif sonuç çok değerli
olurdu, ama yanlış pozitifleri (bileşim/derinlik/parti farkı) ayıklayacak bir tasarım gerekir.

## 6. Literatür: bu tasarımda çalışma var mı?

**Cevap: birebir bulamadım.** Aradığım bileşim: (fare + patojen/kontrol) ∧ (bulk RNA-seq veya repertuvar)
∧ (non-productive/pre-selection IGH) ∧ (rekombinasyon modeli, IGoR/Pgen). Bulduklarım parça parça.
Kanıt düzeyi: **T** = tam metni okudum, **Ö** = özet/arama özeti düzeyinde; tüm künyeleri Europe PMC ile doğruladım.

**A. Yöntem temelleri (seçilimsiz rekombinasyon modeli)**

| Çalışma | Ne yaptı | İlgisi | |
|---|---|---|---|
| Volpe & Kepler 2008, *Immunome Res* 4:3 | İnsan IGH'de ~6.500 productive vs 325 non-productive; non-productive'i seçilimsiz referans aldı | Aynı mantığın erken örneği | Ö |
| Murugan ve ark. 2012, *PNAS* 109:16161 | TCRβ üretim olasılığını non-productive dizilerden çıkardı | Yöntemin kökeni | Ö |
| Elhanati ve ark. 2015, *Phil Trans B* 370:20140243 | İnsan IGH VDJ + SHM olasılıksal çıkarımı | İnsan IGH | Ö |
| Marcou, Mora, Walczak 2018, *Nat Commun* 9:561 (IGoR) | Non-productive dizilerden istatistik; ~5.000 benzersiz TRB çerçeve-dışı dizi yeter, 10⁵'te çok doğru; IGH için naif çerçeve-dışı | Doğrudan araç; dizi sayısı eşiği | T |
| Dokümantasyon (IGoR) | Hazır modeller: insan beş zincir, fare yalnızca TRB | Fare IGH için özel model gerekir | T |

**B. Rekombinasyon istatistiği: genetik mi, çevresel mi?**

| Çalışma | Ne yaptı | İlgisi | |
|---|---|---|---|
| Glanville ve ark. 2011, *PNAS* 108:20066 | İkizler; naif V/D frekansları kalıtsal; çevre, hastalık, lenfosit ablasyonu değiştirmiyor | Çevre etkisinin zayıf olduğu yönünde en güçlü insan verisi | Ö |
| Rubelt ve ark. 2016, *Nat Commun* 7:11112 | İkizlerde naif ve antijen-deneyimli repertuvarlar; kalıtsal bireysel fark | Aynı yönde | Ö |
| Greiff ve ark. 2017, *Cell Rep* 19:1467 | Fare (C57BL/6) pre-B, naif B, plazma hücresi; V kullanımında genetik öngörülebilirlik >%90, antijenin etkisi ≤%40 (plazma hücresi klonal çeşitlilik) | Fare baz çizgisi, halka açık veri (ENA PRJEB18631; açıklaması bu tasarımla birebir örtüşüyor, ArrayExpress E-MTAB-5349 olarak anılıyor) | Ö |
| Collins ve ark. 2015, *Phil Trans B* 370:20140236 | Fare IGH repertuvarı germline odaklı, suşlar arası çok değişken | Suş-özgü referans zorunlu | Ö |
| **Slabodkin ve ark. 2021, *Genome Res* 31:2209** | IGoR modelleri insan ikizlerinde ve **inbred C57BL/6 farelerde** bireyler arası ayırt edilebilir → genetik olmayan modülasyon; modelin zamanla değişimi "gelecek çalışma" | **En yakın çalışma**; hipotezinizin zeminini kuruyor, enfeksiyonu test etmiyor | T |
| Sethna ve ark. 2017, *PNAS* 114:2253 | Fare TCRβ'de E17→6. hafta N insersiyonu artışı; yalnızca çerçeve-dışı dizi | "Dinamik rekombinasyon" için en yakın fare kanıtı, ama gelişimsel, enfeksiyon değil | Ö |
| Feeney 1990, *J Exp Med* 172:1377 | Fetal/neonatal fare IGH'de N bölgesi yok (1/87 DNA, 17/146 RNA) | TdT gelişimsel programı | Ö |
| Russell ve ark. 2022, *eLife* 11:e73475 | İnsan TCRβ'de V(D)J olasılıklarını belirleyen genetik lokuslar (TCRB, MHC; Artemis ve TdT varyantları) | Genetik belirleyiciler | Ö |
| Galson ve ark. 2015, *Front Immunol* 6:531 | İnsan BCR repertuvarında birey içi/arası değişkenlik, zamanla dalgalanma | Zamansal değişkenlik | Ö |
| Barajas-Mora ve ark. 2023, *Nat Immunol* 24:320 | Igk lokusunda enhancer (E34) birincil repertuvarı belirliyor; *S. pneumoniae*'ya karşı koruma | **Ters yön**: rekombinasyon enfeksiyon direncini belirler, enfeksiyon rekombinasyonu değil | Ö |

**C. Çevre, enfeksiyon ve kemik iliği**: Freitas 1991 (*PNAS* 88:5640), Förster & Rajewsky 1990
(*PNAS* 87:4781), Hao & Rajewsky 2001 (*J Exp Med* 194:1151), Yu 1999 (*Nature* 400:682), Lin 1998
(*J Exp Med* 187:79), Ueda 2005 (*J Exp Med* 201:1771), Bockstal 2011 (*J Parasitol Res* 534697);
içerikleri bölüm 5'te.

**D. Bulk RNA-seq'ten repertuvar**

| Çalışma | Ne yaptı | İlgisi | |
|---|---|---|---|
| Blachly ve ark. 2015, *PNAS* 112:4322 | Seçilimsiz RNA-seq'ten tam IGH V-D-J ve SHM (CLL) | Öncül | Ö |
| Bolotin ve ark. 2017, *Nat Biotechnol* 35:908 | MiXCR ile RNA-seq'ten TCR/IG | Araç | Ö |
| Song ve ark. 2021, *Nat Methods* 18:627 (TRUST4) | Bulk BCR doğrulaması ~100 M çift × 150 bp tümör örnekleriyle, MiXCR'dan iyi | Araç; bizim veri çok daha sığ/kısa | Ö |
| Briney ve ark. 2019, *Nature* 566:393 | 10 kişide ~3 milyar IGH dizisi | Yeterli derinliğin ne olduğunu gösteriyor | Ö |

**E. NMD/non-productive transkriptler**: Eberle 2009 (*RNA* 15:1021), Tinguely 2012 (*Mol Cell Biol*
32:107), Srour 2016 (*J Exp Med* 213:109), Lambert, Srour, Delpy 2019 (*BMB Rep* 52:671); içerikleri 4.3'te.

**F. Enfeksiyonda BCR (seçilim bağlamı):** Bashford-Rogers 2019 (*Nature* 574:122) altı bağışıklık-aracılı
hastalıkta kan BCR'sini inceledi (IGHV skewi, izotip); bu, repertuvar farkının seçilimden/izotipten
gelebileceği yönünde.

**G. Halka açık fare BCR verisi.** AIRR Data Commons'ta (VDJServer + iReceptor) fare için yalnızca
**6 çalışma** var ve hiçbiri enfeksiyon-kontrol tasarımında değil; en ilgilileri PRJEB8745 (Collins 2015),
PRJEB18631 (Greiff 2017) ve PRJNA283640 (Cui 2016, SHM modeli). ENA'da başlık taramasıyla fare + repertuvar +
enfeksiyon/bağışıklama temalı birkaç çalışma çıktı (ör. PRJNA714823 influenza-spesifik B hücre BCR,
PRJNA773104 LCMV GC B hücre BCR-seq, PRJNA1209396 viral enfeksiyonda lenf repertuvar dinamiği,
PRJNA682167 *Pneumocystis* BCR). Bunlar yalnızca başlık düzeyinde; çoğu antijen-spesifik veya GC
B hücreleri (güçlü seçilim) olabilir, **uygunlukları ayrıca doğrulanmalı**.

**Greiff 2017 verisi (PRJEB18631) ölçülmüş verimle.** ENA açıklaması: C57BL/6J, 4 kohort (1 muamelesiz,
3 antijenle bağışıklanmış), pre-B, naif B ve plazma hücresi, ek olarak BALB/c ve pet-shop fareleri
(ana makale: *Cell Rep* 2017, V-gen kullanımı, klonal genişleme, çeşitlilik ve boyut incelenmiş;
rekombinasyon modeli parametreleri değil). VDJServer'daki 51 pre-B/naif B repertuvarı için çerçeve-dışı
(≥15 nt, uzunluğu 3'ün katı olmayan) kayıt sayısını sorguladım (eski projenizdeki ADC yöntemiyle).
Sayılar **kayıt**tır, benzersiz birleşim değil; bağışıklama gün 14'te, OVA 100 µg, NP-HEL 100 µg
(ADC'de "HP-HEL"), HBsAg 4 µg:

| Hücre | Kohort | Fare | Çerçeve-dışı kayıt / repertuvar, medyan (min–maks) |
|---|---|---|---|
| Pre-B (kemik iliği) | muamelesiz | 5 | 55.762 (44.420–78.130) |
| Pre-B | OVA | 5 | 68.675 (63.723–102.444) |
| Pre-B | NP-HEL | 5 | 71.877 (60.658–89.175) |
| Pre-B | HBsAg | 4 | 67.710 (26.177–74.547) |
| Naif B (dalak) | muamelesiz | 5 | 139.463 (53.447–279.953) |
| Naif B | OVA | 5 | 70.418 (55.873–85.324) |
| Naif B | NP-HEL | 5 | 72.861 (49.599–118.639) |
| Naif B | HBsAg | 5 | 76.438 (46.509–105.403) |
| Naif B | BALB/c muamelesiz | 4 | 103.250 (62.858–116.565) |
| Naif B | pet-shop | 3 | 20.718 (11.922–23.946) |

Toplam 51 repertuvar, medyan **68.042** çerçeve-dışı kayıt (min 11.922, maks 279.953; toplam 3,9 milyon).
RNA-seq kanda fare başına 6–37 olduğuna göre fark **3 mertebeden fazla**; benzersiz birleşim oranı
%30'a düşse bile C57BL/6J farelerinde tek bir fare IGoR'un ~5.000'lik eşiğini aşar (pet-shop farelerinde
daha düşük: medyan 20.718). Ayrıntı:
[`tables/greiff2017_adc_yield.tsv`](../pilot_igh_yield/tables/greiff2017_adc_yield.tsv).

## 7. Öneri

**Seçenek 1, önerilen: iki katmanlı plan.**

*Katman 1 — mevcut veriyle (1–2 hafta, küme/iş istasyonu):* 244 tam kan kütüphanesinin tamamını
"indir → TRUST4 → sil" akışıyla işle; model içinde kontrol vs hastalık, **tanımlayıcı** IGH repertuvarı:
kütüphane başına IGH contig, izotip dağılımı, V-aile/J kullanımı, CDR3 uzunluğu. Başlık "kan IGH
repertuvarı enfeksiyonla nasıl değişir" olmalı, "rekombinasyon" değil. Zorunlu önlemler:
benzersiz birleşim bazında say; derinlik ve kütüphane başına bilgi miktarını **rarefaction** ile eşitle;
karşılaştırma **model içinde** ve inferans birimi fare olsun; B hücresi/plazma hücresi oranını
(Cd79a, Ms4a1, Jchain, Prdm1, izotip) ortak değişken olarak ekle; negatif kontrol: kontrol-kontrol
rastgele bölme, patojen olmayan HDM, aynı model içinde parti etkisi. Ifnar/Ifngr KO kolu, IFN
bağımlılığını ayırmak için kullanılabilir (n küçük, derinlik düşük: 20–31 M çift).

*Katman 2 — asıl hipotez için doğru veri:*
(a) **İlk gerçek test: Greiff 2017 (PRJEB18631).** Pre-B ve naif B, C57BL/6J, muamelesiz (n=5) vs OVA,
NP-HEL, HBsAg (n=4–5), gün 14. Adımlar: fare IGH IGoR modelini kur ve doğrula (V/D/J FASTA, CDR3 çapaları,
model topolojisi; Slabodkin 2021 bunu yaptı); **bireyler arası gürültü tabanını** ölç (kohort içi fare
çiftleri, teknik replikalar); sonra muamelesiz ↔ bağışıklanmış model uzaklığını bu tabanla karşılaştır.
Hücre tipleri (pre-B, naif B) seçilimden uzak olduğundan sizin sorunuza en yakın hazır test bu;
sınırı: bağışıklama bir enfeksiyon değil, n=4–5 ve tek zaman noktası. Slabodkin'in kullandığı 19 pre-B ve
20 naif B faresi bu verinin dört kohortunun toplamına denk geliyor (5+5+5+4 ve 5+5+5+5); yani fareleri
kohortlara bakmadan "birey" saymış görünüyorlar. Kohort (muamelesiz/bağışıklanmış) bazlı bir model
karşılaştırması yaptıklarını metinde görmedim, ama bunu kesin bir "yapılmadı" olarak değil, kontrol
edilecek bir nokta olarak alın. 2×300 bp okumalarda V-D-J tam, SHM naif/pre-B'de düşük, D çağrısı
RNA-seq'e göre çok daha güvenilir olur.
(b) **Enfeksiyon için uygun halka açık veri taraması:** kontroller, IGH, RNA veya DNA şablon,
UMI, ≥10⁴ dizi/örnek, C57BL/6, seçilimden uzak popülasyon (naif/geçiş B). (c) Uygun veri yoksa
**küçük, tasarlanmış deney**: C57BL/6J dişi, yaş eşleşmeli, kol başına ≥8 fare; kronik bir model ve
zaman serisi (ör. gün 0/7/14/28/56); kemik iliği pre-B/immatür B ve dalak/kan geçiş + naif foliküler B
(FACS), mümkünse RAG2-GFP raportörüyle yeni göç etmiş B hücrelerini ayırma (Yu 1999 aracı);
**gDNA** şablonu (NMD yok, çerçeve-dışı doğal oranda) veya RNA + UMI; kol başına **≥10⁴–10⁵ benzersiz
çerçeve-dışı dizi**; IFNAR/IFNGR KO ve replikalar (teknik/biyolojik). Bu, Slabodkin'in açık bıraktığı
soruyu doğrudan test eder.

**Seçenek 2, yalnızca Singhania verisi, IGoR ile:** yapılabilir ama sonuç düşük çözünürlüklü olur.
Tüm modelleri havuzlayıp "enfekte vs kontrol" (en iyi ihtimalle 700–2.000 vs 1.500–4.400 dizi) bakarsanız
spesifik patojen sorusunu kaybedersiniz ve kollar eşit bilgi içermez. Önermiyorum.

**Karar kriterleri (go / no-go).** Pilotu tüm kohorta genişlettikten sonra:
* kol başına benzersiz çerçeve-dışı dizi **<1.000** → IGoR çalıştırma; yalnızca tanımlayıcı istatistik;
* **1.000–5.000** → IGoR yalnızca V/J seçimi gibi düşük boyutlu parametreler için, belirsizlikle;
* **≥5.000** (tercihen 10⁵) → tam model karşılaştırması.

IGoR'u illa kullanmak gerekmiyor: "hangi V/D/J kullanılıyor" sorusu için TRUST4/MiXCR/IgBLAST çağrıları
yeterli; IGoR'un katkısı üretim olasılığı ve belirsizlik-farkında çıkarım, o da dizi sayısı ister.

## 8. Eski projeden taşınacak dersler

Aynı depodaki eski çalışmanın (insan IGHD, KIARVA/AIRR) commit geçmişinden, bu tasarıma doğrudan uyan dersler:

* Çerçeve-dışı seçimi için **junction ≥15 nt tabanı** (başarısız anotasyonlar sayımı bozuyor) ve
  `productive=false`'a körü körüne güvenmemek.
* **Junction düzeyinde tekilleştirme:** tekrarlar diğer allelin başarısını ölçer, bu allelin değil.
* **Kişi/fare başına eşit ağırlık:** derinlik beş büyüklük mertebesi değişince havuzlanmış yüzde
  derin kütüphaneye göre belirleniyor.
* **D geni çağrısı zor:** kendi D atamanız (`vdj_bias/d_assignment.py`) insan GRCh38 D'lerine göre yazıldı;
  fare için `scripts/fetch_d_germline.py` mantığıyla GRCm39/IMGT D'leri gerekir. Paralog çiftleri
  birlikte okunmalı; çağrılamayanlar açıkça "no_call" olarak ayrı tutulmalı.
* **Orphon** genleri (ör. `/OR`) paydadan çıkarılmalı.
* **Kapsayıcı örnekleme:** uzunluk dağılımından orantılı örnekleme; kısa birleşimler D barındırmaz.
* **Dayanıklı çalıştırma:** konteyner geçici; her kütüphaneden sonra satır ekle, tekrar başlatılabilir olsun.

## 9. Bu raporun sınırları

* **Pilot küçük:** 10 kütüphane, her biri FASTQ başından %12–23'lük dilim. Tam derinliğe yansıtma kaba bir
  büyüklük mertebesi (alt/üst sınır verildi). Tek araç (TRUST4); MiXCR karşılaştırması yapılmadı;
  non-productive çağrılar ikinci bir yöntemle doğrulanmadı; IGoR çalıştırılmadı.
* **Ig ifadesi bir vekil.** HTSeq yalnızca benzersiz hizalanan okumaları sayar (V segmentleri eksik
  sayılır); sabit bölge genleri güvenilir. Vekil, pilotla ρ = 0,87 ama n = 10.
* **Literatür taraması kapsamlı değil.** PubMed/Scopus'a doğrudan erişimim yok; web araması ve Europe
  PMC ile yürüttüm. **Arama aracı üç sorguyu biyogüvenlik kategorisiyle reddetti** (enfeksiyon/inflamasyon
  + B hücre gelişimi + V(D)J/RAG içeren ifadeler). Bunlar zararsız temel immünoloji konuları olduğundan
  yanlış pozitif olduğunu düşünüyorum; sorguları yeniden kurgulayarak aşmaya çalışmadım. Bu yüzden
  "enfeksiyon/inflamasyonun RAG ve V(D)J çıktısına doğrudan etkisi" başlığı daha az taranmış durumda.
  Şu PubMed sorgularını elle çalıştırmanızı öneririm: `("V(D)J recombination" OR "repertoire generation")
  AND (infection OR inflammation) AND (mouse) AND (nonproductive OR "out-of-frame")`;
  `"TdT" AND "B cell precursors" AND (infection OR interferon)`; `"Rag" AND "pre-B" AND "interferon"`.
* **Greiff 2017 verim sayıları ADC kayıt sayısıdır.** Benzersiz birleşim sayısı, VDJServer'ın IgBLAST
  tabanlı `productive`/`junction_length` alanlarının doğruluğu ve replikaların çakışması ayrıca
  denetlenmedi; ham okumalar ENA'da (PRJEB18631) yeniden anotasyon için duruyor. Bağışıklama enfeksiyon
  değildir; gün 14, n=4–5.
* "Birebir çalışma yok" iddiası **kanıt yokluğudur**, yokluğun kanıtı değil.
* Bazı özellikler (ör. Eberle 2009'daki ~40× oranı, Greiff 2017'deki yüzdeler, Collins 2015'teki gen
  sayıları) özet düzeyinde doğrulandı; atıf vermeden önce tam metinden teyit edin.

## Kaynaklar

* Singhania A ve ark. *Nat Commun* 2019;10:2887. GEO GSE119856, PRJNA490475.
* Marcou Q, Mora T, Walczak AM. *Nat Commun* 2018;9:561. IGoR dokümantasyonu: qmarcou.github.io/IGoR
* Slabodkin A ve ark. *Genome Res* 2021;31:2209-2224.
* Song L ve ark. *Nat Methods* 2021;18:627-630. TRUST4: github.com/liulab-dfci/TRUST4
* Bolotin DA ve ark. *Nat Biotechnol* 2017;35:908-911. · Blachly JS ve ark. *PNAS* 2015;112:4322-4327.
* Murugan A ve ark. *PNAS* 2012;109:16161-16166. · Elhanati Y ve ark. *Phil Trans R Soc B* 2015;370:20140243.
* Volpe JM, Kepler TB. *Immunome Res* 2008;4:3.
* Glanville J ve ark. *PNAS* 2011;108:20066-20071. · Rubelt F ve ark. *Nat Commun* 2016;7:11112.
* Greiff V ve ark. *Cell Rep* 2017;19:1467-1478. · Collins AM ve ark. *Phil Trans R Soc B* 2015;370:20140236.
* Sethna Z ve ark. *PNAS* 2017;114:2253-2258. · Feeney AJ. *J Exp Med* 1990;172:1377-1390.
* Russell ML ve ark. *eLife* 2022;11:e73475. · Galson JD ve ark. *Front Immunol* 2015;6:531.
* Barajas-Mora EM ve ark. *Nat Immunol* 2023;24:320-336. · Briney B ve ark. *Nature* 2019;566:393-397.
* Bashford-Rogers RJM ve ark. *Nature* 2019;574:122-126. · Chaudhary N, Wesemann DR. *Front Immunol* 2018;9:462.
* Freitas AA ve ark. *PNAS* 1991;88:5640-5644. · Förster I, Rajewsky K. *PNAS* 1990;87:4781-4784.
* Hao Z, Rajewsky K. *J Exp Med* 2001;194:1151-1164. · Yu W ve ark. *Nature* 1999;400:682-687.
* Han S ve ark. *Science* 1996;274:2094-2097.
* Lin Q, Dong C, Cooper MD. *J Exp Med* 1998;187:79-87. · Ueda Y, Kondo M, Kelsoe G. *J Exp Med* 2005;201:1771-1780.
* Bockstal V, Geurts N, Magez S. *J Parasitol Res* 2011;2011:534697.
* Eberle AB ve ark. *RNA* 2009;15:1021-1028. · Tinguely A ve ark. *Mol Cell Biol* 2012;32:107-117.
* Srour N ve ark. "A plasma cell differentiation quality control ablates B cell clones with biallelic Ig rearrangements and truncated Ig production." *J Exp Med* 2016;213:109-122.
* Lambert JM, Srour N, Delpy L. *BMB Rep* 2019;52:671-678.
