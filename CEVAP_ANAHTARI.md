# Biyoinformatiğe Giriş — Çalışma Kağıdı Cevap Anahtarı

**Ad Soyad:** Semra Kızgındemir

---

## Uygulama 1 — Adli Laboratuvarda Mitokondriyal DNA

### Kalem-Kağıt Isınması

Dizi: 5'-GATCACAGGT-3'

**(a)** Ters tamamlayıcı: **5'-ACCTGTGATC-3'**

**(b)** GC oranı: G=3, C=2 → 5/10 = **%50**

**(c)** 16 569 bp'lik dizide 6 bazlık motif beklentisi: (16569 − 6 + 1) × (1/4)⁶ ≈ 16564 / 4096 ≈ **4,04 kez**

### Görevler

**1.** Kayıt ID: `NC_012920.1`, uzunluk: **16 569 bp**

**2.** A=5124, C=5181, G=2169, T=4094. Evet, 1 adet **N** (belirsiz baz) var, pozisyon **3107**'de. N = "hangi baz olduğu belirlenememiş" anlamında IUPAC belirsizlik kodu. Bu pozisyon tarihsel numaralama uyumu için kasıtlı olarak belirsiz bırakılmıştır.

**3.** GC oranı = (2169 + 5181) / 16568 = **%44,36** (belirsiz baz paydadan çıkarıldı). G (2169) ile C (5181) çok farklı çünkü FASTA'da sadece tek iplik var. Chargaff kuralı (G=C) çift iplikli DNA'nın tamamı için geçerlidir; tek iplikte eşitlik beklenmez.

**4.**

| Enzim | Motif | Gözlenen | Pozisyonlar (1-tabanlı) |
|-------|-------|----------|-------------------------|
| EcoRI | GAATTC | 3 | 4121, 5274, 12640 |
| BamHI | GGATCC | 1 | 14258 |
| HindIII | AAGCTT | 3 | 6203, 11680, 12570 |

Beklenen (eşit olasılık): ~4,04 kez. Gözlenen sayılar beklentiye yakın.

**5.** GATC motifi **23 kez** geçiyor.
- Beklenen (eşit olasılık): 64,71
- Beklenen (gerçek frekanslarla): 51,83
- Gerçek frekans tahmini gözleme daha yakın, ama ikisi de gözlemden çok yüksek. GATC motifi mtDNA'da baskılanmış (under-represented).

**6.**
- Düz metin (1 byte/baz): 16 569 byte = 16,2 KB
- 2-bit kodlama: (16569 × 2 + 7) / 8 = 4 143 byte = 4,0 KB
- Tasarruf: %75

---

## Uygulama 2 — Mitokondri Farklı Konuşur: COX1 Geni

### Kalem-Kağıt Isınması

COX1 ilk 21 baz: ATG TTC GCC GAC CGT TGA CTA

**(a)** mRNA: **AUG UUC GCC GAC CGU UGA CUA**

**(b)** Standart kod (tablo 1): M-F-A-D-R-**\***  → 6. amino asitte dur (TGA = dur)

**(c)** Omurgalı mito kodu (tablo 2): M-F-A-D-R-**W**-L... → protein devam ediyor (TGA = Trp)

### Görevler

**1.** Python dilimi: `genom[5903:7445]`. Gen **1542 baz**, **514 kodon**, 3'ün tam katı.

**2.** mRNA ilk 21 baz: `AUGUUCGCCGACCGUUGACUA` (sadece T→U dönüşümü)

**3.** Standart kodla: **16 dur kodonu (\*)**, ilki **6. amino asitte**. Protein: MFADR*LFST...

**4.** Mito koduyla (tablo 2): **1 dur kodonu** (en sonda). Protein **513 amino asit**. `cds=True` ile çeviri başarılı (hata yok).

**5.**

| Kodon | Sayı | Standart | Mitokondri |
|-------|------|----------|------------|
| TGA | 16 | DUR (*) | Trp (W) |
| ATA | 25 | Ile (I) | Met (M) |
| AGA | 1 | Arg (R) | DUR (*) |
| AGG | 0 | Arg (R) | DUR (*) |

Ekip arkadaşına: "Gen bozuk değil! Mitokondri farklı bir genetik kod kullanıyor. `translate(table=2)` kullanman lazım."

---

## Uygulama 3 — Görme Kaybının Genetik İzi: LHON Mutasyonu

### Kalem-Kağıt Isınması

- Gen içi pozisyon: 11778 − 10760 + 1 = **1019. baz**
- Kodon numarası: (1019−1)/3 + 1 = **340. kodon**
- Kodon içi pozisyon: (1019−1) mod 3 + 1 = **3. harf**

### Görevler

**1.** ND4 uzunluğu: **1378 baz**. 3'ün katı **değil** (kalan: 1). Çünkü mitokondriyal mRNA'lar poli-A kuyruğu alır. Son kodon TA ile bitiyor; poli-A ekleyince TAA (dur kodonu) oluşur.

**2.** ND4 proteini: **459 amino asit** (`to_stop=True` ile)

**3.** m.11778G>A mutasyonu:
- Orijinal baz: G → Mutant baz: A
- Pozisyon **340**'ta amino asit değişimi: **R (Arg) → H (His)**
- Bu bir **anlamlı (missense) mutasyon** → protein işlevi bozulur → LHON hastalığı

**4.** Kodonun 3. harfini değiştirince (wobble pozisyonu) amino asit **değişmez** → **sessiz (synonymous) mutasyon**

**5.** ATP6'da m.8993T>G: Amino asit değişimi **L (Leu) → R (Arg)** → NARP sendromu

---

## Uygulama 4 — LHON Tanı Testi için PCR Primeri

### Kalem-Kağıt Isınması

İleri primer: 5'-TTCTCATAATCGCCCACGGG-3' (20 nt)

**(a)** T=5, C=6, A=3, G=6 (toplam=20)

**(b)** Wallace Tm = 2(5+3) + 4(6+6) = 16 + 48 = **64 °C**

**(c)** GC Tm = 64,9 + 41×(12 − 16,4)/20 = 64,9 + 41×(−4,4)/20 = 64,9 − 9,02 = **55,9 °C**

20 nt için **GC formülü** daha uygundur (Wallace kısa primerler için tasarlanmıştır).

### Görevler

**1.**
- İleri primer: `5'-TTCTCATAATCGCCCACGGG-3'` (11700–11719)
- Geri primer: `5'-AGAGAGTTCTCCCAGTAGGT-3'` (11880–11899 ters tamamlayıcısı)
- Amplikon: **200 bp**
- m.11778 pozisyonu amplikonda: **79. baz** (ortada → iyi tasarım)

**2.**

| Özellik | İleri Primer | Geri Primer |
|---------|-------------|------------|
| GC oranı | %55,0 | %50,0 |
| Wallace Tm | ~64 °C | ~58 °C |
| GC Tm | ~56 °C | ~53 °C |
| NN Tm | ~54,1 °C | ~50,4 °C |

Tm farkı ~3,7 °C → ≤5 °C, primerler uyumlu.

**3.** Her iki primer de genomda **yalnızca 1 kez** eşleşiyor → özgül tasarım.

**4.** GC profili (500 bp pencere, 100 bp adım) çizildi. mtDNA boyunca GC oranı değişkenlik gösteriyor. Amplikon bölgesi orta düzey GC'ye sahip.

---

## Uygulama 5 — Kara Veba'nın İzinde: pPCP1 Plazmiti

### Kalem-Kağıt Isınması

**(a)** Biopython gösterimi: `[6663:7602]`

**(b)** 7602 − 6664 + 1 = **939 bp**

**(c)** 939/3 = **313 kodon** (dur dahil), protein = **312 amino asit**

**(d)** `complement(4815..5888)` → **ters (−) iplik**

### Görevler

**1.** Kimlik: `NC_005816.1`, uzunluk: **9609 bp**, organizma: *Yersinia pestis biovar Microtus str. 91001*, topoloji: **circular**, GC oranı: **%45,3**

**2.** CDS sayısı: **10**. (Tablo notebook'ta.)

**3.** pla geni: 939 bp, 312 amino asit. Kendi çevirimiz kayıttaki `/translation` ile **birebir eşleşiyor**.

**4.** Başlangıç kodonları: ATG (5 gen), GTG (4 gen), ATA (1 gen). GTG/TTG ile başlayan genlerde protein yine M ile başlıyor çünkü ribozom alternatif başlangıç kodonlarını tanır ve her durumda **formil-metionin (fMet)** yerleştirir.

**5.** 10 protein `pPCP1_proteinler.fasta` dosyasına yazıldı. `SeqIO.to_dict` ile geri okunup pla proteinine locus_tag ile erişildi.

**6.** Plazmit haritası çizildi — genler oklarla gösterildi, + iplik yukarıda, − iplik aşağıda.
