# 🧬 Differential Gene Expression Analysis Pipeline

Bu proje, `mock_rep` ve `sars_cov_rep` örnek grupları arasında diferansiyel gen ekspresyonu analizi gerçekleştirmek amacıyla hazırlanmıştır. Aşağıda analiz süreci detaylı adımlarla açıklanmıştır.

---

## 🔄 Analiz İş Akışı

![Workflow](EnformatikWorkflow.png)

Bu görselde veri hazırlığından KEGG analizine kadar tüm analiz adımları özetlenmektedir.

---

## 📊 PCA Görselleştirmesi

Aşağıdaki PCA grafiği, örnekler arasındaki varyansı ve grupların ayrımını göstermektedir:

![PCA](PCA.png)


## 📦 Gerekli Paketler

R ortamında analiz için kullanılan temel kütüphaneler şunlardır: `tidyverse`, `DESeq2`, `gprofiler2`, `ggplot2`.

---

## 📁 Veri Yükleme ve Hazırlama

- `assesment_dataset.tsv` dosyası yüklenir.
- Gen ID'leri (`converted_alias`) alınır ve tekrar edenler filtrelenir.
- Sayımlar tam sayıya yuvarlanır.
- Toplam ifadesi sıfır olan genler filtrelenir.

---

## 🧪 Diferansiyel Gen Ekspresyonu Analizi

- `DESeqDataSetFromMatrix` ile DESeq2 nesnesi oluşturulur.
- Varyans stabilizasyon dönüşümü (VST) uygulanır.
- PCA ile örneklerin dağılımı görselleştirilir.
- `DESeq()` fonksiyonu ile analiz gerçekleştirilir.
- `lfcShrink()` ile log2 fold değişimleri stabilize edilir.

---

## 🎯 Anlamlı Genlerin Filtrelenmesi

- p-adj değeri 0.05'ten küçük olan genler filtrelenir.
- Bu anlamlı genlerin Ensembl ID listesi çıkarılır.

---

## 📈 Fonksiyonel Zenginleştirme Analizi (gProfiler)

gProfiler aracı ile aşağıdaki terim grupları için zenginleştirme analizi yapılır:

- GO:BP (Biyolojik Süreçler)
- GO:MF (Moleküler Fonksiyonlar)
- KEGG
- REAC (Reactome)

Her terim için terim adı, p-değeri ve katkıda bulunan gen sayısı raporlanır.

> **Not:** Analiz, Homo sapiens (insan) türü için yapılmıştır ve `g_SCS` düzeltme yöntemi kullanılmıştır.

---

## 🔄 Entrez ID Dönüşümü

- `converter.tsv` dosyası kullanılarak Ensembl ID’ler Entrez ID’ye dönüştürülür.
- `left_join()` fonksiyonu ile eşleştirme yapılır.

---

## 🧪 KEGG Yol Analizi

1. Yalnızca KEGG kaynaklı terimler filtrelenir.
2. Bu terimler p-değerine göre sıralanır.
3. İlk 5 KEGG yoluna otomatik olarak bağlantılar oluşturulur:

 - https://www.genome.jp/dbget-bin/show_pathway?KEGG:04657
 - https://www.genome.jp/dbget-bin/show_pathway?KEGG:05323
 - https://www.genome.jp/dbget-bin/show_pathway?KEGG:05134
 - https://www.genome.jp/dbget-bin/show_pathway?KEGG:04064
 - https://www.genome.jp/dbget-bin/show_pathway?KEGG:04668


