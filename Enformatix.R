# Gerekli paketler
library(tidyverse)
library(DESeq2)
library(gprofiler2)
library(ggplot2)

# Veri setini yükle
counts_data <- read_tsv("assesment_dataset.tsv")

# Gen ID'lerini çek
gene_ids <- counts_data$converted_alias

# Gen ID'lerde tekrar edenleri filtrele
counts_data <- counts_data[!duplicated(gene_ids), ]
gene_ids <- counts_data$converted_alias

# Gen ID sütununu çıkar, geri kalanı sayım matrisi olarak al
counts_matrix <- counts_data %>% select(-converted_alias)

# DESeq2 için sayıları tam sayıya yuvarla
counts_matrix <- round(counts_matrix)

# Sıfır toplamlı (hiçbir örnekte ifade edilmeyen) genleri filtrele
filter_mask <- rowSums(counts_matrix) > 0
filtered_counts_matrix <- counts_matrix[filter_mask, ]
filtered_gene_ids <- gene_ids[filter_mask]

# Satır isimlerini ata
rownames(filtered_counts_matrix) <- filtered_gene_ids

# Sonuçları yazdır
cat("Toplam gen sayısı:", nrow(counts_matrix), "\n")
cat("Filtrelenmiş gen sayısı:", nrow(filtered_counts_matrix), "\n")

# Sütun isimlerine göre gruplar: mock_rep vs sars_cov_rep
sample_names <- colnames(filtered_counts_matrix)
condition_labels <- ifelse(grepl("mock", sample_names), "mock_rep", "sars_cov_rep")

# Meta-veri oluştur
colData <- data.frame(
  condition = factor(condition_labels)
)
rownames(colData) <- sample_names

# DESeq2 veri nesnesini oluştur
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts_matrix,
  colData = colData,
  design = ~ condition
)

# Diferansiyel ifade analizini çalıştır
dds <- DESeq(dds)

# <--- BURADAN SONRA PCA EKLEMEK MÜMKÜN
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

# Sonuçları al (sars_cov_rep vs mock_rep karşılaştırması)
res <- results(dds, contrast = c("condition", "sars_cov_rep", "mock_rep"))

# Log2 fold change shrink
resLFC <- lfcShrink(dds, coef = "condition_sars_cov_rep_vs_mock_rep", type = "apeglm")


# Sonuçları bir data frame'e dönüştür ve anlamlı olanları filtrele
res_df <- as.data.frame(res) %>%
  rownames_to_column("Geneid") %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

# Anlamlı genlerin sayısını ve ilk birkaç satırını göster
cat("Anlamlı olarak değişen gen sayısı (padj < 0.05):", nrow(res_df), "\n")
head(res_df)

# Anlamlı genlerin Ensembl ID listesini al
significant_genes <- res_df$Geneid

# gProfiler ile zenginleştirme analizi yap
gost_results <- gost(
  query = significant_genes,
  organism = "hsapiens", # İnsan verisi için 'hsapiens'
  sources = c("GO:BP", "GO:MF", "KEGG", "REAC"), # İstenen veritabanları
  user_threshold = 0.05, # Anlamlılık eşiği
  correction_method = "g_SCS"
)

# Sonuçları göster (interaktif bir tablo dönecektir) veya bir data frame olarak al
if (!is.null(gost_results)) {
  enrichment_table <- gost_results$result
  
  # Önemli sütunları seçip göster
  print(enrichment_table %>% select(query, source, term_id, term_name, p_value, intersection_size))
  
  # Sonuçları bir dosyaya da kaydedebilirsiniz
  # write_csv(enrichment_table, "enrichment_results.csv")
} else {
  cat("Zenginleşmiş bir terim bulunamadı.\n")
}

# read_tsv yerine read_csv kullanın
converter <- read_csv("converter.tsv")

# Sütunları kontrol edin
colnames(converter)

# Geçerli verileri filtrele ve sadece gerekli sütunları al
converter_clean <- converter %>%
  filter(!is.na(converted_alias)) %>%
  select(initial_alias, converted_alias)

# Anlamlı genlerle eşleştir
res_df_with_entrez <- res_df %>%
  left_join(converter_clean, by = c("Geneid" = "converted_alias")) %>%
  filter(!is.na(initial_alias))

# Entrez ID'leri al
entrez_ids <- unique(res_df_with_entrez$initial_alias)

cat("Toplam Entrez ID sayısı:", length(entrez_ids), "\n")


# Entrez ID'leri al
entrez_ids <- unique(res_df_with_entrez$initial_alias)

cat("Toplam Entrez ID sayısı:", length(entrez_ids), "\n")

# KEGG terimlerini filtrele
kegg_terms <- enrichment_table %>%
  filter(source == "KEGG") %>%
  arrange(p_value)

# En anlamlı KEGG yollarını yazdır
print(kegg_terms %>% select(term_id, term_name, p_value, intersection_size))

# İlk KEGG yolu için link oluştur
first_kegg_id <- kegg_terms$term_id[1:5]

cat("KEGG Pathway haritası bağlantısı:\n")
cat(paste0("https://www.genome.jp/dbget-bin/show_pathway?", first_kegg_id, "\n"))

