# Gerekli paketler
library(tidyverse)
library(DESeq2)
library(gprofiler2)
library(ggplot2)

# Veri setini y??kle
counts_data <- read_tsv("assesment_dataset.tsv")

# Gen ID'lerini ??ek
gene_ids <- counts_data$converted_alias

# Gen ID'lerde tekrar edenleri filtrele
counts_data <- counts_data[!duplicated(gene_ids), ]
gene_ids <- counts_data$converted_alias

# Gen ID s??tununu ????kar, geri kalan?? say??m matrisi olarak al
counts_matrix <- counts_data %>% select(-converted_alias)

# DESeq2 i??in say??lar?? tam say??ya yuvarla
counts_matrix <- round(counts_matrix)

# S??f??r toplaml?? (hi??bir ??rnekte ifade edilmeyen) genleri filtrele
filter_mask <- rowSums(counts_matrix) > 0
filtered_counts_matrix <- counts_matrix[filter_mask, ]
filtered_gene_ids <- gene_ids[filter_mask]

# Sat??r isimlerini ata
rownames(filtered_counts_matrix) <- filtered_gene_ids

# Sonu??lar?? yazd??r
cat("Toplam gen say??s??:", nrow(counts_matrix), "\n")
cat("Filtrelenmi?? gen say??s??:", nrow(filtered_counts_matrix), "\n")

# S??tun isimlerine g??re gruplar: mock_rep vs sars_cov_rep
sample_names <- colnames(filtered_counts_matrix)
condition_labels <- ifelse(grepl("mock", sample_names), "mock_rep", "sars_cov_rep")

# Meta-veri olu??tur
colData <- data.frame(
  condition = factor(condition_labels)
)
rownames(colData) <- sample_names

# DESeq2 veri nesnesini olu??tur
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts_matrix,
  colData = colData,
  design = ~ condition
)

# Diferansiyel ifade analizini ??al????t??r
dds <- DESeq(dds)

# <--- BURADAN SONRA PCA EKLEMEK M??MK??N
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

# Sonu??lar?? al (sars_cov_rep vs mock_rep kar????la??t??rmas??)
res <- results(dds, contrast = c("condition", "sars_cov_rep", "mock_rep"))

# Log2 fold change shrink
resLFC <- lfcShrink(dds, coef = "condition_sars_cov_rep_vs_mock_rep", type = "apeglm")


# Sonu??lar?? bir data frame'e d??n????t??r ve anlaml?? olanlar?? filtrele
res_df <- as.data.frame(res) %>%
  rownames_to_column("Geneid") %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

# Anlaml?? genlerin say??s??n?? ve ilk birka?? sat??r??n?? g??ster
cat("Anlaml?? olarak de??i??en gen say??s?? (padj < 0.05):", nrow(res_df), "\n")
head(res_df)

# Anlaml?? genlerin Ensembl ID listesini al
significant_genes <- res_df$Geneid

# gProfiler ile zenginle??tirme analizi yap
gost_results <- gost(
  query = significant_genes,
  organism = "hsapiens", # ??nsan verisi i??in 'hsapiens'
  sources = c("GO:BP", "GO:MF", "KEGG", "REAC"), # ??stenen veritabanlar??
  user_threshold = 0.05, # Anlaml??l??k e??i??i
  correction_method = "g_SCS"
)

# Sonu??lar?? g??ster (interaktif bir tablo d??necektir) veya bir data frame olarak al
if (!is.null(gost_results)) {
  enrichment_table <- gost_results$result
  
  # ??nemli s??tunlar?? se??ip g??ster
  print(enrichment_table %>% select(query, source, term_id, term_name, p_value, intersection_size))
  
  # Sonu??lar?? bir dosyaya da kaydedebilirsiniz
  # write_csv(enrichment_table, "enrichment_results.csv")
} else {
  cat("Zenginle??mi?? bir terim bulunamad??.\n")
}

# read_tsv yerine read_csv kullan??n
converter <- read_csv("converter.tsv")

# S??tunlar?? kontrol edin
colnames(converter)

# Ge??erli verileri filtrele ve sadece gerekli s??tunlar?? al
converter_clean <- converter %>%
  filter(!is.na(converted_alias)) %>%
  select(initial_alias, converted_alias)

# Anlaml?? genlerle e??le??tir
res_df_with_entrez <- res_df %>%
  left_join(converter_clean, by = c("Geneid" = "converted_alias")) %>%
  filter(!is.na(initial_alias))

# Entrez ID'leri al
entrez_ids <- unique(res_df_with_entrez$initial_alias)

cat("Toplam Entrez ID say??s??:", length(entrez_ids), "\n")


# Entrez ID'leri al
entrez_ids <- unique(res_df_with_entrez$initial_alias)

cat("Toplam Entrez ID say??s??:", length(entrez_ids), "\n")

# KEGG terimlerini filtrele
kegg_terms <- enrichment_table %>%
  filter(source == "KEGG") %>%
  arrange(p_value)

# En anlaml?? KEGG yollar??n?? yazd??r
print(kegg_terms %>% select(term_id, term_name, p_value, intersection_size))

# ??lk KEGG yolu i??in link olu??tur
first_kegg_id <- kegg_terms$term_id[1:5]

cat("KEGG Pathway haritas?? ba??lant??s??:\n")
cat(paste0("https://www.genome.jp/dbget-bin/show_pathway?", first_kegg_id, "\n"))