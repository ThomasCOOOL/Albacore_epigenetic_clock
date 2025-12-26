#############################################################################################################################################################
##################################################################                                ###########################################################
##################################################################               EWAS             ###########################################################
##################################################################                                ###########################################################
#############################################################################################################################################################
# ===========================================================================================================================================================================
#### 0) Library and Path #### 
# ===========================================================================================================================================================================

#### Library
library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(tibble)
library(purrr)
library(stringr)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(rGREAT)
library(msigdbr)
library(doParallel)
library(WGCNA)
library(iterators)
library(liftOver)
library(tidyr)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(patchwork)

#### Path
root_dir   <- "C:/Users/Download/forSteve3/forSteve"
setwd(root_dir)
meth_file  <- file.path(root_dir, "meth_matrix.csv")
meta_file  <- file.path(root_dir, "Jesstimation_ALB_papier.csv")
cpg_bed_file <- file.path(root_dir, "cpg_bed.tsv")
gff3_file  <- file.path(root_dir, "Thunnus_alalunga.TALA1.1.gff3.gz")
outdir     <- file.path(root_dir, "marginal_rGREAT_results_PRC2_annotated_test")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# CpG window padding (bp on EACH side of CpG midpoint)
pad_bp <- 1

# GREAT domain parameters (basal+extension)
basal_upstream   <- 5000
basal_downstream <- 1000
max_extension    <- 250000

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ===========================================================================================================================================================================
#### 1) Load data & correlations #### 
# ===========================================================================================================================================================================

meth <- fread(meth_file) |> as.data.frame()
stopifnot("cpg_id" %in% names(meth))
meth <- meth[,-1]
rownames(meth) <- meth$cpg_id

meta <- fread(meta_file) |> as.data.frame()
stopifnot(all(c("sample_id","age") %in% names(meta)))
meta$age <- as.numeric(meta$age)

common_samples <- intersect(colnames(meth), meta$sample_id)
if (length(common_samples) < 10) stop("Too few overlapping samples.")
meth <- meth[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$sample_id), ]
stopifnot(identical(colnames(meth), meta$sample_id))

age <- meta$age; names(age) <- meta$sample_id

cor_fun <- function(x) suppressWarnings(cor(x, age, use="pairwise.complete.obs", method="pearson"))
p_fun   <- function(x) suppressWarnings(cor.test(x, age, method="pearson")$p.value)
cors <- apply(meth, 1, cor_fun)
pval <- apply(meth, 1, p_fun)

res <- tibble(
  cpg_id = rownames(meth),
  r      = as.numeric(cors),
  p      = as.numeric(pval)
) |>
  mutate(q = p.adjust(p, method="BH"),
         sign = ifelse(r >= 0, "POS","NEG"))
write_csv(res, file.path(outdir, "CpG_age_marginal_correlations.csv"))

# Histogram
library(ggplot2)

g_hist <- ggplot(res, aes(x = r)) +
  geom_histogram(
    bins = 60,
    color = "black",
    fill = "grey70",
    linewidth = 0.3
  ) +
  labs(
    x = "Pearson correlation coefficient (r)",
    y = "Number of CpG sites"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 14, color = "black"),
    axis.text  = element_text(size = 12, color = "black"),
    axis.line  = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.title = element_blank()
  )

g_hist

ggsave(file.path(outdir, "hist_correlations_final.png"), g_hist, width=6, height=4, dpi=600)

# ===========================================================================================================================================================================
#### 2) CpG → GRanges (+/- pad) #### 
# ===========================================================================================================================================================================

parse_cpg_ids_to_bed <- function(cpg_ids) {
  id0 <- gsub("\\.\\d+$", "", cpg_ids)  
  two_cols <- sub("(.*)[_:\\-](\\d+)$", "\\1\t\\2", id0, perl = TRUE)
  ok <- grepl("\t", two_cols)
  if (!all(ok)) {
    warning("Some CpG IDs could not be parsed: e.g. ", paste(head(id0[!ok]), collapse = ", "))
  }
  mat <- strsplit(two_cols[ok], "\t")
  chr <- vapply(mat, `[`, character(1), 1)
  pos <- as.integer(vapply(mat, `[`, character(1), 2))
  tibble(chr = chr, start = pos, end = pos + 1L, cpg_id = cpg_ids[ok])
}


if (!file.exists(cpg_bed_file)) {
  bedfile <- parse_cpg_ids_to_bed(res$cpg_id)
  write_tsv(bedfile, cpg_bed_file)
}
cpg_bed <- fread(cpg_bed_file) |> as.data.frame()
stopifnot(all(c("chr","start","end","cpg_id") %in% names(cpg_bed)))

gr_all <- makeGRangesFromDataFrame(cpg_bed, keep.extra.columns = TRUE,
                                   seqnames.field = "chr", start.field = "start", end.field = "end",
                                   ignore.strand = TRUE)

if (!is.null(pad_bp) && pad_bp > 0) {
  gr_all <- GenomicRanges::resize(gr_all,
                                  width = pmax(1L, BiocGenerics::width(gr_all) + 2L*pad_bp),
                                  fix   = "center")
}

# ===========================================================================================================================================================================
#### 3) extended_tss from GFF3 with GFF3 gene IDs ####
# ===========================================================================================================================================================================

build_extended_tss_gff_id <- function(gff3_file,
                                      mode = "basalPlusExt",
                                      basal_upstream = 5000,
                                      basal_downstream = 1000,
                                      extension = 50000,
                                      extend_from = "TSS") {
  stopifnot(file.exists(gff3_file))
  txdb <- GenomicFeatures::makeTxDbFromGFF(gff3_file, format = "gff3")
  g_tx <- GenomicFeatures::genes(txdb)
  
  g_all <- rtracklayer::import(gff3_file)
  g_gg  <- g_all[g_all$type == "gene"]
  gff_id <- S4Vectors::mcols(g_gg)$ID
  
  # Map txdb genes to GFF3 genes
  hits <- GenomicRanges::nearest(g_tx, g_gg)
  id_map <- rep(NA_character_, length(g_tx))
  ok <- !is.na(hits)
  id_map[ok] <- as.character(gff_id[hits[ok]])
  
  if (anyNA(id_map)) {
    ov <- GenomicRanges::findOverlaps(g_tx, g_gg, select = "arbitrary")
    ok2 <- is.na(id_map) & !is.na(ov)
    id_map[ok2] <- as.character(gff_id[ov[ok2]])
  }
  
  tss_pos <- ifelse(as.character(GenomicRanges::strand(g_tx)) == "+",
                    GenomicRanges::start(g_tx), GenomicRanges::end(g_tx))
  tss_gr <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(g_tx),
    ranges   = IRanges::IRanges(start = tss_pos, end = tss_pos),
    strand   = GenomicRanges::strand(g_tx)
  )
  S4Vectors::mcols(tss_gr)$gene_id <- ifelse(!is.na(id_map), id_map, names(g_tx))
  
  rGREAT::extendTSS(
    tss_gr,
    mode = mode,
    extend_from = extend_from,
    basal_upstream = basal_upstream,
    basal_downstream = basal_downstream,
    extension = extension
  )
}


build_gene_annot <- function(gff3_file) {
  g_all <- rtracklayer::import(gff3_file)
  g_gg  <- g_all[g_all$type == "gene"]
  # g_gg  <- g_all
  df <- tibble::tibble(
    seqid        = as.character(GenomicRanges::seqnames(g_gg)),
    start        = as.integer(GenomicRanges::start(g_gg)),
    end          = as.integer(GenomicRanges::end(g_gg)),
    strand       = as.character(GenomicRanges::strand(g_gg)),
    gene_id      = as.character(S4Vectors::mcols(g_gg)$locus_tag),
    annotation   = as.character(S4Vectors::mcols(g_gg)$type),
    symbol       = as.character(S4Vectors::mcols(g_gg)$gene %||% S4Vectors::mcols(g_gg)$Name %||% ""),
    description  = as.character(S4Vectors::mcols(g_gg)$description %||% S4Vectors::mcols(g_gg)$product %||% ""),
    gene_biotype = as.character(S4Vectors::mcols(g_gg)$gene_biotype)
  )
  df <- dplyr::filter(df, !is.na(gene_id) & gene_id != "")
  df <- dplyr::mutate(df,
                      tss = dplyr::if_else(strand == "+", start, end)
  )
  df
}


extended_tss <- build_extended_tss_gff_id(
  gff3_file,
  mode = "basalPlusExt",
  basal_upstream = basal_upstream,
  basal_downstream = basal_downstream,
  extension = max_extension
)
extended_tss$gene_id <- gsub(pattern="gene-", replacement = "", extended_tss$gene_id)
saveRDS(extended_tss, file.path(outdir, "tala_extended_tss_GFFID.rds"))

cat("Example extended_tss gene IDs:\n")
print(head(S4Vectors::mcols(extended_tss)$gene_id, 10))
gene_annot <- build_gene_annot(gff3_file)

# ===========================================================================================================================================================================
#### 4) Run msig #### 
# ===========================================================================================================================================================================
hits <- findOverlaps(gr_all, extended_tss,ignore.strand=T, select = "all")

gr_all_nearest <- gr_all[queryHits(hits)]
tss_nearest    <- extended_tss[subjectHits(hits)]
tss_nearest$cpg_id <- gr_all_nearest$cpg_id
tss_nearest$cpg_pos <- tss_nearest$cpg_id
tss_nearest$cpg_pos <- sub("^chr[0-9]+\\.[0-9]+_", "", tss_nearest$cpg_pos)
tss_nearest$cpg_pos <- as.numeric(tss_nearest$cpg_pos)
tss_nearest$diff_dist <- tss_nearest$cpg_pos - tss_nearest$tss_position
tss_nearest$diff_dist <- abs(tss_nearest$diff_dist)
tss_nearest$gene_id <- gsub("gene-","",tss_nearest$gene_id)
tss_nearest$hits <- subjectHits(hits)
max(tss_nearest$diff_dist)

tss_nearest_df <- data.frame(
   cpg_id          = tss_nearest$cpg_id ,
   cpg_pos         = tss_nearest$cpg_pos ,
   gene_id         = tss_nearest$gene_id ,
   tss_position    = tss_nearest$tss_position ,
   tss_range_start = tss_nearest@ranges@start ,
   tss_range_width = tss_nearest@ranges@width ,
   diff_dist       = tss_nearest$diff_dist,
   hit             = tss_nearest$hits
 )
tss_nearest_df$tss_range_end <- tss_nearest_df$tss_range_start + tss_nearest_df$tss_range_width 
tss_nearest_df$tss_range_end <- tss_nearest_df$tss_range_end-1

tss_nearest_df <- tss_nearest_df[,-6]
tss_nearest_df <- tss_nearest_df[c(1:5,8,6,7)]
tss_nearest_df$in_range <- tss_nearest_df$cpg_pos >= tss_nearest_df$tss_range_start &
  tss_nearest_df$cpg_pos <= tss_nearest_df$tss_range_end

# regions_annotated <- tss_nearest_df %>%
#    group_by(cpg_id) %>%
#    slice_min(diff_dist, with_ties = FALSE) %>%
#    ungroup()
 
# regions_annotated$first_hit <- hits_first

# regions_annotated$verif_hit <- regions_annotated$hit == regions_annotated$first_hit

top150_annotated <- res %>%
  mutate(abs_r = abs(r)) %>%
  filter(cpg_id %in% tss_nearest_df$cpg_id) %>%
  arrange(desc(abs_r)) %>%
  filter(p < 0.05)

top150_annotated_NEG <- top150_annotated %>%
  filter(sign == "NEG") %>%
  slice_head(n = 150)

top150_annotated_POS <- top150_annotated %>%
  filter(sign == "POS")%>%
  slice_head(n = 150)


top150_annotated_POS <- merge(tss_nearest_df, top150_annotated_POS,by = "cpg_id", all.x=F, all.y = T)
genes_of_interest_POS <- top150_annotated_POS$gene_id
length(genes_of_interest_POS)
summary(is.na(genes_of_interest_POS))

top150_annotated_NEG <- merge(tss_nearest_df, top150_annotated_NEG,by = "cpg_id", all.x=F, all.y = T)
genes_of_interest_NEG <- top150_annotated_NEG$gene_id
length(genes_of_interest_NEG)
summary(is.na(genes_of_interest_NEG))

genes_background <- tss_nearest_df$gene_id
summary(is.na(genes_background))
length(genes_background)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# Collection biological pathways and gene sets
collection <- msigdbr_collections()

##### Human Phenotype ##### 
msig_human_phenotype <- msigdbr(species = "Homo sapiens", subcollection = "HPO")

# ================= Aging and mortality ================= 
msig_human_phenotype <- msig_human_phenotype |>
  filter(grepl("(^|_)aging(_|$)|(^|_)mortality(_|$)|lifespan|senescence|(^|_)age(_|$)", gs_name, ignore.case = TRUE))

sym_col_human_phenotype <- intersect(c("gene_symbol","human_gene_symbol","hgnc_symbol"), names(msig_human_phenotype))[1]

ga_human_phenotype <- gene_annot %>% dplyr::mutate(symbol_upper = toupper(.data$symbol))
ms_human_phenotype <- msig_human_phenotype   %>% dplyr::mutate(symbol_upper = toupper(.data$gene_symbol))

mapped_human_phenotype <- dplyr::inner_join(ms_human_phenotype, ga_human_phenotype, by = "symbol_upper")
out_mapped_human_phenotype <- split(mapped_human_phenotype$gene_id, mapped_human_phenotype$gs_name)

#### Human Phenotype Aging|Mortality Background
filtered_list_mapped_aging_background <- lapply(out_mapped_human_phenotype, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_aging_background))]))

#### Human Phenotype Top 150POS
filtered_list_mapped_aging_Top_150_POS <- lapply(out_mapped_human_phenotype, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_aging_Top_150_POS))]))

#### Human Phenotype Top 150NEG
filtered_list_mapped_aging_Top_150_NEG <- lapply(out_mapped_human_phenotype, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_aging_Top_150_NEG))]))


# ================= Human Nervous system phenotype =================
msig_human_phenotype_nervous <- msigdbr(species = "Homo sapiens", subcollection = "HPO")
msig_human_phenotype_nervous <- msig_human_phenotype_nervous %>%
  filter(str_detect(as.character(gs_name), fixed("NERVOUS_SYSTEM", ignore_case = TRUE)))

sym_col_human_phenotype_nervous <- intersect(c("gene_symbol","human_gene_symbol","hgnc_symbol"), names(msig_human_phenotype_nervous))[1]

ga_human_phenotype_nervous <- gene_annot %>% dplyr::mutate(symbol_upper = toupper(.data$symbol))
ms_human_phenotype_nervous <- msig_human_phenotype_nervous   %>% dplyr::mutate(symbol_upper = toupper(.data$gene_symbol))

mapped_human_phenotype_nervous <- dplyr::inner_join(ms_human_phenotype_nervous, ga_human_phenotype_nervous, by = "symbol_upper")
out_mapped_human_phenotype_nervous <- split(mapped_human_phenotype_nervous$gene_id, mapped_human_phenotype_nervous$gs_name)

#### Human Nervous system phenotype Background

filtered_list_mapped_nerv_sys_pheno_process_background <- lapply(out_mapped_human_phenotype_nervous, function(x) {
  x[x %in% genes_background]
})

sum(unique(unlist(filtered_list_mapped_nerv_sys_pheno_process_background)) %in% genes_background)

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_nerv_sys_pheno_process_background))]))

#### Human Nervous system phenotype Top 150POS

filtered_list_mapped_nerv_sys_pheno_process_Top_150POS <- lapply(out_mapped_human_phenotype_nervous, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_nerv_sys_pheno_process_Top_150POS))]))

#### Human Nervous system phenotype Top 150NEG

filtered_list_mapped_nerv_sys_pheno_process_Top_150NEG <- lapply(out_mapped_human_phenotype_nervous, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_nerv_sys_pheno_process_Top_150NEG))]))


# length(unique(unlist(filtered_list_mapped_nerv_sys_pheno_process_Top_150NEG)))
# length(unlist(filtered_list_mapped_nerv_sys_pheno_process_Top_150NEG))
# 
# nerv_sys_pheno_top150NEG <- data.frame(
#   Enrichment = rep("Human phenotype", length(unique(unlist(filtered_list_mapped_nerv_sys_pheno_process_Top_150NEG)))),
#   Annotation = rep("Nervous system phenotype", length(unique(unlist(filtered_list_mapped_nerv_sys_pheno_process_Top_150NEG)))),
#   Value = unique(unlist(filtered_list_mapped_nerv_sys_pheno_process_Top_150NEG)),
#   stringsAsFactors = FALSE
# )

# ================= Cancer =================
msig_Cancer <- msigdbr::msigdbr(collection = "C4")
sym_col_Cancer <- intersect(c("gene_symbol","human_gene_symbol","hgnc_symbol"), names(msig_Cancer))[1]

ga_Cancer <- gene_annot %>% dplyr::mutate(symbol_upper = toupper(.data$symbol))
ms_Cancer <- msig_Cancer   %>% dplyr::mutate(symbol_upper = toupper(.data$gene_symbol))

mapped_Cancer <- dplyr::inner_join(ms_Cancer, ga_Cancer, by = "symbol_upper")
out_mapped_Cancer <- split(mapped_Cancer$gene_id, mapped_Cancer$gs_name)

#### Cancer Background
filtered_list_mapped_Cancer_background <- lapply(out_mapped_Cancer, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_Cancer_background))]))

#### Cancer Top 150POS
filtered_list_mapped_Cancer_Top_150POS <- lapply(out_mapped_Cancer, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_Cancer_Top_150POS))]))

#### Cancer Top 150NEG
filtered_list_mapped_Cancer_Top_150NEG <- lapply(out_mapped_Cancer, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_Cancer_Top_150NEG))]))

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### GO Biological process ##### 
# ================= ANATOMICAL_STRUCTURE_DEVELOPMENT https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/ANATOMICAL_STRUCTURE_DEVELOPMENT.html ================= 

ANATOMICAL_STRUCTURE_DEVELOPMENT_tsv <- readr::read_tsv("C:/Users/Download/ANATOMICAL_STRUCTURE_DEVELOPMENT.v2025.1.Hs.tsv")
ANATOMICAL_STRUCTURE_DEVELOPMENT_gene <- unlist(strsplit(ANATOMICAL_STRUCTURE_DEVELOPMENT_tsv$ANATOMICAL_STRUCTURE_DEVELOPMENT[16], split = ","))

msig_GO_Biological_process_anat_str_dev <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "BP")
msig_GO_Biological_process_anat_str_dev <- msig_GO_Biological_process_anat_str_dev[msig_GO_Biological_process_anat_str_dev$gene_symbol %in% ANATOMICAL_STRUCTURE_DEVELOPMENT_gene,]

ga_GO_Biological_process_anat_str_dev <- gene_annot %>% dplyr::mutate(symbol_upper = toupper(.data$symbol))
ms_GO_Biological_process_anat_str_dev <- msig_GO_Biological_process_anat_str_dev %>% dplyr::mutate(symbol_upper = toupper(.data$gene_symbol))

mapped_GO_Biological_process_anat_str_dev <- dplyr::inner_join(ms_GO_Biological_process_anat_str_dev, ga_GO_Biological_process_anat_str_dev, by = "symbol_upper")
out_mapped_GO_Biological_process_anat_str_dev <- split(mapped_GO_Biological_process_anat_str_dev$gene_id, mapped_GO_Biological_process_anat_str_dev$gs_name)

#### ANATOMICAL_STRUCTURE_DEVELOPMENT Background

filtered_list_mapped_anat_str_dev_process_background <- lapply(out_mapped_GO_Biological_process_anat_str_dev, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_anat_str_dev_process_background))]))

#### ANATOMICAL_STRUCTURE_DEVELOPMENT Top 150POS 

filtered_list_mapped_anat_str_dev_process_Top_150POS <- lapply(out_mapped_GO_Biological_process_anat_str_dev, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_anat_str_dev_process_Top_150POS))]))


#### ANATOMICAL_STRUCTURE_DEVELOPMENT Top 150NEG 

filtered_list_mapped_anat_str_dev_process_Top_150NEG <- lapply(out_mapped_GO_Biological_process_anat_str_dev, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_anat_str_dev_process_Top_150NEG))]))

# ================= Developpemental process https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SYSTEM_DEVELOPMENT.html ================= 

msig_GO_Biological_process <- msigdbr::msigdbr(subcollection = "GO:BP")
sym_col_GO_Biological_process <- intersect(c("gene_symbol","human_gene_symbol","hgnc_symbol"), names(msig_GO_Biological_process))[1]

msig_GO_Biological_process <- msig_GO_Biological_process %>%
  filter(str_detect(as.character(gs_name), fixed("SYSTEM_DEVELOPMENT", ignore_case = TRUE)))

ga_GO_Biological_process <- gene_annot %>% dplyr::mutate(symbol_upper = toupper(.data$symbol))
ms_GO_Biological_process <- msig_GO_Biological_process   %>% dplyr::mutate(symbol_upper = toupper(.data$gene_symbol))

mapped_GO_Biological_process <- dplyr::inner_join(ms_GO_Biological_process, ga_GO_Biological_process, by = "symbol_upper")
out_mapped_GO_Biological_process <- split(mapped_GO_Biological_process$gene_id, mapped_GO_Biological_process$gs_name)

#### Developpemental process Background

filtered_list_mapped_GO_Biological_process_background <- lapply(out_mapped_GO_Biological_process, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_GO_Biological_process_background))]))

#### ANATOMICAL_STRUCTURE_DEVELOPMENT Top 150POS 

filtered_list_mapped_GO_Biological_process_Top_150POS <- lapply(out_mapped_GO_Biological_process, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_GO_Biological_process_Top_150POS))]))


#### ANATOMICAL_STRUCTURE_DEVELOPMENT Top 150NEG 

filtered_list_mapped_GO_Biological_process_Top_150NEG <- lapply(out_mapped_GO_Biological_process, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_GO_Biological_process_Top_150NEG))]))



# ================= RNA metabolic process  https://www.gsea-msigdb.org/gsea/msigdb/cards/GOBP_REGULATION_OF_MRNA_METABOLIC_PROCESS ================= 
msig_GO_Biological_process_RNA <- msigdbr::msigdbr(subcollection = "GO:BP")
sym_col_GO_Biological_process_RNA <- intersect(c("gene_symbol","human_gene_symbol","hgnc_symbol"), names(msig_GO_Biological_process))[1]

msig_GO_Biological_process_RNA <- msig_GO_Biological_process_RNA %>%
  filter(str_detect(as.character(gs_name), fixed("REGULATION_OF_MRNA_METABOLIC_PROCESS", ignore_case = TRUE)))

ga_GO_Biological_process_RNA <- gene_annot %>% dplyr::mutate(symbol_upper = toupper(.data$symbol))
ms_GO_Biological_process_RNA <- msig_GO_Biological_process_RNA %>% dplyr::mutate(symbol_upper = toupper(.data$gene_symbol))

mapped_GO_Biological_process_RNA <- dplyr::inner_join(ms_GO_Biological_process_RNA, ga_GO_Biological_process_RNA, by = "symbol_upper")
out_mapped_GO_Biological_process_RNA <- split(mapped_GO_Biological_process_RNA$gene_id, mapped_GO_Biological_process_RNA$gs_name)


#### RNA metabolic process Background

filtered_list_mapped_RNA_metabolic_process_background <- lapply(out_mapped_GO_Biological_process_RNA, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_RNA_metabolic_process_background))]))


#### RNA metabolic process Top 150 POS

filtered_list_mapped_RNA_metabolic_process_Top_150POS <- lapply(out_mapped_GO_Biological_process_RNA, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_RNA_metabolic_process_Top_150POS))]))


#### RNA metabolic process Top 150 NEG

filtered_list_mapped_RNA_metabolic_process_Top_150NEG <- lapply(out_mapped_GO_Biological_process_RNA, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_RNA_metabolic_process_Top_150NEG))]))


# ================= Nervous systeme dev https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/NERVOUS_SYSTEM_DEVELOPMENT.html ================= 
msig_GO_Biological_process_nerv_dev <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "BP")
msig_GO_Biological_process_nerv_dev <- msig_GO_Biological_process_nerv_dev %>%
  filter(str_detect(as.character(gs_name), fixed("NERVOUS_SYSTEM_DEVELOPMENT", ignore_case = TRUE)))

ga_GO_Biological_process_nerv_dev <- gene_annot %>% dplyr::mutate(symbol_upper = toupper(.data$symbol))
ms_GO_Biological_process_nerv_dev <- msig_GO_Biological_process_nerv_dev %>% dplyr::mutate(symbol_upper = toupper(.data$gene_symbol))

mapped_GO_Biological_process_nerv_dev <- dplyr::inner_join(ms_GO_Biological_process_nerv_dev, ga_GO_Biological_process_nerv_dev, by = "symbol_upper")
out_mapped_GO_Biological_process_nerv_dev <- split(mapped_GO_Biological_process_nerv_dev$gene_id, mapped_GO_Biological_process_nerv_dev$gs_name)

#### Nervous systeme dev Background

filtered_list_mapped_DEV_NERV_process_background <- lapply(out_mapped_GO_Biological_process_nerv_dev, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_DEV_NERV_process_background))]))

#### Nervous systeme dev Top 150 pos

filtered_list_mapped_DEV_NERV_process_Top_150_POS <- lapply(out_mapped_GO_Biological_process_nerv_dev, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_DEV_NERV_process_Top_150_POS))]))


#### Nervous systeme dev Top 150 neg

filtered_list_mapped_DEV_NERV_process_Top_150_NEG <- lapply(out_mapped_GO_Biological_process_nerv_dev, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_DEV_NERV_process_Top_150_NEG))]))


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### GO Molecular Function ##### 

msig_GO_Molecular_Function <- msigdbr::msigdbr(subcollection = "GO:MF")
sym_col_GO_Molecular_Function <- intersect(c("gene_symbol","human_gene_symbol","hgnc_symbol"), names(msig_GO_Molecular_Function))[1]

ga_GO_Molecular_Function <- gene_annot %>% dplyr::mutate(symbol_upper = toupper(.data$symbol))
ms_GO_Molecular_Function <- msig_GO_Molecular_Function   %>% dplyr::mutate(symbol_upper = toupper(.data$gene_symbol))
mapped_GO_Molecular_Function <- dplyr::inner_join(ms_GO_Molecular_Function, ga_GO_Molecular_Function, by = "symbol_upper")

out_mapped_GO_Molecular_Function <- split(mapped_GO_Molecular_Function$gene_id, mapped_GO_Molecular_Function$gs_name)

# ================= Transcription factor binding ================= 
filtered_list_TF <- out_mapped_GO_Molecular_Function[grepl("DNA_BINDING_TRANSCRIPTION_FACTOR", names(out_mapped_GO_Molecular_Function))]

#### Transcription factor binding Background

filtered_list_mapped_GO_Biological_function_TF_background <- lapply(filtered_list_TF, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_GO_Biological_function_TF_background))]))

#### Transcription factor binding Top 150 POS

filtered_list_mapped_GO_Biological_function_TF_Top_150_POS <- lapply(filtered_list_TF, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_GO_Biological_function_TF_Top_150_POS))]))

#### Transcription factor binding Top 150 NEG

filtered_list_mapped_GO_Biological_function_TF_Top_150_NEG <- lapply(filtered_list_TF, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_GO_Biological_function_TF_Top_150_NEG))]))



##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### MSigDB pathway ##### 

msig_MSigDB_C2 <- msigdbr::msigdbr(collection = "C2")
msig_MSigDB_perturbation <- msigdbr::msigdbr(subcollection = "CGP")
msig_MSigDB_pathway <- anti_join(msig_MSigDB_C2, msig_MSigDB_perturbation, by = "gene_symbol")

sym_col_MSigDB_pathway <- intersect(c("gene_symbol","human_gene_symbol","hgnc_symbol"), names(msig_MSigDB_pathway))[1]

ga_MSigDB_pathway <- gene_annot %>% dplyr::mutate(symbol_upper = toupper(.data$symbol))
ms_MSigDB_pathway <- msig_MSigDB_pathway   %>% dplyr::mutate(symbol_upper = toupper(.data$gene_symbol))
mapped_MSigDB_pathway <- dplyr::inner_join(ms_MSigDB_pathway, ga_MSigDB_pathway, by = "symbol_upper") 

out_mapped_MSigDB_pathway <- split(mapped_MSigDB_pathway$gene_id, mapped_MSigDB_pathway$gs_name)

# ================= Circadian rhythm ================= 
filtered_list_circadian_pathway <- out_mapped_MSigDB_pathway[grepl("CIRCADIAN", names(out_mapped_MSigDB_pathway))]

#### Developpemental process Background

filtered_list_mapped_circadian_perturbation_background <- lapply(filtered_list_circadian_pathway, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_circadian_perturbation_background))]))

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### MSigDB perturbation #####

msig_MSigDB_perturbation <- msigdbr::msigdbr(subcollection = "CGP")
sym_col_MSigDB_perturbation <- intersect(c("gene_symbol","human_gene_symbol","hgnc_symbol"), names(msig_MSigDB_perturbation))[1]

ga_MSigDB_perturbation <- gene_annot %>% dplyr::mutate(symbol_upper = toupper(.data$symbol))
ms_MSigDB_perturbation <- msig_MSigDB_perturbation   %>% dplyr::mutate(symbol_upper = toupper(.data$gene_symbol))
mapped_MSigDB_perturbation <- dplyr::inner_join(ms_MSigDB_perturbation, ga_MSigDB_perturbation, by = "symbol_upper") 

out_mapped_MSigDB_perturbation <- split(mapped_MSigDB_perturbation$gene_id, mapped_MSigDB_perturbation$gs_name)

# ================= EED ================= 
filtered_list_EED_perturbation <- out_mapped_MSigDB_perturbation[grepl("EED", names(out_mapped_MSigDB_perturbation))]

#### EED Background

filtered_list_mapped_EED_perturbation_background <- lapply(filtered_list_EED_perturbation, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_EED_perturbation_background))]))

#### EED Top 150 POS

filtered_list_mapped_EED_perturbation_Top_150_POS <- lapply(filtered_list_EED_perturbation, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_EED_perturbation_Top_150_POS))]))


#### EED Top 150 NEG

filtered_list_mapped_EED_perturbation_Top_150_NEG <- lapply(filtered_list_EED_perturbation, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_EED_perturbation_Top_150_NEG))]))


# ================= H3K27 ================= 
filtered_list_H3K27_perturbation <- out_mapped_MSigDB_perturbation[grepl("H3K27", names(out_mapped_MSigDB_perturbation))]

#### Developpemental process Background

filtered_list_mapped_H3K27_perturbation_background <- lapply(filtered_list_H3K27_perturbation, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_H3K27_perturbation_background))]))


#### Developpemental process Top 150 POS

filtered_list_mapped_H3K27_perturbation_Top_150_POS <- lapply(filtered_list_H3K27_perturbation, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_H3K27_perturbation_Top_150_POS))]))

#### Developpemental process Top 150 NEG

filtered_list_mapped_H3K27_perturbation_Top_150_NEG <- lapply(filtered_list_H3K27_perturbation, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_H3K27_perturbation_Top_150_NEG))]))


# ================= SUZ12 ================= 
filtered_list_SUZ12_perturbation <- out_mapped_MSigDB_perturbation[grepl("SUZ12", names(out_mapped_MSigDB_perturbation))]

#### SUZ12 Background

filtered_list_mapped_SUZ12_perturbation_background <- lapply(filtered_list_SUZ12_perturbation, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_SUZ12_perturbation_background))]))

#### SUZ12 Top 150 pos

filtered_list_mapped_SUZ12_perturbation_Top_150_POS <- lapply(filtered_list_SUZ12_perturbation, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_SUZ12_perturbation_Top_150_POS))]))


#### SUZ12 Top 150 neg

filtered_list_mapped_SUZ12_perturbation_Top_150_NEG <- lapply(filtered_list_SUZ12_perturbation, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_SUZ12_perturbation_Top_150_NEG))]))



# ================= PRC2 ================= 
filtered_list_PRC2_perturbation <- out_mapped_MSigDB_perturbation[grepl("PRC2", names(out_mapped_MSigDB_perturbation))]

#### PRC2 Background

filtered_list_mapped_PRC2_perturbation_background <- lapply(filtered_list_PRC2_perturbation, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_PRC2_perturbation_background))]))

#### PRC2 Top 150 POS

filtered_list_mapped_PRC2_perturbation_Top_150_POS <- lapply(filtered_list_PRC2_perturbation, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_PRC2_perturbation_Top_150_POS))]))


#### PRC2 Top 150 NEG

filtered_list_mapped_PRC2_perturbation_Top_150_NEG <- lapply(filtered_list_PRC2_perturbation, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_PRC2_perturbation_Top_150_NEG))]))



# ================= Alzheimer's disease ================= 
filtered_list_Alzheimer_perturbation <- out_mapped_MSigDB_perturbation[grepl("ALZHEIMERS", names(out_mapped_MSigDB_perturbation))]

#### Alzheimer's disease Background

filtered_list_mapped_Alzheimer_perturbation_background <- lapply(filtered_list_Alzheimer_perturbation, function(x) {
  x[x %in% genes_background]
})

df_genes_background <- as.data.frame(genes_background)

sum(table(df_genes_background$genes_background[df_genes_background$genes_background %in% unique(unlist(filtered_list_mapped_Alzheimer_perturbation_background))]))

#### Alzheimer's disease Top 150 POS

filtered_list_mapped_Alzheimer_perturbation_Top_150_POS <- lapply(filtered_list_Alzheimer_perturbation, function(x) {
  x[x %in% genes_of_interest_POS]
})

df_genes_of_interest_POS <- as.data.frame(genes_of_interest_POS)

sum(table(df_genes_of_interest_POS$genes_of_interest_POS[df_genes_of_interest_POS$genes_of_interest_POS %in% unique(unlist(filtered_list_mapped_Alzheimer_perturbation_Top_150_POS))]))



#### Alzheimer's disease Top 150 NEG

filtered_list_mapped_Alzheimer_perturbation_Top_150_NEG <- lapply(filtered_list_Alzheimer_perturbation, function(x) {
  x[x %in% genes_of_interest_NEG]
})

df_genes_of_interest_NEG <- as.data.frame(genes_of_interest_NEG)

sum(table(df_genes_of_interest_NEG$genes_of_interest_NEG[df_genes_of_interest_NEG$genes_of_interest_NEG %in% unique(unlist(filtered_list_mapped_Alzheimer_perturbation_Top_150_NEG))]))




# ================= Hypergeometric Test ================= 
N <- 13850  # Total background regions
n <- 239  # Foreground regions
K <- 406  # Background hits for the term (K_θ)
k <- 13   # Foreground hits observed (k_θ)

pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

expected <- n * K / N

fold_enrichment <- k / expected

cat("=== Test hypergéométrique foreground/background ===\n")
cat("N =", N, " | n =", n, " | K =", K, " | k =", k, "\n")
cat("Expected =", round(expected, 3), "\n")
cat("Fold enrichment =", round(fold_enrichment, 3), "\n")
cat("p-value =", format.pval(pval, digits = 3), "\n")

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# ===========================================================================================================================================================================
#### 5) Manhattan plot ####
# ===========================================================================================================================================================================
gr_genes <- GRanges(
  seqnames = gene_annot$seqid,
  ranges   = IRanges(
    start = gene_annot$start,
    end   = gene_annot$end
  ),
  strand  = gene_annot$strand
)

mcols(gr_genes) <- gene_annot %>%
  dplyr::select(
    gene_id,
    annotation,
    symbol,
    description,
    gene_biotype,
    tss
  )

hits_nearest <- distanceToNearest(gr_all, gr_genes)
gr_all_nearest_manhattan <- gr_all[queryHits(hits_nearest)]
gr_genes_nearest    <- gr_genes[subjectHits(hits_nearest)]

gr_genes_nearest$cpg_id <- gr_all_nearest_manhattan$cpg_id

df_gene <- data.frame(
 cpg_id   = gr_all_nearest_manhattan$cpg_id,      
 seqnames = as.character(seqnames(gr_genes_nearest)[queryHits(hits_nearest)]),
 start_cpg    = start(gr_all)[queryHits(hits_nearest)],
 end_cpg      = end(gr_all)[queryHits(hits_nearest)],
 gene_id  = mcols(gr_genes)$gene_id[subjectHits(hits_nearest)],
 start_gene    = as.numeric(gr_genes_nearest@ranges@start)[queryHits(hits_nearest)],
 width_gene      = as.numeric(gr_genes_nearest@ranges@width)[queryHits(hits_nearest)],
 symbol   = as.character((gr_genes_nearest$symbol)[queryHits(hits_nearest)]),
 annotation = as.character((gr_genes_nearest$annotation)[queryHits(hits_nearest)]),
 description = as.character((gr_genes_nearest$description)[queryHits(hits_nearest)])
)

df_gene$end <- df_gene$start_gene + df_gene$width_gene
df_gene$end <- df_gene$end-1
df_gene <- df_gene[,c(1:5,6,11,8,9,10)]
df_gene$distance_cpg_gene <- hits_nearest@elementMetadata@listData$distance

head(df_gene)
length(df_gene$cpg_id)

z_score <- standardScreeningNumericTrait(t(meth),meta$age)
z_score$log10P <- abs(log10(z_score$pvalueStudent))
z_score$cpg_id <- z_score$ID

z_score <- z_score %>%
  tidyr::separate(ID, into = c("CHR", "BP"), sep = "_") %>%
  dplyr::mutate(
    CHR = sub("chr(\\d+)\\.1", "\\1", CHR),
    BP = as.numeric(BP)
  ) 

z_score <- z_score %>%
  mutate(CHR = as.character(CHR),
         CHR = ifelse(CHR == "X", 23,
                      ifelse(CHR == "Y", 24, CHR)),
         CHR = as.numeric(CHR)) %>%
  arrange(CHR, BP)

chrom_lengths <- z_score %>% group_by(CHR) %>% summarize(max_bp = max(BP))
chrom_lengths <- chrom_lengths %>%
  mutate(chr_start = c(0, cumsum(max_bp)[-length(max_bp)]))

z_score <- z_score %>%
  left_join(chrom_lengths, by = "CHR") %>%
  mutate(BPcum = BP + chr_start)

z_score$sign <- ifelse(z_score$cor > 0, "POS", "NEG")

top150_annotated_manhattan <- res %>%
  mutate(abs_r = abs(r)) %>%
  filter(cpg_id %in% df_gene$cpg_id) %>%
  arrange(desc(abs_r)) %>%
  filter(p < 0.05) %>%
  slice_head(n = 800)

top150_annotated_NEG_manhattan <- top150_annotated_manhattan %>%
  filter(sign == "NEG") %>%
  slice_head(n = 150)


top150_annotated_POS_manhattan <- top150_annotated_manhattan %>%
  filter(sign == "POS")%>%
  slice_head(n = 150)

z_score <- z_score %>%
  mutate(logP = -log10(pvalueStudent),
         col = case_when(
           cpg_id %in% top150_annotated_POS_manhattan$cpg_id & sign == "POS" ~ "red",
           cpg_id %in% top150_annotated_NEG_manhattan$cpg_id & sign == "NEG" ~ "blue",
           TRUE ~ "grey70"
         ))

top_zscore <- z_score[z_score$col == "blue" | z_score$col == "red", ]

axisdf <- z_score %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2)

z_score$col[z_score$CHR %% 2 == 0 & z_score$col == "grey70"] <- "gray30"

df_gene$symbol[is.na(df_gene$symbol)] <- df_gene$gene_id[is.na(df_gene$symbol)]

z_score_annotation <- merge(z_score,df_gene,by = "cpg_id")
z_score_annotation$abs_r <- abs(z_score_annotation$cor)

ggplot(z_score_annotation) + 
  geom_histogram(mapping = aes(as.numeric(distance_cpg_gene)))

z_score_POS <- z_score_annotation[z_score_annotation$sign == "POS",]
top_z_score_POS <-  z_score_POS %>%
  arrange(desc(cor)) 
top_z_score_POS$symbol[21:5311] <- NA
top_z_score_POS$col[c(1,7)] <- "chartreuse4"

z_score_NEG <- z_score_annotation[z_score_annotation$sign == "NEG",] 
top_z_score_NEG <- z_score_NEG %>%
  arrange(desc(cor)) 
top_z_score_NEG$symbol[1:2306] <- NA
top_z_score_NEG$col[c(2326)] <- "chartreuse4"

pos_gene <- ggplot(top_z_score_POS, aes(x = BPcum, y = logP, color = col)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = symbol),size = 6, max.overlaps = 50) + 
  scale_color_identity() +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank()
  ) +
  labs(x = NULL, y = expression(-log[10](P)))  +
  geom_hline(yintercept = 10.867884, color = "purple", linetype = "dashed") 

neg_gene <- ggplot(top_z_score_NEG, aes(x = BPcum, y = logP, color = col)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = symbol),size = 6, max.overlaps = 50) + 
  scale_color_identity() +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    title = NULL
  ) +
  labs(x = "Chromosome", y = expression(-log[10](P))) + scale_y_reverse() +
  geom_hline(yintercept = 7.878529, color = "purple", linetype = "dashed") 

manhattan_plot <- ggarrange(pos_gene, neg_gene, common.legend = TRUE,nrow = 2)
manhattan_plot

##### Plot per CpGs

cpg_top1 <- meth[rownames(meth) %in% c("chr17.1_14176962"),]
cpg_top1 <- t(cpg_top1)
cpg_top1 <- as.data.frame(cpg_top1)
cpg_top1$age <- meta$age

plot_six6 <- ggplot(cpg_top1, aes(x = chr17.1_14176962, y = age)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "chr17.1_14176962 (egapx_017075 = six6)",
    y = "Age (years)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", size = 1),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.x = element_text(color = "chartreuse4", face = "bold")
  )
  

cpg_six3a <- meth[rownames(meth) %in% c("chr14.1_22159909"),]
cpg_six3a <- t(cpg_six3a)
cpg_six3a <- as.data.frame(cpg_six3a)
cpg_six3a$age <- meta$age

plot_six3a <- ggplot(cpg_six3a, aes(x = chr14.1_22159909, y = age)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x= "chr14.1_22159909 (six3a)",
    y = "Age (years)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", size = 1),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.x = element_text(color = "chartreuse4", face = "bold")
  )

cpg_top1_NEG <- meth[rownames(meth) %in% c("chr3.1_12985832"),]
cpg_top1_NEG <- t(cpg_top1_NEG)
cpg_top1_NEG <- as.data.frame(cpg_top1_NEG)
cpg_top1_NEG$age <- meta$age

plot_cpgNEG <- ggplot(cpg_top1_NEG, aes(x = chr3.1_12985832, y = age)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "chr3.1_12985832 (lmx1bb)",
    y = "Age (years)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", size = 1),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.x = element_text(color = "chartreuse4", face = "bold")
  )
  
  

manhattan_plot_final <- ggdraw() +
  # Figure A : en haut, toute la largeur, 2/3 de la hauteur
  draw_plot(manhattan_plot, x = 0, y = 1/3, width = 1, height = 2/3) +
  
  # Figure B : en bas à gauche (1/3 largeur)
  draw_plot(plot_six6, x = 0, y = 0, width = 1/3, height = 1/3) +
  
  # Figure C : au centre
  draw_plot(plot_six3a, x = 1/3, y = 0, width = 1/3, height = 1/3) +
  
  # Figure D : à droite
  draw_plot(plot_cpgNEG, x = 2/3, y = 0, width = 1/3, height = 1/3) +
  
  # Labels
  draw_plot_label(
    label = c("A", "B", "C", "D"),
    x = c(0,    0,      1/3,      2/3),
    y = c(1,    1/3,    1/3,      1/3),
    size = 15
  )
manhattan_plot_final
ggsave(file.path("C:/Users/Download/forSteve3/forSteve/", "manhattan_plot_publi_final_250kb_top150CpGs.png"), manhattan_plot_final, width=24, height=12, dpi=600)

# ===========================================================================================================================================================================
#### 6 Enrichment plot ####
final_data <- read.csv2("D:/UCLA_Matteo/Analyse_Thomas/MammalMethyl/EWAS/Summary_table_EWAS_tss_corrigé_150CpGs_multiple_annotation.csv")
final_data <- final_data[c(-8),]

final_data$Annotation <- factor(final_data$Annotation, levels = rev(unique(final_data$Annotation)))

final_data$Tissue <- "Muscle"

final_data$p.value_POS <- as.numeric(final_data$p.value_POS)
final_data$p.value_NEG <- as.numeric(final_data$p.value_NEG)

final_data$logP_p.value_POS <- -log10(final_data$p.value_POS)
final_data$logP_p.value_NEG <- -log10(final_data$p.value_NEG)

final_data$p.value_POS_exp <- format(final_data$p.value_POS, scientific = TRUE)
final_data$p.value_NEG_exp <- format(final_data$p.value_NEG, scientific = TRUE)

final_data$LabelPOS <- paste0(final_data$p.value_POS_exp,"\n(", final_data$Top150cpgsPOS,")")

p1 <- ggplot(final_data, aes(x = Annotation, y = Background, fill = Enrichment_analysis)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "", y = "Annotation size",fill = "Enrichtement analysis:") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.position = "bottom"
  ) + scale_fill_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) 

p2_250cpgsPOS <- ggplot(final_data, aes(x = Tissue, y = Annotation, fill = logP_p.value_POS)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = expression(-log[10](p))
  ) +   labs(x = "", y = "") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.position = "right"
  ) + geom_text(aes(label = LabelPOS), size = 5) 

final_data$LabelNEG <- paste0(final_data$p.value_NEG_exp,"\n(", final_data$Top150cpgsNEG ,")")

p2_250cpgsNEG <- ggplot(final_data, aes(x = Tissue, y = Annotation, fill = logP_p.value_NEG)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = expression(-log[10](p))
  ) +   labs(x = "", y = "") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.position = "right"
  ) + geom_text(aes(label = LabelNEG), size = 5) 


final_data$logP_p.value_NEG[6] <- 0
final_data$logP_p.value_NEG[12] <- 0

maxLogP <- max(c(final_data$logP_p.value_POS, final_data$logP_p.value_NEG), na.rm = TRUE)

heat_pal <- scale_fill_gradient(
  low = "white",
  high = "red",
  limits = c(0, maxLogP),
  name = expression(-log[10](p))
)

p2_250cpgsPOS <- ggplot(final_data, aes(x = Tissue, y = Annotation, fill = logP_p.value_POS)) +
  geom_tile(color = "grey90", linewidth = 0.3) +
  heat_pal +
  geom_text(
    data = subset(final_data, p.value_POS < 0.05),
    aes(label = LabelPOS),
    size = 4.8, fontface = "bold", color = "black"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",     # <<< pas de légende ici
    plot.margin = margin(0,10,0,10)
  )


p2_250cpgsNEG <- ggplot(final_data, aes(x = Tissue, y = Annotation, fill = logP_p.value_NEG)) +
  geom_tile(color = "grey90", linewidth = 0.3) +
  heat_pal +
  geom_text(
    data = subset(final_data, p.value_NEG < 0.05),
    aes(label = LabelNEG),
    size = 4.8, fontface = "bold"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 15),
    legend.position = "right",
    plot.margin = margin(0,10,0,10)
  ) +
  xlab("Muscle")   # <<< un seul label clair et propre


final_plot <- p1 + p2_250cpgsPOS + p2_250cpgsNEG +
  plot_layout(widths = c(4.5, 1.2, 1.2)) &
  theme(plot.background = element_rect(fill="white", color=NA))

final_plot

ggsave(file.path("C:/Users/Download/forSteve3/forSteve/", "enrichment_plot_top150POS_top150NEG_250kb_TSS_correct_multi_annotation.png"), final_plot, width=16, height=8, dpi=600)
# ===========================================================================================================================================================================






