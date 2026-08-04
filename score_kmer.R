#!/usr/bin/env Rscript
# =============================================================================
# Score biologique — VERSION K-MER (projet Chromoanagenesis-AML)
# -----------------------------------------------------------------------------
# Portage du "score biologique" (repo Score-biologique, basé sur Arriba) vers
# une version qui fonctionne à partir de fichiers PAUVRES en information : les
# fichiers k-mer. On n'a PLUS les sorties Arriba : on RECONSTRUIT les infos
# manquantes (strand -> direction -> type/classe chimérique) et on remplace les
# composantes Arriba par des équivalents k-mer.
#
# ---------------------------------------------------------------------------
# ENTRÉES
#   1. Fichier(s) k-mer patho (TSV, sans en-tête, 3 colonnes) :
#        <gene5'_chr5'_bp5'_gene3'_chr3'_bp3'_index|index...>  <patient>  <comptage>
#      ex :  ACTB_chr7_5570155_PTMAP5_chr13_82264248_813|2628   BA2816R   7.19
#      (le 1er gène est en 5', le 2e en 3' ; une ligne = un k-mer)
#   2. Blacklist normaux (kmer_score/kmer_normaux_33.tsv) : MÊME format que (1),
#      comptages k-mer chez ~33 échantillons normaux.
#   3. Deux BED : left = gène 5', right = gène 3' :
#        chr  start  end  <gene_chr_bp_index...>  .  strand
#      -> apportent le STRAND (+/-) de chaque gène.
#   4. Annotation GENCODE (hg38) : téléchargée UNE FOIS puis mise en CACHE
#      (--gtf pour fournir un fichier local et éviter le réseau).
#
# ---------------------------------------------------------------------------
# RECONSTRUCTIONS
#   strand -> direction :
#       gène 5' / +  -> downstream      gène 5' / -  -> upstream
#       gène 3' / +  -> upstream        gène 3' / -  -> downstream
#   direction + chr + distance -> type (classe chimérique) :
#       chr différents           -> Translocation (Classe 1)
#       même chr, dir. égales    -> Inversion     (Classe 4)
#       même chr, dir. opposées  -> Délétion (Classe 2) ou Duplication (Classe 3)
#                                   selon l'ordre des breakpoints ; Read-through
#                                   (Classe 2) si délétion colinéaire même brin < seuil kb.
#
# ---------------------------------------------------------------------------
# SCORE  =  Σ(fraction × poids) / Σ(poids)     (max 11 pts par défaut)
#   spécificité   (poids 4) : 2/3·(1 − présence_normaux) + 1/3·expr_patho
#                             présence_normaux = moyenne(freq_norm, val_norm)
#                             (absence des normaux dominante ; expr patho module)
#   type          (poids 3) : Transloc/Inversion 1 · Délét/Duplic 0.5 · RT 0.25 · sinon 0
#   reading_frame (poids 2) : moyenne des 2 breakpoints ; par breakpoint :
#                             bordure exacte 1 · <5nt de la bordure 0.9 · CDS 0.5 · reste 0
#   WHO           (poids 2) : fusion d'intérêt WHO -> 1 · sinon 0
#   (confidence & reads Arriba : supprimées ; l'expression patho vit dans spéc.)
#
# Agrégation des comptages = MOYENNE par fusion (tous k-mers × échantillons
# positifs ; les absents ne diluent pas). Vaut pour patho ET normaux.
#
# Priorités : P1 ≥ 65 % · P2 ≥ 40 % · P3 ≥ 20 % · sinon NP.
# Un poids à 0 retire la composante du calcul ET du dénominateur.
#
# ---------------------------------------------------------------------------
# USAGE
#   Rscript score_kmer.R \
#       --dir-kmer  kmer_score/patho \
#       --blacklist kmer_score/kmer_normaux_33.tsv \
#       --bed-left  kmer_score/left.bed \
#       --bed-right kmer_score/right.bed \
#       [--gtf annot_cache/gencode.v45.annotation.gtf.gz] \
#       --dir-out   analyse_fusions_kmer
#   Rscript score_kmer.R --help
# =============================================================================

# locale UTF-8 : rendu correct des accents et du "²" (sinon locale C -> "..")
for (loc in c("C.UTF-8", "en_US.UTF-8", "fr_FR.UTF-8", "C.utf8")) {
  if (suppressWarnings(Sys.setlocale("LC_CTYPE", loc)) != "") break
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})
# device cairo si dispo : rendu correct des accents et du "²" dans les figures
if (isTRUE(capabilities("cairo"))) options(bitmapType = "cairo")
# figures en PDF vectoriel ; cairo_pdf gère l'UTF-8 (accents, "²"), sinon pdf standard
PDF_DEV <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf

# ── CONFIG PAR DÉFAUT ────────────────────────────────────────────────────────
opt <- list(
  dir_kmer   = "kmer_score/patho",                 # dossier des *.tsv patho
  blacklist  = "kmer_score/kmer_normaux_33.tsv",   # blacklist normaux
  bed_left   = "kmer_score/left.bed",              # BED gène 5'
  bed_right  = "kmer_score/right.bed",             # BED gène 3'
  gtf        = "",                                  # GTF local (sinon téléchargé)
  gtf_url    = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz",
  cache_dir  = "annot_cache",                       # cache de l'annotation
  dir_out    = "analyse_fusions_kmer",
  n_top      = 30,                                  # figures
  rt_kb      = 300,                                 # seuil read-through (kb)
  bp_tol     = 5,                                   # tolérance bordure d'exon (nt)
  w_type = 3, w_frame = 2, w_spec = 4, w_who = 2   # spéc. > type (surchargeable en CLI)
)

# ── PARSER CLI ───────────────────────────────────────────────────────────────
alias <- c(type = "w_type", frame = "w_frame", spec = "w_spec", who = "w_who")
string_opts <- c("dir_kmer", "blacklist", "bed_left", "bed_right",
                 "gtf", "gtf_url", "cache_dir", "dir_out")
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || "-h" %in% args) {
  cat("Score biologique k-mer.\n",
      "Chemins  : --dir-kmer --blacklist --bed-left --bed-right --gtf --gtf-url",
      "--cache-dir --dir-out\n",
      "Poids    : --type (4) --frame (2) --spec (2) --who (2)  [0 = retirée]\n",
      "Autres   : --n-top (30) --rt-kb (300) --bp-tol (5)\n",
      "Voir l'entête du script pour le détail.\n")
  quit(status = 0)
}
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (grepl("^--", a)) {
    a2 <- sub("^--", "", a)
    if (grepl("=", a2)) {
      kv <- strsplit(a2, "=", fixed = TRUE)[[1]]; key <- kv[1]; val <- kv[2]
    } else {
      key <- a2
      if (i < length(args) && !grepl("^--", args[i + 1])) { val <- args[i + 1]; i <- i + 1 }
      else val <- "TRUE"
    }
    key <- gsub("-", "_", key)
    if (key %in% names(alias)) key <- alias[[key]]
    if (key %in% names(opt))
      opt[[key]] <- if (key %in% string_opts) val else suppressWarnings(as.numeric(val))
    else warning("Option inconnue ignorée : --", key)
  }
  i <- i + 1
}

DIR_OUT <- opt$dir_out
DIR_FIG <- file.path(DIR_OUT, "figures")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)
theme_set(theme_bw(base_size = 12))

cat("=== CONFIGURATION (score k-mer) ===\n")
cat("  k-mer patho :", opt$dir_kmer, "\n")
cat("  blacklist   :", opt$blacklist, "\n")
cat("  BED 5' / 3' :", opt$bed_left, "|", opt$bed_right, "\n")
cat("  GTF local   :", ifelse(nzchar(opt$gtf), opt$gtf, "(téléchargé + caché)"), "\n")
cat("  sortie      :", normalizePath(opt$dir_out, mustWork = FALSE), "\n")
cat(sprintf("  Poids : type=%s frame=%s spec=%s who=%s | RT<%skb | tol=%snt\n\n",
            opt$w_type, opt$w_frame, opt$w_spec, opt$w_who, opt$rt_kb, opt$bp_tol))

# ── Fusions WHO d'intérêt + alias HGNC (repris du score original) ─────────────
WHO_FUSIONS_INTEREST <- c(
  "RUNX1--RUNX1T1", "CBFB--MYH11", "PML--RARA", "KMT2A--MLLT3",
  "DEK--NUP214", "NUP98--NSD1", "BCR--ABL1", "NRIP1--MIR99AHG"
)
GENE_ALIASES <- c(
  MLL = "KMT2A", HRX = "KMT2A", ALL1 = "KMT2A", MLL1 = "KMT2A", TRX1 = "KMT2A",
  AF9 = "MLLT3", LTG9 = "MLLT3", AF6 = "MLLT4", AFDN = "MLLT4", ENL = "MLLT1",
  ETO = "RUNX1T1", MTG8 = "RUNX1T1", CBFA2T1 = "RUNX1T1", BCR1 = "BCR"
)
resolve_alias <- function(g) ifelse(g %in% names(GENE_ALIASES), GENE_ALIASES[g], g)

# min-max robuste : renvoie 0.5 si tout est égal (évite NaN / tout-à-zéro)
minmax <- function(x) {
  x <- as.numeric(x); r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || diff(r) == 0) return(rep(0.5, length(x)))
  (x - r[1]) / diff(r)
}
# variante pour les comptages normaux : si tout est égal (ex. aucune fusion
# blacklistée), renvoie 0 -> "absent des normaux" = spécifique (et non 0.5 neutre)
minmax0 <- function(x) {
  x <- as.numeric(x); r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

# ── Parsing des identifiants ─────────────────────────────────────────────────
# fusion_id : gene5_chr5_bp5_gene3_chr3_bp3[_index|index...]
FUSION_RE <- "^(.+?)_(chr[0-9A-Za-z]+)_([0-9]+)_(.+?)_(chr[0-9A-Za-z]+)_([0-9]+)(?:_([0-9|]+))?$"
parse_fusion_id <- function(id) {
  m <- str_match(id, FUSION_RE)
  tibble(gene5 = m[, 2], chr5 = m[, 3], bp5 = suppressWarnings(as.numeric(m[, 4])),
         gene3 = m[, 5], chr3 = m[, 6], bp3 = suppressWarnings(as.numeric(m[, 7])),
         index = m[, 8])
}
# nom BED : gene_chr_bp_index... -> gene, chr, bp
BED_RE <- "^(.+?)_(chr[0-9A-Za-z]+)_([0-9]+)"
parse_bed_name <- function(nm) {
  m <- str_match(nm, BED_RE)
  tibble(gene = m[, 2], chr = m[, 3], bp = suppressWarnings(as.numeric(m[, 4])))
}

# ── 1. LECTURE DES FICHIERS K-MER (patho) ────────────────────────────────────
read_kmer_tsv <- function(path) {
  read_tsv(path, col_names = c("fusion_id", "sample", "count"),
           col_types = cols(fusion_id = col_character(), sample = col_character(),
                            count = col_double()), progress = FALSE)
}
kmer_files <- list.files(opt$dir_kmer, pattern = "\\.tsv$", full.names = TRUE)
if (length(kmer_files) == 0) {
  # tolère aussi un fichier unique passé en --dir-kmer
  if (file.exists(opt$dir_kmer) && !dir.exists(opt$dir_kmer)) kmer_files <- opt$dir_kmer
  else stop("Aucun fichier .tsv k-mer trouvé.\n  --dir-kmer = '", opt$dir_kmer,
            "'\n  cherché (absolu) : ", normalizePath(opt$dir_kmer, mustWork = FALSE),
            "\n  dossier courant  : ", getwd(),
            "\n  -> vérifie le chemin (absolu conseillé) : le fichier/dossier n'existe pas ici.")
}
cat(length(kmer_files), "fichier(s) k-mer patho\n")
patho_raw <- map_dfr(kmer_files, read_kmer_tsv) %>% filter(!is.na(fusion_id))

# identité d'une fusion = tout SAUF l'index k-mer et le patient
patho <- bind_cols(patho_raw, parse_fusion_id(patho_raw$fusion_id)) %>%
  filter(!is.na(bp5), !is.na(bp3)) %>%
  mutate(fusion_core = paste(gene5, chr5, bp5, gene3, chr3, bp3, sep = "_"))
if (nrow(patho) == 0) stop("Aucune ligne k-mer parsable (format fusion_id ?).")

# clé de matching blacklist : gene5_bp5_gene3_bp3 (avec breakpoints, SANS chr ni index)
mk_bl_key <- function(g5, b5, g3, b3) paste(g5, b5, g3, b3, sep = "_")

# ── 2. AGRÉGATION MOYENNE PAR FUSION (patho) ─────────────────────────────────
# comptage fusion = moyenne des comptages bruts sur tous ses k-mers × patients positifs
patho_agg <- patho %>%
  group_by(fusion_core) %>%
  summarise(
    gene5 = first(gene5), chr5 = first(chr5), bp5 = first(bp5),
    gene3 = first(gene3), chr3 = first(chr3), bp3 = first(bp3),
    mean_count_patho = mean(count, na.rm = TRUE),
    n_patients_pos   = n_distinct(sample),
    n_kmers          = n_distinct(coalesce(index, "")),
    .groups = "drop") %>%
  mutate(bl_key = mk_bl_key(gene5, bp5, gene3, bp3))
cat(nrow(patho_agg), "fusions uniques (agrégation MOYENNE)\n")

# ── 3. BLACKLIST NORMAUX -> fréquence + valeur ───────────────────────────────
if (file.exists(opt$blacklist)) {
  bl_raw <- read_kmer_tsv(opt$blacklist) %>% filter(!is.na(fusion_id))
  bl <- bind_cols(bl_raw, parse_fusion_id(bl_raw$fusion_id)) %>%
    filter(!is.na(bp5), !is.na(bp3)) %>%
    mutate(bl_key = mk_bl_key(gene5, bp5, gene3, bp3)) %>%
    group_by(bl_key) %>%
    summarise(mean_count_normaux = mean(count, na.rm = TRUE),
              n_normaux_pos = n_distinct(sample), .groups = "drop")
  n_normaux_total <- suppressWarnings(as.integer(str_match(basename(opt$blacklist), "([0-9]+)")[, 2]))
  if (is.na(n_normaux_total) || n_normaux_total == 0) {
    n_normaux_total <- n_distinct(bl_raw$sample)
  }
  n_normaux_total <- max(n_normaux_total, 1L)   # évite une division par zéro
  cat(nrow(bl), "fusions dans la blacklist (", n_normaux_total, "normaux)\n")
} else {
  warning("Blacklist introuvable : ", opt$blacklist, " -> spécificité = absence supposée")
  bl <- tibble(bl_key = character(), mean_count_normaux = double(), n_normaux_pos = integer())
  n_normaux_total <- 1L
}

# ── 4. BED -> STRAND de chaque gène ──────────────────────────────────────────
read_bed_strand <- function(path) {
  if (!file.exists(path)) { warning("BED introuvable : ", path); return(tibble()) }
  b <- read_tsv(path, col_names = c("chr", "start", "end", "name", "score", "strand"),
                col_types = cols(.default = col_character()), progress = FALSE)
  bind_cols(b["strand"], parse_bed_name(b$name)) %>%
    filter(!is.na(bp)) %>% distinct(gene, chr, bp, .keep_all = TRUE) %>%
    transmute(gene, chr, bp, strand)
}
strand_left  <- read_bed_strand(opt$bed_left)   # gène 5'
strand_right <- read_bed_strand(opt$bed_right)  # gène 3'

lookup_strand <- function(tbl, gene, chr, bp) {
  if (nrow(tbl) == 0) return(rep(NA_character_, length(gene)))
  key <- paste(gene, chr, bp, sep = "_")
  ref <- setNames(tbl$strand, paste(tbl$gene, tbl$chr, tbl$bp, sep = "_"))
  unname(ref[key])
}

# ── 5. STRAND -> DIRECTION -> TYPE / CLASSE ──────────────────────────────────
# direction du partenaire : gène 5'/+ -> downstream, 5'/- -> upstream ;
#                           gène 3'/+ -> upstream,   3'/- -> downstream.
direction_5p <- function(strand) case_when(strand == "+" ~ "downstream",
                                           strand == "-" ~ "upstream", TRUE ~ NA_character_)
direction_3p <- function(strand) case_when(strand == "+" ~ "upstream",
                                           strand == "-" ~ "downstream", TRUE ~ NA_character_)

# type Arriba reconstruit à partir de (chr, direction, position, brin)
reconstruct_type <- function(chr5, bp5, dir5, str5, chr3, bp3, dir3, str3, rt_kb) {
  n <- length(chr5); out <- rep(NA_character_, n)
  for (k in seq_len(n)) {
    if (chr5[k] != chr3[k]) { out[k] <- "Translocation"; next }        # chr différents
    if (is.na(dir5[k]) || is.na(dir3[k])) { out[k] <- "Indéterminé"; next }
    if (dir5[k] == dir3[k]) { out[k] <- "Inversion"; next }            # dir. égales
    # directions opposées -> délétion ou duplication selon l'ordre des breakpoints
    del <- (dir5[k] == "downstream" && bp5[k] <  bp3[k]) ||
           (dir5[k] == "upstream"   && bp5[k] >  bp3[k])
    if (del) {
      colinear <- !is.na(str5[k]) && !is.na(str3[k]) && str5[k] == str3[k]
      dist_kb  <- abs(bp5[k] - bp3[k]) / 1000
      out[k] <- if (colinear && dist_kb < rt_kb) "Read-through" else "Délétion"
    } else out[k] <- "Duplication"
  }
  out
}
classe_chimerique <- function(type_base) case_when(
  type_base == "Translocation" ~ "Classe 1",
  type_base %in% c("Délétion", "Read-through") ~ "Classe 2",  # read-through = délétion particulière
  type_base == "Duplication"   ~ "Classe 3", type_base == "Inversion" ~ "Classe 4",
  TRUE ~ NA_character_)
TYPE_COLORS <- c(Translocation = "#d62728", Inversion = "#9467bd", "Délétion" = "#ff7f0e",
                 Duplication = "#1f77b4", "Read-through" = "#2ca02c",
                 "Indéterminé" = "grey65")

# ── 6. ANNOTATION GENCODE (téléchargement + cache) -> contexte des breakpoints ─
get_gtf_path <- function() {
  if (nzchar(opt$gtf) && file.exists(opt$gtf)) return(opt$gtf)
  dir.create(opt$cache_dir, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(opt$cache_dir, basename(opt$gtf_url))
  if (file.exists(dest) && file.size(dest) > 0) { cat("GTF (cache) :", dest, "\n"); return(dest) }
  cat("Téléchargement GENCODE ->", dest, "(une seule fois)...\n")
  ok <- tryCatch({ download.file(opt$gtf_url, dest, mode = "wb", quiet = TRUE); TRUE },
                 error = function(e) { warning("Échec téléchargement GTF : ", conditionMessage(e)); FALSE })
  if (ok && file.exists(dest) && file.size(dest) > 0) dest else NA_character_
}

# Charge, pour l'ensemble de gènes voulu, les bordures d'exons et les intervalles CDS.
load_annotation <- function(genes) {
  genes <- unique(genes[!is.na(genes)])
  empty <- list(exon = tibble(gene = character(), chr = character(), pos = double()),
                cds  = tibble(gene = character(), chr = character(), start = double(), end = double()))
  gtf <- get_gtf_path()
  if (is.na(gtf) || !file.exists(gtf)) { warning("Pas d'annotation -> reading_frame neutre (0.5)"); return(NULL) }
  gf <- tempfile(); writeLines(genes, gf)
  # pré-filtre shell : exon/CDS + gènes voulus (zcat -f gère .gz et le plein texte)
  cmd <- sprintf("zcat -f %s | awk -F'\t' '$3==\"exon\"||$3==\"CDS\"' | grep -F -f %s",
                 shQuote(gtf), shQuote(gf))
  raw <- tryCatch(suppressWarnings(readLines(pipe(cmd))), error = function(e) character(0))
  if (length(raw) == 0) { warning("Annotation vide après filtre -> reading_frame neutre"); return(NULL) }
  parts <- str_split_fixed(raw, "\t", 9)
  gn <- str_match(parts[, 9], 'gene_name "([^"]+)"')[, 2]
  df <- tibble(chr = parts[, 1], feature = parts[, 3],
               start = suppressWarnings(as.numeric(parts[, 4])),
               end = suppressWarnings(as.numeric(parts[, 5])), gene = gn) %>%
    filter(gene %in% genes)                              # match EXACT du gene_name
  exons <- df %>% filter(feature == "exon") %>%
    transmute(gene, chr, pos = start) %>%
    bind_rows(df %>% filter(feature == "exon") %>% transmute(gene, chr, pos = end)) %>%
    distinct()
  cds <- df %>% filter(feature == "CDS") %>% transmute(gene, chr, start, end) %>% distinct()
  list(exon = exons, cds = cds)
}

# contexte d'un breakpoint, nuancé selon la distance à la bordure d'exon :
#   "exon_boundary" = match PARFAIT (dist 0) ; "exon_near" = à <= tol nt d'une
#   bordure ; sinon "CDS" (codant) ou "other" (intron/UTR/intergénique).
bp_context <- function(gene, chr, bp, annot, tol) {
  n <- length(bp); out <- rep("other", n)
  if (is.null(annot)) return(rep(NA_character_, n))       # -> fraction neutre plus loin
  norm_chr <- function(x) sub("^chr", "", x)             # GENCODE "chr22" vs Ensembl "22"
  for (k in seq_len(n)) {
    g <- gene[k]; c <- norm_chr(chr[k]); p <- bp[k]
    if (is.na(g) || is.na(p)) { out[k] <- NA_character_; next }
    ex <- annot$exon %>% filter(gene == g, norm_chr(chr) == c)
    if (nrow(ex) > 0) {
      d <- min(abs(ex$pos - p))
      if (d == 0)   { out[k] <- "exon_boundary"; next }   # bordure exacte
      if (d <= tol) { out[k] <- "exon_near";     next }   # à <= tol nt de la bordure
    }
    cd <- annot$cds %>% filter(gene == g, norm_chr(chr) == c)
    if (nrow(cd) > 0 && any(cd$start <= p & p <= cd$end)) { out[k] <- "CDS"; next }
  }
  out
}
# fractions : bordure exacte 1.0 · proche (<= tol) 0.9 · CDS 0.5 · reste 0 · NA neutre
ctx_frac <- function(ctx) case_when(
  ctx == "exon_boundary" ~ 1.0, ctx == "exon_near" ~ 0.9,
  ctx == "CDS" ~ 0.5, ctx == "other" ~ 0.0, TRUE ~ 0.5)

# ── 7. ASSEMBLAGE DES CARACTÉRISTIQUES ───────────────────────────────────────
feat <- patho_agg %>%
  left_join(bl, by = "bl_key") %>%
  mutate(
    mean_count_normaux = coalesce(mean_count_normaux, 0),
    n_normaux_pos      = coalesce(n_normaux_pos, 0L),
    g5n = resolve_alias(gene5), g3n = resolve_alias(gene3),
    strand5 = lookup_strand(strand_left,  gene5, chr5, bp5),
    strand3 = lookup_strand(strand_right, gene3, chr3, bp3),
    direction5 = direction_5p(strand5),
    direction3 = direction_3p(strand3),
    type_base  = reconstruct_type(chr5, bp5, direction5, strand5,
                                  chr3, bp3, direction3, strand3, opt$rt_kb),
    classe_chimerique = classe_chimerique(type_base),
    distance_bp  = if_else(chr5 == chr3, abs(bp5 - bp3), NA_real_),
    is_who = paste(g5n, g3n, sep = "--") %in% WHO_FUSIONS_INTEREST |
             paste(g3n, g5n, sep = "--") %in% WHO_FUSIONS_INTEREST
  )

# reading_frame : contexte des 2 breakpoints via GENCODE
annot <- load_annotation(c(feat$gene5, feat$gene3))
feat <- feat %>%
  mutate(ctx5 = bp_context(gene5, chr5, bp5, annot, opt$bp_tol),
         ctx3 = bp_context(gene3, chr3, bp3, annot, opt$bp_tol),
         frac_frame = (ctx_frac(ctx5) + ctx_frac(ctx3)) / 2)

# ── 8. SCORE PONDÉRÉ ─────────────────────────────────────────────────────────
MAX_SCORE <- with(opt, w_type + w_frame + w_spec + w_who)
if (MAX_SCORE <= 0) stop("Tous les poids sont nuls : score impossible.")

feat <- feat %>%
  mutate(
    # composante type
    frac_type = case_when(
      type_base %in% c("Translocation", "Inversion") ~ 1.0,
      type_base %in% c("Délétion", "Duplication")    ~ 0.5,
      type_base == "Read-through"                    ~ 0.25,
      TRUE ~ 0),
    # composante spécificité (contraste patho vs normaux, absence dominante)
    expr_patho       = minmax(log1p(mean_count_patho)),
    val_norm         = minmax0(log1p(mean_count_normaux)),
    freq_norm        = pmin(n_normaux_pos / n_normaux_total, 1),
    presence_normaux = (freq_norm + val_norm) / 2,
    frac_spec = (2/3) * (1 - presence_normaux) + (1/3) * expr_patho,
    # composante WHO
    frac_who = if_else(is_who, 1.0, 0),
    # score
    score_type  = frac_type  * opt$w_type,
    score_frame = frac_frame * opt$w_frame,
    score_spec  = frac_spec  * opt$w_spec,
    score_who   = frac_who   * opt$w_who,
    score_total = score_type + score_frame + score_spec + score_who,
    score_norm  = score_total / MAX_SCORE,
    priorite = factor(case_when(
      score_norm >= 0.65 ~ "P1", score_norm >= 0.40 ~ "P2",
      score_norm >= 0.20 ~ "P3", TRUE ~ "NP"), levels = c("P1", "P2", "P3", "NP"))
  ) %>%
  mutate(gene_pair = paste(g5n, g3n, sep = "--"))

# ── 8b. UNE FUSION = UNE PAIRE DE GÈNES (on garde la variante au meilleur score) ─
# Une même paire de gènes peut apparaître avec plusieurs breakpoints (plusieurs
# k-mers) : on ne la compte qu'UNE fois, en retenant la variante au score le plus
# haut. feat_all conserve toutes les variantes pour la traçabilité.
feat_all <- feat
feat <- feat_all %>%
  group_by(gene_pair) %>%
  slice_max(order_by = score_norm, n = 1, with_ties = FALSE) %>%
  ungroup()
cat(nrow(feat), "fusions après déduplication par paire de gènes (sur",
    nrow(feat_all), "variantes de breakpoints)\n")

# ── 8c. STATS PAR PATIENT (pour les cartes de priorisation) ──────────────────
# expression d'une fusion chez un patient = moyenne de ses k-mers chez ce patient.
# Concentration = max / somme = part de l'expression portée par le patient
# principal (100 % = tout chez un seul patient ; petit = étalé sur beaucoup).
per_pat <- patho %>%
  mutate(gene_pair = paste(resolve_alias(gene5), resolve_alias(gene3), sep = "--")) %>%
  group_by(gene_pair, sample) %>%
  summarise(expr = mean(count, na.rm = TRUE), .groups = "drop")
foc <- per_pat %>%
  group_by(gene_pair) %>%
  summarise(max_patient = max(expr, na.rm = TRUE),
            sum_patient = sum(expr, na.rm = TRUE),
            n_pat_pos   = sum(expr > 0), .groups = "drop") %>%
  mutate(concentration = ifelse(sum_patient > 0, max_patient / sum_patient, 0),
         repartition = factor(case_when(
           n_pat_pos <= 1 ~ "1 patient (unique)",
           n_pat_pos <= 3 ~ "2–3 patients (restreinte)",
           TRUE           ~ "≥ 4 patients (diffuse)"),
           levels = c("1 patient (unique)", "2–3 patients (restreinte)", "≥ 4 patients (diffuse)")))
feat <- feat %>% left_join(foc, by = "gene_pair")

# ── 9. SORTIE TABLE ──────────────────────────────────────────────────────────
out_cols <- c(
  "gene_pair",
  "fusion_core", "gene5", "chr5", "bp5", "gene3", "chr3", "bp3",
  "strand5", "strand3", "direction5", "direction3",
  "type_base", "classe_chimerique", "distance_bp", "is_who",
  "mean_count_patho", "n_patients_pos", "n_kmers", "expr_patho",
  "n_normaux_pos", "freq_norm", "mean_count_normaux", "val_norm", "presence_normaux",
  "max_patient", "n_pat_pos", "concentration", "repartition",
  "ctx5", "ctx3", "frac_frame",
  "frac_type", "frac_spec", "frac_who",
  "score_type", "score_frame", "score_spec", "score_who",
  "score_total", "score_norm", "priorite")
res <- feat %>% select(any_of(out_cols)) %>% arrange(desc(score_norm), gene_pair)
write_tsv(res, file.path(DIR_OUT, "fusions_score_kmer.tsv"))

# toutes les variantes de breakpoints (avant déduplication) — traçabilité
res_all <- feat_all %>% select(any_of(out_cols)) %>%
  arrange(gene_pair, desc(score_norm), fusion_core)
write_tsv(res_all, file.path(DIR_OUT, "fusions_score_kmer_all_variants.tsv"))

# fusions "chromo-spécifiques" : absentes des normaux (présence = 0)
res_spec <- res %>% filter(presence_normaux == 0)
write_tsv(res_spec, file.path(DIR_OUT, "fusions_specifiques_kmer.tsv"))

# ── 10. FIGURES ──────────────────────────────────────────────────────────────
PRIO_COLORS <- c(P1 = "#d62728", P2 = "#ff7f0e", P3 = "#9467bd", NP = "#bdbdbd")
N_TOP <- opt$n_top
fig_df <- feat %>% mutate(fusion_label = paste(gene5, gene3, sep = "--"),
                          type_base = factor(type_base, levels = names(TYPE_COLORS)))

# 10a. Classement par score biologique
p_rank <- fig_df %>% slice_max(score_norm, n = N_TOP, with_ties = FALSE) %>%
  mutate(fusion_ord = reorder(fusion_label, score_norm),
         etiq = paste0(coalesce(classe_chimerique, "—"), " · ",
                       coalesce(recode(ctx5, exon_boundary = "exon", exon_near = "exon≈",
                                       CDS = "cds", other = "—"), "—"))) %>%
  ggplot(aes(score_norm, fusion_ord, fill = priorite)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = c(0.65, 0.40, 0.20), linetype = c("dashed", "dotted", "dotdash"),
             color = "grey50", linewidth = 0.3) +
  geom_text(aes(label = etiq), hjust = -0.05, size = 2.2, color = "grey25") +
  scale_fill_manual(values = PRIO_COLORS, drop = FALSE, name = "Priorité") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.18),
                     breaks = seq(0, 1, 0.2)) +
  labs(title = paste0("Top ", N_TOP, " fusions par score biologique (k-mer)"),
       subtitle = paste0("Score / ", MAX_SCORE, " pts  |  seuils P1 65% · P2 40% · P3 20%"),
       x = "Score normalisé", y = NULL) +
  theme(axis.text.y = element_text(size = 7), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_FIG, "score_classement.pdf"), p_rank, width = 12, height = 9, device = PDF_DEV)

# 10b. Répartition des types chimériques reconstruits
p_type <- fig_df %>% count(type_base, classe_chimerique, name = "n") %>%
  ggplot(aes(n, fct_reorder(type_base, n, sum), fill = classe_chimerique)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "Classe chimérique", na.value = "grey70") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Types chimériques reconstruits (strand + direction + distance)",
       x = "Nombre de fusions", y = NULL) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_FIG, "repartition_types.pdf"), p_type, width = 9, height = 5, device = PDF_DEV)

# 10c. Décomposition du score (composantes actives)
comp_def <- tibble::tribble(
  ~col, ~label, ~weight,
  "score_type",  "Type",         opt$w_type,
  "score_frame", "Reading frame", opt$w_frame,
  "score_spec",  "Spécificité",  opt$w_spec,
  "score_who",   "WHO",          opt$w_who) %>% filter(weight > 0)
dec_df <- fig_df %>% slice_max(score_norm, n = N_TOP, with_ties = FALSE) %>%
  mutate(fusion_ord = reorder(fusion_label, score_norm)) %>%
  pivot_longer(all_of(comp_def$col), names_to = "comp", values_to = "val") %>%
  mutate(comp = factor(setNames(comp_def$label, comp_def$col)[comp], levels = comp_def$label))
p_dec <- ggplot(dec_df, aes(val, fusion_ord, fill = comp)) +
  geom_col(width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(title = "Décomposition du score biologique (k-mer)",
       x = paste0("Points cumulés (max = ", MAX_SCORE, ")"), y = NULL) +
  theme(axis.text.y = element_text(size = 7), legend.position = "bottom")
ggsave(file.path(DIR_FIG, "score_decomposition.pdf"), p_dec, width = 11, height = 9, device = PDF_DEV)

# ── 10d/10e. CARTES DE PRIORISATION (score × expression max par patient) ──────
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
REPART_COLORS <- c("1 patient (unique)" = "#d62728",
                   "2–3 patients (restreinte)" = "#ff7f0e",
                   "≥ 4 patients (diffuse)" = "#1f77b4")
map_df <- fig_df %>% filter(!is.na(max_patient), max_patient > 0)
# labels : fusions notables (zone P1 OU expression dans le top 10 %)
lab_df <- if (nrow(map_df) > 0) map_df %>%
  filter(score_norm >= 0.65 |
         max_patient >= quantile(max_patient, 0.90, na.rm = TRUE)) else map_df
add_labels <- function(p) {
  if (nrow(lab_df) == 0) return(p)
  if (has_repel)
    p + ggrepel::geom_text_repel(data = lab_df, aes(label = fusion_label), size = 2.4,
          color = "grey20", max.overlaps = 30, min.segment.length = 0, seed = 1)
  else
    p + geom_text(data = lab_df, aes(label = fusion_label), size = 2.2,
                  color = "grey20", vjust = -0.7)
}
prioris_map <- function(colvar, cscale) {
  p <- ggplot(map_df, aes(x = score_norm, y = max_patient)) +
    annotate("rect", xmin = 0.65, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "#d62728", alpha = 0.05) +
    geom_vline(xintercept = 0.65, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = 0.40, linetype = "dotted", color = "grey55", linewidth = 0.3) +
    geom_point(aes(size = concentration, color = .data[[colvar]]), alpha = 0.85) +
    scale_y_log10() +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_size_area(max_size = 12, limits = c(0, 1),
                    breaks = c(0.25, 0.5, 0.75, 1), labels = percent_format(accuracy = 1),
                    name = "Concentration\n(part du patient principal)") +
    cscale +
    labs(x = "Score biologique", y = "Comptage k-mer max chez un patient (log)") +
    theme(plot.title = element_text(face = "bold"))
  add_labels(p)
}
if (nrow(map_df) > 0) {
  p_map_rep <- prioris_map("repartition",
      scale_color_manual(values = REPART_COLORS, name = "Répartition entre patients", na.value = "grey70")) +
    labs(title = "Carte de priorisation des fusions chromo-spécifiques",
         subtitle = "Score biologique × expression max par patient · taille = concentration sur le patient principal")
  ggsave(file.path(DIR_FIG, "carte_priorisation_repartition.pdf"), p_map_rep,
         width = 12, height = 9, device = PDF_DEV)

  p_map_type <- prioris_map("type_base",
      scale_color_manual(values = TYPE_COLORS, drop = FALSE, name = "Type chimérique", na.value = "grey65")) +
    labs(title = "Carte de priorisation — colorée par type chimérique",
         subtitle = "Score biologique × expression max par patient · taille = concentration sur le patient principal")
  ggsave(file.path(DIR_FIG, "carte_priorisation_type.pdf"), p_map_type,
         width = 12, height = 9, device = PDF_DEV)
}

# ── 10f. CHARGE PAR PATIENT (fusions chromo-spécifiques portées par patient) ──
spec_pairs <- feat %>% filter(presence_normaux == 0) %>%
  transmute(gene_pair, type_base = factor(type_base, levels = names(TYPE_COLORS)))
burden <- per_pat %>% filter(expr > 0) %>%
  inner_join(spec_pairs, by = "gene_pair") %>%
  count(sample, type_base, name = "n")
if (nrow(burden) > 0) {
  # code anonyme par patient (P01 = plus forte charge) : la figure affiche les
  # codes (souvent nombreux), les vrais noms sont dans la table de correspondance.
  tot <- burden %>% group_by(sample) %>% summarise(t = sum(n), .groups = "drop") %>%
    arrange(desc(t)) %>% mutate(patient_id = sprintf("P%02d", row_number()))
  burden <- burden %>% left_join(tot, by = "sample")

  # table de correspondance code -> patient (+ charge totale et détail par type)
  corr <- burden %>%
    pivot_wider(id_cols = c(patient_id, sample, t), names_from = type_base,
                values_from = n, values_fill = 0) %>%
    rename(n_fusions_total = t) %>% arrange(desc(n_fusions_total))
  write_tsv(corr, file.path(DIR_OUT, "charge_par_patient_correspondance.tsv"))

  p_burden <- burden %>%
    mutate(patient_id = factor(patient_id, levels = rev(tot$patient_id))) %>%
    ggplot(aes(n, patient_id, fill = type_base)) +
    geom_col() +
    scale_fill_manual(values = TYPE_COLORS, drop = TRUE, name = "Type chimérique") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Charge de fusions chromo-spécifiques par patient",
         subtitle = paste0("Codes anonymes (P01 = charge la plus forte) · noms réels dans ",
                           "charge_par_patient_correspondance.tsv"),
         x = "Nombre de fusions", y = NULL) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 6))
  ggsave(file.path(DIR_FIG, "charge_par_patient.pdf"), p_burden, width = 10, height = 8, device = PDF_DEV)
}

# ── 11. RÉSUMÉ CONSOLE ───────────────────────────────────────────────────────
cat("\n=== RÉSUMÉ ===\n")
cat("Fusions totales      :", nrow(res), "\n")
cat("Spécifiques (0 normal):", nrow(res_spec), "\n")
cat("\nTypes reconstruits :\n"); print(count(res, type_base))
cat("\nPriorités :\n"); print(count(res, priorite))
cat("\nScore /", MAX_SCORE, "pts | Tables :", DIR_OUT, "| Figures :", DIR_FIG, "\n")
