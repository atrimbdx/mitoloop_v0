#!/usr/bin/env Rscript
# run_morfee.R : annotation MORFEE d'un VCF annoté par ANNOVAR (hg38)
# Usage : Rscript run_morfee.R in.hg38_multianno.vcf MORFEE_DATA.RData out.tsv "CHAMP1,CHAMP2,..."
# Sortie : TSV sans en-tête CHROM POS REF ALT puis un champ MORFEE par colonne ("." si absent)
#
# Version MitoLoop (même sortie, plus rapide), sans modifier MORFEE :
#  1. seuls les variants que MORFEE annote sont transmis (Func.ensGene UTR5, UTR3
#     ou exonic : MORFEE passe tous les autres) ;
#  2. l'annotation GENCODE (environ 3 millions de lignes, où MORFEE cherche chaque
#     transcrit par grep) est réduite aux transcrits cités par ANNOVAR pour ces
#     variants : grep y trouve les mêmes lignes ;
#  3. les variants sont répartis sur MORFEE_CORES processus (défaut THREADS, sinon
#     4 ; parallel::mclapply, mémoire partagée), chaque variant étant annoté
#     indépendamment des autres par MORFEE.

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(MORFEE)
  library(parallel)
})

args   <- commandArgs(trailingOnly = TRUE)
vcf_in <- args[1]; rdata <- args[2]; out <- args[3]
fields <- strsplit(args[4], ",")[[1]]
cores  <- as.integer(Sys.getenv("MORFEE_CORES", Sys.getenv("THREADS", "4")))
if (is.na(cores) || cores < 1) cores <- 1

t0 <- Sys.time()
load(rdata)                       # charge l'objet MORFEE_DATA
vcf <- readVcf(vcf_in)
message(sprintf("MORFEE : données chargées en %.0f s, %d variants",
                as.numeric(difftime(Sys.time(), t0, units = "secs")), nrow(vcf)))

# 1. Variants que MORFEE annote (même test que morfee.annotation)
inf  <- info(vcf)
keep <- vapply(seq_len(nrow(vcf)), function(i)
  any(as.character(inf[i, "Func.ensGene"]) %in% c("UTR5", "UTR3", "exonic")), logical(1))
vcf  <- vcf[keep]
message(sprintf("MORFEE : %d variants en UTR5, UTR3 ou exonic", nrow(vcf)))

flat <- function(x) {
  if (is(x, "List") || is.list(x))
    x <- vapply(as.list(x), function(e) paste(e[!is.na(e)], collapse = "&"), "")
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA
  gsub("[;= \t]", "_", x)          # caractères interdits dans un champ INFO
}

if (nrow(vcf) == 0) {
  file.create(out)
  quit(save = "no")
}

# 2. Transcrits cités par ANNOVAR (ensGene) : annotation GENCODE réduite à eux
inf  <- info(vcf)
txt  <- unlist(lapply(c("GeneDetail.ensGene", "AAChange.ensGene"), function(f)
  if (f %in% colnames(inf)) as.character(unlist(inf[[f]])) else character(0)))
enst <- unique(unlist(regmatches(txt, gregexpr("ENST[0-9]+", txt))))
ann_id <- sub("\\..*$", "", MORFEE_DATA[["GENCODE_ANNOT"]]$transcript_id)
MORFEE_DATA[["GENCODE_ANNOT"]] <- MORFEE_DATA[["GENCODE_ANNOT"]][which(ann_id %in% enst), ]
md_id <- sub("\\..*$", "", MORFEE_DATA[["GENCODE_METAD"]]$V1)
# la première ligne est gardée en tête : morfee.annotation la duplique (GRCh38)
MORFEE_DATA[["GENCODE_METAD"]] <- MORFEE_DATA[["GENCODE_METAD"]][unique(c(1, which(md_id %in% enst))), ]
message(sprintf("MORFEE : %d transcrits, %d lignes GENCODE gardées", length(enst),
                nrow(MORFEE_DATA[["GENCODE_ANNOT"]])))

# 3. Annotation, variants répartis entre processus
chunks <- split(seq_len(nrow(vcf)), cut(seq_len(nrow(vcf)), min(cores, nrow(vcf)), labels = FALSE))
run1 <- function(ix) tryCatch(morfee.annotation(vcf[ix], MORFEE_DATA),
                              error = function(e) { message("MORFEE : ", conditionMessage(e)); NULL })
t1 <- Sys.time()
res_all <- if (length(chunks) > 1) mclapply(chunks, run1, mc.cores = length(chunks)) else lapply(chunks, run1)
message(sprintf("MORFEE : annotation en %.0f s sur %d processus",
                as.numeric(difftime(Sys.time(), t1, units = "secs")), length(chunks)))

tabs <- list()
for (res in res_all) {
  for (o in res) {
    if (is.null(o) || length(o) == 0) next
    df <- data.frame(CHROM = as.character(seqnames(o)), POS = start(o),
                     REF = as.character(ref(o)), ALT = as.character(unlist(alt(o))),
                     stringsAsFactors = FALSE)
    inf_o <- info(o)
    for (f in fields) df[[f]] <- if (f %in% colnames(inf_o)) flat(inf_o[[f]]) else NA_character_
    tabs[[length(tabs) + 1]] <- df
  }
}
all <- do.call(rbind, tabs)

if (is.null(all) || nrow(all) == 0) {
  file.create(out)
  quit(save = "no")
}

key <- paste(all$CHROM, all$POS, all$REF, all$ALT, sep = "\t")
agg <- sapply(fields, function(f)
  tapply(all[[f]], key, function(x) {
    u <- unique(x[!is.na(x)])
    if (length(u)) paste(u, collapse = "&") else "."
  }))
keys <- names(tapply(key, key, length))
agg  <- matrix(agg, ncol = length(fields), dimnames = list(keys, fields))
keep <- apply(agg, 1, function(r) any(r != "."))
res_df <- data.frame(key = rownames(agg)[keep], agg[keep, , drop = FALSE],
                     stringsAsFactors = FALSE, check.names = FALSE)

write.table(res_df, out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
