#' Read a restricted diploid biallelic SNP VCF into population counts
#'
#' @param file Path to a local VCF or gzip-compressed VCF. BCF is unsupported.
#' @param populations Named population vector; names must match VCF sample IDs
#'   exactly. The vector is aligned by sample ID, not input order.
#' @param pass_only Retain only records with FILTER exactly `PASS`. If FALSE,
#'   failed and unfiltered records are also read. An audit records exclusions.
#' @param block_size Optional positive integer physical block width in base pairs.
#'   Labels combine chromosome with `floor((POS - 1) / block_size)`. This does
#'   not estimate linkage or guarantee independent blocks.
#' @param chunk_size Number of variant lines read at a time.
#' @return A `snpflow_counts` object with locus metadata and an `import` audit.
#'   The parser streams genotypes into counts but retains all count columns in
#'   memory. It is intended for moderate files, not out-of-core biobank analysis.
#' @details Only single-base A/C/G/T REF and ALT records with distinct alleles are
#'   retained. Multiallelic records, indels, symbolic alleles, and other bases are
#'   skipped and counted. Unsupported ploidy and malformed retained records cause
#'   an error. GT must be the first FORMAT field. Phased and unphased diploid
#'   calls are accepted. Entire and partial missing calls contribute no alleles;
#'   they are never converted to reference homozygotes. Depth, genotype quality,
#'   likelihoods, INFO, and per-sample FT are not used: apply suitable quality
#'   filtering before import. Duplicate chromosome-position records are rejected.
#' @export
#' @examples
#' file <- system.file("extdata", "tiny.vcf", package = "snpflowR")
#' x <- read_snp_vcf(file, c(a1 = "A", a2 = "A", b1 = "B", b2 = "B"))
#' x$import
read_snp_vcf <- function(file, populations, pass_only = TRUE, block_size = NULL,
                         chunk_size = 1000L) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !file.exists(file)) {
    stop("file must identify an existing local VCF.", call. = FALSE)
  }
  .flag(pass_only, "pass_only")
  .whole_number(chunk_size, "chunk_size", 1)
  if (!is.null(block_size)) .whole_number(block_size, "block_size", 1)
  con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, "rt") else base::file(file, "rt")
  on.exit(close(con), add = TRUE)
  first <- readLines(con, n = 1L, warn = FALSE)
  if (length(first) != 1L || !grepl("^##fileformat=VCFv4\\.[0-5]$", first)) {
    stop("Expected a VCFv4.0 through VCFv4.5 fileformat header.", call. = FALSE)
  }
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("VCF has no #CHROM header.", call. = FALSE)
    if (startsWith(line, "##")) next
    header <- strsplit(line, "\t", fixed = TRUE)[[1]]
    required <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    if (length(header) < 11L || !identical(header[1:9], required)) {
      stop("Invalid VCF header or fewer than two samples.", call. = FALSE)
    }
    break
  }
  samples <- .labels(header[-seq_len(9)], length(header) - 9L, "VCF sample IDs")
  .labels(names(populations), length(populations), "Population-vector sample IDs")
  if (!setequal(names(populations), samples) || anyNA(populations) ||
      any(!nzchar(trimws(as.character(populations))))) {
    stop("populations must name every VCF sample exactly once with a nonempty population label.", call. = FALSE)
  }
  populations <- as.character(populations[match(samples, names(populations))])
  pops <- unique(populations)
  if (length(pops) < 2L) stop("At least two populations are required.", call. = FALSE)
  audit <- c(records = 0L, retained = 0L, filtered = 0L, non_biallelic_snp = 0L,
             missing_genotypes = 0L, partial_missing_genotypes = 0L)
  alt_list <- called_list <- metadata <- list()
  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    for (line in lines) {
      audit["records"] <- audit["records"] + 1L
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) != length(header) || any(!nzchar(fields))) {
        stop("Malformed VCF record ", audit["records"], ": wrong field count or empty field.", call. = FALSE)
      }
      if (pass_only && fields[7] != "PASS") {
        audit["filtered"] <- audit["filtered"] + 1L
        next
      }
      if (!(fields[4] %in% c("A", "C", "G", "T")) ||
          !(fields[5] %in% c("A", "C", "G", "T")) || fields[4] == fields[5]) {
        audit["non_biallelic_snp"] <- audit["non_biallelic_snp"] + 1L
        next
      }
      pos <- suppressWarnings(as.numeric(fields[2]))
      if (!grepl("^[0-9]+$", fields[2]) || !is.finite(pos) || pos < 1 || pos > 2^53 ||
          fields[1] == ".") stop("Invalid chromosome or position.", call. = FALSE)
      fmt <- strsplit(fields[9], ":", fixed = TRUE)[[1]]
      if (fmt[1] != "GT" || anyDuplicated(fmt)) {
        stop("GT must be the first FORMAT field and FORMAT keys must be unique.", call. = FALSE)
      }
      sample_fields <- strsplit(fields[-seq_len(9)], ":", fixed = TRUE)
      if (any(lengths(sample_fields) > length(fmt))) stop("Too many genotype subfields.", call. = FALSE)
      gt <- vapply(sample_fields, `[`, character(1), 1L)
      if (any(!grepl("^([01.](/|\\|)[01.]|\\.)$", gt))) {
        stop("Unsupported or malformed genotype at ", fields[1], ":", fields[2],
             "; expected diploid biallelic GT (0/1 or 0|1) or missing.", call. = FALSE)
      }
      absent <- grepl(".", gt, fixed = TRUE)
      partial <- absent & !(gt %in% c(".", "./.", ".|."))
      audit["missing_genotypes"] <- audit["missing_genotypes"] + sum(absent)
      audit["partial_missing_genotypes"] <- audit["partial_missing_genotypes"] + sum(partial)
      dosage <- rep(NA_real_, length(gt))
      dosage[!absent] <- as.numeric(substr(gt[!absent], 1, 1)) + as.numeric(substr(gt[!absent], 3, 3))
      a <- vapply(pops, function(p) sum(dosage[populations == p], na.rm = TRUE), numeric(1))
      n <- vapply(pops, function(p) 2 * sum(!is.na(dosage[populations == p])), numeric(1))
      audit["retained"] <- audit["retained"] + 1L
      k <- audit["retained"]
      alt_list[[k]] <- a
      called_list[[k]] <- n
      metadata[[k]] <- data.frame(locus = paste(fields[1], fields[2], fields[4], fields[5], sep = ":"),
        chromosome = fields[1], position = pos, ref = fields[4], alt = fields[5],
        id = fields[3], stringsAsFactors = FALSE)
    }
  }
  if (!length(metadata)) stop("No usable biallelic SNP records remain.", call. = FALSE)
  loci <- do.call(rbind, metadata)
  if (anyDuplicated(paste(loci$chromosome, loci$position, sep = ":"))) {
    stop("Duplicate chromosome-position records; normalize and deduplicate upstream.", call. = FALSE)
  }
  alt <- do.call(cbind, alt_list)
  called <- do.call(cbind, called_list)
  colnames(alt) <- colnames(called) <- loci$locus
  blocks <- if (is.null(block_size)) NULL else paste(loci$chromosome, floor((loci$position - 1) / block_size), sep = ":")
  ans <- snp_counts(alt, called, loci, blocks)
  ans$import <- list(format = first, audit = audit, pass_only = pass_only,
    block_size = block_size, quality_filtering = "site FILTER only; genotype quality must be filtered upstream")
  ans$history <- list(list(operation = "read_snp_vcf", import = ans$import))
  ans$sample_sizes <- table(factor(populations, levels = pops))
  ans
}
