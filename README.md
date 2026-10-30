# snpflowR

**Directional genetic signals from SNPs, with explicit uncertainty and failure diagnostics.**

[![R-CMD-check](https://github.com/codewithPauline/snpflowR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/codewithPauline/snpflowR/actions/workflows/R-CMD-check.yaml)
[![Comparator-audit](https://github.com/codewithPauline/snpflowR/actions/workflows/comparator-audit.yaml/badge.svg)](https://github.com/codewithPauline/snpflowR/actions/workflows/comparator-audit.yaml)
[![Uncertainty-engine](https://github.com/codewithPauline/snpflowR/actions/workflows/uncertainty-engine.yaml/badge.svg)](https://github.com/codewithPauline/snpflowR/actions/workflows/uncertainty-engine.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE.md)

**Status: research development (`0.0.0.9000`). Not a CRAN release or a validated
tool for establishing contemporary migration.** See the [research plan](RESEARCH_PLAN.md)
and [validation record](validation/README.md) before interpreting results.

`snpflowR` estimates how each population differs from a hypothetical migrant
pool constructed from a pair of populations. It implements the published
[Sundqvist et al. (2016)](https://doi.org/10.1002/ece3.2096) construction for
biallelic SNPs and reports paired bootstrap intervals for the directional
contrast. The original migrant-pool method is **not a new invention of this
package**. The proposed contribution is an auditable SNP workflow and a rigorous
assessment of when its directional signals can and cannot be trusted.

## What is implemented

| Component | Current behavior |
|---|---|
| Input | Population allele counts; diploid hard-call matrices; a restricted local SNP VCF reader |
| Differentiation | Directional and symmetric plug-in multilocus GST |
| Relative signal | Descriptive migration transform, with singularities exposed |
| Uncertainty | Paired locus/block percentile intervals; nominal, with calibration not established |
| Reporting | Explicit `unresolved` labels, exclusion counts, interval diagnostics, provenance |
| QC | Minimum allele counts, optional pooled MAF filter, auditable retention decisions |
| Validation | Independent arithmetic, pinned comparator sources, invariance and parser tests, seeded simulations |

There is no Weir-Cockerham FST estimator, Jost's D estimator, demographic
likelihood, genotype-likelihood input, haploid/polyploid parser, or fitted
absolute migration rate in this version. The bundled Wright-Fisher simulator
has independent loci; it does not simulate realistic linkage.

## Install the development version

```r
# Install remotes from CRAN first if needed:
# install.packages("remotes")
remotes::install_github("codewithPauline/snpflowR", build_vignettes = TRUE)
library(snpflowR)
```

## A reproducible first analysis

```r
# Independent SNPs under a specified neutral simulation, not a real dataset.
m <- matrix(0, 2, 2)
m[1, 2] <- 0.02  # A contributes 2% of B's gene pool per generation
sim <- simulate_snp_flow(n_loci = 300, migration = m, seed = 42)

fit <- estimate_flow(
  sim$counts,
  n_boot = 999,
  resample = "locus",  # explicit independence assumption
  seed = 43
)
fit$pairs
plot(fit)
```

This example returns **unresolved even though migration was simulated** in the
recorded development environment. Known movement does not guarantee an
identifiable allele-frequency signal. Neither a large score nor a narrow
interval is proof of migration.

For genotype matrices, rows are individuals and columns are loci:

```r
g <- rbind(c(0, 1, NA), c(1, 2, 0), c(2, 1, 1), c(1, 0, 2))
counts <- genotype_counts(g, populations = c("A", "A", "B", "B"))
estimate_flow(counts)  # point estimates only; no supported direction labels
```

For a quality-filtered diploid SNP VCF:

```r
vcf <- system.file("extdata", "tiny.vcf", package = "snpflowR")
populations <- c(a1 = "A", a2 = "A", b1 = "B", b2 = "B")
counts <- read_snp_vcf(vcf, populations, block_size = 100000)
counts$import
```

The tiny fixture demonstrates import only. It contains too few loci and blocks
for inference. For real data, choose blocks using chromosome structure and LD
information, justify their size, and assess sensitivity. Physical windows do
not themselves establish block independence. The reader checks site FILTER,
but **does not filter genotype depth, genotype quality, likelihoods, or sample
FT**. Filter these upstream. It counts skipped variant records and partial
missing calls; missing genotypes are never made into reference homozygotes.

## Read the output correctly

- Matrix entry `[A, B]` compares A with the A/B migrant pool. A smaller
  directional GST corresponds to a larger descriptive A-to-B migration transform.
- `delta_gst = GST[B, pool] - GST[A, pool]`. Positive favors A as source under
  the model; negative favors B. The interval is for this genetic contrast.
- `direction = "unresolved"` means reporting safeguards are not satisfied or
  the interval includes the specified threshold. It does **not** establish
  symmetric migration or absence of movement.
- `bootstrap_supported_signal` means a bootstrap-supported genetic contrast,
  not a validated demographic conclusion, calibrated p-value, or probability
  that migration occurred.
- `provenance$calibration_status = "not_established"` records that nominal
  intervals have not demonstrated reliable coverage for the genetic contrast.
- Scaling is relative to the analyzed population set. Zero differentiation
  produces an infinite transform; no epsilon is inserted. Relative scaling
  is unavailable in that case.
- Loci fixed for opposite alleles have no shared migrant pool. They contribute
  to symmetric GST but are excluded from directional estimates and counted in
  `n_no_shared`. This is an explicit departure from the original paper's
  maximal-differentiation convention; it needs sensitivity analysis.

Unequal population sizes, sampling imbalance, founder events, ancestral
structure, unsampled populations, selection, SNP ascertainment, sequencing
error, and linkage can produce or obscure apparent directions. Locus/block
resampling conditions on the observed individuals and does not remove those
problems. Do not interpret the island-model transformation as an estimate of
absolute `Nm`.

**The initial 900-replicate pilot found substantial confounding:** in the
specified no-migration scenarios, supported signals occurred in 24/100
unequal-sample replicates, 100/100 unequal-effective-size replicates, and
73/100 founder-history replicates. These scenario-specific results are
development evidence, not universal performance rates. All raw outcomes are
retained in [validation/results](validation/results).

**The point-estimator comparison now covers 110 arithmetic fixtures** against
pinned diveRsity and corridoR sources. The 105 finite, nonsingular fixtures
agree to an absolute tolerance of 1e-12, and two wholly undefined cases also
agree. Three boundary cases document deliberate normalization differences and
numerical instability in one comparator; none were discarded. See the
[complete comparator audit](validation/COMPARATOR_AUDIT.md). This checks
computation, not biological reliability. divMigrate's individual bootstrap and
snpflowR's locus/block contrast bootstrap are different procedures; their
support labels are not interchangeable.

**The prospective sampling study is complete: 15,125 datasets across 17
scenarios. None of the four resampling methods met the prespecified overall
calibration screen.** In a fixed-panel null with samples of 10:80 individuals,
the block method reported false genetic-direction support in 83.2% of datasets.
For balanced samples with four identical SNP copies per block, false support
was 33.5% under locus resampling and 6.1% under block resampling. These are
scenario-specific results; changing resampling alone did not resolve bias.
All outcomes, Monte Carlo intervals, and literal-reference checks are preserved
in the [sampling-uncertainty report](validation/uncertainty/README.md).
The individual and crossed alternatives remain research prototypes; the
package's numerical estimator and thresholds have not been changed.

## Scientific and software development

- [Methods and worked examples](vignettes/snpflowR.Rmd)
- [Research questions, benchmark design, and release gates](RESEARCH_PLAN.md)
- [Validation scripts, results, and limitations](validation/README.md)
- [Pinned comparator results and boundary cases](validation/COMPARATOR_AUDIT.md)
- [Prospective sampling-uncertainty study and calibration failures](validation/uncertainty/README.md)
- [Developer instructions](DEVELOPMENT.md)
- [Contributing](CONTRIBUTING.md)
- [Report a bug](https://github.com/codewithPauline/snpflowR/issues)

Publication will require broader uncertainty-method comparisons and calibration,
expanded independent benchmarks, carefully sourced empirical data, and external
scientific review. CRAN submission is a separate software-distribution review. Neither
acceptance nor scientific novelty is claimed at this stage.

## Credit and citation

Package author and maintainer: **Pauline Owusu-Ansah**.
Use `citation("snpflowR")` for the software and underlying method citations.

Sundqvist L, Keenan K, Zackrisson M, Prodohl P, Kleinhans D (2016).
Directional genetic differentiation and relative migration.
*Ecology and Evolution* 6:3461-3475. <https://doi.org/10.1002/ece3.2096>.

Related software includes
[`diveRsity::divMigrate`](https://github.com/kkeenan02/diveRsity) and
[`corridoR::directional_gst`](https://paulmaier.github.io/corridoR/reference/directional_gst.html).
Point-estimate agreement and documented exceptions are recorded in the comparator
audit. General method equivalence or equivalent uncertainty estimates are not claimed.
