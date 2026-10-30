# Pinned point-estimator comparison

Status: completed exploratory software audit. This is not a demographic
validation study, an independent expert review, or evidence of calibrated
confidence intervals. The 110 arithmetic fixtures are separate from the
1,100 simulated datasets in the [pilot record](README.md).

## What was executed

| Comparator | Version in DESCRIPTION | Immutable source revision |
|---|---|---|
| diveRsity | 1.9.91 | [`22943bc95c80f46418578e7de205d4841976d140`](https://github.com/kkeenan02/diveRsity/tree/22943bc95c80f46418578e7de205d4841976d140) |
| corridoR | 0.2.0 | [`9e2e0cb3aa4e8947ce4a0827703dc3bc5d333d43`](https://github.com/paulmaier/corridoR/tree/9e2e0cb3aa4e8947ce4a0827703dc3bc5d333d43) |

The audit executes the **unmodified estimator source functions**, loaded into
separate R environments. It does not install or validate the complete external
packages, their landscape dependencies, network plots, or other statistics.

- For diveRsity, the actual `rgp` GENEPOP reader, compiled `myTab` and `pwHt`
  helpers, and `divMigrate(stat = "gst", boots = 0, plot_network = FALSE)` run.
  The public function returns relative scores. Raw GST is additionally recovered
  from the upstream heterozygosity helper using `(mean(Ht)-mean(Hs))/mean(Ht)`;
  this adapter is not an additional public divMigrate output.
- For corridoR, the actual `new_popcounts`, conversion methods, and
  `directional_gst` functions run on the same chromosome counts.
- Every external file is verified against its pinned Git blob ID before any
  source is executed. The [manifest](results/comparator-source-manifest.json)
  also retains SHA-256 checksums. Upstream sources remain in an ignored cache;
  they are not copied into the snpflowR package or relicensed as its own work.
  The respective upstream licenses are GPL (>= 2) and MIT.

## Inputs and comparison rule

Ten deterministic fixtures cover three populations with distinct frequencies,
unequal sample sizes, missing population-locus calls, added shared fixed loci,
added opposite fixed loci, all-opposite fixation, all-shared fixation, an
identical pair within three populations, one fixed population with zero pool
diversity, and two identical populations. An additional 100 cases use frozen
seeds 810001-810100, three populations, and 50 loci. Random allele counts are
restricted to the interior (both alleles observed) by design.

All observed population-locus totals contain at least two chromosome copies;
the audit does not compare single-copy or haploid input policies.
These are numerical test inputs, not demographic simulations or estimates of
how commonly boundary conditions occur in real data. Genotypes are constructed
by deterministic chromosome pairing, not a Hardy-Weinberg sampling model.
Counts round-trip through snpflowR's genotype input and the real diveRsity
GENEPOP parser. Every population-locus frequency is checked before comparison.
Input construction, seeds, and file fingerprints are retained in the scripts
and [fixture manifest](results/comparator-fixtures.csv).

All 648 off-diagonal entries per implementation are retained, totaling 1,944
rows of [raw output](results/comparator-values.csv). Matrices are compared in
the same source-row/destination-column orientation, without transposition to
force agreement. Absolute tolerance is 1e-12. NA and NaN both represent an
undefined estimate for agreement testing; their distinct output types are
recorded. Infinite and finite values are not treated as equivalent.

## Results in the recorded R 4.3.3 environment

There are 107 fixtures with agreement for both GST and relative scores:
105 have finite directional estimates, and two are wholly undefined
(all-opposite or all-shared fixation). This includes missing data, unequal
samples, and added opposite fixation. Among finite estimates over **all** fixtures,
the largest absolute GST difference is 4.31e-16 and the largest relative-score
difference is 1.94e-14 (rounded upward). These are floating-point differences,
not biological effect sizes or accuracy estimates.

| Comparator against snpflowR | GST disagreements | Relative-score disagreements | Explanation |
|---|---:|---:|---|
| diveRsity | 0/648 | 4/648 | Normalization policy when a pair has zero GST |
| corridoR | 1/648 | 9/648 | Normalization policy and floating-point boundary behavior |

See [summary](results/comparator-summary.csv),
[every discrepant entry](results/comparator-discrepancies.csv), and
[captured warnings](results/comparator-warnings.csv). The observed differences
were investigated after the first audit run; this is an exploratory diagnostic
exercise, not a preregistered equivalence test. No fixture was removed and the
comparison tolerance was not increased to conceal a disagreement.

### Identical pair within three populations

The A/B pair has exactly equal frequencies and nonzero diversity, so its GST
is zero and `(1/GST-1)/4` is infinite. snpflowR retains that singularity and
withholds the entire relative-score matrix. diveRsity removes the infinite
transforms, then normalizes the four remaining finite directed entries.
Those four availability differences are a deliberate policy distinction.

In the pinned corridoR source, roundoff in heterozygosity subtraction instead
gave A/B GST approximately -2.73e-16. Transforming that tiny negative number
produced relative scores near -1.61e14 for the identical pair. For two identical
populations alone, the same effect produced normalized scores of 1 in both
directions. Exact numerical artifacts can depend on BLAS/compiler/platform;
these values describe the recorded run, not a guaranteed result everywhere.

### Fixed population with no diversity in its comparison pool

A has only the reference allele at every locus; B remains polymorphic. The
shared pool also has only the reference allele. For A versus that pool,
both Ht and Hs are zero, and GST is **undefined (0/0)**, not a defined zero.
snpflowR reports NA and leaves the pair unresolved. diveRsity also returns
an undefined value. The pinned corridoR source produced negative infinity
for this entry in the recorded run, followed by a negative relative score.
The reverse direction remains defined; its availability is retained explicitly.

### No shared alleles

Both pinned implementations effectively omit opposite-fixed loci from the
multilocus directional calculation, consistent with snpflowR's documented
exclusion policy. Adding such loci leaves the directional estimate unchanged
in these fixtures. An entirely opposite-fixed comparison is undefined.
Agreement between implementations does not settle the methodological choice:
the original paper discusses a maximal-differentiation convention. Sensitivity
to a fully specified alternative aggregation remains research work.

## Uncertainty methods are not equivalent

The pinned divMigrate source resamples **individuals within populations** and
compares separate 95% intervals for the normalized directional scores using
an interval non-overlap rule. snpflowR currently resamples **loci or blocks**,
conditions on sampled individuals, and forms an interval for the **paired GST
contrast**; its default also adjusts across unordered population pairs.
corridoR's audited point-estimator function does not supply bootstrap intervals.

Agreement of point estimates cannot justify substituting one procedure's
support label for another's. This audit does not run or validate the external
bootstrap routines, compare interval coverage, or establish which approach is
better. Individual/hierarchical resampling and a defined coverage target remain
priorities in the [research plan](../RESEARCH_PLAN.md).

## Reproduce and prevent regressions

From the repository root, with snpflowR installed, Python 3, Rcpp, and an R
C++ toolchain available:

```sh
R CMD INSTALL .
Rscript validation/run_comparator_audit.R
```

The R script invokes `fetch_comparators.py`, verifies all cached/downloaded
files, and fails on an unrecognized discrepancy. No external code is patched.
Network access is needed on the first run; a verified cache works offline.
Runtime package users and ordinary `R CMD check` do not need Rcpp, Python, or
either external package. The dedicated Comparator-audit workflow repeats the
source comparison; the normal five-platform checks run offline regression
tests against small retained reference outputs from both comparators.

CI requires agreement in the ordinary domain and checks snpflowR's boundary
behavior independently. It does not require a particular upstream roundoff
artifact to persist across platforms. All boundary outputs remain available
in the workflow artifact, including when their numerical details change.

Primary source inspection:
[divMigrate.R](https://github.com/kkeenan02/diveRsity/blob/22943bc95c80f46418578e7de205d4841976d140/R/divMigrate.R),
[pwHt.cpp](https://github.com/kkeenan02/diveRsity/blob/22943bc95c80f46418578e7de205d4841976d140/src/pwHt.cpp),
[genetic_stats.R](https://github.com/paulmaier/corridoR/blob/9e2e0cb3aa4e8947ce4a0827703dc3bc5d333d43/R/genetic_stats.R).
Underlying method: [Sundqvist et al. (2016)](https://doi.org/10.1002/ece3.2096).
