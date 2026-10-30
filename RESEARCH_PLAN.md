# Research plan and release gates

Owner: Pauline Owusu-Ansah. Established 2026-10-28.
Status: active research development; neither a preregistration nor a submission.

## Research question

Under what sampling designs and demographic conditions can a SNP-based
migrant-pool contrast usefully indicate asymmetric gene flow, and when should
the analysis withhold a biological directional conclusion?

The starting point is the published Sundqvist et al. (2016) method. A new name
and R implementation alone do not establish scientific novelty. Candidate
contributions must be demonstrated: SNP-specific diagnostics, better handling
of sampling/linkage uncertainty, and a transparent map of failure conditions.
If the evidence supports a limitations/benchmark paper rather than an improved
estimator, report that outcome honestly.

## Separate three questions

1. **Computational correctness:** does the implementation calculate the defined
   statistic and resample the specified units correctly?
2. **Statistical calibration:** do intervals cover a precisely defined genetic
   estimand, and how often do null genetic contrasts receive directional labels?
3. **Demographic validity:** when does the sign indicate the actual direction of
   migration rather than drift, size imbalance, ancestry, or ascertainment?

Passing the first question cannot establish the second or third. A nonzero
genetic contrast under zero migration can be a real contrast with a misleading
demographic interpretation. The initial pilot already exhibits this behavior.

## Work packages

### 1. Reproducible foundation

- Implement count/genotype inputs and a carefully limited VCF parser.
- Define orientation, aggregation, missingness, fixed loci, and numerical boundaries.
- Test independent arithmetic, allele/population relabeling, and whole-block pairing.
- Preserve package metadata, method credit, source provenance, seeds, and outcomes.
- Keep development claims distinct from a stable release.

### 2. Estimator and comparator audit

Completed: a [pinned source comparison](validation/COMPARATOR_AUDIT.md) on
110 arithmetic fixtures, with parser round trips, agreement in the ordinary
domain, and preserved numerical/normalization differences at the boundaries.
This covers point estimates only. Bootstrap procedures differ and require a
separate study; no claim of full-package equivalence is made.

- Reproduce paper examples where the supported input domain permits it.
- Compare pinned `diveRsity::divMigrate` and `corridoR::directional_gst` versions
  on shared valid SNP fixtures and document every discrepancy.
- Examine the consequences of excluding no-shared-allele loci versus other
  explicit conventions; do not silently force agreement.
- Evaluate finite-sample bias, zero counts, sample-size balancing/rarefaction,
  and individual or hierarchical bootstrap alternatives using held-out seeds.
- Consider a correctly named symmetric FST estimator and a separately justified
  Jost's D implementation only with independent references and tests.

### 3. Calibration and robustness study

Completed first prospective sampling evaluation: [17 scenarios, 15,125 datasets](validation/uncertainty/README.md),
with a public protocol frozen before evaluation. None of the locus, block,
individual, or crossed procedures met the overall screen. The study measures
coverage for known genetic contrasts; it does not calibrate demographic claims.
Unequal sampling exposes strong finite-sample bias, while resampling units affect
linkage and individual-dependency uncertainty. These failures block promotion of
the research prototypes or a general calibration claim.

The next development cycle must investigate bias correction and balanced-sampling
strategies, with separate attention to zero-differentiation boundaries. Any method
selected using these outcomes needs a new protocol and fresh evaluation seeds.

Freeze a versioned protocol before the confirmatory study. Exploratory pilot
results must remain labeled exploratory. Use a development grid for design and
distinct, unreused simulation seeds for evaluation. No selection of favorable
replicates, retrospective threshold tuning on the evaluation set, or replacing
failed runs without an audit.

Vary migration (none, symmetric, one-way, asymmetric, changing through time),
effective sizes, split times, founder/bottleneck history, ghost populations,
sample sizes, neutral SNP counts, missingness, ascertainment, genotype error,
recombination/LD, block lengths, and population numbers. Include realistic
scenarios generated independently of the implementation. Forward migration
fractions, counts of migrants, and backward coalescent rates must be distinguished.

Prespecify at least 1,000 replicates per primary null condition, with more when
the desired Monte Carlo precision requires it. Report binomial uncertainty;
1,000 replicates do not justify an exact 5% calibration claim. Select alternative
condition replicate counts by a stated precision target, not convenience.

Primary summaries:

- False biological direction calls under zero and symmetric migration.
- Correct, wrong, and unresolved proportions over **all** replicates.
- Sign recovery among finite estimates, separately from supported-call accuracy.
- Coverage and width for a declared genetic estimand, using an appropriate
  independent high-precision target or analytical truth.
- Effects of sample balancing, QC, private alleles, block size, and ancestry.
- Runtime and peak memory with input dimensions, machine, and versions recorded.

Current safeguards (20 loci, 20 units, nominal tails) are engineering defaults,
not empirically established biological thresholds. Do not present them as
calibration. Distinguishing founders from ongoing asymmetric migration may
require additional demographic information or a different model.

### 4. Empirical demonstration

Select at least one openly licensed, traceable dataset with population metadata,
sampling dates, quality-control information, and an independent biological
reason to assess connectivity. Preserve accessions, provenance, licenses, and
checksums. Use an appropriate demographic comparator where practical. Without
ground-truth movement, an empirical example cannot establish predictive accuracy.
No empirical dataset or agreement with field movement is currently claimed.

### 5. External review and distribution

Arrange independent population-genetics and statistical review of the equations,
estimand, simulations, and interpretation. Record review comments and responses.
Seek user testing with clear supported input sizes and informative error reports.

Only after that work: freeze an API, create a stable version and release notes,
archive an immutable release (with a real DOI if obtained), prepare a manuscript
with limitations and comparator results, and consider a journal matched to the
demonstrated contribution. Never fabricate a DOI, preregistration, acceptance,
review, or award. CRAN distribution and journal review are separate processes.

## Release gates

| Gate | Evidence needed before a stable scientific release |
|---|---|
| Defined method | Written estimand, assumptions, orientation, all edge cases, method attribution |
| Correct implementation | Independent numeric checks, fixtures, parser tests, invariants, no regressions |
| Comparator agreement | Pinned external software comparisons with documented differences |
| Calibration | Frozen confirmatory design, adequate replication, uncertainty, held-out evaluation |
| Confounding | Explicit failure regions and validated scope; no unsupported migration claims |
| Empirical use | Reproducible, licensed demonstration with honest limits on ground truth |
| Review | Independent expert feedback addressed and logged |
| Packaging | Current R release/devel checks, Windows/macOS/Linux results, docs and examples |
| Submission | Name/Archive/Bioconductor checks, maintainer details, license audit, CRAN policies rechecked |

The implementation, pilots, comparator audit, and prospective sampling study
provide evidence for the release gates. The sampling study found calibration
failures rather than establishing a calibrated method. Empirical demonstration
and independent review also remain open. The full release criteria are unmet.

## Primary references and policy sources

- Sundqvist et al. (2016): <https://doi.org/10.1002/ece3.2096>
- `corridoR` documentation: <https://paulmaier.github.io/corridoR/reference/directional_gst.html>
- `diveRsity` source: <https://github.com/kkeenan02/diveRsity>
- msprime demographic conventions: <https://tskit.dev/msprime/docs/stable/demography.html>
- VCF specification: <https://samtools.github.io/hts-specs/VCFv4.5.pdf>
- Writing R Extensions: <https://cran.r-project.org/doc/manuals/r-release/R-exts.html>
- CRAN policies: <https://cran.r-project.org/web/packages/policies.html>
- CRAN submission checklist: <https://cran.r-project.org/web/packages/submission_checklist.html>

Sources reviewed 2026-10-28; recheck distribution policies before submission.
