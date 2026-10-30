# Sampling uncertainty and directional support: protocol v1

Owner: Pauline Owusu-Ansah. Established 2026-10-29 (America/New_York).
This protocol and the tested research implementation will be committed before
running the evaluation seeds. This is an internal prospective evaluation, not
an external preregistration, independent validation, or a manuscript submission.
The machine-readable configuration is `protocol-v1.json`.

## Aim and estimands

Evaluate the nominal 95% intervals for the **genetic contrast**, and the
frequency of directional labels when that contrast is truly zero. A zero
genetic contrast is not a claim of zero or symmetric demographic migration.
No migration process is simulated in this experiment. The previous demographic
pilots already show why genetic and biological interpretation must be separated.

For latent population allele frequencies p and q at a locus, let
`r = sqrt(p*q)/(sqrt(p*q)+sqrt((1-p)*(1-q)))`, the hypothetical pool frequency.
For direction A, define `u_A=(p-r)^2/2`, `v_A=(p+r)*(1-(p+r)/2)`; use q for B.
The contrast is `sum(u_B)/sum(v_B) - sum(u_A)/sum(v_A)`.

Two targets are retained for every replicate:

1. **Realized-panel target:** apply the statistic to the exact latent allele
   frequencies at the sampled loci. This conditions on the genomic panel and
   measures uncertainty from sampling individuals.
2. **Genomic-distribution target:** replace each sum by an expectation over the
   specified four equally probable frequency types. This target averages across
   genomic blocks as well as individuals. It is available analytically from four
   types; no estimated high-precision simulation truth is substituted.

The primary coverage target is the realized panel for a fixed-panel scenario
and the genomic distribution for a random-panel scenario. Secondary panel
coverage is reported for all methods in every random-panel scenario. The true
contrast is computed independently using two-allele heterozygosities, without
calling snpflowR's estimator. Null targets are analytically exactly zero.

## Data generation

Each dataset has two populations and 80 independent genomic blocks. At each
block there is either one SNP or four exactly duplicated SNPs (a deliberate
complete-LD stress test, not a realistic recombination model). Equal block sizes
allow an exact collapsed implementation without changing the resampling law.

The four p types are 0.1, 0.2, 0.3, 0.4. In differentiated null scenarios,
q=1-p, making both directional GST values equal. Identical-population null
scenarios use q=p and test the nonregular boundary where true GST is zero.
Alternatives use q=0.35, 0.4, 0.45, 0.5; the reversed alternative swaps p and q.
Fixed panels contain 20 blocks of each type. Random panels independently draw
80 types with replacement from the uniform four-type distribution.

Diploid genotypes are Binomial(2, frequency). Individual counts are 40:40 or
10:80. Selected scenarios set 20% of whole diploid calls missing independently
of genotype; all four copies of a linked block share the missing call.
In row-structure scenarios each individual receives an independent latent
shift of +0.08 or -0.08 with equal probability, shared across all its blocks.
Genotypes are conditionally independent given that shift and block frequencies.
This induces correlation across loci within individuals; marginal allele
frequencies remain exactly p and q. Marginal genotypes need not satisfy HWE.
This is a controlled dependency model, not a fitted admixture history.

There is no MAF filter, post hoc locus selection, replacement of unfavorable
replicates, genotype error, empirical dataset, selection, family sampling, ghost
population, or demographic model. The scope is explicitly narrower than a full
population-genetic validation study.

## Methods and reporting

All methods use the same plug-in point estimate and 999 percentile bootstrap
draws (R quantile type 7, nominal 95%). Population order, alleles, and directions
are paired within every replicate. The four methods are:

| Method | Resampled units | Fixed units |
|---|---|---|
| locus | SNP columns, independently | Sampled individuals |
| block | Complete genomic blocks | Sampled individuals |
| individual | Whole genotype rows within each population | Realized genomic panel |
| crossed | Whole rows within populations and whole blocks, independently | None of those sampled units |

Individual and crossed resampling are **research prototypes**, not new public
package APIs. Crossed resampling is not a nested hierarchical bootstrap. Rows
and blocks are crossed factors; multiplying their independent bootstrap weights
may overstate interaction variance. Results for bootstrap means in the cited
literature do not prove validity for this nonlinear migrant-pool statistic.

With one SNP per block, locus and block methods use identical draws. For four
identical SNP copies, independent locus draws aggregate to Multinomial(320,1/80)
block weights, while block draws use Multinomial(80,1/80). The crossed method
uses the same row weights as the individual method and the same block weights
as the block method, reducing Monte Carlo noise in paired method comparisons.

Point-estimate locus exclusions match snpflowR: fewer than two observed allele
copies, no shared alleles, and shared fixation contribute zero directional
components, and an all-zero denominator remains undefined. All exclusions are
counted. Thresholds are not fitted to the evaluation outcomes.

Supported genetic direction requires finite nonsingular point estimates,
at least 20 contributing blocks and 20 usable SNPs, at least 95% finite draws,
at least five expected draws per nominal tail, and a percentile interval wholly
above or below zero with a concordant point estimate. There is one population
pair, so Bonferroni adjustment has no effect. For the one-factor methods this
matches the package safeguards on these data; all scenarios have 80 blocks.
For individual/crossed prototypes the same genomic-information gate is used,
without adding an unvalidated minimum of 20 sampled individuals that would
automatically suppress every 10:80 result. Interval coverage is evaluated
separately from the support gate. None of these are validated biological rules.

## Replication, seeds, and failure accounting

Twelve primary null scenarios have 1,000 datasets each. Five alternatives have
625 each: **15,125 datasets and 60,500 method-result rows**. Near a 5% rate,
1,000 replicates give Monte Carlo SE about 0.69 percentage points; this is not
precision sufficient to establish an exact 5% error rate. With 625 replicates,
the worst-case Monte Carlo SE for an alternative proportion is 2 percentage
points. Report Wilson 95% intervals and Monte Carlo SE, not only percentages.

Evaluation dataset seed = 60000000 + 100000 * scenario_id + replicate.
Bootstrap seed = 90000000 + 100000 * scenario_id + replicate.
Development/smoke checks use the separate 50000000 range, with no overlap.
R RNG is Mersenne-Twister/Inversion/Rejection. There is no empirical tuning of
the interval level, number of replicates, thresholds, or scenario list after
examining evaluation outcomes. Any amendment requires a separately logged
protocol version and new evaluation seeds.

Every attempted dataset and all four method outcomes are retained. Execution
errors are recorded and cause an incomplete-study status; they are not silently
replaced. Undefined intervals count as failures for unconditional coverage and
are separately reported; conditional coverage among defined intervals is also
shown. Early stopping based on significance or favorable performance is forbidden.
Each completed scenario is checkpointed for deterministic continuation.

## Performance and decision criteria

Report coverage for both targets, undefined interval frequency, mean/median
width, point-estimate bias and empirical SD, correct/wrong/unresolved support,
and false genetic-direction labels under each null. Use all planned datasets
as the denominator. Keep the two identical-population boundary nulls separate
from the differentiated nulls. Record runtime, peak resident memory, software
versions, seeds, source checksums, and run completion status.

For a provisional per-scenario screen, the Wilson interval for coverage must
lie within 92.5%-97.5%, and the upper Wilson bound for false directional support
must be <=7.5% under each null. These are prespecified **practical screening
tolerances**, not a formal equivalence test or proof of calibration. Very wide
intervals are not rewarded as calibrated merely for exceeding 95% coverage.
No method will be labeled generally calibrated on the basis of selected favorable
scenarios; no automatic change of package defaults follows this experiment.
Any improved procedure motivated by these outcomes must face new held-out seeds.

## References

- Morris, White & Crowther (2019), simulation design and Monte Carlo uncertainty:
  https://doi.org/10.1002/sim.8086
- Owen (2007), crossed row/column bootstrap:
  https://doi.org/10.1214/07-AOAS122 ; https://arxiv.org/abs/0712.1111
- Owen & Eckles (2012), product-weight bootstrap and conservative variance:
  https://doi.org/10.1214/12-AOAS547
- Sundqvist et al. (2016), migrant-pool directional differentiation:
  https://doi.org/10.1002/ece3.2096
