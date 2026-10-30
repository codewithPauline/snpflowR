# Comparator reference outputs

`comparator-counts.csv` holds the first five deterministic SNP fixtures in
`validation/comparator_fixtures.R`. `comparator-expected.csv` holds the actual
outputs from both pinned external estimator sources, extracted from
`validation/results/comparator-values.csv`. These numeric outputs are not copies
of upstream code. They test orientation, unequal counts, missing population-loci,
and the shared/opposite fixation conventions without runtime downloads.

Sources: diveRsity 1.9.91 at `22943bc95c80f46418578e7de205d4841976d140` and
corridoR 0.2.0 at `9e2e0cb3aa4e8947ce4a0827703dc3bc5d333d43`.
See `validation/COMPARATOR_AUDIT.md` in the source repository for the execution
scope, original licenses, checksums, and all boundary differences. No complete
external-package or bootstrap equivalence is implied.

Reference files are reviewed fixtures, not automatically rewritten by tests.
