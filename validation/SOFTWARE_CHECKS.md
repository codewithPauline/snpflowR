# Software verification record

## Latest milestone: pinned comparator audit

Completed 2026-10-29 (America/New_York). Research version: `0.0.0.9000`.
Tested source and workflow revision:
`71c0c553f852e26cfacbc4b82686ba98a218a0fa`.

**All six GitHub jobs passed.** Each of the five package-check jobs reported
`Status: OK` and **146 passing assertions**, with no failures, test warnings,
or skips. The separate comparator job verified eight upstream source blobs and
completed the 110-fixture audit, retaining the documented boundary differences.

| Check | Result | Evidence |
|---|---|---|
| Linux, R release | Passed | Package-check workflow below |
| Linux, R devel | Passed | Package-check workflow below |
| Linux, R oldrel-1 | Passed | Package-check workflow below |
| Windows, R release | Passed | Package-check workflow below |
| macOS, R release | Passed | Package-check workflow below |
| Pinned comparator sources, Linux | Passed | Comparator workflow below |

- [Package-check run 40981479667](https://github.com/codewithPauline/snpflowR/actions/runs/40981479667)
- [Comparator run 40981479706](https://github.com/codewithPauline/snpflowR/actions/runs/40981479706)
- [Retained job conclusions and log excerpts](results/comparator-ci-checks.json)
- [Full comparator methods, outputs, and exceptions](COMPARATOR_AUDIT.md)

The fresh local R 4.3.3 `--as-cran` check also completed with **0 errors,
0 warnings, and 1 note**. The note concerns a new submission and the development
version's `9000` component. All 146 assertions, examples, rebuilt vignettes,
and PDF/HTML manual checks passed. The [complete log](results/R-CMD-check-comparators-R4.3.3.log)
was recorded after publishing the new workflow, so its README link was checked
at an existing URL.

```sh
# From the checkout's parent directory:
R CMD build snpflowR-repo
_R_CHECK_SYSTEM_CLOCK_=FALSE R_RD4PDF='times,hyper' \
  R CMD check --as-cran --output=comparator-verified-check snpflowR_0.0.0.9000.tar.gz
```

The local clock-network probe was disabled; the manual used installed Times
fonts. CRAN incoming URL checks remained enabled. GitHub's package-check action
uses its standard incoming/clock settings and checks without building the PDF
manual. Package functions, scientific thresholds, and the estimator were not
changed by this milestone. The added 42 assertions protect reference agreement
and boundary behavior; they do not calibrate demographic conclusions.

## Initial foundation milestone (historical record)

Date: 2026-10-28. Package version: `0.0.0.9000` (research development).
Source and workflow revision: `5069e9824ffa336398ad2777a037031e5c32b555`.

## Local package check

Ubuntu 24.04.3, R 4.3.3: **0 errors, 0 warnings, 1 note**.
All **104 test assertions passed**, with no skipped tests or test warnings.
Examples, vignette rebuilding, documentation, the PDF manual, and HTML manual
validation passed. The [full check log](results/R-CMD-check-R4.3.3.log) is retained.

Commands, executed from the checkout's parent directory:

```sh
R CMD build snpflowR-repo
_R_CHECK_SYSTEM_CLOCK_=FALSE R_RD4PDF='times,hyper' \
  R CMD check --as-cran --output=verified-check snpflowR_0.0.0.9000.tar.gz
```

The clock-network probe was disabled for the execution environment. The PDF
manual used the installed Times font configuration instead of Inconsolata;
manual generation and validation still ran. CRAN incoming checks were enabled.
The remaining note identifies a new submission and the development version's
large `9000` component. A release version and fresh submission checks are needed
before any CRAN submission. This record is not CRAN acceptance.

## Continuous integration

**All five jobs passed** for the recorded source revision. Each job reported
`Status: OK` and 104 passing assertions, with no failures, test warnings, or skips.

| Platform | R channel | Result |
|---|---|---|
| Linux | release | Passed |
| Linux | oldrel-1 | Passed |
| Linux | devel | Passed |
| Windows | release | Passed |
| macOS | release | Passed |

[GitHub Actions run 40786850652](https://github.com/codewithPauline/snpflowR/actions/runs/40786850652)
checks Linux with R release, oldrel-1, and devel, plus Windows and macOS with R
release. These jobs build the vignette and run `R CMD check --as-cran --no-manual`.
The action's standard configuration disables CRAN incoming and clock-network
checks. Warnings cause failure; tests are not skipped.

The first workflow run failed before checking the package because the action's
`error-on` input was parsed as the R function `warning`, not a character string.
The recorded revision corrects the expression to `"warning"`. It also enforces
LF line endings for text files and generated validation CSVs. It changes no
estimator, simulation model, reporting threshold, or scientific result.

## Scope of verification

Software checks establish that the declared calculations and interfaces run
under the tested environments. They do not establish demographic identifiability,
interval calibration, scientific novelty, or publication readiness. See the
[1,100-dataset pilot record](README.md) and [research plan](../RESEARCH_PLAN.md)
for the observed limitations and remaining scientific work.
