# Complete uncertainty study tables

Percentages with Wilson 95% Monte Carlo intervals. These describe the specified simulations, not universal performance.

## Primary-target interval coverage

| Scenario | Locus | Block | Individual | Crossed |
|---|---:|---:|---:|---:|
| null_fixed_balanced | 93.50 [91.80, 94.87] | 93.50 [91.80, 94.87] | 97.20 [95.98, 98.06] | 100.00 [99.62, 100.00] |
| null_fixed_unequal | 16.80 [14.61, 19.24] | 16.80 [14.61, 19.24] | 0.10 [0.02, 0.56] | 1.90 [1.22, 2.95] |
| null_random_balanced | 93.70 [92.02, 95.05] | 93.70 [92.02, 95.05] | 98.00 [96.93, 98.70] | 100.00 [99.62, 100.00] |
| null_random_unequal | 19.20 [16.88, 21.76] | 19.20 [16.88, 21.76] | 0.50 [0.21, 1.17] | 2.70 [1.86, 3.90] |
| null_identical_balanced | 91.70 [89.83, 93.25] | 91.70 [89.83, 93.25] | 99.40 [98.70, 99.72] | 100.00 [99.62, 100.00] |
| null_identical_unequal | 51.50 [48.40, 54.59] | 51.50 [48.40, 54.59] | 17.20 [14.99, 19.66] | 59.90 [56.83, 62.89] |
| null_linked_balanced | 66.50 [63.52, 69.36] | 93.90 [92.24, 95.22] | 97.60 [96.45, 98.38] | 99.90 [99.44, 99.98] |
| null_linked_unequal | 1.70 [1.06, 2.71] | 17.70 [15.46, 20.19] | 0.50 [0.21, 1.17] | 2.80 [1.94, 4.02] |
| null_missing_balanced | 92.10 [90.26, 93.62] | 92.10 [90.26, 93.62] | 98.50 [97.54, 99.09] | 99.90 [99.44, 99.98] |
| null_missing_unequal | 10.10 [8.38, 12.12] | 10.10 [8.38, 12.12] | 0.20 [0.05, 0.73] | 0.50 [0.21, 1.17] |
| null_row_structure_balanced | 62.60 [59.56, 65.55] | 62.60 [59.56, 65.55] | 93.80 [92.13, 95.13] | 97.90 [96.81, 98.62] |
| null_row_structure_unequal | 32.30 [29.47, 35.26] | 32.30 [29.47, 35.26] | 9.00 [7.38, 10.93] | 25.20 [22.61, 27.98] |
| alternative_fixed | 78.24 [74.84, 81.30] | 78.24 [74.84, 81.30] | 2.56 [1.58, 4.12] | 48.16 [44.27, 52.08] |
| alternative_random | 69.12 [65.39, 72.62] | 69.12 [65.39, 72.62] | 7.20 [5.42, 9.50] | 48.00 [44.11, 51.92] |
| alternative_linked | 25.28 [22.03, 28.83] | 67.36 [63.59, 70.92] | 5.44 [3.92, 7.51] | 45.60 [41.73, 49.52] |
| alternative_reversed | 75.36 [71.83, 78.58] | 75.36 [71.83, 78.58] | 3.36 [2.21, 5.08] | 45.60 [41.73, 49.52] |
| alternative_row_structure | 60.16 [56.27, 63.92] | 60.16 [56.27, 63.92] | 43.84 [40.00, 47.76] | 69.60 [65.88, 73.08] |

Fixed panels use the realized-panel target; random panels use the genomic-distribution target. See summary.csv for secondary panel coverage.

## False genetic-direction support under the null

| Scenario | Locus | Block | Individual | Crossed |
|---|---:|---:|---:|---:|
| null_fixed_balanced | 6.50 [5.13, 8.20] | 6.50 [5.13, 8.20] | 2.80 [1.94, 4.02] | 0.00 [0.00, 0.38] |
| null_fixed_unequal | 83.20 [80.76, 85.39] | 83.20 [80.76, 85.39] | 99.90 [99.44, 99.98] | 98.10 [97.05, 98.78] |
| null_random_balanced | 6.30 [4.95, 7.98] | 6.30 [4.95, 7.98] | 2.00 [1.30, 3.07] | 0.00 [0.00, 0.38] |
| null_random_unequal | 80.80 [78.24, 83.12] | 80.80 [78.24, 83.12] | 99.40 [98.70, 99.72] | 97.30 [96.10, 98.14] |
| null_identical_balanced | 8.30 [6.75, 10.17] | 8.30 [6.75, 10.17] | 0.60 [0.28, 1.30] | 0.00 [0.00, 0.38] |
| null_identical_unequal | 48.50 [45.41, 51.60] | 48.50 [45.41, 51.60] | 82.80 [80.34, 85.01] | 40.10 [37.11, 43.17] |
| null_linked_balanced | 33.50 [30.64, 36.48] | 6.10 [4.78, 7.76] | 2.40 [1.62, 3.55] | 0.10 [0.02, 0.56] |
| null_linked_unequal | 98.30 [97.29, 98.94] | 82.30 [79.81, 84.54] | 99.50 [98.83, 99.79] | 97.20 [95.98, 98.06] |
| null_missing_balanced | 7.90 [6.38, 9.74] | 7.90 [6.38, 9.74] | 1.50 [0.91, 2.46] | 0.10 [0.02, 0.56] |
| null_missing_unequal | 89.90 [87.88, 91.62] | 89.90 [87.88, 91.62] | 99.80 [99.27, 99.95] | 99.50 [98.83, 99.79] |
| null_row_structure_balanced | 37.40 [34.45, 40.44] | 37.40 [34.45, 40.44] | 6.20 [4.87, 7.87] | 2.10 [1.38, 3.19] |
| null_row_structure_unequal | 67.70 [64.74, 70.53] | 67.70 [64.74, 70.53] | 91.00 [89.07, 92.62] | 74.80 [72.02, 77.39] |

Each null has 1,000 datasets. The null is zero genetic contrast; demographic migration is not simulated.

## Alternative direction counts

Cells are correct / wrong / unresolved among 625 datasets.

| Scenario | Locus | Block | Individual | Crossed |
|---|---:|---:|---:|---:|
| alternative_fixed | 625 / 0 / 0 | 625 / 0 / 0 | 625 / 0 / 0 | 625 / 0 / 0 |
| alternative_random | 625 / 0 / 0 | 625 / 0 / 0 | 625 / 0 / 0 | 625 / 0 / 0 |
| alternative_linked | 625 / 0 / 0 | 625 / 0 / 0 | 625 / 0 / 0 | 625 / 0 / 0 |
| alternative_reversed | 625 / 0 / 0 | 625 / 0 / 0 | 625 / 0 / 0 | 625 / 0 / 0 |
| alternative_row_structure | 625 / 0 / 0 | 625 / 0 / 0 | 625 / 0 / 0 | 625 / 0 / 0 |

Perfect sign recovery in these alternatives coexists with poor interval coverage. It does not validate migration inference or general performance.
