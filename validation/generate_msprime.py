"""Independent coalescent benchmark; run from repository root.

Requires validation/requirements-msprime.txt. No snpflowR simulation code is used.
Ancestry blocks are independent replicates; recombination is simulated within them.
"""
import argparse
import csv
import json
import platform
from pathlib import Path

import msprime
import numpy as np
import tskit


def reference_gst(p, q):
    """Direct two-allele heterozygosity arithmetic, independent of R code."""
    a = np.column_stack([p, 1 - p])
    b = np.column_stack([q, 1 - q])
    pool = np.sqrt(a * b)
    usable = (pool.sum(axis=1) > 0) & ~((p == q) & ((p == 0) | (p == 1)))
    a, pool = a[usable], pool[usable]
    pool = pool / pool.sum(axis=1, keepdims=True)
    hs = (2 - (a*a).sum(axis=1) - (pool*pool).sum(axis=1)) / 2
    ht = 1 - (((a+pool)/2)**2).sum(axis=1)
    return float((ht-hs).sum() / ht.sum()) if ht.sum() else float("nan")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--output", type=Path, default=Path("validation/generated"))
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("replicates must be positive")
    args.output.mkdir(parents=True, exist_ok=True)
    results = Path("validation/results")
    results.mkdir(exist_ok=True)
    scenarios = [("coalescent_null", 0, 0), ("coalescent_symmetric", .001, .001),
                 ("coalescent_A_to_B", .002, 0), ("coalescent_B_to_A", 0, .002)]
    manifest = []
    for s, (name, m_ab, m_ba) in enumerate(scenarios):
        demography = msprime.Demography()
        for pop in ["A", "B", "ancestor"]:
            demography.add_population(name=pop, initial_size=1000)
        # msprime rates describe backward lineage movement. Reverse labels
        # to encode the forward fraction of the destination replaced by migrants.
        demography.set_migration_rate(source="B", dest="A", rate=m_ab)
        demography.set_migration_rate(source="A", dest="B", rate=m_ba)
        assert demography.migration_matrix[1, 0] == m_ab
        assert demography.migration_matrix[0, 1] == m_ba
        demography.add_population_split(time=1000, derived=["A", "B"], ancestral="ancestor")
        demography.validate()
        for rep in range(1, args.replicates + 1):
            seed = 3_000_000 + s*100_000 + rep
            ancestry = msprime.sim_ancestry(
                samples={"A": 40, "B": 40}, ploidy=2, demography=demography,
                sequence_length=10000, recombination_rate=1e-8,
                num_replicates=30, model="hudson", random_seed=seed,
            )
            records = []
            for block, ts in enumerate(ancestry):
                ts = msprime.sim_mutations(ts, rate=1e-8,
                    model=msprime.BinaryMutationModel(), random_seed=seed + (block+1)*1000)
                sample_pop = np.array([ts.node(u).population for u in ts.samples()])
                assert (sample_pop == 0).sum() == 80 and (sample_pop == 1).sum() == 80
                for variant in ts.variants():
                    gt = variant.genotypes
                    if not np.all((gt == 0) | (gt == 1)):
                        raise ValueError("Unexpected missing or multiallelic genotype")
                    if gt.min() == gt.max():
                        continue
                    a, b = int(gt[sample_pop == 0].sum()), int(gt[sample_pop == 1].sum())
                    records.append([f"b{block}_s{variant.site.id}", str(block), a, 80, b, 80])
            if not records:
                raise ValueError("No variants: benchmark must record rather than replace this replicate")
            file = args.output / f"{name}_{rep:03d}.csv"
            with file.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["locus", "block", "alt_a", "called_a", "alt_b", "called_b"])
                writer.writerows(records)
            p = np.array([r[2] / r[3] for r in records])
            q = np.array([r[4] / r[5] for r in records])
            g_ab, g_ba = reference_gst(p, q), reference_gst(q, p)
            manifest.append(dict(scenario=name, replicate=rep, ancestry_seed=seed,
                file=str(file), truth="A -> B" if m_ab > m_ba else "B -> A" if m_ba > m_ab else "null",
                m_ab=m_ab, m_ba=m_ba, n_loci=len(records), n_blocks=30,
                reference_gst_ab=g_ab, reference_gst_ba=g_ba))
        print(f"{name}: generated {args.replicates} independent datasets", flush=True)
    with (results / "msprime-manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest[0].keys())
        writer.writeheader()
        writer.writerows(manifest)
    protocol = dict(msprime=msprime.__version__, tskit=tskit.__version__,
        numpy=np.__version__, python=platform.python_version(), replicates=args.replicates,
        blocks=30, block_length=10000, recombination_rate=1e-8, mutation_rate=1e-8,
        effective_sizes=[1000, 1000], ancestral_size=1000, split_generations=1000,
        n_ind=[40, 40], model="Hudson coalescent; binary mutation; segregating sites retained",
        rate_convention="forward destination-replacement fraction approximated by backward lineage rate",
        limitations="Limited balanced-sample grid; no validation of demographic identifiability or CI coverage")
    (results / "msprime-protocol.json").write_text(json.dumps(protocol, indent=2)+"\n")


if __name__ == "__main__":
    main()
