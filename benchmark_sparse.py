"""Benchmark: sparse vs dense vectorize_counts_and_tree.

Measures runtime and peak memory for varying tree sizes and sparsity levels,
then saves plots to benchmark_results/.
"""

import time
import tracemalloc
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import os

from skbio import TreeNode
from skbio.diversity._util import vectorize_counts_and_tree


def create_star_tree(n_tips):
    """Create a simple star tree with n_tips tips."""
    newick = "({})root;".format(",".join(f"t{i}:1.0" for i in range(n_tips)))
    return TreeNode.read([newick])


def measure(counts, taxa, tree):
    """Return (time_seconds, peak_memory_MB) for one call."""
    tracemalloc.start()
    t0 = time.perf_counter()
    result, _, _ = vectorize_counts_and_tree(counts, taxa, tree)
    elapsed = time.perf_counter() - t0
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return elapsed, peak / 1024 / 1024, result


def run_scaling_benchmark():
    """Benchmark across tree sizes at fixed sparsity."""
    n_samples = 50
    tip_sizes = [500, 1000, 2000, 5000]
    sparsity = 0.05

    results = {
        "tips": [],
        "dense_time": [],
        "sparse_time": [],
        "dense_mem": [],
        "sparse_mem": [],
    }

    print("=== Scaling benchmark (sparsity=5%) ===")
    for n_tips in tip_sizes:
        print(f"  n_tips={n_tips} ...", end=" ", flush=True)
        tree = create_star_tree(n_tips)
        taxa = np.array([f"t{i}" for i in range(n_tips)])

        rng = np.random.default_rng(42)
        counts_dense = (rng.random((n_samples, n_tips)) < sparsity).astype(np.float64)
        counts_sparse = sp.csr_matrix(counts_dense)

        dt, dm, res_d = measure(counts_dense, taxa, tree)
        st, sm, res_s = measure(counts_sparse, taxa, tree)
        np.testing.assert_allclose(res_d, res_s)

        results["tips"].append(n_tips)
        results["dense_time"].append(dt)
        results["sparse_time"].append(st)
        results["dense_mem"].append(dm)
        results["sparse_mem"].append(sm)

        print(f"dense={dt:.3f}s/{dm:.1f}MB  sparse={st:.3f}s/{sm:.1f}MB")

    return results


def run_sparsity_benchmark():
    """Benchmark across sparsity levels at fixed tree size."""
    n_tips = 2000
    n_samples = 50
    sparsity_levels = [0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001]

    results = {
        "sparsity": [],
        "dense_time": [],
        "sparse_time": [],
        "dense_mem": [],
        "sparse_mem": [],
    }

    print("\n=== Sparsity benchmark (n_tips=2000) ===")
    for sparsity in sparsity_levels:
        print(f"  sparsity={sparsity} ...", end=" ", flush=True)
        tree = create_star_tree(n_tips)
        taxa = np.array([f"t{i}" for i in range(n_tips)])

        rng = np.random.default_rng(42)
        counts_dense = (rng.random((n_samples, n_tips)) < sparsity).astype(np.float64)
        counts_sparse = sp.csr_matrix(counts_dense)

        dt, dm, res_d = measure(counts_dense, taxa, tree)
        st, sm, res_s = measure(counts_sparse, taxa, tree)
        np.testing.assert_allclose(res_d, res_s)

        results["sparsity"].append(sparsity)
        results["dense_time"].append(dt)
        results["sparse_time"].append(st)
        results["dense_mem"].append(dm)
        results["sparse_mem"].append(sm)

        nnz_pct = (counts_sparse.nnz / (n_samples * n_tips)) * 100
        print(
            f"dense={dt:.3f}s/{dm:.1f}MB  "
            f"sparse={st:.3f}s/{sm:.1f}MB  nnz={nnz_pct:.1f}%"
        )

    return results


def plot_results(scaling, sparsity_res, outdir="benchmark_results"):
    os.makedirs(outdir, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Top-left: Runtime vs tree size
    ax = axes[0, 0]
    ax.plot(scaling["tips"], scaling["dense_time"], "o-", label="Dense", linewidth=2)
    ax.plot(scaling["tips"], scaling["sparse_time"], "s-", label="Sparse", linewidth=2)
    ax.set_xlabel("Number of tree tips")
    ax.set_ylabel("Runtime (seconds)")
    ax.set_title("Runtime vs Tree Size (sparsity=5%)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Top-right: Memory vs tree size
    ax = axes[0, 1]
    ax.plot(scaling["tips"], scaling["dense_mem"], "o-", label="Dense", linewidth=2)
    ax.plot(scaling["tips"], scaling["sparse_mem"], "s-", label="Sparse", linewidth=2)
    ax.set_xlabel("Number of tree tips")
    ax.set_ylabel("Peak memory (MB)")
    ax.set_title("Peak Memory vs Tree Size (sparsity=5%)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Bottom-left: Runtime vs sparsity
    ax = axes[1, 0]
    ax.plot(
        sparsity_res["sparsity"],
        sparsity_res["dense_time"],
        "o-",
        label="Dense",
        linewidth=2,
    )
    ax.plot(
        sparsity_res["sparsity"],
        sparsity_res["sparse_time"],
        "s-",
        label="Sparse",
        linewidth=2,
    )
    ax.set_xlabel("Sparsity (fraction of non-zero entries)")
    ax.set_ylabel("Runtime (seconds)")
    ax.set_title("Runtime vs Sparsity (n_tips=2000)")
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Bottom-right: Memory vs sparsity
    ax = axes[1, 1]
    ax.plot(
        sparsity_res["sparsity"],
        sparsity_res["dense_mem"],
        "o-",
        label="Dense",
        linewidth=2,
    )
    ax.plot(
        sparsity_res["sparsity"],
        sparsity_res["sparse_mem"],
        "s-",
        label="Sparse",
        linewidth=2,
    )
    ax.set_xlabel("Sparsity (fraction of non-zero entries)")
    ax.set_ylabel("Peak memory (MB)")
    ax.set_title("Peak Memory vs Sparsity (n_tips=2000)")
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Sparse vs Dense: vectorize_counts_and_tree Performance",
        fontsize=14,
        fontweight="bold",
        y=1.01,
    )
    fig.tight_layout()
    fig.savefig(
        os.path.join(outdir, "benchmark_comparison.png"), dpi=150, bbox_inches="tight"
    )
    plt.close(fig)

    print(f"\nPlot saved to {outdir}/benchmark_comparison.png")


if __name__ == "__main__":
    scaling = run_scaling_benchmark()
    sparsity_res = run_sparsity_benchmark()
    plot_results(scaling, sparsity_res)
