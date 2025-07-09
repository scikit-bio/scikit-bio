scikit-bio roadmap
==================
This is a dynamic document tracking current development directions. It may be used for inspiration and encouragement and for coordination of efforts.

Type annotations
----------------
Progress has been made in adding type annotations to the scikit-bio codebase. Future development is encouraged to use type annotations by default, and existing modules will be updated on an individual basis. The benefits of type annotations include code readability, improved IDE autocompletion, easier bug detection, and checking correctness in downstream applications.

Interoperability
----------------
We are systematically adding type annotations to the scikit-bio codebase to improve code readability, IDE autocompletion, and bug detection in both development and downstream applications. Recent progress includes type annotations for distance statistics, alignment modules, and other core functionality shipped in v0.7.0. New code development incorporates type annotations by default, while existing modules are being updated incrementally. This ongoing effort will enhance the developer experience and facilitate more robust integration with downstream tools like QIIME.

GPU support
-----------
Modern omics datasets often exceed memory limitations and benefit from hardware acceleration. To address this, we are adopting the Python array API standard, which enables seamless interoperability with GPU-accelerated libraries like CuPy and JAX, as well as distributed computing frameworks like Dask. Initial implementation has focused on compositional transformation algorithms, with successful testing on PyTorch GPU objects. This foundation will be extended to additional computational bottlenecks including UniFrac, PERMANOVA, mantel, and permdisp, enabling significant performance improvements for large-scale analyses.

Support for free-threaded Python
--------------------------------
CPython 3.13 introduces a free-threaded build that removes the Global Interpreter Lock, enabling true multi-threading for CPU-bound tasks. This presents significant opportunities for scikit-bio's computationally intensive operations. However, adoption depends on ecosystem readiness - key dependencies including statsmodels and h5py do not yet provide free-threaded builds. We will monitor dependency support and evaluate implementation once the required packages become available, with particular focus on algorithms that would benefit most from parallel execution.

Benchmarking
------------
Performance monitoring is critical for preventing regression in computationally intensive bioinformatics workflows. While scikit-bio previously had lightweight benchmarks using asv, these were not maintained and became outdated. We have established a dedicated benchmarking repository (https://github.com/scikit-bio/scikit-bio-benchmarks) to systematically track performance across releases, starting with version 0.6.1 and continuing for all future releases. The current framework includes basic benchmarks with plans to expand coverage to key algorithms and data structures. This infrastructure will help ensure that performance improvements and potential regressions are identified early in the development process.
