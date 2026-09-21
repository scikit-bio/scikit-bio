Computation and Performance
===========================

Biological datasets continue to grow in size and complexity. scikit-bio prioritizes
computational efficiency and scalability to very large datasets, alongside
numerical correctness and reproducibility. Broad support for computer architectures
and environments is also a central goal. The aim is to help researchers carry out
reliable analyses on available resources, whether a laptop or a supercomputer cluster.

Most analyses can start with the default settings, which execute efficient, CPU-based
numerical computing using the vectorized array operations of `NumPy
<https://numpy.org/>`_ and, if necessary, `Cython <https://cython.org/>`_ for
loop-intensive calculations. For more demanding work, scikit-bio offers alternative
implementations of selected analyses, involving parallelization, `Numba
<https://numba.pydata.org/>`_ engines, support for GPU computation, and support for
alternative array backends. This guide explains these options and how to choose among
them.


.. _array_backends:

Array backends
--------------

Multiple array libraries are available or emerging in the Python scientific computing
ecosystem to support diverse computational requirements. The `Python array API standard
<https://data-apis.org/array-api/latest/>`_ defines a common API with which many array
libraries comply. An increasing number of scikit-bio functions are implemented using
array-API-compatible code to provide **native support** for multiple array backends.
This allows supported computations to use the corresponding backend without converting
the input array to NumPy on the CPU. Support and any limitations depend on the function
and backend, and are noted in each function's documentation. The output will also use
the same format when applicable.

We will demonstrate this using the :func:`~skbio.stats.composition.clr` function, which
performs the centered log-ratio (CLR) transformation of compositional data. A NumPy array
is created and supplied to the function, which returns a new NumPy array containing the
result::

    import numpy as np
    from skbio.stats.composition import clr

    rng = np.random.default_rng(42)
    arr = rng.integers(1, 1000, size=(1000, 10000))
    result = clr(arr)

Let's now use `JAX <https://docs.jax.dev/en/latest/>`_ as the array backend. One just
needs to cast the NumPy array into a JAX array before the function call. This simple
change can significantly improve performance, depending on the workload and computing
environment::

    import jax.numpy as jnp

    jarr = jnp.asarray(arr)
    result = clr(jarr)

While any compliant array backend may work automatically, scikit-bio focuses on
validated support for five common backends. A table of backend and device support is
provided on the documentation page of each function that supports array backends.
For example, the CLR function supports the following combinations:

+---------+---------+---------+
| Backend | CPU     | GPU     |
+=========+=========+=========+
| NumPy   | |check| | n/a     |
+---------+---------+---------+
| CuPy    | n/a     | |check| |
+---------+---------+---------+
| PyTorch | |check| | |check| |
+---------+---------+---------+
| JAX     | |check| | |check| |
+---------+---------+---------+
| Dask    | |check| | n/a     |
+---------+---------+---------+

scikit-bio is not dependent on any of these array libraries except for NumPy. To
utilize a particular array backend that best fits your task and computational resources,
you will need to *install that library separately*.

See below for :ref:`GPU computing <gpu_computing>` using array backends.


.. _compute_engines:

Compute engines
---------------

Some scikit-bio functions offer more than one implementation to perform the same task.
These tasks often involve loop-intensive computations that benefit from compiled code.
Each implementation is referred to as a **compute engine**. Such functions have an
``engine`` parameter letting you choose between compute engines. Three options are
currently available:

- ``cython``: The default engine implemented in `Cython <https://cython.org/>`_,
  which is translated into C code and compiled into binaries when scikit-bio is built. It
  often offers C-level performance in heavy computing tasks. Cython is the legacy
  engine of many scikit-bio functions. Pure Python implementations are also referred
  to as ``cython`` for simplicity.

- ``numba``: An alternative engine implemented in `Numba <https://numba.pydata.org/>`_,
  which stays as source code at deployment and is only compiled before execution
  (just-in-time (JIT) compilation). This process adds overhead to the first call of
  the function. Once compiled, the binaries are cached for reuse in subsequent calls,
  avoiding that overhead. Numba is not a required dependency of scikit-bio. To utilize
  the Numba engine, you will need to :install:`install Numba <#numba>` separately.

- ``fast``: Let the function choose the faster available engine for you. This selection
  is determined based on the development team's benchmarks on representative datasets,
  though it does not guarantee faster execution for every dataset. Each function's
  documentation page explains the choice and its impact.

.. note::
    The scikit-bio development team is currently expanding Numba engines for valuable
    functions. They are often more efficient than the legacy Cython engines.

Both Cython and Numba engines of each function perform the same task and the results
are usually consistent, although numerical identity is not guaranteed due to the
different numerical and stochastic behaviors between the engines.

Without explicitly specifying an engine in the function call (i.e., ``engine=None``),
the choice of an engine will be controlled by the :ref:`global configuration
<configuration>` option ``compute_engine``, which defaults to ``cython``. If you want
to apply one choice (e.g., ``numba``) throughout a workflow, you just need to set the
global option once prior to execution::

    from skbio import set_config

    set_config('compute_engine', 'numba')


.. _gpu_computing:

GPU computing
-------------

GPUs are often advantageous over CPUs in large-scale analyses as they enable massively
parallel calculations. Support for GPU computing in scikit-bio is provided via two
mechanisms:

1\. :ref:`Array backends <array_backends>` provide backend-specific GPU computing
support out-of-the-box. Let's continue with the CLR example to demonstrate this. We
will use `CuPy <https://cupy.dev/>`_, a library featuring CUDA GPU-accelerated array
computing. By converting a NumPy array into a CuPy array, the data is moved to the GPU,
and the subsequent CLR transformation will take place on the GPU automatically::

    import cupy as cp

    carr = cp.asarray(arr)
    result = clr(carr)

For backends that support both CPU and GPU, one may need to specify which device to
store the data on using the ``device`` parameter. The following example uses `PyTorch
<https://pytorch.org/>`_, a common deep learning library supporting both CPU- and
GPU-resident arrays (tensors)::

    import torch

    assert torch.cuda.is_available()
    tensor = torch.tensor(arr, device='cuda')
    result = clr(tensor)

2\. The :ref:`Numba engines <compute_engines>` of some functions are capable of GPU
computing through the extensions `numba-cuda <https://nvidia.github.io/numba-cuda/>`_
(for CUDA GPUs) and `numba-hip <https://github.com/ROCm/numba-hip>`_ (for ROCm GPUs).
Refer to the :install:`installation instructions <#numba>`. With a compatible GPU
extension installed, these Numba engines can automatically use GPU computation for
supported CuPy and PyTorch arrays residing on the GPU. The function will fall back to
the array-API implementation for other array formats or if the GPU kernel is
unavailable.

We will demonstrate this using the :func:`~skbio.stats.distance.permanova` function,
which performs the Permutational Multivariate Analysis of Variance (PERMANOVA) on a
distance matrix and a grouping vector. By default, the input distance matrix is a
CPU-resident NumPy array and the compute engine is Cython (see compute_engines_)::

    from skbio.stats.distance import permanova, randdm

    dm = randdm(10000)
    grouping = [0] * 5000 + [1] * 5000
    res = permanova(dm, grouping)

By specifying ``engine='numba'``, the computation will be performed by a Numba engine,
but still using the CPU::

    res = permanova(dm, grouping, engine='numba')

Now we will move the distance matrix data to the GPU using CuPy. This is supported by
scikit-bio's :class:`~skbio.stats.distance.DistanceMatrix` class::

    import cupy as cp
    from skbio.stats.distance import DistanceMatrix

    cdm = DistanceMatrix(cp.asarray(dm.data), dm.ids)

Then call the ``permanova`` function with ``engine='numba'``. With a compatible GPU,
the Numba kernel and some additional array operations run on the GPU, and the distance
matrix avoids a CPU round-trip. This can provide significant performance gains,
depending on the dataset and computing environment::

    res = permanova(cdm, grouping, engine='numba')


.. _parallelization:

Parallelization
---------------

Most computer systems have multiple CPU cores that can work simultaneously.
**Parallelization** divides a calculation among the CPU cores. Some scikit-bio
functions do this directly. Others delegate to upstream libraries that perform parts of
the work in parallel (e.g., NumPy parallelizes some linear algebra calculations using
BLAS). These happen automatically, without changing any of your code. By default,
parallel functions in scikit-bio utilize all available CPU cores when possible.

Sometimes you may want to limit how many resources an analysis uses, particularly on a
shared computer or when running several analyses at once. CPU work is organized into
**threads**, which are separate streams of work within a program. Currently there
is no per-function parameter to control the number of threads used. However, you can
still control this behavior using the approaches described below.

Parallel :ref:`Cython engines <compute_engines>` use `OpenMP <https://www.openmp.org/>`_
to manage threads. You can control the number of OpenMP threads used by a Python program
by setting the ``OMP_NUM_THREADS`` environment variable. For example, the following
command sets the default OpenMP thread count to four for ``script.py``::

    OMP_NUM_THREADS=4 python script.py

Alternatively, this can be configured inside a Python script. Note that this must be
set before the ``import`` statement of whichever functionality you will be using::

    import os
    os.environ['OMP_NUM_THREADS'] = '4'

    from skbio import some_function

More granular control can be achieved using `threadpoolctl
<https://github.com/joblib/threadpoolctl>`_, which needs to be installed separately.
This lets you limit the number of threads for individual code blocks, including those
used by upstream libraries such as NumPy and SciPy::

    from skbio import some_function
    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=4):
        some_function()

See the :install:`installation instructions <#parallelization>` for building scikit-bio
with or without OpenMP support.

Parallel :ref:`Numba engines <compute_engines>` manage threads through Numba itself.
One may specify the number of threads by setting the ``NUMBA_NUM_THREADS`` environment
variable::

    NUMBA_NUM_THREADS=4 python script.py

This can be combined with ``OMP_NUM_THREADS`` if your program combines Cython and Numba
engines::

    OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4 python script.py

Within Python, ``numba.set_num_threads`` can configure Numba's maximum thread count
before Numba is imported. See the `Numba thread controls
<https://numba.readthedocs.io/en/stable/user/threading-layer.html#setting-the-number-of-threads>`_
for details.

Note that the approaches discussed above do not provide a universal limit for GPU
execution.


.. _binary_acceleration:

Optional binary acceleration
----------------------------

`scikit-bio-binaries <https://github.com/scikit-bio/scikit-bio-binaries>`_ is a
separate package written in C++, providing optimized implementations for selected
analyses such as PERMANOVA and PCoA. scikit-bio can use these automatically when the
package is installed and the calculation's inputs are supported. See the
:install:`installation instructions <#binary-acceleration>` for setup. Each function's
Notes describe where it can be used.

This option is only available when the Cython engine is selected (including when
``fast`` selects Cython), and an eligible calculation may use the external
implementation. Selecting Numba bypasses this package.


.. |check| unicode:: U+2713
   :trim:
