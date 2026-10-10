r"""Numba kernels over a BPTree (:mod:`skbio.tree.bp.numba`)
==========================================================

.. currentmodule:: skbio.tree.bp.numba

The navigation operations of :class:`~skbio.tree.BPTree` as functions that can
be called from your own compiled code: a Numba ``@njit`` function on the CPU,
or a ``@cuda.jit`` kernel on a GPU (``numba_cuda_mlir.cuda`` or ``numba.cuda`` on
NVIDIA, ``numba.hip`` on AMD).

A method such as ``bp.parent(i)``, called from Python, costs far more in
interpreter overhead than the navigation itself, and a batch method (e.g.
:meth:`~skbio.tree.BPTree.lca_batch`) spreads that overhead only over
independent queries, answered in one call. Neither expresses a computation
whose steps depend on each other, e.g. the lowest common ancestor of a *set* of
nodes, or a walk to the root. Inside a compiled loop, a navigation function
costs only the operation itself, and Python is not involved until the loop
returns.

Functions
---------

.. autosummary::
   :toctree: .

   cpu
   gpu

Classes
-------

.. autosummary::
   :toctree: .
   :template: namedtuple.rst

   Navigation

Usage
-----
The functions are compiled per target, and a tree's arrays are exported per
tree, so a kernel is written once and runs on any tree:

1. :func:`cpu` (or :func:`gpu`) returns the navigation functions as a named
   tuple, conventionally ``nav``. Your kernel calls them as ``nav.parent(T, i)``.
2. :meth:`BPTree.numba_arrays() <skbio.tree.BPTree.numba_arrays>` exports a
   tree's arrays as ``T``, which is passed to the kernel and on to each
   function. Node attributes are separate arrays:
   :attr:`~skbio.tree.BPTree.lengths` and :attr:`~skbio.tree.BPTree.edges`.
3. Nodes are positions in the parentheses (their opening parenthesis), as for
   the ``BPTree`` methods. :meth:`~skbio.tree.BPTree.tip_positions` and
   :meth:`~skbio.tree.BPTree.positions` give the positions of the tips and of
   named nodes, and :meth:`~skbio.tree.BPTree.name` maps a position back.

The lowest common ancestor of each sample's tips, where sample ``s`` owns
``tips[offsets[s]:offsets[s + 1]]``:

.. code-block:: python

    import numpy as np
    from numba import njit, prange
    from skbio.tree import BPTree
    from skbio.tree.bp import numba as bpn

    nav = bpn.cpu()

    @njit(parallel=True)
    def sample_lcas(T, tips, offsets):
        out = np.empty(offsets.size - 1, np.intp)
        for s in prange(out.size):
            a = tips[offsets[s]]
            for k in range(offsets[s] + 1, offsets[s + 1]):
                a = nav.lca(T, a, tips[k])
            out[s] = a
        return out

    bp = BPTree.read(["((a:1,b:2)c:3,(d:4,e:5)f:6)root;"])
    tips = bp.positions(["a", "b", "b", "d", "e"])
    offsets = np.array([0, 2, 5])  # sample 0 owns a and b; sample 1, b, d and e
    lcas = sample_lcas(bp.numba_arrays(), tips, offsets)
    [bp.name(i) for i in lcas]  # ['c', 'root']

The same on a GPU: ``T`` is uploaded once and cached on the tree, and the
kernel's other arrays must be on the device too. ``cuda.to_device`` copies a
NumPy array there, and ``cuda.as_cuda_array`` wraps a CuPy array or a PyTorch
tensor on the GPU without a copy (but see the warnings below on PyTorch). With
the deprecated numba-cuda instead of numba-cuda-mlir, import ``cuda`` from
``numba``.

.. code-block:: python

    from numba_cuda_mlir import cuda

    gnav = bpn.gpu(cuda)

    @cuda.jit
    def sample_lcas_gpu(T, tips, offsets, out):
        s = cuda.grid(1)
        if s < out.shape[0]:
            a = tips[offsets[s]]
            for k in range(offsets[s] + 1, offsets[s + 1]):
                a = gnav.lca(T, a, tips[k])
            out[s] = a

    n = offsets.size - 1
    out = cuda.device_array(n, dtype=np.intp)
    threads = 128
    blocks = (n + threads - 1) // threads
    sample_lcas_gpu[blocks, threads](
        bp.numba_arrays(gpu=cuda), cuda.to_device(tips), cuda.to_device(offsets), out
    )
    lcas = out.copy_to_host()

Conventions
-----------
Each function takes ``T`` first, then the arguments of the ``BPTree`` method
of the same name, and returns what that method returns, except that:

- ``minselect`` returns -1 where the method returns None;
- ``is_tip`` and ``is_ancestor`` return a boolean;
- ``count`` takes its ``tips`` argument positionally: ``count(T, i, tips)``.

The -1 is because a compiled function returns one type: an integer that is
sometimes None would need a check wherever it is used, and cannot be the
return value of a GPU device function at all. -1 is never a position.

So the sentinels are those of the methods:

- ``parent`` of the root is -1;
- ``first_child``, ``last_child``, ``next_sibling`` and ``previous_sibling``
  return 0 when there is no such node. 0 is the root's position, which is
  never a child or a sibling, so the value is unambiguous, but it is a
  position;
- ``level_ancestor`` returns -1 for ``d <= 0``, and the root for a ``d``
  beyond the node's depth: it stops at the root rather than signal that it
  overshot;
- ``minselect`` returns -1 when there is no ``q``-th minimum, including for
  ``q < 1``.

Warnings
--------
The functions do not check their arguments: a position outside
``[0, T.size)`` reads outside the tree's arrays. On a CPU that returns an
arbitrary result; on a GPU it can fault and leave the device unusable for the
rest of the process. The batch methods of ``BPTree`` check their queries.

Test a result for its sentinel before passing it on: a sentinel given to
another function does not fail, it answers for a different node.

- -1 indexes from the end, as in NumPy, so ``nav.depth(T, nav.parent(T, 0))``
  returns the excess at the last position (0), with no error. Numba's
  ``boundscheck=True`` does not catch it either, since -1 is in range once
  wrapped.
- 0 is the root, so ``nav.first_child(T, nav.first_child(T, tip))`` returns
  the root's first child rather than "none".

.. code-block:: python

    p = nav.parent(T, i)
    if p != -1:          # i is not the root
        ...
    c = nav.first_child(T, i)
    if c != 0:           # i is not a tip
        ...

A kernel runs on a stream of the Numba GPU module, not on one of the library
that made its arrays. CuPy tells Numba which stream produced an array, so the
kernel waits for it, but PyTorch does not: a kernel can read a tensor before
the work queued on PyTorch's stream has written it. Before launching a kernel
on PyTorch tensors, wait for that stream, e.g. with
``torch.cuda.current_stream().synchronize()``. ``T`` itself needs no wait.

``T`` is opaque: pass it through, but do not rely on its fields, which are
the tree's index and may change between versions. The exception is
``T.size``, the number of parentheses.

The first call of a kernel compiles it, together with the functions it
calls, which takes far longer than running it, the more so for ``lca``,
``rmq`` and the operations built on them. Numba compiles a kernel per
argument type, not per tree, so later trees reuse it.

Numba is an optional dependency: this module imports without it, and
:func:`cpu` and :func:`gpu` raise ``ImportError``.

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from collections import namedtuple

from . import _bp_numba

__all__ = ["Navigation", "cpu", "gpu"]


# The public functions: every navigation operation with a BPTree method of the
# same name. The index operations they are built from (rank, select, excess,
# open, enclose, fwdsearch, bwdsearch) stay private.
Navigation = namedtuple(
    "Navigation",
    [
        "root",
        "close",
        "depth",
        "parent",
        "is_tip",
        "first_child",
        "last_child",
        "next_sibling",
        "previous_sibling",
        "preorder_rank",
        "preorder_select",
        "postorder_rank",
        "postorder_select",
        "is_ancestor",
        "count",
        "level_ancestor",
        "level_next",
        "lca",
        "deepest_node",
        "height",
        "rmq",
        "rMq",
        "mincount",
        "minselect",
    ],
)
Navigation.__doc__ = """The navigation functions of a BPTree, for Numba code.

Each field is a compiled function ``f(T, ...)`` mirroring the BPTree method of
the same name; see :mod:`skbio.tree.bp.numba` for the conventions.
"""

# By target: "cpu", or the Numba GPU module's name.
_navigation = {}


def _require_numba():
    if not _bp_numba.NUMBA_AVAILABLE:
        raise ImportError(
            "The BPTree Numba functions require Numba, which is not installed."
        )


def _select(primitives):
    return Navigation(*(getattr(primitives, name) for name in Navigation._fields))


def cpu():
    """The navigation functions compiled for the CPU.

    Returns
    -------
    Navigation
        Functions callable from Numba ``@njit`` code (and, one call at a time,
        from Python), e.g. ``nav.parent(T, i)``. The same object is returned
        on every call.

    Raises
    ------
    ImportError
        If Numba is not installed.

    See Also
    --------
    gpu
    skbio.tree.BPTree.numba_arrays

    """
    _require_numba()
    if "cpu" not in _navigation:
        _navigation["cpu"] = _select(_bp_numba.CPU)
    return _navigation["cpu"]


def gpu(module):
    """The navigation functions compiled as device functions of a GPU.

    Parameters
    ----------
    module : module
        The Numba GPU module: ``numba_cuda_mlir.cuda`` or ``numba.cuda`` on
        NVIDIA, ``numba.hip`` on AMD.

    Returns
    -------
    Navigation
        Device functions, callable from the kernels of ``module`` (e.g.
        ``@cuda.jit``), e.g. ``nav.parent(T, i)``. Each is compiled with the
        first kernel that calls it. The same object is returned on every call
        with the same module.

    Raises
    ------
    ImportError
        If Numba is not installed.

    See Also
    --------
    cpu
    skbio.tree.BPTree.numba_arrays

    """
    _require_numba()
    name = module.__name__
    if name not in _navigation:
        _navigation[name] = _select(_bp_numba.gpu_primitives(module))
    return _navigation[name]
