
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import numpy as np
import pandas as pd
from scipy.sparse import issparse

from skbio.stats.distance import DistanceMatrix
from skbio.table._tabular import _create_table, _create_table_1d
from ._ordination_results import OrdinationResults


def graphembed(
    graph,
    method="sketching",
    dimensions=128,
    decay=0.3,
    nbiter=5,
    symetric=True,
    output_format=None,
    **kwargs
):
    r"""Perform graph embedding using graphembed_rs.

    This function wraps the ``graphembed_rs`` package, which must be installed
    separately. It calculates low-dimensional vector representations of the nodes
    in a graph, preserving high-order proximity or using recursive sketching.

    Parameters
    ----------
    graph : DistanceMatrix, DataFrame, ndarray, or scipy.sparse
        The adjacency matrix of the graph. Positive values represent edge weights.
        If a DistanceMatrix is provided, distances are inversely related to affinity,
        so care must be taken (e.g., converting to affinity/similarity first is
        recommended for optimal embedding).
    method : str, optional
        The embedding algorithm to use. Available options are ``"sketching"``
        (recursive sketching) or ``"hope"`` (High-Order Proximity preserved Embedding).
        Default is ``"sketching"``.
    dimensions : int, optional
        The number of dimensions for the resulting embedding. Default is 128.
    decay : float, optional
        The decay factor for sketching. Used only when method is ``"sketching"``.
        Default is 0.3.
    nbiter : int, optional
        The number of iterations/hops to consider. Default is 5.
    symetric : bool, optional
        Whether the graph is symmetric. Used only when method is ``"sketching"``.
        Default is True.
    output_format : optional
        Standard table parameters. See :ref:`table_params` for details.
    **kwargs : dict
        Additional arguments passed to the underlying `graphembed_rs` functions.

    Returns
    -------
    OrdinationResults
        Object that stores the embedding results.

    Raises
    ------
    ImportError
        If `graphembed_rs` is not installed.

    Notes
    -----
    This function relies on writing the graph to a temporary file via an edge list
    format before invoking the Rust-based execution engine.

    See Also
    --------
    skbio.stats.ordination.OrdinationResults

    """
    try:
        import graphembed_rs as ge
    except ImportError:
        raise ImportError(
            "The graphembed function requires the graphembed_rs package. "
            "Please install it using 'pip install graphembed_rs'."
        )


    if isinstance(graph, DistanceMatrix):
        ids = graph.ids
        adj = graph.data
    elif isinstance(graph, pd.DataFrame):
        ids = graph.index.astype(str).tolist()
        adj = graph.values
    else:
        num_nodes = graph.shape[0]
        ids = [str(i) for i in range(num_nodes)]
        adj = graph


    with tempfile.NamedTemporaryFile(mode="w+", suffix=".txt", delete=False) as tf:
        tmp_filepath = tf.name

        if issparse(adj):
            adj_coo = adj.tocoo()
            for r, c, v in zip(adj_coo.row, adj_coo.col, adj_coo.data):
                if v != 0:
                    tf.write(f"{ids[r]}\t{ids[c]}\t{v}\n")
        else:
            # Write dense matrix
            n = adj.shape[0]
            for r in range(n):
                for c in range(n):
                    v = adj[r, c]
                    if v != 0:
                        tf.write(f"{ids[r]}\t{ids[c]}\t{v}\n")

    try:
        if method == "sketching":
            with tempfile.TemporaryDirectory() as td:
                out_path = os.path.join(td, "embedded")
                ge.embed_sketching(
                    tmp_filepath,
                    decay=decay,
                    dim=dimensions,
                    nbiter=nbiter,
                    symetric=symetric,
                    output=out_path,
                    **kwargs
                )

                out_file = out_path if os.path.isfile(out_path) else f"{out_path}.txt"
                if not os.path.exists(out_file):
                    if os.path.exists(out_path):
                        out_file = out_path
                    else:
                        raise RuntimeError(
                            "graphembed_rs did not produce an output file."
                        )

                res_df = pd.read_csv(out_file, sep=r"\s+", header=None, index_col=0)
                res_df.index = res_df.index.astype(str)
                res_df = res_df.reindex(ids).fillna(0)
                coordinates = res_df.values

        elif method == "hope":
            with tempfile.TemporaryDirectory() as td:
                res = ge.embed_hope_rank(
                    tmp_filepath,
                    target_rank=dimensions,
                    nbiter=nbiter,
                    **kwargs
                )

                if res is None:
                    raise RuntimeError("HOPE embedding didn't return in memory.")

                if isinstance(res, dict):
                    coords = [res.get(node_id, np.zeros(dimensions)) for node_id in ids]
                    coordinates = np.array(coords)
                elif isinstance(res, np.ndarray):
                    coordinates = res
                else:

                    coordinates = np.array(res)

        else:
            raise ValueError(f"Unknown graphembed method: {method}")

    finally:
        if os.path.exists(tmp_filepath):
            os.remove(tmp_filepath)

    long_method_name = f"Graph Embedding via {method.capitalize()}"
    axis_labels = [f"Dim{i}" for i in range(1, coordinates.shape[1] + 1)]

    pseudo_eigvals = np.zeros(coordinates.shape[1])
    pseudo_prop = np.zeros(coordinates.shape[1])

    return OrdinationResults(
        short_method_name="GraphEmbed",
        long_method_name=long_method_name,
        eigvals=_create_table_1d(
            pseudo_eigvals, index=axis_labels, backend=output_format
        ),
        samples=_create_table(
            coordinates,
            index=ids,
            columns=axis_labels,
            backend=output_format,
        ),
        proportion_explained=_create_table_1d(
            pseudo_prop, index=axis_labels, backend=output_format
        ),
    )
