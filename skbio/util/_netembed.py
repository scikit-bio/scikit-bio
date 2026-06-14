# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------
"""scikit-bio utility mixin for graph-/network-embedding.

This module exposes a :class:`NetEmbedMixin` that mirrors the pattern used in
:pyfile:`skbio/util/_plotting.py`.  The mixin provides **lazy access** to the
optional Rust extension package *graphembed-rs*, keeping the dependency strictly
optional at import time.
"""

from __future__ import annotations

import importlib
import os
from types import ModuleType
from typing import Any, List


class NetEmbedMixin:
    """Mixin granting on-demand access to *graphembed_rs* network embeddings."""

    def _get_ge(self) -> None:
        """Import *graphembed_rs* once and cache the modules on ``self``.

        Raises
        ------
        ImportError
            If *graphembed_rs* is unavailable.
        """
        err_msg = (
            "Network embedding requires the optional dependency "
            "`graphembed_rs` (>= 0.1.2).  Install with\n\n"
            "    pip install graphembed_rs\n\n"
            "or build from source with `maturin develop`."
        )

        if hasattr(self, "_ge"):
            if self._ge is None:
                raise ImportError(err_msg)
            return

        try:
            os.environ.setdefault("RUST_LOG", "info")
            ge: ModuleType = importlib.import_module("graphembed_rs.graphembed_rs")
            ge_utils: ModuleType = importlib.import_module("graphembed_rs.load_utils")
        except ImportError as e:
            self._ge = None
            self._ge_utils = None
            raise ImportError(err_msg) from e
        else:
            self._ge = ge
            self._ge_utils = ge_utils

    def embed_hope(
        self,
        edge_list: str,
        *,
        target_rank: int = 128,
        nbiter: int = 4,
        output: str | None = "embedding_output",
    ) -> str:
        """Embed *edge_list* using HOPE; returns the path to *output*."""
        self._get_ge()
        return self._ge.embed_hope_rank(
            edge_list,
            target_rank=target_rank,
            nbiter=nbiter,
            output=output,
        )

    def embed_sketching(
        self,
        edge_list: str,
        *,
        decay: float = 0.3,
        dim: int = 128,
        nbiter: int = 5,
        symetric: bool = True,
        output: str | None = "embedding_output",
    ) -> str:
        """Embed *edge_list* with NodeSketch; returns the path to *output*."""
        self._get_ge()
        return self._ge.embed_sketching(
            edge_list,
            decay=decay,
            dim=dim,
            nbiter=nbiter,
            symetric=symetric,
            output=output,
        )

    def validate_sketching(
        self,
        edge_list: str,
        *,
        decay: float = 0.3,
        dim: int = 128,
        nbiter: int = 3,
        nbpass: int = 1,
        skip_frac: float = 0.2,
        symetric: bool = True,
        centric: bool = True,
    ) -> List[float]:
        """Return AUC scores (list of floats) for *nbpass* validation rounds."""
        self._get_ge()
        return self._ge.validate_sketching(
            edge_list,
            decay=decay,
            dim=dim,
            nbiter=nbiter,
            nbpass=nbpass,
            skip_frac=skip_frac,
            symetric=symetric,
            centric=centric,
        )

    def load_embedding(self, bson_file: str) -> Any:
        """Load an embedding (NumPy array) produced by the above helpers."""
        self._get_ge()
        return self._ge_utils.load_embedding_bson(bson_file)

    @property
    def ge(self) -> ModuleType:
        """The underlying :pyobj:`graphembed_rs.graphembed_rs` module."""
        self._get_ge()
        return self._ge

    @property
    def ge_utils(self) -> ModuleType:
        """The underlying :pyobj:`graphembed_rs.load_utils` helpers."""
        self._get_ge()
        return self._ge_utils
