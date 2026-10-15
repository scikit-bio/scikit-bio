# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from collections.abc import Hashable
from collections.abc import Sequence as PySequence
from typing import Any, TypeAlias

from numpy.typing import NDArray

from ._sequence import Sequence


# ------------------------------------------------
# SequenceLike : a finite, ordered sequence of hashable symbols
# ------------------------------------------------

SequenceLike: TypeAlias = (
    Sequence | str | bytes | bytearray | PySequence[Hashable] | NDArray[Any]
)
