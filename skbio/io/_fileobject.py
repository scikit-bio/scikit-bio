# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import h5py


def is_binary_file(file):
    return isinstance(
        file, (io.BufferedReader, io.BufferedWriter, io.BufferedRandom, h5py.Group)
    )


# Everything beyond this point will be some kind of hack needed to make
# everything work. It's not pretty and it doesn't make great sense much
# of the time. I am very sorry to the poor soul who has to read beyond.


class FlushDestructorMixin:
    def __del__(self):
        # By default, the destructor calls close(), which flushes and closes
        # the underlying buffer. Override to only flush.
        if not self.closed:
            self.flush()


class SaneTextIOWrapper(FlushDestructorMixin, io.TextIOWrapper):
    pass


class WrappedBufferedRandom(FlushDestructorMixin, io.BufferedRandom):
    pass


class CompressedMixin(FlushDestructorMixin):
    """Act as a bridge between worlds."""

    def __init__(self, before_file, *args, **kwargs):
        self.streamable = kwargs.pop("streamable", True)
        self._before_file = before_file
        super(CompressedMixin, self).__init__(*args, **kwargs)

    @property
    def closed(self):
        return self.raw.closed or self._before_file.closed

    def close(self):
        super(CompressedMixin, self).close()

        # The above will not usually close before_file. We want the
        # decompression to be transparent, so we don't want users to deal with
        # this edge case. Instead we can just close the original now that we
        # are being closed.
        self._before_file.close()


class CompressedBufferedReader(CompressedMixin, io.BufferedReader):
    pass


class CompressedBufferedWriter(CompressedMixin, io.BufferedWriter):
    pass


class IterableStringReaderIO(io.StringIO):
    def __init__(self, iterable, newline):
        self._iterable = iterable
        super(IterableStringReaderIO, self).__init__("".join(iterable), newline=newline)


class IterableStringWriterIO(IterableStringReaderIO):
    def close(self):
        if not self.closed:
            backup = self.tell()
            self.seek(0)
            for line in self:
                self._iterable.append(line)
            self.seek(backup)
        super(IterableStringWriterIO, self).close()
