from skbio.io import create_format
from skbio.metadata._metadata import Metadata

sk_metadata = create_format("sk_metadata")


@sk_metadata.reader(Metadata)
def _sk_metadata_read(fh):
    return Metadata.load(fh)
