# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.

from skbio import DNA, RNA
import numpy as np

num_bases = 1000000
size = int(num_bases / 4)
short_len = 100

dna_template_bytes = [ord(x) for x in 'ACGT']
dna_template_bytes_gapped = [ord(x) for x in 'AC-.']
rna_template_bytes = [ord(x) for x in 'ACGU']

dna_bytes = np.array(dna_template_bytes * size, dtype=np.uint8)
dna_bytes_short = dna_bytes[:short_len]
dna_bytes_gapped = np.array(dna_template_bytes_gapped * size, dtype=np.uint8)
rna_bytes = np.array(rna_template_bytes * size, dtype=np.uint8)

dna_seq = DNA(dna_bytes)
dna_seq_short = DNA(dna_bytes_short)
dna_gapped = DNA(dna_bytes_gapped)
rna_seq = RNA(rna_bytes)

motif_1 = "GGTGCAAGCCGGTGGAAACA"
motif_1_regex = '(' + motif_1 + ')'


def consume_iterator(iterator):
    for _ in iterator:
        pass


class BenchmarkSuite:

    def time_object_creation(self):
        DNA(dna_bytes, validate=False)

    def time_object_creation_validate(self):
        DNA(dna_bytes)

    def time_reverse_complement(self):
        dna_seq.reverse_complement()

    def time_degap_all(self):
        dna_seq.degap()

    def time_translate(self):
        rna_seq.translate()

    def time_search_for_motif(self):
        consume_iterator(dna_seq.find_with_regex(motif_1_regex))

    def time_kmer_count_5(self):
        dna_seq_short.kmer_frequencies(5)

    def time_kmer_count_25(self):
        dna_seq_short.kmer_frequencies(25)

    def time_gc_content(self):
        dna_seq.gc_content()

    def time_search_for_motif_in_gapped(self):
        consume_iterator(
            dna_seq.find_with_regex(motif_1_regex, ignore=dna_seq.gaps()))
