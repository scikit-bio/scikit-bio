from io import StringIO
import sys
from skbio import io, TabularMSA, DNA
import skbio
from skbio.io.format.fasta import _fasta_sniffer

# fs = '\n'.join([r"@seq1 description 1", r"AACACCAAACTTCTCCACC", r"ACGTGAGCTACAAAAG", r"+seq1 description 1", r"''''Y^T]']C^CABCACC", r"'^LB^CCYT\T\Y\WF", r"@seq2 description 2", r"TATGTATATATAACATATACATATATACATACATA", r"+", r"]KZ[PY]_[YY^'''AC^\\'BT''C'\AT''BBB"])
# fh = StringIO(fs)
# fl = [">seq1 Turkey\n", "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT\n", ">seq2 Salmo gair\n", "AAGCCTTGGCAGTGCAGGGTGAGCCGTGG\n", "CCGGGCACGGTAT\n", ">seq3 H. Sapiens\n", "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA\n", ">seq4 Chimp\n", "AAACCCTTGCCG\n", "TTACGCTTAAAC\n", "CGAGGCCGGGAC\n", "ACTCAT\n", ">seq5 Gorilla\n", "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA\n"]
# msa = TabularMSA.read(fh, constructor=DNA, variant='sanger')
# print(msa)
#msa = TabularMSA.read(fl, constructor=DNA)
#print(msa)

for r in skbio.read(sys.stdin, format="fastq"):
    print(r)