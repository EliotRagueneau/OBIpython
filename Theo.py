from typing import *
from eliot_main import get_genetic_code, write_csv, read_csv


def read_flat_file(filename: str) -> str:
    """Load a file in memory by returning a string

            This function is written by Théo Gauvrit.

            Args:
                filename: file to open

            Returns:
                string of the whole file (with \n)
        """
    # Voir with statements sur Sam et Max
    pass


def reversed_complement(dna_seq: str) -> str:
    """Return the reversed complement of the given DNA sequence

            This function is written by .

            Args:
                dna_seq: DNA sequence to be reversed.

            Returns:
                reversed complement DNA sequence
        """
    pass


def find_orf(seq: str, threshold: int, code_table_id: int) -> [dict]:
    """Give a list of all ORF in the sequence if they are grater than the threshold

            This function is written by .

            Args:
                seq: Sequence to analyse
                threshold: Minimum size of the ORF in the list
                code_table_id: NCBI identifier of the translation table used on this sequence

            Returns:
                list of ORF
                    start: start position (in bp)
                    stop: stop position (in bp)
                    length: ORF length (in bp)
                    protein: translated protein sequence if available.
        """
    pass
