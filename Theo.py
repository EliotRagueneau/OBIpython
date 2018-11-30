from typing import *
from eliot_main import get_genetic_code, write_csv, read_csv


def read_flat_file(filename):
    """Load a file in memory by returning a string

            This function is written by Théo Gauvrit.

            Args:
                filename: file to open

            Returns:
                string of the whole file (with \n)
        """
       
    fichier=open(filename,"r")
    txt = fichier.read()
    fichier.close()
    return txt # ← 'LOCUS       NM_000518 ' ...


def reversed_complement(dna_seq):
    """Return the reversed complement of the given DNA sequence

            This function is written by Théo GAUVRIT.

            Args:
                dna_seq: DNA sequence to be reversed.

            Returns:
                reversed complement DNA sequence
        """
    complement_seq=[]
    dna_seq=dna_seq.upper()
    complement_dict={"A":"T","C":"G","G":"C","T":"A"}
    
    for nucleotide in dna_seq[::-1]:
        complement_seq.append(complement_dict[nucleotide])
        
    return "".join(complement_seq)


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
