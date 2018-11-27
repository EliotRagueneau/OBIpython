from typing import *
from eliot_main import get_genetic_code, write_csv, read_csv
import random

list_of_orf = []  # Liste d'ORF sur la quelle tu vas pouvoir essayer tes fonctions
for _ in range(30):
    orf = {"start": random.randint(1, 1000)}
    orf["stop"] = orf["start"] + 3 * random.randint(100, 1000)
    orf["length"] = orf["stop"] - orf["start"]
    orf["protein"] = "M" * orf["length"]
    list_of_orf.append(orf)


def get_lengths(orf_list: list) -> [int]:
    """Give the list of orf lengths from a list of orf

        This function is written by Mélissa Sadouki.

        Args:
            orf_list: list of ORF.

        Returns:
            list of lengths of ORF
        """
    pass


def get_longest_orf(orf_list: List[dict]) -> dict:
    """Give the longest orf from a list of orf

            This function is written by Mélissa Sadouki.

            Args:
                orf_list: list of ORF.

            Returns:
                longest ORF
        """
    pass


def get_top_longest_orf(orf_list: List[dict], value: float) -> [dict]:
    """Return the value% top longest orfs from a list of orf

            This function is written by Mélissa Sadouki.

            Args:
                orf_list: list of ORF.
                value: 0 > float > 1 : % of the top longest orf to show

            Returns:
                list of top ORF
        """
    pass
