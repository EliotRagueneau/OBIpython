from typing import *
from eliot_main import get_genetic_code, write_csv, read_csv
from Theo import read_flat_file


def get_features(txt: str) -> str:
    """Extract features lines from flat text and return them

            This function is written by Kévin Merchadou.

            Args:
                txt: flat text with features to extract

            Returns:
                string of features
        """
    # Voir package builtin re : peut être utile ici
    print("Théo est une grosse merde !")


def get_genes(features: str) -> [dict]:
    """Extract gene and CDS data from their section inside features table

            This function is written by Kévin Merchadou.

            Args:
                features: text with gene to extract

            Returns:
                list of ORF (either as ORF or as dict)
        """
    # Voir package builtin re : peut être utile ici
    pass


def read_gen_bank(filename: str) -> Dict[str, Union[str, List[dict]]]:
    """Parse a GenBank file

            This function is written by Kévin Merchadou.

            Args:
                filename: .gb file to parse

            Returns:
                dictionary of :

                    features:
                            description: entry title (genbank descriptor DEFINITION)
                            type: myBio keywords only - dna, rna, or protein
                            data: sequence data only if available otherwise set to ‘xxx’. When the sequence is too large
                                  the entry does not contain data.
                            ID: Identifier (locus)
                            length: sequence length
                            gbtype: molecule type as described in a genbank entry.
                            organism: organism²
                            codeTableID: NCBI genetic code table identifier

                    list of gene = ORF (dict)
                            start: start position (in bp)
                            stop: stop position (in bp)
                            frame: frame (1, 2, 3,- 1, -2, or -3)
                            length: gene length (in bp)
                            name: gene name if available. By default, ‘unknown’.
                            protein: translated protein sequence if available. By default, ‘xxx’.
                            product: product name if available. By default, ‘unknown’.
        """
    pass
