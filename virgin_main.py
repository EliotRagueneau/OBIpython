from typing import *  # Permet de préciser le type du contenu des listes et dictionnaires


def find_orf(seq: str, threshold: int, code_table_id: int) -> [dict]:
    """Give a list of all ORF in the sequence if they are grater than the threshold

            This function is written by .

            Args:
                seq: Sequence to analyse
                threshold: Minimum size of the ORF in the list
                code_table_id: NCBI identifier of the translation table used on this sequence

            Returns:
                list of ORF
        """
    pass


def get_genetic_code(ncbi_id: int) -> Tuple[dict, dict]:
    """Give the initiation codons and the translation table by its id

            This function is written by .

            Args:
                ncbi_id: ncbi translation table identifier

            Returns:
                transl_table: translating dictionary
                start_table: list of start codons
        """
    pass


def get_lengths(orf_list: list) -> [int]:
    """Give the list of orf lengths from a list of orf

        This function is written by .

        Args:
            orf_list: list of ORF.

        Returns:
            list of lengths of ORF
        """
    pass


def get_longest_orf(orf_list: List[dict]) -> dict:
    """Give the longest orf from a list of orf

            This function is written by .

            Args:
                orf_list: list of ORF.

            Returns:
                longest ORF
        """
    pass


def get_top_longest_orf(orf_list: List[dict], value: float) -> [dict]:
    """Return the value% top longest orfs from a list of orf

            This function is written by .

            Args:
                orf_list: list of ORF.
                value: 0 > float > 1 : % of the top longest orf to show

            Returns:
                list of top ORF
        """
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


def read_csv(filename: str, separator: str = ";") -> [dict]:
    """Read a csv file delimited by separator

            This function is written by .

            Args:
                filename: .csv file to read
                separator: separator between data in the csv file

            Returns:
                list of dictionary of features:
                    [
                     {colonne 1: element 1, colonne 2: element 2}
                     {colonne 1: element 3, colonne 2: element 4}
                     {colonne 1: element 5, colonne 2: element 6}
                     {colonne 1: element 7, colonne 2: element 8}
                    ]
        """
    pass


def write_csv(filename: str, data: List[dict], separator: str = ";"):
    """Write a csv file delimited by separator from a list of dictionary

                This function is written by .

                Args:
                    filename: .csv file to read
                    data: list of dict where each dict is a line in the product file,
                                            and their keys the first row of the file
                    separator: separator between data in the csv file

                Returns:
                    None
            """
    pass


def read_flat_file(filename: str) -> str:
    """Load a file in memory by returning a string

            This function is written by .

            Args:
                filename: file to open

            Returns:
                string of the whole file (with \n)
        """
    # Voir with statements sur Sam et Max
    pass


def get_features(txt: str) -> str:
    """Extract features lines from flat text and return them

            This function is written by .

            Args:
                txt: flat text with features to extract

            Returns:
                string of features
        """
    # Voir package builtin re : peut être utile ici
    pass


def get_genes(features: str) -> [dict]:
    """Extract gene and CDS data from their section inside features table

            This function is written by .

            Args:
                features: text with gene to extract

            Returns:
                list of ORF (either as ORF or as dict)
        """
    # Voir package builtin re : peut être utile ici
    pass


def read_gen_bank(filename: str) -> Dict[str: str, str: List[dict]]:
    """Parse a GenBank file

            This function is written by .

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
                            organism: organism
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


if __name__ == '__main__':
    pass
