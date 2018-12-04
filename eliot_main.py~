from typing import *


class ORF:
    """Open Reading Frame defined by start and stop nucleotide"""

    def __init__(self, source: str, start: int, stop: int, code_table_id: int, protein: str, strand: str):
        ##########################
        # À corriger / Améliorer #
        ##########################
        self.source = source
        self.strand = strand
        if self.strand == "+":
            self.start = start
            self.stop = stop
        else:
            self.start = len(source) - start
            self.stop = len(source) - stop
        self.code_table = code_table_id
        self.protein = protein
        self.length = abs(self.stop - self.start)
        self.product = None  # A rajouter
        self.name = None  # A rajouter

    def __len__(self):
        return self.length

    def __repr__(self):
        return "{} : {}..{}".format(self.strand, self.start, self.stop)

    def as_dict(self):
        return {"strand": self.strand,
                "start": self.start,
                "stop": self.stop,
                "transl_table": self.code_table,
                "protein": self.protein
                }


def find_orf(seq: str, threshold: int, code_table_id: int) -> [ORF]:
    """Give a list of all ORF in the sequence if they are grater than the threshold

            This function is written by Eliot Ragueneau.

            Args:
                seq: Sequence source
                threshold: Minimum size of the ORF in the list
                code_table_id: NCBI identifier of the translation table used on this sequence

            Returns:
                list of ORF
        """
    ##########################
    # À corriger / Améliorer #
    ##########################
    transl_table, start_table = get_genetic_code(code_table_id)

    strands = {"+": seq, "-": reversed_complement(seq)}
    orf_list = []
    for strand in strands:
        dict_init = {}
        for i in range(len(strands[strand])):
            if strands[strand][i:i + 3] in start_table:
                dict_init[i] = None
        list_stop = []  # Seulement pour n'avoir que les ORFs maximum
        for init in dict_init:
            prot = ""
            for i in range(init + 3, len(strands[strand]), 3):
                codon = strands[strand][i - 3: i]
                if len(codon) == 3:
                    aa = transl_table[codon]
                    if aa == "*":
                        if i - init > threshold and i not in list_stop:  # Seulement pour n'avoir que les ORFs maximum
                            list_stop.append(i)
                            orf_list.append(ORF(seq,
                                                start=init,
                                                stop=i,
                                                code_table_id=code_table_id,
                                                protein=prot,
                                                strand=strand))
                        break
                    prot += aa
                else:
                    break

    return orf_list


def get_genetic_code(ncbi_id: int) -> Tuple[dict, dict]:
    """Give the initiation codons and the translation table by its id

            This function is written by Eliot Ragueneau.

            Args:
                ncbi_id: ncbi translation table identifier

            Returns:
                transl_table: translating dictionary
                start_table: list of start codons
        """
    with open("Translation_tables/{}.txt".format(ncbi_id), "r") as file:
        lines = [line.strip() for line in file.readlines()]
        transl_table = {lines[2][i] + lines[3][i] + lines[4][i]: lines[0][i]
                        for i in range(len(lines[0]))}
        start_table = {lines[2][i] + lines[3][i] + lines[4][i]: lines[1][i]
                       for i in range(len(lines[0])) if lines[1][i] == "M"}
        return transl_table, start_table


def get_lengths(orf_list: list) -> [int]:
    """Give the list of orf lengths from a list of orf

        This function is written by Eliot Ragueneau.

        Args:
            orf_list: list of ORF.

        Returns:
            list of lengths of ORF
        """
    return [len(orf) for orf in orf_list]


def get_longest_orf(orf_list: List[ORF]) -> ORF:
    """Give the longest orf from a list of orf

            This function is written by Eliot Ragueneau.

            Args:
                orf_list: list of ORF.

            Returns:
                longest ORF
        """
    return max(orf_list, key=lambda orf: len(orf))


def get_top_longest_orf(orf_list: List[ORF], value: float) -> [ORF]:
    """Return the value% top longest orfs from a list of orf

            This function is written by Eliot Ragueneau.

            Args:
                orf_list: list of ORF.
                value: 0 > float > 1 : % of the top longest orf to show

            Returns:
                list of top ORF
        """
    orf_list.sort(key=lambda orf: len(orf))
    return orf_list[int(- value * len(orf_list)) - 1:]


def reversed_complement(dna_seq: str):
    """Return the reversed complement of the given DNA sequence

            This function is written by Eliot Ragueneau.

            Args:
                dna_seq: DNA sequence to be reversed.

            Returns:
                reversed complement DNA sequence
        """
    complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([complement_dict[i] for i in dna_seq[::-1]])


def read_csv(filename: str, separator: str = ";") -> [dict]:
    with open(filename, 'r') as file:
        list_lines = [line.strip() for line in file.readlines()]
        keys = list_lines[0].split(separator)
        return [{key: element for key in keys for element in line.split(separator)} for line in list_lines[1:]]


def write_csv(filename: str, data: list, separator: str = ";"):
    if isinstance(data[0], ORF):
        data = [orf.as_dict() for orf in data]
    filename += ".csv" if filename[-4:] != ".csv" else ""
    with open(filename, "w") as file:
        file.write("{}\n".format(separator.join(data[0])))  # Write only the keys
        for dictionary in data:
            file.write("{}\n".format(separator.join([str(dictionary[key]) for key in dictionary])))


class GenBank:
    """GenBank data structure (More interesting than dictionaries)"""

    def __init__(self, filename: str):
        """Parse a GenBank file

                    This function is written by .

                    Args:
                        filename: .gb file to parse

                    Returns:
                        list of ORF (either as ORF or as dict)
                """
        text = self.read_flat_file(filename)
        features = self.get_features(text)
        self.genes = self.get_genes(features)

        self.definition = None
        self.type = None  # dna, rna, or protein
        self.data = None  # only if available otherwise set to ‘xxx’. if seq too large, the entry does not contain data.
        self.id = None
        self.length = None
        self.gbtype = None
        self.organism = None
        self.code_table_id = None

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def get_genes(features: str) -> [ORF]:
        """Extract gene and CDS data from their section inside features table

                This function is written by .

                Args:
                    features: text with gene to extract

                Returns:
                    list of ORF (either as ORF or as dict)
            """
        # Voir package builtin re : peut être utile ici
        pass

    @staticmethod
    def read(filename):
        """Parse a GenBank file

                This function is written by .

                Args:
                    filename: .gb file to parse

                Returns:
                    list of ORF (either as ORF or as dict)
            """
        return GenBank(filename)


def read_fasta(filename: str) -> str:
    dna = ""
    with open(filename, 'r') as fasta:
        for fasta_line in fasta:
            if fasta_line[0] != ">":
                dna += fasta_line.strip()


if __name__ == '__main__':
    dna = read_fasta("influenza.fasta")
    list_of_orf = find_orf(dna[:50001], 300, 11)
    list_of_orf.sort(key=lambda orf: orf.start, reverse=True)
    print(get_lengths(list_of_orf))
    print(get_longest_orf(list_of_orf))
    print(get_top_longest_orf(list_of_orf, 0.2))
