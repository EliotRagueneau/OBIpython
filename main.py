class ORF:
    def __init__(self, source: str, start: int, stop: int, code_table: dict, protein: str, strand: str):
        self.source = source
        self.strand = strand
        if self.strand == "+":
            self.start = start
            self.stop = stop
        else:
            self.start = len(source) - start
            self.stop = len(source) - stop
        self.code_table = code_table
        self.protein = protein
        self.length = abs(self.stop - self.start)

    def __len__(self):
        return self.length

    def __repr__(self):
        return "{} : {}..{}".format(self.strand, self.start, self.stop)

    def translate(self):
        protein_seq = ""
        for i in range(self.start, self.stop, 3):
            protein_seq += self.code_table[self.source[i: i + 3]]


def find_orf(seq: str, threshold: int, code_table_id: int) -> list:
    """Give a list of all ORF in the sequence if they are grater than the threshold

            This function is written by Eliot Ragueneau.

            Args:
                seq: Sequence source
                threshold: Minimum size of the ORF in the list
                code_table_id: NCBI identifier of the translation table used on this sequence

            Returns:
                list of ORF
        """
    transl_table, start_table = get_genetic_code(code_table_id)

    strands = {"+": seq, "-": reversed_complementary(seq)}
    for strand in strands:
        dict_init = {}
        for i in range(len(strands[strand])):
            if strands[strand][i:i + 3] in start_table:
                dict_init[i] = None

        orf_list = []
        for init in dict_init:
            prot = ""
            for i in range(init + 3, len(strands[strand]), 3):
                codon = strands[strand][i - 3: i]
                if len(codon) == 3:
                    aa = transl_table[codon]
                    if aa == "*":
                        if i - init > threshold:
                            orf_list.append(ORF(sequence,
                                                start=init,
                                                stop=i,
                                                code_table=transl_table,
                                                protein=prot,
                                                strand=strand))
                        break
                    prot += aa
                else:
                    break

    return orf_list


def get_genetic_code(ncbi_id: int) -> tuple:
    """Give the initiation codons and the translation table by its id

            This function is written by Eliot Ragueneau.

            Args:
                ncbi_id: identifiant ncbi de la table de traduction à utiliser

            Returns:
                transl_table: dictionnaire de traduction
                start_table: liste des codons d'inititation
        """
    with open("Translation_tables/{}.txt".format(ncbi_id), "r") as file:
        lines = [line.strip() for line in file.readlines()]
        transl_table = {lines[2][i] + lines[3][i] + lines[4][i]: lines[0][i]
                        for i in range(len(lines[0]))}
        start_table = {lines[2][i] + lines[3][i] + lines[4][i]: lines[1][i]
                       for i in range(len(lines[0])) if lines[1][i] == "M"}
        return transl_table, start_table


def get_lengths(orf_list: list):
    """Give the list of orf lengths from a list of orf

        This function is written by .

        Args:
            orf_list: list of ORF.

        Returns:
            list of lengths of ORF
        """
    return [len(orf) for orf in orf_list]


def get_longest_orf(orf_list: list):
    """Give the longest orf from a list of orf

            This function is written by Eliot.

            Args:
                orf_list: list of ORF.

            Returns:
                longest ORF
        """
    return max(orf_list, key=lambda orf: len(orf))


def get_top_longest_orf(orf_list: list, value: float):
    """Give the value% top longest orfs from a list of orf

            This function is written by Eliot.

            Args:
                orf_list: list of ORF.
                value: 0 > float > 1 : % of the top longest orf to show

            Returns:
                list of top ORF
        """
    orf_list.sort(key=lambda orf: len(orf))
    return orf_list[int(- value * len(orf_list)) - 1:]


def reversed_complementary(sequence):
    complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([complement_dict[i] for i in sequence[::-1]])


sequence = ""
with open("influenza.fasta", 'r') as fasta:
    for line in fasta:
        if line[0] != ">":
            sequence += line.strip()

orf_list = find_orf(sequence, 300, 11)
orf_list.sort(key=lambda orf: orf.start)
print(orf_list)
