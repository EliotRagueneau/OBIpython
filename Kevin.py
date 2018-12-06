from typing import *

from Theo import read_flat_file


def get_features(txt: str) -> str:
    """Extract features lines from flat text and return them
            This function is written by Kévin Merchadou.
            Args:
                txt: flat text with features to extract
            Returns:
                string of features
        """
    lines = txt.split("\n")
    features = ""
    for i in range(len(lines)):
        if lines[i] == "FEATURES             Location/Qualifiers":
            features = "".join(lines[i + 1:])
            break
    return features


def get_genes(features: str) -> [dict]:
    """Extract gene and CDS data from their section inside features table

            This function is written by Kévin Merchadou.

            Args:
                features: text with gene to extract

            Returns:
                list of ORF (either as ORF or as dict)
        """
    list_cds = []
    list_genes = []
    bloc = "".join(features.split("          "))
    bloc_cds = bloc.split("     ")
    for ligne in range(len(bloc_cds)):
        if bloc_cds[ligne][0:3] == "CDS":
            list_cds.append(ligne)
    for ligne in range(len(list_cds)):
        cds = {"start": None, "stop": None, "frame": None, "length": None, "name": "unknown", "protein": "xxx",
               "product": "unknown"}
        bloc_cleaned = bloc_cds[list_cds[ligne]].split("/")
        if bloc_cleaned[0][6:17] == "complement(":
            start = int((bloc_cleaned[0][17:].strip()).split("..")[0])
            inter_stop = (bloc_cleaned[0][17:].strip()).split("..")[1]
            stop = int(inter_stop[:len(inter_stop) - 1])
            cds["start"] = start
            cds["stop"] = stop
            cds["frame"] = 1 - (start % 3)
            cds["length"] = stop - start
        else:
            start = int((bloc_cleaned[0].strip()).split("..")[0][6:])
            stop = int((bloc_cleaned[0].strip()).split("..")[1])
            cds["start"] = start
            cds["stop"] = stop
            cds["frame"] = 1 + (start % 3)
            cds["length"] = stop - start
        for ligne2 in range(len(bloc_cleaned)):
            if bloc_cleaned[ligne2][:5] == "gene=":
                cds["name"] = (bloc_cleaned[ligne2][5:].split('"'))[1]
            elif bloc_cleaned[ligne2][:8] == "product=":
                cds["product"] = (bloc_cleaned[ligne2][8:].split('"'))[1]
            elif bloc_cleaned[ligne2][:12] == "translation=":
                sequence = (bloc_cleaned[ligne2][12:].split('"')[1]).split(" ")
                protein = ""
                for aa in range(len(sequence)):
                    protein = protein + sequence[aa]
                    cds["protein"] = protein
        list_genes.append(cds)
    return list_genes
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
                                data: sequence data only if available otherwise set to ‘xxx’. When the sequence is too
                                      large the entry does not contain data.
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
    features = {"Description": None, "ID": None, "length": None, "organism": None, "type": None, "genes": None,
                "sequence": None}
    genes = get_genes(get_features(read_flat_file(filename)))
    string = read_flat_file(filename)
    lignes = string.split("\n")
    for i in range(len(lignes)):
        mots = lignes[i].split(" ")
        for j in range(len(mots)):
            if mots[j] == "DEFINITION":
                features["Description"] = " ".join(mots[j + 2:])
            elif mots[j] == "VERSION":
                features["ID"] = mots[j + 5]
            elif mots[j] == "bp":
                features["length"] = mots[j - 1] + " bp"
            elif mots[j] == "ORGANISM":
                features["organism"] = " ".join(mots[j + 2:])
            elif mots[j] in ["DNA", "RNA", "Protein"]:
                features["type"] = mots[j]
            elif mots[j] == "ORIGIN":
                features["sequence"] = "\n".join(lignes[i + 1:-3][10:])

    features["genes"] = genes
    return lignes
