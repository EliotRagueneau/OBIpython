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
    conc=txt.split("\n")
    for i in range(len(conc)):
        if conc[i]=="FEATURES             Location/Qualifiers":
            features = "".join(conc[i+1:])
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
    list_CDS=[]
    list_genes=[]
    bloc=("").join(features.split("          "))
    bloc_CDS=bloc.split("     ")
    for ligne in range(len(bloc_CDS)):
        if bloc_CDS[ligne][0:3]=="CDS":
            list_CDS.append(ligne)
    for ligne in range(len(list_CDS)):
        CDS={"start":None,"stop":None,"frame":None,"length":None,"name":"unknown","protein":"xxx","product":"unknown"}
        bloc_cleaned=bloc_CDS[list_CDS[ligne]].split("/")
        if bloc_cleaned[0][6:17]=="complement(":
            start=int((bloc_cleaned[0][17:].strip()).split("..")[0])
            inter_stop=(bloc_cleaned[0][17:].strip()).split("..")[1]
            stop=int(inter_stop[:len(inter_stop)-1])
            CDS["start"]=start
            CDS["stop"]=stop
            CDS["frame"]=1-(start%3)
            CDS["length"]=stop-start
        else:
            start=int((bloc_cleaned[0].strip()).split("..")[0][6:])
            stop=int((bloc_cleaned[0].strip()).split("..")[1])
            CDS["start"]=start
            CDS["stop"]=stop
            CDS["frame"]=1+(start%3)
            CDS["length"]=stop-start
        for ligne2 in range(len(bloc_cleaned)):
            if bloc_cleaned[ligne2][:5]=="gene=":
                CDS["name"]=(bloc_cleaned[ligne2][5:].split('"'))[1]
            elif bloc_cleaned[ligne2][:8]=="product=":
                CDS["product"]=(bloc_cleaned[ligne2][8:].split('"'))[1]
            elif bloc_cleaned[ligne2][:12]=="translation=":
                sequence=(bloc_cleaned[ligne2][12:].split('"')[1]).split(" ")
                protein=""
                for aa in range(len(sequence)):
                    protein=protein+sequence[aa]
                    CDS["protein"]=protein
        list_genes.append(CDS)
    return list_genes
    pass


def read_gen_bank(filename: str) -> Dict[str, Union[str, List[dict]]]:
    Features={"Description":None,"ID":None,"length":None,"organism":None,"type":None,"genes":None,"sequence":None}
    Genes=get_genes(get_features(read_flat_file(filename)))
    String=read_flat_file(filename)
    ligne=String.split("\n")
    for i in range(len(ligne)):
        mot=ligne[i].split(" ")
        for j in range(len(mot)):
            if mot[j]=="DEFINITION":
                Features["Description"]=" ".join(mot[j+2:])    
            elif mot[j]=="VERSION":
                Features["ID"]=mot[j+5]
            elif mot[j]=="bp":
                Features["length"]=mot[j-1]+" bp"
            elif mot[j]=="ORGANISM":
                Features["organism"]=" ".join(mot[j+2:])
            elif mot[j] in ["DNA","RNA","Protein","rRNA","mRNA","cRNA","tRNA"]:
                Features["type"]=mot[j]
            elif mot[j]=="ORIGIN":
                Features["sequence"]=("\n").join(ligne[i+1:-3])

    Features["genes"]=Genes
    return Features
    
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
                 a=get_genes(get_features(read_flat_file("sequence.gb")))
b=read_gen_bank(read_flat_file("sequence.gb"))
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
