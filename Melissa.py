from typing import *
from eliot_main import get_genetic_code, write_csv, read_csv
import random


def get_lengths(list_of_orf: list) -> [int]:
    """Give the list of orf lengths from a list of orf
        This function is written by Melissa Sadouki.
        Args:
            list_of_orf: list of ORF.
        Returns:
            list of lengths of ORF
        """

    for i in range (len(list_of_orf)):
        length=[]
        length_list.append(list_of_orf[i]["length"])
    return length_list


def get_longest_orf(list_of_orf: List[dict]) -> dict:
    """Give the longest orf from a list of orf
            This function is written by Melissa Sadouki.
            Args:
                list_of_orf: list of ORF.
            Returns:
                longest ORF
        """
    max_length=list_of_orf[0]["length"]
    for i in range (len(list_of_orf)) :
        if list_of_orf[i]["length"] > max_length :
            max_length = list_of_orf[i]["length"]
            max_orf = list_of_orf[i]
    return max_orf


def get_top_longest_orf(list_of_orf: List[dict], value: float) -> [dict]:
    """Return the value% top longest orfs from a list of orf
            This function is written by Melissa Sadouki.
            Args:
                list_of_orf: list of ORF.
                value: 0 > float > 1 : % of the top longest orf to show
            Returns:
                list of top ORF
        """
    
    key = lambda orf : orf ["length"]
    sorted_length = sorted(list_of_orf, key= key)

    print ("Entrez le pourcentage, entre 0 et 1, des séquences les plus longues, que vous souhaitez afficher.")
    value = float(input())
    if value > 0 and value < 1 :
        valuePourcent=value*100
        pourcentage = (len(sorted_length)*value)
    else :
        print ("Veuillez entrer le pourcentage avec des nombres compris entre 0 et 1.")
    
    print ("Les", valuePourcent,"\% de séquences ORF les plus longues sont les suivantes :")
    for i in range (pourcentage) :
        print (sorted_length[i])

###Main###

if __name__=='__main__':
    list_of_orf = []  # Liste d'ORF sur laquelle tu vas pouvoir essayer tes fonctions
    for _ in range(30):
        orf = {"start": random.randint(1, 1000)}
        orf["stop"] = orf["start"] + 3 * random.randint(100, 1000)
        orf["length"] = orf["stop"] - orf["start"]
        orf["protein"] = "M" * orf["length"]
        list_of_orf.append(orf)

    print(get_lengths(list_of_orf))
    print (get_longest_orf(list_of_orf))
    print (get_top_longest_orf(list_of_orf, 0.1))