from typing import *
from eliot_main import get_genetic_code, write_csv, read_csv, read_fasta


def read_flat_file(filename):
    """Load a file in memory by returning a string

            This function is written by Theo Gauvrit.

            Args:
                filename: file to open

            Returns:
                string of the whole file (with \n)
        """
       
    fichier=open(filename,"r")
    txt = fichier.read()
    fichier.close()
    return txt #  'LOCUS       NM_000518 ' ...


def reversed_complement(dna_seq):
    """Return the reversed complement of the given DNA sequence

            This function is written by Theo GAUVRIT.

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


def find_orf(seq: str, threshold:int,code_table_id: int):
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
                    frame:1,2,3 ou -1,-2,-3
        """
    transl_table, start_table = get_genetic_code(code_table_id)

    ORF_lidi=[]
    #ORF_lidi[e]={"start":0,"stop":0,"length":0,"protein":0,"frame":0}
    start_lidi=[]
    listestop=[]
    stop_lidi=[]
    orftrans=[]
    rever=reversed_complement(seq)
 
    
    for z in transl_table.keys():
        if transl_table[z]=='*':
            listestop.append(z)
            
    for a in range(3):
        for pos in range(a,len(seq)):
            if seq[pos:pos+3] in start_table:
                start_lidi.append({"pos":pos,"frame":a+1})
            if rever[pos:pos+3] in start_table:
                start_lidi.append({"pos":len(seq)-pos,"frame":-a-1})
                
            if seq[pos:pos+3] in listestop:
                stop_lidi.append({"pos":pos,"frame":a+1})
            if rever[pos:pos+3] in listestop:
                stop_lidi.append({"pos":len(seq)-pos,"frame":-a-1})


    for i in range(len(stop_lidi)):
        for j in range(i,len(stop_lidi)):
           if stop_lidi[j]["frame"]==stop_lidi[i]["frame"] and stop_lidi[j]["pos"]>stop_lidi[i]["pos"]:
               
               for k in range(len(start_lidi)):
                  
                   
                  if stop_lidi[j]["frame"]>0:
                      
                      print(2135465464)
                      break
                  if start_lidi[k]["frame"]==stop_lidi[i]["frame"] and ((stop_lidi[j]["pos"])-(start_lidi[k]["pos"]))>threshold:
                        print(3165464646464)
                         
                        for l in range(start_lidi[k]["pos"],stop_lidi[j]["pos"]):
                              print(123123)
                              print("ouiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii")
     #                           orftrans.append(transl_table[seq[l:l+3]])
     #                       ORF_lidi.append({"start":k["pos"],"stop":j["pos"],"length":len(orftrans)*3,"protein":orftrans,"frame":k["frame"]})
                           
                   # else:
                   #     if k["frame"]==i["frame"] and (k["pos"]-j["pos"])<(-threshold):
                   #         for l in range(k["pos"]-3,j["pos"]):
                   #             orftrans.append(transl_table[seq[l:l+3]])
                   #         ORF_lidi.append({"start":k["pos"],"stop":j["pos"],"length":len(orftrans)*3,"protein":orftrans,"frame":k["frame"]})
                           
                   
                  
               
               
               
               

    
    
                

                
    return listestop,start_lidi,stop_lidi,ORF_lidi


if __name__=="__main__":
    
   seqentry=read_fasta("influenza.fasta")
   a,b,c,d=find_orf(seqentry,300,11)
  
   
