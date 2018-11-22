def getGeneticCode(NCBI_ID):
    with open("Translation_tables/{}.txt".format(NCBI_ID), "r") as file:
        lines = [line.strip() for line in file.readlines()]
        transl_table = {lines[2][i] + lines[3][i] + lines[4][i]: lines[0][i]
                        for i in range(len(lines[0]))}
        start_table = {lines[2][i] + lines[3][i] + lines[4][i]: lines[1][i]
                       for i in range(len(lines[0])) if lines[1][i] == "M"}
        print(transl_table)
        print(start_table)
        for codon_start in start_table:

            print(codon_start)


getGeneticCode(11)
