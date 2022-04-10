### ECS 129 ###

# filename = AA_seq.csv
def build_AA_table(filename):
    infile = open(filename, 'r')
    AA_table = []

    for line in infile:
        line = line.strip()
        (code, name) = line.split(",")
        AA_table.append((code,name))
    # AA table is a tuple list of AA's with their three letter codon followed by their 3 letter AA name
    return AA_table


def readfile(filename):
    infile = open(filename, 'r')
    DNA = infile.readline()
    return DNA

# find complmentary strand

def clean_DNA(DNA):
    L = []
    desired = "ATGC"
    for i in DNA:
        if i in desired:
            L.append(i)
        else:
            continue
    seperator = ""
    clean_DNA = seperator.join(L)
    return(clean_DNA)

def complement(DNA):
    L = []
    for i in DNA:
        if i == 'G':
            L.append('C')
        elif i == 'C':
            L.append("G")
        elif i == 'A':
            L.append("T")
        elif i == 'T':
            L.append("A")
        else:
            continue

    separator = ""

    clean_complement_DNA = separator.join(L)
    return clean_complement_DNA
'''else:
    L = []
    five = input("Input 3' -> 5'")
    L.append("5'-")
    for i in five:
        if i == 'G':
            L.append('C')
        elif i == 'C':
            L.append("G")
        elif i == 'A':
            L.append("T")
        elif i == 'T':
            L.append("A")
        else:
            continue
            # print("error")
    L.append("-3'")
    separator = " "
    print(separator.join(L))'''

def RNAise(DNA):
    RNA = DNA.replace("T","U")
    return RNA
# recursively go through and chunk it

def AUG_finder(DNA):
    positions = []
    pos = 0
    while pos+3 <= len(DNA):
        pos += 1
        if DNA[pos] == 'A' and DNA[pos+1] == 'U' and DNA[pos+2] == 'G':
            #print("Check", DNA[pos:], pos)
            positions.append(pos)
        else:
            continue
    return positions
# pos reset to 0 after stop codon along with genotype
# to keep the index all we need is a simple equation, pos - (len genotype)*3 ****(pos + 3 - (len(genotype)*3)))
'''def translate(DNA, positions_list, AA_table, pos = 0, genotype = "", translated_list = []):
    if pos == 0:
        pos = positions_list.pop(0)
    chopped = DNA[pos:]
    for i in AA_table:
        if i[0] == chopped[:3]:
            if i[1] == "X":
                genotype = genotype + (i[1])
                #print("LENGTH", (pos + 3 - (len(genotype)*3)))

                #translated_dict[pos - len(genotype)*3] =
                #translated_list.append(genotype)
                pos = 0
                genotype = ""
                #print("Working here")
                #print("Finished",finished)
                return translate(DNA, positions_list, AA_table, pos, genotype, translated_list)
            else:
                genotype = genotype+(i[1])
                pos += 3
                return translate(DNA, positions_list, AA_table, pos, genotype, translated_list)
        else:
            continue
    return translated_list
    #print(genotype)
    #print(translated_list)'''
# pos reset to 0 after stop codon along with genotype
# to keep the index all we need is a simple equation, pos - (len genotype)*3 ****(pos + 3 - (len(genotype)*3)))
# The key value function of a dictionary prevents the genes without stop codons from being added
def translate_dict(DNA, positions_list, AA_table, pos = 0, genotype = "", translated_dict = {}):
    if pos == 0:
        if positions_list == []:
            return translated_dict
        else:
            pos = positions_list.pop(0)
    chopped = DNA[pos:]
    for i in AA_table:
        if i[0] == chopped[:3]:
            if i[1] == "X":
                translated_dict["3' -> 5'  "+str(pos - len(genotype)*3)] = genotype
                #print("DICT", translated_dict)
                pos = 0
                genotype = ""
                return translate_dict(DNA, positions_list, AA_table, pos, genotype, translated_dict)
            else:
                genotype = genotype+(i[1])
                pos += 3
                return translate_dict(DNA, positions_list, AA_table, pos, genotype, translated_dict)
        else:
            continue
    return translated_dict

def translate_dict_template(DNA, positions_list, AA_table, pos = 0, genotype = "", translated_dict = {}):
    if pos == 0:
        if positions_list == []:
            return translated_dict
        else:
            pos = positions_list.pop(0)
    chopped = DNA[pos:]
    for i in AA_table:
        if i[0] == chopped[:3]:
            if i[1] == "X":
                translated_dict["5' -> 3'  "+str(pos - len(genotype)*3)] = genotype
                #print("DICT", translated_dict)
                pos = 0
                genotype = ""
                return translate_dict_template(DNA, positions_list, AA_table, pos, genotype, translated_dict)
            else:
                genotype = genotype+(i[1])
                pos += 3
                return translate_dict_template(DNA, positions_list, AA_table, pos, genotype, translated_dict)
        else:
            continue
    return translated_dict

# Awaiting Patrice's repsonse
def input_check(DNA):
    desired = "GCTA"
    if len(DNA) == 0:
        print("Empty")
        return False
    for i in DNA:
        if i not in desired:
            print("Failure")
            return False
        else:
            continue
    print("Exited")
    return DNA



# print(AA_table[0])

def main_2():
    #DNA = "TAGGCGAATGAATCTATTCTTACTGTATCGAAGAATGGCCTCGCGGAGGCATGTGTCATGCTAGCGTGCGGGGTACTCTAGTTATCCATATGGTCCACAGGACACTCGTTGCTTTCGGATTTGCCCTTTATGCGCCGGTTTTCAGCCACGCTTATGCTCAGCATCGTTATAACCAGACCGATACTAGATCTATAAAGTCC"
    #input_check(DNA)
    DNA = readfile("DNA.txt")
    #print("DNA", DNA)
    template_DNA = clean_DNA(DNA)
    complement_RNA = RNAise(template_DNA)
    template_RNA = complement(template_DNA)
    template_RNA = RNAise(template_RNA)
    template_RNA = template_RNA[::-1]
    AA_table = [('UUU', 'F'), ('UUC', 'F'), ('UUA', 'L'), ('UUG', 'L'), ('UCU', 'S'), ('UCC', 'S'), ('UCA', 'S'), ('UCG', 'S'), ('UAU', 'Y'), ('UAC', 'Y'), ('UAA', 'X'), ('UAG', 'X'), ('UGU', 'C'), ('UGC', 'C'), ('UGA', 'X'), ('UGG', 'W'), ('CUU', 'L'), ('CUC', 'L'), ('CUA', 'L'), ('CUG', 'L'), ('CCU', 'P'), ('CCC', 'P'), ('CCA', 'P'), ('CCG', 'P'), ('CAU', 'H'), ('CAC', 'H'), ('CAA', 'Q'), ('CAG', 'Q'), ('CGU', 'R'), ('CGC', 'R'), ('CGA', 'R'), ('CGG', 'R'), ('AUU', 'I'), ('AUC', 'I'), ('AUA', 'I'), ('AUG', 'M'), ('ACU', 'T'), ('ACC', 'T'), ('ACA', 'T'), ('ACG', 'T'), ('AAU', 'N'), ('AAC', 'N'), ('AAA', 'K'), ('AAG', 'K'), ('AGU', 'S'), ('AGC', 'S'), ('AGA', 'R'), ('AGG', 'R'), ('GUU', 'V'), ('GUC', 'V'), ('GUA', 'V'), ('GUG', 'V'), ('GCU', 'A'), ('GCC', 'A'), ('GCA', 'A'), ('GCG', 'A'), ('GAU', 'D'), ('GAC', 'D'), ('GAA', 'E'), ('GAG', 'E'), ('GGU', 'G'), ('GGC', 'G'), ('GGA', 'G'), ('GGG', 'G')]
    #AA_table = build_AA_table('AA_seq.csv')
    #print(AA_table)
    #print(template_RNA)
    #print(complement_RNA)
    # position on the DNA is same as complement's position, meaning template is backwards
    positions = AUG_finder(template_RNA)
    positions_backwards = AUG_finder(complement_RNA)
    #print("Positions read from template DNA", positions)
    #print("Positions read from complement DNA", positions_backwards)
    protein_3to5 = translate_dict(complement_RNA, positions, AA_table)
    protein_5to3 = translate_dict_template(template_RNA, positions_backwards, AA_table)
    print("3to5",protein_3to5)
    #print("Protein_5to3", protein_5to3)
    #print("Protein_3to5", protein_3to5)
    return


#main()
'''DDNA = 'ATG'
print("original DNA: 5'-", DDNA)
print("Complementry: 3'-", complement(DDNA))'''

'''complement_RNA = RNAise(DDNA)

template_RNA = complement(DDNA)
template_RNA = RNAise(template_RNA)
template_RNA = template_RNA[::-1]'''

#print("Two strings going from 5'-3' for translation")
#print(complement_RNA)
#print(template_RNA)
#print("Backwards", DDNA[::-1])
'''DNA = readfile("DNA.txt")
print("complement:", complement(DNA))
print("original  :", clean_DNA(DNA))'''

main_2()
