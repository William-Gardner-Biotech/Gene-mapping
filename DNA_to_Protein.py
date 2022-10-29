# ECS 129
# William Gardner
# Final Version
# 3/2/2022

# important packages for hashing

import os
import random
import hashlib

# Test case for running any length DNA sequences through

def random_DNA(length):
    return ''.join(random.choice('CGTA') for x in range(length))

# Simple prompt for the terminal interface that gets the filename from the same folder as this program

def acquire_file():
    choice = input("Please input the name of your .txt DNA file or PRESS ENTER TO TERMINATE PROGRAM:\n")
    if choice == "":
        return 0
    elif ".txt" in choice:
        try:
            open(choice, 'r')
            return choice
        except:
            print(
                "Filename invalid, check your input and ensure your .txt file is in the same folder as this program\n")
            return acquire_file()
    else:
        print("Invalid input\n")
        return acquire_file()

# Takes the string of DNA and makes a variable

def readfile(filename):
    infile = open(filename, 'r')
    DNA = infile.readline()
    return DNA


# Checks input file to make sure rest of function works and it is valid input

def clean_DNA(DNA):
    L = []
    DNA = DNA.upper()
    c_DNA = 0
    # prevents empty input files
    if len(DNA) == 0:
        return c_DNA
    desired = "ATGC"
    x = 0
    # iterates over DNA and can find errors
    for i in DNA:
        x += 1
        if i in desired:
            L.append(i)
        else:
            # Allows user to proofread DNA and resubmit
            print("The character", i, "@ position", x,
                  "is not a valid input, Please make sure the .txt file fits the parameters and try again.\n")
            return c_DNA
    seperator = ""
    c_DNA = seperator.join(L)
    return c_DNA

# Forms a complement DNA that is 3'to5'

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

# Simple, for clarity in main function

def RNAise(DNA):
    RNA = DNA.replace("T", "U")
    return RNA

# Iterate over whole DNA sequence to find all AUG regardless of reading frame

def AUG_finder(DNA):
    positions = []
    pos = 0
    # No need to check ends of DNA when only 3 remain because there will be no stop codon and no translation
    while pos + 3 < len(DNA):
        # 3 wide window for checking for the start codon
        if DNA[pos] == 'A' and DNA[pos + 1] == 'U' and DNA[pos + 2] == 'G':
            positions.append(pos)
            pos += 1
        else:
            pos += 1
            continue
    return positions

###### Not ACTIVE Recusion does not work
# pos reset to 0 after stop codon along with genotype
# to keep the index all we need is a simple equation, pos - (len genotype)*3 ****(pos + 3 - (len(genotype)*3)))
# The key value function of a dictionary prevents the genes without stop codons from being added

def translate_dict(DNA, positions_list, AA_table, pos=0, genotype="", translated_dict={}):
    if pos == 0:
        if positions_list == []:
            return translated_dict
        else:
            pos = positions_list.pop(0)
    chopped = DNA[pos:]
    for i in AA_table:
        if i[0] == chopped[:3]:
            if i[1] == "X":
                translated_dict["3' -> 5'  " + str(pos - len(genotype) * 3)] = genotype
                # print("DICT", translated_dict)
                pos = 0
                genotype = ""
                return translate_dict(DNA, positions_list, AA_table, pos, genotype, translated_dict)
            else:
                genotype = genotype + (i[1])
                pos += 3
                return translate_dict(DNA, positions_list, AA_table, pos, genotype, translated_dict)
        else:
            continue
    return translated_dict


# Function that builds the complete protein sequence over multiple iterations until positions list is empty

def gene_list(mRNA, positions, AA_table, translated_dict={}):
    for x in range(len(positions)):
        # pop removes the position to be sequenced next
        pos = positions.pop()
        # Boolean removes proteins with no stop codon
        if translator(mRNA, AA_table, pos) == None:
            continue
        else:
            # forms dictionary value for the completed sequence (postion, sequence) Also position is relative to template
            translated_dict["5'->3' " + str(pos)] = translator(mRNA, AA_table, pos)
    # print("NO WAY", translated_dict)
    return translated_dict


# Second Identical function to keep complement mRNA seperate so we don't lose mix up 5'to3'

def gene_list_template(mRNA, positions, AA_table, translated_dict={}):
    for x in range(len(positions)):
        # pop removes the position to be sequenced next
        pos = positions.pop()
        if translator(mRNA, AA_table, pos) == None:
            continue
        else:
            translated_dict["3'->5' " + str(pos)] = translator(mRNA, AA_table, pos)
    # print("NO WAY", translated_dict)
    return translated_dict


# Ultimately like the ribosome, builds protein sequences using concatenation

def translator(mRNA, AA_table, pos):
    genotype = ""
    while len(mRNA[pos:]) > 2:
        if match(mRNA[pos:pos + 3], AA_table) == 'X':
            # stop codon
            return genotype
        else:
            # continue to add AA to the sequence
            genotype = genotype + (match(mRNA[pos:pos + 3], AA_table))
        pos += 3


# Translator() helper function that pairs the codon with AA table and returns it to be concatenated

def match(codon, AA_table):
    for i in AA_table:
        if i[0] == codon:
            return i[1]
        else:
            continue


# Merge 2 dictionaries

def merge(dict1, dict2):
    for key in dict2:
        dict1[key] = dict2[key]
    return dict1


# A list of tuples that checks second position always as first position is position of start codon

def longest_gene(proteins, longest=[(0, "")]):
    for key in proteins:
        if len(proteins[key]) > len(longest[0][1]):
            longest = [(key, proteins[key])]
        # Appends equal length proteins to the list
        # Above Boolean recreates whole new list so multiple equal length genes don't linger
        elif len(proteins[key]) == len(longest[0][1]):
            longest.append((key, proteins[key]))
        else:
            continue
    return longest


# For clarity in rest of functions

def hash(obj):
    obj = hashlib.sha256(str(obj).encode('utf-8'))
    return obj.hexdigest()


# Create or read the old hash table .csv and create the hash table to use within the program

def hash_open_file():
    # We initiate a new dictionary every time we open the file, read the old and rebuild it
    hash_table = {}
    if os.path.exists("Protein_database.csv"):
        infile = open('Protein_database.csv', 'r')
        for line in infile:
            line = line.strip()
            (key, protein_code) = line.split(',')
            hash_table[key] = protein_code
        # print("File found and hash table created!")
        return hash_table
    else:
        # print("NO HASH TABLE FOUND, CREATED NEW TABLE")
        return hash_table


# Builds a database where keys are sha256 codes and only access is through original DNA strand

def enter_values(hash_table, DNA, protein_dict):
    secured = hash(DNA)
    protein_dict = str(protein_dict)
    protein_dict = protein_dict.replace(",", ".")
    hash_table[secured] = protein_dict
    return hash_table


## Lookup a hash table value using the DNA key

def access(hash_table, DNA):
    DNA = hash(DNA)
    locked_code = hash_table[(DNA)]
    return locked_code


## Check for duplicate DNA before translation and prevents collision

def search_secure(DNA, hash_table):
    if hash(DNA) in hash_table:
        return (access(hash_table, DNA), True)
    else:
        return (0, False)


## Write stored hash table on a .csv format

def hash_out_file(hash_table, hashed_DNA):
    outfile = open('Protein_database.csv', 'w')
    for i in hash_table.keys():
        if hash_table[i] == "{}":
            continue
        else:
            outfile.write(i + ',' + hash_table[i] + '\n')
    outfile.close()


def main():
    #DNA = random_DNA(150000)
    file_name = acquire_file()
    if file_name == 0:
        print("Application Terminated")
        return
    DNA = readfile(file_name)
    template_DNA = clean_DNA(DNA)
    if template_DNA == 0:
        print("Application Terminated")
        return
    # initiate hash table
    secure_DNA_dictionary = hash_open_file()
    # if input matches a hash key function exits and saves time
    if (search_secure(template_DNA, secure_DNA_dictionary)[1]) == True:
        return print("\nDNA already sequenced, Stored Sequence:", search_secure(template_DNA, secure_DNA_dictionary)[0])
    complement_RNA = RNAise(template_DNA)
    # template takes three steps
    template_RNA = complement(template_DNA)
    template_RNA = RNAise(template_RNA)
    template_RNA = template_RNA[::-1]
    # Constructed by hand
    AA_table = [('UUU', 'F'), ('UUC', 'F'), ('UUA', 'L'), ('UUG', 'L'), ('UCU', 'S'), ('UCC', 'S'), ('UCA', 'S'),
                ('UCG', 'S'), ('UAU', 'Y'), ('UAC', 'Y'), ('UAA', 'X'), ('UAG', 'X'), ('UGU', 'C'), ('UGC', 'C'),
                ('UGA', 'X'), ('UGG', 'W'), ('CUU', 'L'), ('CUC', 'L'), ('CUA', 'L'), ('CUG', 'L'), ('CCU', 'P'),
                ('CCC', 'P'), ('CCA', 'P'), ('CCG', 'P'), ('CAU', 'H'), ('CAC', 'H'), ('CAA', 'Q'), ('CAG', 'Q'),
                ('CGU', 'R'), ('CGC', 'R'), ('CGA', 'R'), ('CGG', 'R'), ('AUU', 'I'), ('AUC', 'I'), ('AUA', 'I'),
                ('AUG', 'M'), ('ACU', 'T'), ('ACC', 'T'), ('ACA', 'T'), ('ACG', 'T'), ('AAU', 'N'), ('AAC', 'N'),
                ('AAA', 'K'), ('AAG', 'K'), ('AGU', 'S'), ('AGC', 'S'), ('AGA', 'R'), ('AGG', 'R'), ('GUU', 'V'),
                ('GUC', 'V'), ('GUA', 'V'), ('GUG', 'V'), ('GCU', 'A'), ('GCC', 'A'), ('GCA', 'A'), ('GCG', 'A'),
                ('GAU', 'D'), ('GAC', 'D'), ('GAA', 'E'), ('GAG', 'E'), ('GGU', 'G'), ('GGC', 'G'), ('GGA', 'G'),
                ('GGG', 'G')]
    positions = AUG_finder(template_RNA)
    positions_backwards = AUG_finder(complement_RNA)
    protein_3to5 = gene_list_template(template_RNA, positions, AA_table)
    protein_5to3 = gene_list(complement_RNA, positions_backwards, AA_table)
    protein_dict = merge(protein_5to3, protein_3to5)
    print("\nComplete Gene Map:", protein_dict)
    print("\nNumber of Genes in the given DNA strand:", len(protein_dict))
    print("\nLongest Gene:", longest_gene(protein_dict))
    print("\nGene Density:", len(protein_dict) / (len(template_DNA)))
    enter_values(secure_DNA_dictionary, DNA, protein_dict)
    hash_out_file(secure_DNA_dictionary, hash(DNA))
    print(
        "\nDNA successfully hashed and stored for your privacy! To access the protein sequence at any time, input the same DNA as you entered during this run.")
    return

main()