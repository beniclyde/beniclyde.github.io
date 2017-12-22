"""This module is for various mutations of sequence data.

Functions include random RNA generator, mutation programs that modify the RNA strings
and mutations finder programs that attempt to identify the mutation.

Copyright Benjamin Clyde, 2017
"""
import random
import gzip
import pickle
from ipywidgets import interact
from IPython.display import HTML, display
from DNA import rna_codon_table
from DNA import translate

rna_codon_list_no_stop = [k for k in rna_codon_table if rna_codon_table[k] != 'STOP']

def random_rna(length):
    """
    Function to take a given length and return a random RNA sequence
    Function ensures the RNA sequence starts with a start codon and ends with a stop codon
    """
    rna_alphabet = list('AUGC')
    rna = [random.choice(rna_codon_list_no_stop) for i in range(int((length-6)/3))]
    rna = ''.join(rna)
    l = len(rna)
    for i in range(3, l, 3):
        if rna[i:i+3] == "UAG" or rna[i:i+3] == "UAA" or rna[i:i+3] == "UGA":
            rna = rna.replace(rna[i], random.choice(rna_alphabet))
    rna = "AUG" + rna + "UGA"
    return rna

def snp(sequence):
    """
    Function takes a sequence string and introduces a single nucleotide polymorphism (SNP)
    Randomly finds an index in the string and replaces the value at the index with a new
    letter from the given alphabet. If the value is the same as the previous value, the
    function calls itself again to make sure the sequence changes.
    """
    rna_alphabet = list('AUGC')
    l = len(sequence)
    i = random.randint(0, l)
    mutation = random.choice(rna_alphabet)
    while sequence[i] == mutation:
        mutation = random.choice(rna_alphabet)

    new_seq = sequence[:i] + mutation + sequence[i+1:]
    if new_seq[i] == sequence[i]:
        snp(new_seq)
    return new_seq

def find_snp(rna, mutated_rna):
    """
    Function compares the old sequence to the mutated sequence to find the SNP
    Returns the position, the previous nucleotide and the new nucleotide
    """
    for i in range(0, len(rna)):
        if not rna[i] == mutated_rna[i]:
            result = ("At position {}, {} mutated to {}").format(i, rna[i], mutated_rna[i])
    if translate(rna) == translate(mutated_rna):
        return result + ", silent mutation"
    else:
        return result

def find_codon_diff(protein, mutated_protein):
    """
    Function to find the effect the SNP had on the protein.
    First compares the length of the proteins to see if the SNP caused a premature
    stop codon, otherwise, returns the new codon.
    """
    if len(protein) == len(mutated_protein):
        for i in range(0, len(protein)):
            if not protein[i] == mutated_protein[i]:
                return "Amino acid {} at position {} changed to {}".format(protein[i], i, mutated_protein[i])
    else:
        return (protein[len(mutated_protein)+1], len(mutated_protein) + 1, "STOP")

def insertion(rna):
    """
    Function to take an RNA string and randomly choose an index for the string and then
    randomly add 1-9 nucleotides sequentially and return the new RNA string
    """
    rna_alphabet = list('AUGC')
    l = len(rna)
    index = random.randint(0, l)
    nuc = random.randint(1, 10)
    insert = [random.choice(rna_alphabet) for i in range(nuc)]
    insert = ''.join(insert)
    new_seq = rna[:index] + insert + rna[index:]
    return new_seq

def find_insertion(rna, inserted_rna):
    """
    Function to find the inserted sequence into an RNA seq
    First takes the difference in length between the two sequences to find the amount
    of nucleotides that were inserted, then compare the sequences until it can determine
    where the sequences differed.
    """
    len_diff = len(inserted_rna) - len(rna)
    for i in range(0, len(rna)):
        if rna[i] == inserted_rna[i]:
            continue
        else:
            try:
                if rna[i] == inserted_rna[i+len_diff] and rna[i:len(rna)] == \
                inserted_rna[i+len_diff:len(inserted_rna)]:
                    return "At position {}, nucleotide(s) {} were inserted".\
                            format(i, inserted_rna[i:i+len_diff])
                else:
                    raise ValueError("Too many nucleotides added to RNA")
            except:
                raise ValueError("Something went wrong")

def deletion(rna):
    """
    Function to take an RNA string and randomly choose an index for the string and then
    randomly delete 1-9 nucleotides sequentially and return the new RNA string
    """
    l = len(rna)
    index = random.randint(0, l)
    nuc = random.randint(1, 10)
    new_seq = rna[:index] + rna[index+nuc:]
    return new_seq

def find_deletion(rna, deleted_rna):
    """
    Function to find the deleted sequence from an RNA seq
    First takes the difference in length between the two sequences to find the amount
    of nucleotides that were inserted, then compare the sequences until it can determine
    where the sequences differed.
    """
    len_diff = len(rna) - len(deleted_rna)
    for i in range(0, len(rna)):
        if rna[i] == deleted_rna[i]:
            continue
        else:
            try:
                if deleted_rna[i] == rna[i+len_diff] and rna[i+len_diff:len(rna)] == \
                deleted_rna[i:len(deleted_rna)]:
                    return "At position {}, nucleotide(s) {} were deleted".\
                            format(i, rna[i:i+len_diff])
                else:
                    raise ValueError("Too many nucleotides deleted RNA")
            except:
                raise ValueError("Something went wrong")

def mutator(rna):
    """
    Function to take an RNA string and randomly apply one of the mutation functions.
    Probably should adjust the random to closely match the frequency the mutations
    occur in nature, but couldn't find good numbers.
    """
    muta = random.randint(0, 2)
    if muta == 0:
        return snp(rna)
    if muta == 1:
        return insertion(rna)
    if muta == 2:
        return deletion(rna)

def find_mutation(rna, mutated_rna):
    """
    Function to take an RNA string and the function mutated string to find the mutation
    regardless of type. Limited to finding one mutation at a time, will look into
    multiple mutations in a later version
    """
    if len(rna) == len(mutated_rna):
        return find_snp(rna, mutated_rna)
    if len(rna) > len(mutated_rna):
        return find_deletion(rna, mutated_rna)
    if len(rna) < len(mutated_rna):
        return find_insertion(rna, mutated_rna)

def data_gen(number):
    """
    Given an input number, generates an RNA string, mutates the string, and
    then pickles the data for later use.
    """
    results = []
    for i in range(0, number):
        length = random.randint(60, 82)
        start_rna = random_rna(length)
        start_prot = translate(start_rna)
        mut_rna = mutator(start_rna)
        mut_prot = translate(mut_rna)
        result = find_mutation(start_rna, mut_rna)
        record = (start_rna, start_prot, mut_rna, mut_prot, result)
        results.append(record)
        with gzip.open("rna_mutation_data.pickle.gz", 'wb') as fo:
            pickle.dump(record, fo)
    return results
