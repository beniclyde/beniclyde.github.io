"""This module is for various manipulations of sequence data.

Functions include transcription, translation and other common uses in genetics.

Copyright Benjamin Clyde, 2017
"""

def nucleotide_count(sequence):
    """
    Counts nucleotides in a sequence.
    Arguments: sequence string
    Returns: list of tuples of the nucleotide and the count
    """
    if not set(sequence) == set("AGTC"): #error check for invalid characters
        raise TypeError("Invalid nucleotide in sequence")
    nuc_count = []
    for i in set(sequence):
        nuc_count.append((i, sequence.count(i)))
    return nuc_count

def transcribe(dna):
    """
    Transcribes a DNA sequence by replacing thymine (T) with uracil (U)
    """
    return dna.replace('T', 'U')


def complement(sequence):
    """
    Finds the complement of the sequence by list comprehension using dictionary of complement values
    """
    basecomplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'U': 'A', 'N': 'N'}
    bases = list(sequence)
    comp = [basecomplement[base] for base in bases]
    return ''.join(comp)

def reverse(sequence):
    """
    Reverses a sequence, used to get the 5' orientation of the complement strand of a sequence
    """
    return sequence[::-1]

def gc_content(sequence):
    """
    Calculates the GC content of a sequence, which is important to know in molecular biology
    """
    gc_percent = (sequence.count('G') + sequence.count('C')) / len(sequence)
    return gc_percent * 100


def parse_fasta(fasta):
    """
    Parses a fasta file into a dict with the key as the ID tag, and the value as the sequence
    """
    sequences = {}
    for line in fasta:
        if line.startswith(">"):
            try:
                sequences[key] = sequence
            except NameError:
                pass
            key = line[1:].strip()
            sequence = ""

        else:
            sequence += line.strip()
    return sequences


def codons(sequence, rf=0):
    """
    Finds the list of codons in a sequence in a given reading frame.

    In a sequence there are 6 possible reading frames, 3 on the forward strand, 3 on the complement.
    Arguments: sequence for splitting in to codons
    rf: reading frame. 0 = 1st forward frame; 1 = 2nd forward frame; 2 = 3rd forward frame
        3 = 1st complement frame; 4 = 2nd complement frame; 5 = 3rd complement frame
    Returns: list of codons in sequence for specified reading frame
    """
    if rf == 0:
        stoppoint = len(sequence) - (len(sequence) % 3)
        return [sequence[i:i+3] for i in range(rf, stoppoint, 3)]
    if rf == 1:
        stoppoint = len(sequence[rf:]) - (len(sequence[rf:]) % 3)
        return [sequence[i:i+3] for i in range(rf, stoppoint, 3)]
    if rf == 2:
        stoppoint = len(sequence[2:]) - (len(sequence[2:]) % 3)
        return [sequence[i:i+3] for i in range(rf, stoppoint, 3)]
    if rf == 3 or rf == 4 or rf == 5:
        new_s = transcribe(complement(sequence))
        stoppoint = len(new_s[rf-3:]) - (len(new_s[rf-3:]) % 3)
        return [new_s[i:i+3] for i in range(rf-3, stoppoint, 3)]
    else:
        raise ValueError("Reading frame rf must be a value from 0-5")

rna_codon_table = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
                   "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
                   "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
                   "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
                   "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
                   "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                   "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
                   "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                   "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
                   "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                   "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
                   "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
                   "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
                   "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                   "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
                   "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}


def translate(mrna):
    """
    Translates an mRNA sequence into a protein sequence using the mRNA codon look up dictionary
    """
    protein = ""
    for i in range(0, len(mrna) - (len(mrna) % 3), 3):
        symbol = rna_codon_table[mrna[i:i+3]]
        if symbol == "STOP":
            break
        else:
            protein += symbol
    return protein


def six_frames(mrna):
    """
    Translate sequence into protein in all 6 reading frames,
    prints the mrna with complement and the translated 6 frames aligned
    with the nucleotides
    """
    for rf in reversed(range(3)):
        protein = " "*rf
        for i in codons(mrna, rf=rf):
            if not rna_codon_table[i] == "STOP":
                protein += ((" ") + rna_codon_table[i] + " ")
            else:
                protein += ' * '
        print(protein)
    print(mrna)
    print(transcribe(complement(mrna)))
    for rf in range(3, 6):
        protein = " "*(rf-3)
        for i in codons(mrna, rf=rf):
            if not rna_codon_table[i] == "STOP":
                protein += ((" ") + rna_codon_table[i] + " ")
            else:
                protein += ' * '
        print(protein)

def find_start_codons(sequence, codon='AUG'):
    """
    Returns a list of indexes for all start codons in a sequence
    """
    return [i for i in range(len(sequence)) if sequence.startswith(codon, i)]

def orf_protein(fasta):
    """
    Function to take the input of of a fasta file, parse it, and then return the set
    of all possible proteins from open reading frames. To be considered, the sequence
    has to start with a start codon and then goes until it finds the first stop codon
    in the reading frame. Sequences that do not hit a stop codon are not included
    """
    parsed = parse_fasta(fasta)
    parsed_list = list(parsed.values())
    results = []
    for v in parsed_list:
        comp = complement(v)
        reverse_comp = reverse(comp)
        rna1 = transcribe(v)
        rna2 = transcribe(reverse_comp)

        indices_forward = find_start_codons(rna1)
        indices_reverse = find_start_codons(rna2)
        length = len(rna1)

        for i in indices_forward:
            found_stop = False
            protein = ""
            for j in range(i, length, 3):
                if len(rna1[j:j+3]) == 3:
                    codon = rna_codon_table[rna1[j:j+3]]
                    if codon == "STOP":
                        found_stop = True
                        break
                    if not codon:
                        break
                    protein += codon
            if found_stop:
                results.append(protein)
        #return results
        for i in indices_reverse:
            found_stop = False
            protein = ""
            for j in range(i, length, 3):
                if len(rna1[j:j+3]) == 3:
                    codon = rna_codon_table[rna2[j:j+3]]
                    if codon == "STOP":
                        found_stop = True
                        break
                    if not codon:
                        break
                    protein += codon
            if found_stop:
                results.append(protein)
    return set(results)

dna_codon_table = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
                   "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
                   "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
                   "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
                   "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
                   "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                   "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
                   "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                   "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
                   "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                   "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
                   "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
                   "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
                   "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                   "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
                   "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

class DNA(str):
    """
    Class for representing DNA sequences as strings
    """
    basecomplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'U': 'A'}
    dna_codon_table = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
                       "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
                       "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
                       "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
                       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
                       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
                       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
                       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
                       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
                       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
                       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
                       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

    def __init__(self, seq):
        """Create DNA object to string seq."""
        if set(seq.upper()).issubset({'A', 'T', 'C', 'G', 'N'}):
            self.seq = seq.upper()
        else:
            raise ValueError("DNA has illegal character")

    def transcribe(self):
        """
        Transcribes a DNA sequence by replacing thymine (T) with uracil (U)
        """
        return self.seq.replace('T', 'U')

    def complement(self):
        """
        Finds the complement of the sequence by list comprehension using dict of complement values
        """
        basecomplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'U': 'A', 'N': 'N'}
        bases = list(self)
        comp = [basecomplement[base] for base in bases]
        return ''.join(comp)

    def reverse(self):
        """
        Reverses a sequence, used to get the 5' orientation of the complement strand of a sequence
        """
        return self[::-1]

    def gc_content(self):
        """
        Calculates the GC content of a sequence, which is important to know in molecular biology
        """
        gc_percent = (self.count('G') + self.count('C')) / len(self)
        return gc_percent * 100

    def codons(self, rf=0):
        """
        Finds the list of codons in a sequence in a given reading frame.

        In a sequence there are 6 possible reading frames, 3 on forward strand, 3 on complement.
        Arguments: sequence for splitting in to codons
        rf: reading frame. 0 = 1st forward frame; 1 = 2nd forward frame; 2 = 3rd forward frame
            3 = 1st complement frame; 4 = 2nd complement frame; 5 = 3rd complement frame
        Returns: list of codons in sequence for specified reading frame
        """
        if rf == 0:
            stoppoint = len(self) - (len(self) % 3)
            return [self[i:i+3] for i in range(rf, stoppoint, 3)]
        if rf == 1:
            stoppoint = len(self[rf:]) - (len(self[rf:]) % 3)
            return [self[i:i+3] for i in range(rf, stoppoint, 3)]
        if rf == 2:
            stoppoint = len(self[2:]) - (len(self[2:]) % 3)
            return [self[i:i+3] for i in range(rf, stoppoint, 3)]
        if rf == 3 or rf == 4 or rf == 5:
            new_s = complement(self)
            stoppoint = len(new_s[rf-3:]) - (len(new_s[rf-3:]) % 3)
            return [new_s[i:i+3] for i in range(rf-3, stoppoint, 3)]
        else:
            raise ValueError("Reading frame rf must be a value from 0-5")

    def translate(self):
        """
        Translates an mRNA sequence into a protein sequence using the mRNA codon look up dictionary
        """
        protein = ""
        for i in range(0, len(self) - (len(self) % 3), 3):
            symbol = dna_codon_table[self[i:i+3]]
            if symbol == "STOP":
                break
            else:
                protein += symbol
        return protein

    def six_frames(self):
        """
        Translate sequence into protein in all 6 reading frames,
        prints the mrna with complement and the translated 6 frames aligned
        with the nucleotides
        """
        for i in reversed(range(3)):
            protein = " "*i
            for j in self.codons(self, rf=i):
                if not dna_codon_table[j] == "STOP":
                    protein += ((" ") + dna_codon_table[j] + " ")
                else:
                    protein += ' * '
            print(protein)
        print(self)
        print(complement(self))
        for i in range(3, 6):
            protein = " "*(i-3)
            for j in self.codons(self, rf=i):
                if not dna_codon_table[j] == "STOP":
                    protein += ((" ") + dna_codon_table[j] + " ")
                else:
                    protein += ' * '
            print(protein)

    def __str__(self):
        return "DNA = {} \n RNA = {} \n Protein = {}".format(self.seq, self.transcribe(), self.translate())

