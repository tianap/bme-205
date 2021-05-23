#!/usr/bin/env python
########################################################################
# File: problem15.py
#  executable: problem15.py < input.txt > output.out
# Purpose: finds all DNA sequences that code for a given peptide
#   stderr: errors and status
#   stdout:
#
# Author: Tiana Pereira
# History:      10/27/2020 Created
#
########################################################################

import sys


########################################################################
# Global functions
########################################################################

def reverseComplement(seq):
    """
    Finds the reverse complement for a given DNA sequence
    :param seq: a DNA sequence
    :return: the reverse complement of seq
    """
    basePairs = {
        'A': 'T',
        'G': 'C',
        'T': 'A',
        'C': 'G'
    }
    rc = ''
    for i in seq:
        rc += basePairs[i]
    return rc[::-1]


class Peptide:
    """
    Class containing a DNA sequence, a peptide sequence, and finds the DNA subsequences coding for the peptide sequence.

    3 Attributes:
        self.dnaSeq
        self.peptide
        self.peptideSeqs

    4 methods:
        self.__init__()
        self.translation()
        self.findPeptideSeqs()
        self.printPeptideSeqs()

    NOTE: There is one other variable, dnaCodonTable, that is not included in the constructor since it is the same across all Peptide
    objects. This table stores all of the codons and the amino acids they code for, as well as the stops.
    """

    def __init__(self, dnaSeq, peptide):
        """
        Constructor for the Peptide class. Takes two parameters:
        :param dnaSeq: a DNA sequence
        :param peptide: a peptide sequence

        Contains the following attributes:
            self.dnaSeq  # (str) the DNA sequence provided by the user
            self.peptide  # (str) the peptide sequence provided by the user
            self.peptideSeqs  # (list) stores the DNA subsequences that code for the self.peptide
        """
        self.dnaSeq = dnaSeq
        self.peptide = peptide
        self.peptideSeqs = []

    '''
    The dictionary below is taken from: https://pythonforbiologists.com/dictionaries
    '''
    dnaCodonTable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

    def translation(self, seq):
        """
        Finds the subsequences of a given seq and determines whether they code for the peptide. Also considers the
        reverse complements of seq.
        :param seq: a DNA sequence
        :return: updates self.peptideSeqs

        """

        thisPeptide = ''
        if len(seq) % 3 == 0 and len(self.peptide) * 3 == len(seq):  # if the sequence is a valid length
            # Find the amino acids for the forward sequence
            for i in range(0, len(seq), 3):
                thisPeptide += self.dnaCodonTable[seq[i:i + 3]]
            # If this amino acid sequence is the same as the peptide
            if thisPeptide == self.peptide:
                self.peptideSeqs.append(seq)

            thisPeptide = ''
            rc = reverseComplement(seq)
            for i in range(0, len(seq), 3):
                thisPeptide += self.dnaCodonTable[rc[i:i + 3]]
            if thisPeptide == self.peptide:
                self.peptideSeqs.append(seq)

    def findPeptideSeqs(self):
        """
        Slices the DNA sequences into subsequences that potentially code for the peptide. Calls on self.translation() to
        translate each subsesquence into a peptide.
                """

        for i in range(0, len(self.dnaSeq)):
            seq = self.dnaSeq[i:(i + 3 * len(self.peptide))]
            self.translation(seq)

    def printPeptideSeqs(self):
        """
        Prints the DNA sequences coding for the peptide by calling on self.findPeptideSeqs().
        """
        self.findPeptideSeqs()
        for seq in self.peptideSeqs:
            print(seq)


def main():
    """
    Reads from STDIN and creates an object Peptide with the given peptide and DNA sequence.
    :return:
    """
    peptide = ''
    sequence = ''
    lineCount = 0

    for line in sys.stdin:
        if lineCount == 0:
            sequence = line.strip()
            lineCount += 1
        else:
            peptide = line.strip()

    myPeptide = Peptide(sequence, peptide)
    myPeptide.printPeptideSeqs()


if __name__ == "__main__":
    main()
