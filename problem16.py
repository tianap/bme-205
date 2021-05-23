#!/usr/bin/env python
########################################################################
# File: problem16.py
#  executable: problem16.py < input.txt > output.out
# Purpose: Generates a theoretical spectrum given an amino acid string Peptide

#
# Author: Tiana Pereira
# History:      10/30/2020 Created
#
########################################################################

import sys


class Cyclospectrum:
    """
    Class containing the cyclospectrum generated from a amino acid string, peptide.

    4 Attributes:
        self.peptide
        self.tSpectrum
        self.subPeptides
        self.n

    4 Methods:
        self.__init__()
        self.findSubpeptides()
        self.findMasses()
        self.printCyclospectrum()

    NOTE: There is also one extra variable, self.aaMass, that is not contained in the constructor class since it is the same
    across all Cyclopeptide objects.

    """
    aaMass = {
        'G': 57,
        'A': 71,
        'S': 87,
        'P': 97,
        'V': 99,
        'T': 101,
        'C': 103,
        'I': 113,
        'L': 113,
        'N': 114,
        'D': 115,
        'K': 128,
        'Q': 128,
        'E': 129,
        'M': 131,
        'H': 137,
        'F': 147,
        'R': 156,
        'Y': 163,
        'W': 186
    }

    def __init__(self, peptide):
        """
        Constructs the Cyclospectrum class. Takes one argument, peptide, a sequence of amino acids.

        :param peptide: peptide sequence provided by user

        Attributes:
            self.peptide  # (str) stores the peptide sequence provided by the user
            self.tSpectrum  # (list) stores the spectrum values, starting with the integer mass 0
            self.subPeptides  # (list) stores the subPeptides, including the full peptide sequence
            self.n  # (int) the length of the peptide
        """
        self.peptide = peptide
        self.tSpectrum = [0]
        self.subPeptides = [self.peptide]
        self.n = len(peptide)

    def findSubpeptides(self):
        """
        Finds all subpeptides of a cyclic peptide sequence.

        :return: updates self.subPeptides
        """
        for i in range(0, self.n):
            for j in range(1, self.n):
                if i + j > self.n:
                    self.subPeptides.append(self.peptide[i:] + self.peptide[:i + j - self.n])
                else:
                    self.subPeptides.append(self.peptide[i:i + j])

    def findMasses(self):
        """
        Calls on self.findSubpeptides and finds the masses for each subpeptide, then updates self.tSpectrum.

        :return: updates self.tSpectrum
        """
        self.findSubpeptides()
        for seq in self.subPeptides:
            mass = 0
            for aa in seq:
                mass += self.aaMass[aa]
            self.tSpectrum.append(mass)

    def printCyclospectrum(self):
        """
        Prints the theoretical spectrum for self.peptide.
        """
        self.findMasses()
        print(' '.join(str(i) for i in sorted(self.tSpectrum)))


def main():
    """
    Reads from STDIN and creates a Cyclospectrum object using the read peptide sequence. Then calls the function
    printCyclospectrum() to print the theoretical spectrum of the peptide.
    :return:
    """
    seq = ''
    for line in sys.stdin:
        seq = line.strip()

    myPeptide = Cyclospectrum(seq)
    myPeptide.printCyclospectrum()


if __name__ == "__main__":
    main()
