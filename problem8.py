#!/usr/bin/env python
import sys


class Composition:
    """
    Represents the k-mer composition of a given string.
    """

    def __init__(self, k, text):
        """

        """
        self.k = k
        self.text = text
        self.kmerComp = []

    def findComposition(self):
        """
        Loops through self.text and generates all substrings of length k.

        :return: updates self.kmerComp
        """

        for i in range(0, len(self.text) - self.k + 1):
            self.kmerComp.append(self.text[i:i + self.k])

    def printComposition(self):
        """
        Calls on self.findComposition and prints the found kmers within self.kmerComp.

        :return:  prints the kmer composition
        """
        self.findComposition()
        for kmer in self.kmerComp:
            print(kmer)


def main():
    """
    Reads in STDIN and calls on class Composition() using the integer and string from STDIN as arguments.

    Calls on method printComposition() to print the kmer compositions.
    """

    lineCount = 0
    k = 0
    text = ''
    for line in sys.stdin:
        if lineCount == 0:  # the first line
            k = int(line.strip())
        else:
            text = line.strip()
        lineCount += 1

    myComposition = Composition(k, text)
    myComposition.printComposition()


if __name__ == "__main__":
    main()
