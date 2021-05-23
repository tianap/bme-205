#!/usr/bin/env python
import sys


class GenomePath:
    def __init__(self, patterns):
        """
        Initializes the class GenomePath with a given parameter patterns, a list of k-mers.

        :param patterns: a list of kmers
        """
        self.genomePathString = ''
        self.patterns = patterns

    def findString(self):
        """
        Loops through the k-mers and assembles the genome string.

        Called upon by self.printString()
        :return:
        """
        # Set the start of the string to be the first k-mer
        self.genomePathString += self.patterns[0]

        # For every other k-mer after the first one, add the last character
        for pattern in self.patterns[1:]:
            self.genomePathString += pattern[:-1]

    def printString(self):
        """
        Calls on self.findString() then prints the updated self.genomePathString
        """
        self.findString()
        print(self.genomePathString)

def main():
    """
    Reads from STDIN and creates an object myGenomePath using the arguments from STDIN.
    """
    patterns = []
    for line in sys.stdin:
        patterns.append(line.strip())

    myGenomePath = GenomePath(patterns)
    myGenomePath.printString()

