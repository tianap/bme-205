#!/usr/bin/env python
import sys


class deBruijn:
    """
    Object that contains the de Bruijn graph of a string. Takes two parameters, a string 'sequence' and an integer 'k'

    Attributes:
        self.sequence
        self.k
        self.deBruijnDict

    Methods:
        self.findNodes()
        self.printDeBruijn()
    """
    def __init__(self, sequence, k):
        """
        Initializes the deBruijn object. Takes two parameters, sequence and k.
        :param sequence: a string DNA sequence
        :param k: an integer k
        """
        self.sequence = sequence  # a DNA string, passed by th user.
        self.deBruijnDict = dict()  # dictionary storing all (k-1) sized nodes. Each key is a unique (k-1)mer of
                                    # self.sequence, and each value is a list of nodes whose prefix is the same as the
                                    # suffix of the key. In this way, each edge is the shared sequence of a prefix and
                                    # suffix for two nodes.
        self.k = k  # an integer k, passed by the user.

    def findNodes(self):
        """
        Finds the (k-1)-mers, aka nodes, of self.sequence.

        For a given (k-1)mer in self.sequence, the next (k-1)mer shifted right by one character is appended to the list
        of connecting nodes stored as a value for the given (k-1)mer.

        ex.
        self.deBruijnDict = {
            kmer1 = [kmer2],
            kmer2 = [kmer3, kmer4]
        }
        :return: updates self.deBruijnDict
        """
        for i in range(0, len(self.sequence) - self.k + 1):
            node = self.sequence[i:i + self.k - 1]
            nextNode = self.sequence[i + 1:i + self.k]
            if node not in self.deBruijnDict.keys():
                self.deBruijnDict[node] = []
            self.deBruijnDict[node].append(nextNode)

    def printDeBruijn(self):
        """
        Calls on self.findNodes() to find the nodes for self.deBruijnDict, then prints out the sorted nodes

        :return: prints the resulting adjacency list
        """
        self.findNodes()
        for key in sorted(self.deBruijnDict.keys()):
            if len(self.deBruijnDict[key]) != 0:
                print(key + ' -> ' + ','.join(str(i) for i in sorted(self.deBruijnDict[key])))


def main():
    """
    Creates an object deBruijn using the parameters from STDIN.
    :return:
    """
    k = 0
    sequence = ''
    lineCount = 0

    for line in sys.stdin:
        if lineCount == 0:  # the first line
            k = int(line.strip())
        else:
            sequence = line.strip()
        lineCount += 1

    myDeBruijn = deBruijn(sequence, k)
    myDeBruijn.printDeBruijn()


if __name__ == "__main__":
    main()
