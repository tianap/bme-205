#!/usr/bin/env python
import sys


def prefix(pattern):
    return pattern[:-1]


def suffix(pattern):
    return pattern[1:]


class deBruijn:
    """
    Contains the de Bruijn representation for a given set of kmers.
    """
    def __init__(self, patterns):
        """
        Initializes the deBruijn object with a parameter 'patterns,' a given list of k-mers.
        :param patterns: a list of k-mers

        Attributes:
            self.patterns  # list of kmers
            self.edges  # list of tuples representing the edges of a deBruijn graph
            self.deBruijnDict  # dictionary storing all (k-1)-mers

        """
        self.patterns = patterns
        self.edges = []
        self.deBruijnDict = dict()

    def findNodes(self):
        """
        Finds the edges (the prefix and suffix) for each k-mer in self.patterns, appends to self.edges.

        Called on by self.findConnections()
        :return:
        """
        for pattern in self.patterns:
            self.edges.append((prefix(pattern), suffix(pattern)))

    def findConnections(self):
        """
        Finds the connections across all edges in self.edges.

        Called upon by self.printDeBruijn()
        :return: updates self.deBruijnDict
        """
        self.findNodes()
        edgeSize = len(self.edges)
        for i in range(0, edgeSize):
            node = self.edges[i][0]
            if node not in self.deBruijnDict.keys():
                self.deBruijnDict[node] = []
                for edge in self.edges:
                    if edge[0] == node:
                        self.deBruijnDict[node].append(edge[1])

    def printDeBruijn(self):
        """
        Calls on self.findConnections() then prints the keys and values of self.deBruijnDict

        """
        self.findConnections()
        for key in sorted(self.deBruijnDict.keys()):
            if len(self.deBruijnDict[key]) != 0:
                print(key + ' -> ' + ','.join(str(i) for i in sorted(self.deBruijnDict[key])))


def main():
    """
    Reads from STDIN and uses these arguments for the class deBruijn.
        """
    patterns = []

    for line in sys.stdin:
        patterns.append(line.strip())

    myDeBruijn = deBruijn(patterns)
    myDeBruijn.printDeBruijn()


if __name__ == "__main__":
    main()