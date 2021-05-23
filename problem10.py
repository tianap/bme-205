#!/usr/bin/env python
import sys

"""
Global Methods:
    suffix(pattern)
    prefix(pattern)
"""


def suffix(pattern):
    """
    Returns the string containing the last k-1 characters of pattern, aka the suffix.
    :param pattern: a string
    :return: the suffix of pattern
    """
    return pattern[1:len(pattern)]


def prefix(pattern):
    """
    Returns the string containing the first k-1 characters of pattern, aka the prefix.
    :param pattern: a string
    :return: the prefix of pattern
    """
    return pattern[:len(pattern) - 1]


class OverlapGraph:
    """
    Contains the adjacency list from a given list of kmers
    """

    def __init__(self, kmers):
        """
        Intializes the OverlapGraph object. Takes in one parameter, a list of kmers.

        :param kmers:
        """
        self.kmers = kmers
        self.adjacencyList = []

    def assembleGraph(self):
        """
        Generates the adjacency list self.adjacencyLIst from the list of k-mers given by the user.

        Called upon by self.printAdjList()
        :return: updates self.adjacencyList
        """
        for kmer1 in self.kmers:
            for kmer2 in self.kmers:
                if suffix(kmer1) == prefix(kmer2):
                    self.adjacencyList.append((kmer1, kmer2))

    def printAdjList(self):
        """
        Calls on method self.assembleGraph() and prints the resulting self.adjacencyList

        :return: Prints the adjacency list
        """
        self.assembleGraph()
        self.adjacencyList.sort(key=lambda x: x[0])
        for i in self.adjacencyList:
            print(i[0] + " -> " + i[1])


def main():
    """
    Reads from STDIN and creates an OverlapGraph object using the patterns passed as the argument.

    """
    pattern_list = []
    for line in sys.stdin:
        pattern_list.append(line.strip())

    myOverlapGraph = OverlapGraph(pattern_list)
    myOverlapGraph.printAdjList()


if __name__ == "__main__":
    main()
