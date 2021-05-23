#!/usr/bin/env python
import sys
import random

sys.setrecursionlimit(10 ** 6)
"""
Imported modules:
    sys  # used for reading in files and setting the recursion limit to be higher
    random  # used to randomly select items from lists
"""


class Graph:
    """
    The class Graph takes a set of DNA sequences of size k and constructs the nodes and edges needed to generate the
    sequence composition.

    It contains five separate functions in addition to def __init__(). The only function called in main() is the
    method printString(), which calls on the other methods (more information provided in the method docstrings).

    Attributes:
        self.patterns
        self.edges
        self.nodes
        self.stringPath
        self.sequence

    Methods:
        findEdges()
        findNodes()
        getStart()
        findPath(node v)
        printString()

    NOTE: most of these methods are used within problem13.py, I just decided not to import them.
    """

    def __init__(self, patterns):
        """
        The initialization method for the class Graph. Takes a list of DNA sequences of size k, 'patterns', which will
        then be used to create the attributes self.edges and self.nodes. Assumes the strings contained within the list
        of patterns are all the same size.

        Attributes:
            self.patterns  # list set to the list of kmers given by a user
            self.edges  # list of tuples representing edges for the graph, generated by self.findEdges()
                ex. [(a,b),(c,d),(d,a),(b,c)]
            self.nodes  # dictionary in which each key is a pattern and each value is a tuple containing two lists: the
                        # first containing the incoming nodes, the second containing the outgoing nodes
                ex.
                    self.nodes = {
                        'a': ([d],[b])
                        'b': ([a],[c])
                        'c': ([b],[d])
                        'd': ([c],[a])
                    }
            self.stringPath  # list of patterns in the Eulerian path, generated by self.findPath()
            self.sequence  # the outputted sequence composition of all kmers within patterns

        :param patterns: the list of kmers given by a user
        """
        self.patterns = patterns
        self.edges = []
        self.nodes = dict()
        self.stringPath = []
        self.sequence = ''

    def findEdges(self):
        """
        Finds all the edges between the kmers within self.patterns by comparing the (k-1) length suffixes and prefixes
        of each kmer. If the suffix of one pattern is equal to the prefix of another pattern, then a tuple containing
        both patterns is added to the self.edges.

        Called upon by self.findNodes().

        :return: updates self.edges
        """
        for pattern in self.patterns:
            for otherPattern in self.patterns:
                if pattern[1:] == otherPattern[:-1]:
                    self.edges.append((pattern, otherPattern))

    def findNodes(self):
        """
        Uses self.edges generated in self.findEdges() to generate the dictionary self.nodes. Refer to the __init__
        method for more information.

        For every edge in self.pattern, each pattern contained in edge is added to the self.node. If a pattern shows up
        in the first position of the tuple (node1), then we update the pattern in the second position (node2) to have
        an incoming node from node1. We also update the value of node1 in self.nodes to include node2 in the outgoing
        node list.

        Called upon by self.getStart()

        :return: updates self.nodes
        """
        self.findEdges()  # call on self.findEdges()
        for edge in self.edges:
            node1 = edge[0]
            node2 = edge[1]
            # Add an outgoing edge to node1 (edge[0])
            if node1 not in self.nodes.keys():
                self.nodes[node1] = ([], [])  # (inEdges, outEdges)
            self.nodes[node1][1].append(node2)  # add the outgoing edge
            # Add an incoming edge to node2 (edge[0])
            if node2 not in self.nodes.keys():
                self.nodes[node2] = ([], [])
            self.nodes[node2][0].append(node1)

    def getStart(self):
        """
        Finds the start node of the Eulerian path. The start node of an Eulerian path is a node whose number of outgoing
        edges - number of incoming edges is equal to 1. This method uses this rule by going through each of the nodes in
        self.nodes and calculating the difference between the length of the outgoing node list minus the length of the
        incoming node list. If this is equal to 1, that node is set to the startNode. This node is then returned to be
        used in self.findPath().

        If there is no node that meets this condition, then a node is chosen at random (this is the case of an Eulerian
        cycle).

        :return: the node in which |outgoing nodes| - |incoming nodes| = 1
        """
        self.findNodes()
        startNode = random.choices(list(self.nodes.keys()))[0]
        for node in self.nodes.keys():
            # Checks for start node
            if len(self.nodes[node][1]) - len(self.nodes[node][0]) == 1:  # if outDegrees - inDegrees = 1
                startNode = node
        return startNode

    def findPath(self, v):
        """
        This method uses a recursive implementation to find the Eulerian cycles within a graph. It begins by setting
        the outgoing nodes to a list called 'vEdges.' In the case of its first run, the parameter v is provided by the
        graph's starting node found in self.getStart().

        Then while there are still outgoing edges of node v, we randomly select one of the outgoing nodes and
        remove this edge from self.edges. We then call findPath() again, but this time with the selected node. Once we
        have reached all the edges of a node, we insert this node to the front of list representing the Eulerian path,
        self.stringPath.

        NOTE: I remembered using a recursive algorithm for finding Eulerian paths in my old algorithms class
        and based my code off of the content on this page: https://cp-algorithms.com/graph/euler_path.html

        :param v: node V
        :return: updates self.stringPath

        """
        vEdges = self.nodes[v][1]  # outgoing edges

        while len(vEdges) > 0:
            newV = random.choices(vEdges)[0]
            self.nodes[v][1].remove(newV)
            self.findPath(newV)
        self.stringPath.insert(0, v)

    def printString(self):
        """
        Using the patterns within self.stringPath, prints the generated sequence.

        Calls on self.findPath() using the node outputted from self.getStart().

        NOTE: this is the only function that needs to be called in main(), since all of the other functions are called
        on each other.
        """
        self.findPath(self.getStart())
        for i in range(0, len(self.stringPath)):
            if i == 0:
                self.sequence += self.stringPath[i]
            else:
                self.sequence += self.stringPath[i][-1]
        print(self.sequence)


def main():
    """
    The main() method initializes the list of patterns that is generated by the class Graph and its methods. It reads
    sys.stdin and appends all kmers contained within an input file. We then initialize the object Graph using these
    patterns, then call upon the method printString() contained within the Graph class.
    :return:
    """
    patterns = []
    lineCount = 0

    # Read the input file
    for line in sys.stdin:
        if lineCount > 0:
            patterns.append(line.strip())
        lineCount += 1

    myGraph = Graph(patterns)
    myGraph.printString()


if __name__ == "__main__":
    main()
