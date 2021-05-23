#!/usr/bin/env python

########################################################################
# File: problem18.py
#  executable: py problem18.py < input.txt > output.out
# Purpose: Given a directed acyclic graph (DAG), returns the length of the longest path and the longest path.
#
# Classes:
#   Node
#   Edge
#   DAG
#
# Author: Tiana Pereira
# History:      11/3/2020 Created
#
########################################################################

import sys, random, copy, re


class Node:
    """
    A class to represent a node in a graph.

    Methods
    -------
    def __init__()
    def addInEdge(i)
    def addOutEdge(i)
    """

    def __init__(self, name, inComing, outGoing):
        """
        Constructs all the necessary attributes for the Node object.

        Parameters
        ----------
        :param name: name of the node
        :param inComing: the names of the incoming nodes
        :param outGoing: the names of the outgoing nodes

        Attributes
        ----------
        name : str
            name of the node, given by user
        inEdges : list
            the names of the incoming nodes, assigned to inComing parameter
        outEdges : list
            the names of the outgoing nodes, assigned to outGoing parameter
        maxPredecessor : int
            preceding node whose score is the max of all preceding nodes, initialized to -1.
        score : int
            stores the score of the node (score(maxPredecessor) + weight(maxPredecessor -> node)), initialized to -1000
        """

        self.name = name
        self.inEdges = inComing
        self.outEdges = outGoing
        self.maxPredecessor = -1
        self.score = -1000

    def addInEdge(self, i):
        """
        Adds a preceding node to present an edge into the node.
        :param i: a new node to add
        :return: updates self.inEdges
        """

        self.inEdges.append(i)

    def addOutEdge(self, i):
        """
        Adds a succeding node to present an edge out of the node.
        :return: updates self.outEdges
        :param i: a new node to add
        :return: updates self.outEdges
        """

        self.outEdges.append(i)


class Edge:
    """
    A class to represent an edge in a graph.

    """

    def __init__(self, name, weight):
        """
        Constructs all the necessary attributes for the Edge object.

        Parameters
        ----------
        :param name: name of the edge
        :param weight: the weight of the edge

        Attributes
        ----------
        name : tuple
            name of the edge in the form of a tuple, given by user
        weight : int
            the weight of the edge
        """

        self.name = name
        self.weight = weight


class DAG:
    """
    A class to represent a directed acyclic graph (DAG). Contains the methods that construct the longest path and its
    score.

    Methods
    -------
    parseInput()
    topOrdering()
    maxPredecessor()
    longestPath()
    outputPath()

    """

    def __init__(self, inputList):
        """
        Constructs all the necessary attributes for the DAG object.

        Parameters
        ----------
        :param inputList: list of strings contains the start node, end node, and all edges with their weights

        Attributes
        ----------
        input : list
            assigned to inputList
        nodes : dict
            contains all of the nodes in the DAG (keys are node names, values are references to Node objects)
        edges : dict
            contains all of the edges in the DAG (keys are edge names, values are references to Edge objects)
        topOrder : list
            the topological order of the DAG
        nodesCopy : dict
            a deep copy of self.nodes, used to generate the topological order of the DAG
        startNode : str
            the start node provided by inputList, initialized to -1
        endNode : str
            the end node provided by inputList, initialized to -1
        path : list
            the longest path
        """

        self.input = inputList
        self.nodes = dict()
        self.edges = dict()
        self.topOrder = []
        self.nodesCopy = dict()
        self.startNode = -1
        self.endNode = -1
        self.path = []

    def parseInput(self):
        """
        Parses the start node, end node, and edges out of self.input. Generates Node and Edge objects, and updates
        self.nodes and self.edges to include these references.

        :return: updates self.startNode, self.endNode, self.nodes, self.edges
        """

        # Iterate through list version of input
        for line in range(0, len(self.input)):
            # Add the start node
            if line == 0:  # the first line
                self.startNode = self.input[line]

            # Add the end node
            elif line == 1:
                self.endNode = self.input[line]

            # Add an edge
            else:
                # Create an Edge() object using the values provided
                values = re.findall('[0-9]+', self.input[line])
                name = (values[0], values[1])  # name is represented by a tuple
                newEdge = Edge(name, int(values[2]))

                # Add new edge to this graph's edges
                self.edges[name] = newEdge  # key = name (type tuple)

                # Add/update incoming node, depending on whether it appears in self.nodes
                if values[0] not in self.nodes.keys():  # if the node has not been added yet
                    newNode = Node(values[0], [], [values[1]])
                    self.nodes[values[0]] = newNode
                else:
                    self.nodes[values[0]].addOutEdge(values[1])

                # Add/update outgoing node, depending on whether it appears in self.nodes
                if values[1] not in self.nodes.keys():
                    newNode = Node(values[1], [values[0]], [])
                    self.nodes[values[1]] = newNode
                else:
                    self.nodes[values[1]].addInEdge(values[0])

        # Create a deep copy of self.nodes (used for generating the topological ordering)
        self.nodesCopy = copy.deepcopy(self.nodes)

    def topOrdering(self):
        """
        Calls on parseInput(), then generates the topological ordering of the DAG. Based off of the pseudocode for
        TopologicalOrdering(Graph) from the textbook/lecture.

        :return: updates self.topOrder to contain the topological ordering of the DAG.
        """

        self.parseInput()
        candidates = []  # all nodes with no incoming edges

        # Add all nodes with no input edges to the candidates list
        for node in self.nodes.keys():
            if len(self.nodes[node].inEdges) == 0:
                candidates.append(node)

        # Loop through all candidate nodes
        while len(candidates) > 0:
            a = random.choice(candidates)  # arbitrary node from candidates
            self.topOrder.append(a)  # add to list
            candidates.remove(a)  # remove from candidates

            # for outgoing edge from node a to node b
            for b in self.nodesCopy[a].outEdges:
                # Remove a from incoming edges of b, modifying the Node object
                self.nodesCopy[b].inEdges.remove(a)

                # If there are no other incoming edges, add to candidates
                if len(self.nodesCopy[b].inEdges) == 0:
                    candidates.append(b)

    def maxPredecessor(self, a):
        """
        Finds the maximum predecessor of a given node, updates the maxPredecessor attribute of the Node object

        NOTE: By default, the score of each node is initialized to -1000 when the Node object is constructed.

        :param a: a node
        :return: the max predecessor node, the score of node a
        """

        # Initialize the maxNode to be -1
        maxNode = -1

        # Loop through incoming edges of a
        for predecessor in self.nodes[a].inEdges:
            # If the score of the predecessor + edge weight is greater than its current node
            if self.nodes[predecessor].score + self.edges[(predecessor, a)].weight > self.nodes[a].score:
                maxNode = predecessor  # update the max predecessor
                self.nodes[a].score = self.nodes[predecessor].score + self.edges[(predecessor, a)].weight  # update score
        # updates maxPredecessor attribute if node a
        self.nodes[a].maxPredecessor = maxNode
        return maxNode, self.nodes[a].score

    def longestPath(self):
        """
        Calls upon self.topOrdering() to generate the topological order, then calls on self.maxPredecessor to generate
        the score of the longest path. Based on the LongestPath(Graph, source, sink) algorithm in the textbook/lecture.

        :return: updates self.path
        """

        self.topOrdering()  # generate the topological order
        self.nodes[self.startNode].score = 0  # initialize the score of the start node to be 0

        # Slice the self.topOrder to start with self.startNode and end with self.endNode
        for a in self.topOrder[self.topOrder.index(self.startNode):self.topOrder.index(self.endNode) + 1]:
            if a != self.startNode:  # skip the source
                # Updates the maxPredecessor of node a
                self.maxPredecessor(a)

    def outputPath(self):
        """
        Calls on self.longestPath() to find the longest path, then outputs the output path by inserting the max
        predecessor starting from self.endNode to self.startNode.
        :return: prints the score of the longest path and the longest path itself.
        """
        # Begin by calling longestPath() to find the path with the highest scores
        self.longestPath()

        # Insert the last node the last node into the path
        thisNode = self.endNode
        self.path.insert(0, self.endNode)

        # Iterate through the max predecessors for each node
        while thisNode != self.startNode:
            # Add the max predecessor for each node
            self.path.insert(0, self.nodes[thisNode].maxPredecessor)

            # Reset thisNode to the max predecessor node of thisNode
            thisNode = self.nodes[thisNode].maxPredecessor

        # Print the output
        print(self.nodes[self.endNode].score)
        print(*self.path, sep='->')


def main():
    """
    Reads from stdin and creates a DAG object with this input, then calls on outputPath() to print the score of the
    longest path and the longest path itself.
    """
    inputList = []
    for line in sys.stdin:
        inputList.append(line.strip())
    myDAG = DAG(inputList)
    myDAG.outputPath()


if __name__ == "__main__":
    main()
