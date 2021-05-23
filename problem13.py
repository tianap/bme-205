#!/usr/bin/env python
import sys
import random
import re

"""
Imported modules
    sys  # used to read in files from stdin
    random  # used to randomly select items fom lists
    re  # used for parsing through the input file in main()
"""


class Eulerian:
    """
    The Eulerian class contains all of the nodes and edges of a given list of edges. Its methods are able to generate
    the Eulerian path across all nodes.

    It contains 4 methods in addition to the initialization method def __init__. The only method that is called on in
    main() is self.printPath(), which calls on other methods in order to construct the Eulerian path.

    Attributes:
        self.order
        self.edges
        self.nodes

    Methods:
          self.findNodes()
          self.getStart()
          self.findPath()
          self.printPath()

    """
    def __init__(self, edgeList):
        """
        The initialization method for the class Eulerian. Takes in the parameter edgeList from a user.

        :param edgeList: a list of tuples representing edges

        Attributes:
            self.order  # a list that keeps track of the Eulerian path
            self.edges  # the list of edges within the graph, assigned to the parameter edgeList
            self.nodes  # dictionary in which each key is a pattern and each value is a tuple containing two lists: the
                        # first containing the incoming nodes, the second containing the outgoing nodes
                ex.
                    self.nodes = {
                        'a': ([d],[b])
                        'b': ([a],[c])
                        'c': ([b],[d])
                        'd': ([c],[a])
                    }
        """
        self.order = []
        self.edges = edgeList  # list of tuples representing edges
        self.nodes = dict()  # Node: ([inComing],[outGoing])

    def findNodes(self):
        """
        Uses self.edges generated in self.findEdges() to generate the dictionary self.nodes. Refer to the __init__
        method for more information.

        For every edge in self.pattern, each pattern contained in edge is added to the self.node. If a pattern shows up
        in the first position of the tuple (node1), then we update the pattern in the second position (node2) to have
        an incoming node from node1. We also update the value of node1 in self.nodes to include node2 in the outgoing
        node list.

        Called upon in main().
        :return: updates self.nodes
        """
        for edge in self.edges:
            # Add an outgoing edge to node edge[0]
            node1 = edge[0]
            node2 = edge[1]
            if node1 not in self.nodes.keys():
                self.nodes[node1] = ([], [])  # (inEdges, outEdges)
            self.nodes[node1][1].append(node2)  # add the outgoing edge

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
        self.order.

        NOTE: I remembered using a recursive algorithm for finding Eulerian paths in my old algorithms class
        and based my code off of the content on this page: https://cp-algorithms.com/graph/euler_path.html

        :param v: node V
        :return: updates self.order

        """
        vEdges = self.nodes[v][1]  # outgoing edges
        while len(vEdges) > 0:
            newV = random.choices(vEdges)[0]
            self.nodes[v][1].remove(newV)
            self.findPath(newV)
        self.order.insert(0, v)

    def printPath(self):
        """
        Prints the Eulerian path, with each node separated by '->'. Calls upon self.findPath() along with the node
        returned from self.getStart() as its parameter.

        Called upon by main().

        :return: prints the Eulerian path
        """
        self.findPath(self.getStart())
        print(*self.order, sep='->')


def main():
    """
    The main methods reads in an adjacency list representing edges from file from stdin. It then intializes an Eulerian
    object with these edges, represented as a list of tuples. We then find the nodes for this object using the method
    findNodes(), then call printPath() to print the Eulerian path.

    :return:
    """
    edges = []
    for line in sys.stdin:
        lineList = re.split(', | -> |,', line.strip())
        for i in lineList[1:]:
            if i.isdigit():
                edges.append((int(lineList[0]), int(i)))

    eulerPath = Eulerian(edges)
    eulerPath.findNodes()
    eulerPath.printPath()


if __name__ == "__main__":
    main()
