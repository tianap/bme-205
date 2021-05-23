#!/usr/bin/env python

########################################################################
# File: problem21.py
#  executable: py problem21.py < input.txt > output.out
# Purpose: Given a string x, the alphabet sigma, states, a transition matrix, and emission matrix, returns the
# unconditional probability of all possible paths.
#
# Classes:
#   Node
#   HMM
#
# Author: Tiana Pereira
# History:      11/20/2020 Created
#
########################################################################
import sys, math


class Node:
    """
    Class representing a node in an HMM.

    Attributes:
        symbol : str
            the symbol observed at position i in string x in HMM class
        state : str
            the state for this node
        position : int
            the position in string x in HMM class
        maxPred : tuple
            stores the tuple (k,i) of the preceding node with the max score, initialized to -1 until calculated
        score : float
            the calculated score of thenode
    """

    def __init__(self, x, k, i):
        """
        Constructs the necessary attributes for the Node class
        :param x: symbol observed as position i of string x in HMM class
        :param k: state of this node
        :param i: position in string x in HMM class
        """

        self.symbol = x
        self.state = k
        self.position = i
        self.maxPred = -1
        self.score = -1


class HMM:
    """
    Class representing an HMM and its hidden path, states, and transition matrix.

    Attributes:
        input : list
            input list taken from sys.stdin
        x : str
            an observed string x
        sigma : list
            list of symbols contained in x
        hiddenPath : str
            hidden path of the HMM
        states : list
            states of the HMM
        emission : list
            list of lists representing emission matrix of the HMM
        transition : list
            lists of lists representing transition matrix of the HMM
        accessionDict : dict
            dictionary storing the accession values of a char symbol corresponding to an index for the emission matrix
        nodeDict : dict
            dictionary storing all the pointers to Node objects


    Methods:
        parseInput():
            Calculates the probability of the hidden path occurring, prints said probability.
        buildAccessDict():
            Builds the accession dictionary
        setNodes():
            Creates Node objects for each node in the HMM
        findHiddenPath():
            Generates the hidden path by following the backtracking pointers
        weight():
            Calculates the weight of an edge between two nodes
        findMax():
            Finds the maximum preceding node
        viterbi():
            Implements the Viterbi alogorithm to find the maximum probability Pr(x, pi)

    """

    def __init__(self, userInput):
        """
        Constructs necessary attributes for the HMM class.

        Parameters:
        :param userInput: list of lines read in from sys.stdin

        """

        self.input = userInput
        self.x = 0  # string x containing alphabet sigma
        self.sigma = []
        self.hiddenPath = []
        self.states = []
        self.emission = []
        self.transition = []
        self.accessionDict = dict()
        self.nodeDict = dict()

    def parseInput(self):
        """
        Parses the input from STDIN and assigns respective values to x, sigma, states, emission, and transition.
        """

        for i in range(0, len(self.input)):
            if self.input[i].strip() != '--------':
                if i == 0:
                    self.x = self.input[i].strip()
                elif i == 2:
                    self.sigma = self.input[i].strip().split()
                elif i == 4:
                    self.states = self.input[i].strip().split()
                elif i == 7:
                    for j in range(0, len(self.states)):
                        self.transition.append([float(i) for i in self.input[i].strip().split()[1:]])
                        print(i, self.input[i].strip())
                        i += 1
                    i += 2
                    for j in range(0, len(self.states)):
                        self.emission.append([float(i) for i in self.input[i].strip().split()[1:]])
                        i += 1

    def buildAccessDict(self):
        """
        Calls on self.parseInput(), then builds accession dictionary by looping through self.sigma and self.states and
        keys and values corresponding to their symbol (char) and integer position (i).

        Called upon by self.findProb()
        :return: updates self.accessionDict
        """

        self.parseInput()
        for i in range(0, len(self.sigma)):
            self.accessionDict[self.sigma[i]] = i

        for i in range(0, len(self.states)):
            self.accessionDict[self.states[i]] = i

    def setNodes(self):
        """
        Calls on self.buildAccessDict() then creates Node objects that are added to self.nodeDict.

        :return: updates self.nodeDict
        """
        self.buildAccessDict()
        for k in self.states:  # for each state
            for i in range(0, len(self.x)):  # for each position i
                newNode = Node(self.x[i], k, i)
                self.nodeDict[(k, i)] = newNode

    def findHiddenPath(self):
        """
        Calls on self.viterbi() to go through the Viterbi algorithm found in the book. Then generates the hidden path by
        using backtracking pointers.

        Implementation:
            self.viterbi() calls on a variety of other methods that determines the maximum-scoring preceding node, which
            is stored in the maxPred attribute of the Node object. This attribute will store the the tuple (k,i) that
            represents the state of the node at position i. By design, the values in self.nodeDict point to Node objects.
            Therefore, when generating the hidden path string, this algorithm goes through backtracking pointers
            contained in each of the nodes.

        """
        self.viterbi()

        # Find the end state
        endNodes = [key for key in self.nodeDict.keys() if key[1] == len(self.x) - 1]
        maxProb = 0
        lastNode = 0
        for node in endNodes:
            if self.nodeDict[node].score > maxProb:
                maxProb = self.nodeDict[node].score
                lastNode = node
        # Add the last state of the hidden path
        self.hiddenPath.insert(0, self.nodeDict[lastNode].maxPred[0])

        lastNode = self.nodeDict[lastNode].maxPred

        for i in range(len(self.x), 0, -1):
            if lastNode == -1:  # Time to end the loop
                break

            # Insert the state of this node in the first position
            self.hiddenPath.insert(0, lastNode[0])

            # Reset lastNode to be the maximum-scoring preceding node
            lastNode = self.nodeDict[lastNode].maxPred

        # Print the formatted hidden path
        print(''.join(self.hiddenPath))

    def weight(self, l, k, x):
        """
        Calculates the weight of an edge between nodes of states k (with symbol x) and l.
        :param l: preceding state l
        :param k: current state k
        :param x: observed symbol at node k

        :return weight (emission*transmission)
        """
        prevState = self.accessionDict[l]
        state = self.accessionDict[k]
        symbol = self.accessionDict[x]
        emission = self.emission[state][symbol]
        transition = self.transition[prevState][state]
        return emission * transition

    def findMax(self, k, i):
        """
        Finds the maximum-scoring node preceding a given node of state k at position i. Based off of the definitions of
        score in the textbook.
        :param k: a state k
        :param i: position i in self.x
        :return maxScore, prevState: the maximum score from the preceding node and the state of such preceding node
        """
        maxScore = 0
        prevState = 0

        # Loop through the nodes at position i-1
        for l in self.states:
            prevScore = self.nodeDict[(l, i - 1)].score * self.weight(l, k, self.x[i])
            if prevScore > maxScore:  # Replace maximum score with current score, store the state l
                maxScore = prevScore
                prevState = (l, i - 1)

        return maxScore, prevState

    def viterbi(self):
        """
        Calls on self.setNodes(), then implements the Viterbi algorithm from the textbook.
        """

        self.setNodes()
        for i in range(0, len(self.x)):  # from source to sink
            if i == 0:
                for k in self.states:
                    access = (self.accessionDict[k], self.accessionDict[self.x[i]])
                    self.nodeDict[(k, i)].score = self.emission[access[0]][
                        access[1]]  # score(source) = 1, weight of source to states A and B equal
            else:
                for k in self.states:
                    maxL = self.findMax(k, i)
                    self.nodeDict[(k, i)].score = maxL[0]
                    self.nodeDict[(k, i)].maxPred = maxL[1]


def main():
    """
    Reads from STDIN, then passes the input as an argument for the HMM class.
    """

    myInput = []
    for line in sys.stdin:
        myInput.append(line)
    myHMM = HMM(myInput)
    myHMM.findHiddenPath()


if __name__ == "__main__":
    main()
