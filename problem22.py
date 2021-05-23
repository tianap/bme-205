#!/usr/bin/env python

########################################################################
# File: problem22.py
#  executable: py problem22.py < input.txt > output.out
# Purpose: Given a string x, the alphabet sigma, states, a transition matrix, and emission matrix, returns the
# conditional probability that the HMM emits the string x
# Classes:
#   Node
#   HMM
#
# Author: Tiana Pereira
# History:      11/23/2020 Created
#
########################################################################
import sys


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
        forward : float
            the calculated forward score of the node, initialized to -1; represents the sum of all edge weights and score
            of preceding nodes
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
        self.forward = -1


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
        prob : float
            the conditional probability that the HMM emits the string x


    Methods:
        parseInput():
            Calculates the probability of the hidden path occurring, prints said probability.
        buildAccessDict():
            Builds the accession dictionary
        setNodes():
            Creates Node objects for each node in the HMM
        findProb():
            Calculates the probability self.prob
        weight():
            Calculates the weight of an edge between two nodes
        findForward():
            Finds the forward score of a node
        viterbi():
            Implements the Viterbi alogorithm to find the maximum probability Pr(x, pi)

    """

    def __init__(self, userInput):
        """
        Constructs necessary attributes for the HMM class.

        Parameters:
        :param states: list of states of the HMM
        :param transitionMatrix: 2x2 list of lists representing the transition matrix

        Assumptions:
            This class and construction assumes that there are only two states involved in the path, and that the
            transition matrix is a 2x2 matrix.
        """

        self.input = userInput
        self.x = ''  # string x containing alphabet sigma
        self.sigma = []
        self.states = []
        self.emission = []
        self.transition = []
        self.accessionDict = dict()
        self.nodeDict = dict()
        self.prob = 0

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
                        i += 1
                    i += 2
                    for j in range(0, len(self.states)):
                        self.emission.append([float(i) for i in self.input[i].strip().split()[1:]])
                        i += 1

    def buildAccessDict(self):
        """
        Builds accession dictionary by looping through self.sigma and self.states and keys and values
        corresponding to their symbol (char) and integer position (i).

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

    def findProb(self):
        """
        Calls on self.viterbi(), then finds the probability self.prob that the HMM emits x by summing through forward
        scores, prints self.prob.
        """

        self.viterbi()
        endNodes = [key for key in self.nodeDict.keys() if key[1] == len(self.x) - 1]
        for node in endNodes:
            self.prob += self.nodeDict[node].forward
        print(self.prob)

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

    def findForward(self, k, i):
        """
        Calculates the forward score for a given node to be stored in the forward attribute of the Node object.
        :param k:
        :param i:
        :return:
        """
        forward = 0

        for l in self.states:
            forward += self.nodeDict[(l, i - 1)].forward * self.weight(l, k, self.x[i])

        return forward

    def viterbi(self):
        """
        Calls on self.setNodes(), then implements the Viterbi algorithm from the textbook.
        """

        self.setNodes()
        for i in range(0, len(self.x)):  # from source to sink
            if i == 0:
                for k in self.states:
                    access = (self.accessionDict[k], self.accessionDict[self.x[i]])
                    self.nodeDict[(k, i)].forward = self.emission[access[0]][access[1]] * (1 / len(self.states))
            else:
                for k in self.states:
                    self.nodeDict[(k, i)].forward = self.findForward(k, i)


def main():
    """
    Reads from STDIN, then passes the input as an argument for the HMM class.
    """

    myInput = []
    for line in sys.stdin:
        myInput.append(line)
    myHMM = HMM(myInput)
    myHMM.findProb()


if __name__ == "__main__":
    main()
