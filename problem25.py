#!/usr/bin/env python

########################################################################
# File: problem25.py
#  executable: py problem25.py < input.txt > output.out
# Purpose: Given a string x, the alphabet sigma, states, a transition matrix, and emission matrix, returns the
# conditional probability Pr(pi_i=k|x) that the HMM was in statek at step i (for each state k and for each step i).
#
# Classes:
#   Node
#   HMM
#
# Author: Tiana Pereira
# History:      12/10/2020 Created
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
            the calculated forward score of the node, initialized to 0; represents the sum of all edge weights and score
            of preceding nodes
        backward : float
            the calculated backward score of the node, initialized to 0; represents the sum of all edge weights and scores
            of succeeding nodes
        prob : float
            the calculated probability that the emitted string x was in state k at time i, calculated by the
            forward-backward algorithm; initialized to 0
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
        self.forward = 0
        self.backward = 0
        self.prob = 0


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
        findTotalProb():
            Calculates the probability self.prob to be used in the forward-backward algorithm
        findNodeProb():
            Calculates the probability that the emitted string was in state k at time i, updates the prob attribute of
            the node object
        weight():
            Calculates the weight of an edge between two nodes
        findForward():
            Finds the forward score of a node
        findBackward():
            Finds the backward score of a node
        viterbi():
            Implements the Viterbi alogorithm to find the maximum probability Pr(x, pi)
        softDecoding():
            Calls on self.findTotalProb() and self.findNodeProb() to calculate all probabilities for each node, then
            prints to STDOUT.

    """

    def __init__(self, userInput):
        """
        Constructs necessary attributes for the HMM class.

        Parameters:
            userInput : list of lines from input from main()

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

    def findTotalProb(self):
        """
        Calls on self.viterbi(), then finds the probability self.prob that the HMM emits x by summing through forward
        scores, prints self.prob.
        """

        self.viterbi()
        endNodes = [key for key in self.nodeDict.keys() if key[1] == len(self.x) - 1]
        for node in endNodes:
            self.prob += self.nodeDict[node].forward

    def findNodeProb(self):
        """
        Loops through the nodes in reverse from the sink to the source, finds the backward sum of succeeding nodes for
        each node. Then calculates the probability of the nodes (the probability that the emitted string was in state k
        at time i), based off of the forward-backward algorithm.

        NOTE: When finding the backwards sum of the last nodes before the sink, the scores of backwards are set to 1.
        This is because the sink node (which is not included in the node dictionaries) has a probability of 1, and the
        weight between each of the nodes connecting to the sink have a weight of 1 (since there is no other node to go
        to).

        If the node is not in the last position, then it calls self.findBackward(k, i).

        :return: updates the prob attribute of each node in self.nodeDict
        """
        # Find the backwards probabilities, starting at the last position of the emitted string x:
        for i in reversed(range(len(self.x))):
            # For all states
            for k in self.states:
                # The last nodes before the sink, backward = 1
                if i == len(self.x) - 1:
                    self.nodeDict[(k, i)].backward = 1
                else:  # calculate the backwards using self.findBackward()
                    self.nodeDict[(k, i)].backward = self.findBackward(k, i)
                # Calculate the probability
                self.nodeDict[(k, i)].prob = round(
                    (self.nodeDict[(k, i)].forward * self.nodeDict[(k, i)].backward) / self.prob, 4)

    def weight(self, l, k, x):
        """
        Calculates the weight of an edge between nodes of states k (with symbol x) and another state l.

        :param l: other state l
        :param k: current state k
        :param x: observed symbol at node k

        :return weight (emission*transmission)
        """

        otherState = self.accessionDict[l]
        state = self.accessionDict[k]
        symbol = self.accessionDict[x]
        emission = self.emission[state][symbol]
        transition = self.transition[otherState][state]
        return emission * transition

    def findForward(self, k, i):
        """
        Calculates the forward score for a given node to be stored in the forward attribute of the Node object.
        :param k: the state of the node
        :param i: the position of the node
        :return: updates the forward attribute of the Node object
        """
        forward = 0

        for l in self.states:
            forward += self.nodeDict[(l, i - 1)].forward * self.weight(l, k, self.x[i])

        return forward

    def findBackward(self, k, i):
        """
        Calculates the backward score for a given node by summing over the backward scores for succeeding nodes times
        the weight between the current node and each of the succeeding nodes.

        NOTE: the weight of the nodes is calculated between the current node and its next node, so the emitted symbol is
        in the i+1 position.

        :param k: state of the current node
        :param i: position of the current node
        :return: updates the backward attribute of the Node object

        """

        # Initialize backward to 0
        backward = 0

        # For each succeeding node
        for l in self.states:
            backward += self.nodeDict[(l, i + 1)].backward * self.weight(k, l, self.x[i + 1])
        return backward

    def viterbi(self):
        """
        Calls on self.setNodes(), then implements the Viterbi algorithm from the textbook. Forward scores are calculated
        here (backward scores are calculated in self.findNodeProb() and self.findBackward().
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

    def softDecoding(self):
        """
        Calls on self.findTotalProb() and self.findNodeProb() for the probabilities of each node (using the
        forward-backward algorithm). Then prints the formatted probabilities
        :return:
        """
        self.findTotalProb()
        self.findNodeProb()
        print(' '.join(self.states))
        for i in range(len(self.x)):
            row = []
            for k in self.states:
                row.append(str(self.nodeDict[(k, i)].prob))
            print(' '.join(row))


def main():
    """
    Reads from STDIN, then passes the input as an argument for the HMM class.
    """

    myInput = []
    for line in sys.stdin:
        myInput.append(line)
    myHMM = HMM(myInput)
    myHMM.softDecoding()


if __name__ == "__main__":
    main()
