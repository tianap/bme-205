#!/usr/bin/env python

########################################################################
# File: problem23.py
#  executable: py problem23.py < input.txt > output.out
# Purpose: Given a string x, the alphabet sigma, states, and a hidden path, calculates the transition and emission
#  matrices
# Classes:
#   HMM
#
# Author: Tiana Pereira
# History:      12/3/2020 Created
#
########################################################################
import sys


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


    Methods:
        parseInput():
            Calculates the probability of the hidden path occurring, prints said probability.
        buildAccessDict():
            Builds the accession dictionary

    """

    def __init__(self, userInput):
        """
        Constructs necessary attributes for the HMM class.

        Parameters:
        :param userInput: list of all lines from STDIN

        Assumptions:
            This class and construction assumes input from STDIN is formatted the same as problem23 in Rosalind.
        """

        self.input = userInput
        self.x = ''  # string x containing alphabet sigma
        self.hiddenPath = ''
        self.sigma = []
        self.states = []
        self.emission = []
        self.transition = []
        self.accessionDict = dict()

    def parseInput(self):
        """
        Parses the input from STDIN and assigns respective values to x, sigma, states, emission, and transition.
        """

        for i in range(0, len(self.input)):
            if self.input[i].strip() != '--------':
                if i == 0:  # string x
                    self.x = self.input[i].strip()
                elif i == 2:  # sigma
                    self.sigma = self.input[i].strip().split()
                elif i == 4:  # hiddenPath
                    self.hiddenPath = self.input[i].strip()
                else:
                    self.states = self.input[i].strip().split()

    def buildAccessDict(self):
        """
        Builds accession dictionary by looping through self.sigma and self.states and keys and values
        corresponding to their symbol (char) and integer position (i).

        :return: updates self.accessionDict
        """
        self.parseInput()
        for i in range(0, len(self.sigma)):
            self.accessionDict[self.sigma[i]] = i

        for i in range(0, len(self.states)):
            self.accessionDict[self.states[i]] = i

    def buildMatrices(self):
        """
        Builds transition and emission matrices with the number of rows corresponding to the number of states.

        :return: builds self.transition, self.emission
        """
        for i in range(0, len(self.states)):
            self.transition.append([])
            self.emission.append([])

    def findParameters(self):
        """
        Calls on self.buildAccessDict(), self.buildMatrices(), self.findTransitions(), self.findEmissions(), to find
        the probabilities in the transition and emission matrices.

        This is called upon by self.printParameters()
        """
        self.buildAccessDict()
        self.buildMatrices()
        self.findTransitions()
        self.findEmissions()

    def findTransitions(self):
        """
        Calculates the transition probabilities based on the hidden path by finding the number of times a transition
        between two states, l and k, occurs divided by the total number of transitions from state l.

        Assumptions:
            - Since we are only finding the number of transitions FROM a particular state, just need to count the number
            of times that state appears between the first and second-to-last position of the hidden Path.

        :return: updates self.transition with probabilities
        """

        for l in self.states:
            total = self.hiddenPath[:-1].count(l)
            for k in self.states:
                count = 0
                for i in range(0, len(self.hiddenPath) - 1):
                    trans = self.hiddenPath[i:i + 2]
                    if trans[0] == l and trans[1] == k:
                        count += 1
                # Store the probability
                try:
                    prob = count / total
                except ZeroDivisionError:
                    if count == 0 and total == 0:
                        prob = 1 / len(self.states)
                    else:
                        prob = 0
                self.transition[self.accessionDict[l]].append("{0}".format(round(prob,3)))

    def findEmissions(self):
        """
        Calculates the emission probabilities based on the hidden path and string x by the following formula for a given
        symbol from sigma and a state k:

        Pr(symbol|k) / # times k appears in hidden path

        :return: updates self.emission with probabilities
        """
        for k in self.states:
            total = self.hiddenPath.count(k)
            for symbol in self.sigma:
                count = 0
                for i in range(0, len(self.x)):
                    if self.x[i] == symbol and self.hiddenPath[i] == k:
                        count += 1

                try:
                    prob = count / total
                except ZeroDivisionError:
                    if count == 0 and total == 0:
                        prob = 1 / len(self.states)
                    else:
                        prob = 0
                self.emission[self.accessionDict[k]].append("{0}".format(round(prob, 3)))

    def printParameters(self):
        """
        Calls on self.findParameters(), then prints tab-separated lines that display matrices (according to the Rosalind
        format). NOTE: This is the only method of the HMM class that is directly called on in main().

        :return: prints the transition and emission matrices
        """
        self.findParameters()
        print(' '.join(self.states))
        i = 0
        for state in self.states:
            print(state + ' ' + ' '.join(self.transition[i]))
            i += 1
        print('--------')
        print(' ' + ' '.join(self.sigma))
        i = 0
        for state in self.states:
            print(state + ' ' + ' '.join(self.emission[i]))
            i += 1


def main():
    """
    Reads from STDIN, then passes the input as an argument for the HMM class.
    """

    myInput = []
    for line in sys.stdin:
        myInput.append(line)
    myHMM = HMM(myInput)
    myHMM.printParameters()


if __name__ == "__main__":
    main()
