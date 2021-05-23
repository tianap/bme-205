#!/usr/bin/env python

########################################################################
# File: problem20.py
#  executable: py problem20.py < input.txt > output.out
# Purpose: Given a string x, the alphabet sigma that x was constructed, a hidden path, the states, and the emission
# matrix, returns the conditional probability that string x is emitted by the HMM.
#
# Classes:
#   HMM
#
# Author: Tiana Pereira
# History:      11/17/2020 Created
#
########################################################################
import sys


class HMM:
    """
    Class representing an HMM and its hidden path, states, and transition matrix.

    Attributes:
        x : str
            an observed string x
        sigma : list
            a list of the symbols in x
        hidden : str
            hidden path of the HMM
        states : list
            states of the HMM
        emission : list of list, representing 2x2 emission matrix
            transition matrix of the HMM
        accessionDict : dict
            dictionary storing the accession values of a char symbol corresponding to an index for the emission matrix
        self.prob : float
            conditional probability that string x will be emitted by HMM

    Methods:
        buildAccessionDict():
            Creates a dictionary that stores the integer value corresponding to the characters in sigma and states
        findProb():
            Calculates the conditional probability of the string x occurring, prints said probability.


    """

    def __init__(self, x, sigma, hiddenPath, states, emissionMatrix):
        """
        Constructs necessary attributes for the HMM class.

        Parameters:
        :param x: a string x
        :param sigma: list of the symbols in x
        :param hiddenPath: string representing the hidden path
        :param states: list of states of the HMM
        :param emissionMatrix: 2x2 list of lists representing the emission matrix

        Assumptions:
            This class and construction assumes that there are only two states involved in the path, and that the
            emission matrix is a 2x2 matrix.
        """

        self.x = x  # string x containing alphabet sigma
        self.sigma = sigma
        self.hidden = hiddenPath
        self.states = states
        self.emission = emissionMatrix
        self.accessionDict = dict()
        self.prob = 1

    def buildAccessionDict(self):
        """
        Builds accession dictionary by looping through self.sigma and self.states and keys and values
        corresponding to their symbol (char) and integer position (i).

        Called upon by self.findProb()
        :return: updates self.accessionDict
        """
        for i in range(0, len(self.sigma)):
            self.accessionDict[self.sigma[i]] = i

        for i in range(0, len(self.states)):
            self.accessionDict[self.states[i]] = i

    def findProb(self):
        """
        Calls on self.buildAccessionDict then calculates the conditional probability that string x will be emitted by the
        HMM given the hidden path.

        Implementation:

            for each position in hidden path and emitted symbols
                multiply the probability by the emission of x_i according to pi_i

        NOTE: Based off of the key:value structure of self.accessionDict, emission values are accessed according to the
        integer values that correspond to the symbol in each sequence.

        :return: prints self.prob
        """
        # Build accession dictionary
        self.buildAccessionDict()

        # Loop through the hidden path
        for i in range(0, len(self.hidden)):
            # Multiply the probability by the emission value of the ith value of x according to the ith value of pi
            self.prob *= self.emission[self.accessionDict[self.hidden[i]]][self.accessionDict[self.x[i]]]
        print(self.prob)


def main():
    """
    Reads from STDIN and formats the input to be passed as arguments for the HMM class.

    """

    input = []
    lineCount = 0
    emissionMatrix = []
    for line in sys.stdin:
        if lineCount == 0:
            input.append(line.strip())  # string x (input[0])
        elif lineCount == 2:
            input.append(line.strip().split())  # alphabet sigma (input[1])
        elif lineCount == 4:
            input.append(line.strip())  # hidden path (input[2])
        elif lineCount == 6:
            input.append(line.strip().split())  # states (input[3])
        elif lineCount == 9:
            emissionMatrix.append([float(i) for i in line.strip().split()[1:]])  # emissionMatrix row 1
        elif lineCount == 10:
            emissionMatrix.append([float(i) for i in line.strip().split()[1:]])  # emissionMatrix row 2
        lineCount += 1
    
    # Create an HMM object using the generated matrix and the first and second elements of the input list.
    myHMM = HMM(input[0], input[1], input[2], input[3], emissionMatrix)
    # Call on findProb() to print the calculated probability
    myHMM.findProb()


if __name__ == "__main__":
    main()
