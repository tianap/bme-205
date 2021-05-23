#!/usr/bin/env python

########################################################################
# File: problem19.py
#  executable: py problem19.py < input.txt > output.out
# Purpose: Given a hidden path, states, and transition matrix of an HMM, returns the probability of the path.
#
# Classes:
#   HMM
#
# Author: Tiana Pereira
# History:      11/15/2020 Created
#
########################################################################
import sys


class HMM:
    """
    Class representing an HMM and its hidden path, states, and transition matrix.

    Attributes:
        hidden : str
            hidden path of the HMM
        states : list
            states of the HMM
        transition : list of list, representing 2x2 matrix
            transition matrix of the HMM
        prob : float
            probability of the hidden path

    Methods:
        findProb():
            Calculates the probability of the hidden path occurring, prints said probability.

    """

    def __init__(self, hiddenPath, states, transitionMatrix):
        """
        Constructs necessary attributes for the HMM class.

        Parameters:
        :param hiddenPath: string representing the hidden path
        :param states: list of states of the HMM
        :param transitionMatrix: 2x2 list of lists representing the transition matrix

        Assumptions:
            This class and construction assumes that there are only two states involved in the path, and that the
            transition matrix is a 2x2 matrix.
        """

        self.hidden = hiddenPath
        self.states = states
        self.transition = transitionMatrix
        self.prob = 1 / len(states)  # assuming the initial probabilities are equal, pr(A) = pr(B)

    def findProb(self):
        """
        Calculates the probability of the hidden path occurring, given the transition matrix of the states.

        NOTE: Assumes the initial probabilities are equal (i.e. Pr(A) = Pr(B)). Therefore, the starting probability
        before calculation is set to 0.5.

        :return: prints self.prob
        """
        # Loop through the hidden path
        for i in range(0, len(self.hidden) - 1):
            # Look at subsequences of len(2)
            seq = self.hidden[i:i + 2]
            matrixPos = []

            # Find the matrix position given the sequence of len(2)
            # NOTE: The integers 0 and 1 correspond to the positions in self.transition
            # ex. [0,0] corresponds to 'AA' with Pr(AA) = self.transition[0][0]
            for base in seq:
                if base == 'A':
                    matrixPos.append(0)
                elif base == 'B':
                    matrixPos.append(1)
            # Multiply the probabilities
            self.prob *= self.transition[matrixPos[0]][matrixPos[1]]
        print(self.prob)


def main():
    """
    Reads from STDIN and formats the input to be passed as arguments for the HMM class.

    """
    input = []
    lineCount = 0
    matrix = []
    for line in sys.stdin:
        if lineCount == 0:
            input.append(line.strip())
        elif lineCount == 2:
            input.append(line.strip().split())
        elif lineCount == 5:
            matrix.append([float(i) for i in line.strip().split()[1:]])
        elif lineCount == 6:
            matrix.append([float(i) for i in line.strip().split()[1:]])
        lineCount += 1

    # Create an HMM object using the generated matrix and the first and second elements of the input list.
    myHMM = HMM(input[0], input[1], matrix)
    # Call on findProb() to print the calculated probability
    myHMM.findProb()


if __name__ == "__main__":
    main()
