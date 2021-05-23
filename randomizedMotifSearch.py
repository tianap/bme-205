#!/usr/bin/env python3

########################################################################
# File:randomizedMotifSearch.py
#  executable: randomizedMotifSearch.py
# Purpose: given a fasta file, returns the consensus sequence along with its associated profile score (the sum of
# entropies across each position in the final profile.
#   stderr: errors and status
#   stdout:
#
# Author: Tiana Pereira
# History:      10/12/2020 Created
#
########################################################################

########################################################################
# Import random and math modules, as well as FastAreader from fastaReader
########################################################################
import random, math
from fastaReader import FastAreader


########################################################################
# CommandLine
########################################################################
class CommandLine():
    """
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    available within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    """

    def __init__(self, inOpts=None):
        """
        CommandLine constructor.

        Implements a parser to interpret the command line argv string using argparse.

        There are three arguments that are passed to the class Consensus():
            -i iterations (int)
            -k motif length (int)
            -p pseudocount (float)
        """

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='python randomizedMotifSearch.py -i int --maxMotif int --cutoff int < input.fa > output.out'
        )
        # Be sure to go over the argument information again
        self.parser.add_argument('-i', type=int, action='store',
                                 help='number of iterations (int)')
        self.parser.add_argument('-k', type=int, action='store',
                                 help='motif length (int)')
        self.parser.add_argument('-p', type=int, action='store',
                                 help='pseudocount (float)')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class Usage(Exception):
    """
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    """

    def __init__(self, msg):
        self.msg = msg


class Consensus:
    """
    Contains all attributes and methods related to finding the Consensus sequence of a given set of fasta sequences.

    Methods:
        randomMotifs()
        count()
        findProfile()
        motifs()
        entropyScore()

    """

    def __init__(self, i, k, p, Dna):
        """
        Constructor for Consensus class
        :param i: the number of iterations to run Randomized Motif Search (int)
        :param k: motif length (int)
        :param p: psuedocount (float)
        :param Dna: list of sequences taken from fasta file (list)

        There are also 2 additional attributes:
        1. self.bestScore - used to store the best entropy score for a set of motifs; intentionally set to a high,
        impossible number since we want to look for the smallest entropy score.
        2. self.bestMotifs - used to store the motifs whose entropy scores are the best
        """

        self.i = i
        self.k = k
        self.p = p
        self.Dna = Dna
        self.bestScore = 1000000
        self.bestMotifs = []

    def randomMotifs(self):
        """
        Given a list of sequences, randomly generates motifs of length k.
        (Implements random function randint())

        :return: random list of motifs
        """
        randomMotifs = []
        for seq in self.Dna:  # Loop through each sequence passed by user
            start = random.randint(0, len(seq) - self.k)  # generate random start
            randomMotifs.append(seq[start:start + self.k])
        return randomMotifs

    def count(self, motifs):
        """
        Given a list of motifs, calculates the counts of nucleotides in positions across the list of motifs.
        All counts are added along with the pseudocount, self.p, which is given by the user via the command line.
        This method is called on by findProfile(), and the returned dictionary of counts is used to find the
        distribution profile for the given set of motifs.

        :param motifs: a list of randomly selected motifs (generated by randomMotifs())
        :return: a dictionary storing the counts of all nucleotides in their respective position.

        NOTE: this dictionary is meant to simulate a matrix where each row (i) is a nucleotide and each column (j) is a
        position. Thus, countDict(i,j) = the number of nucleotides i at position j.

        Example:
                pos_1 pos_2 ... pos_k
            A:    3     2   ...   2
            T:    1     2   ...   2
            C:    2     2   ...   2
            G:    2     2   ...   3

        """

        countDict = {  # keeps the counts of A,C,T,G across all sequences in self.Dna
            'A': [],
            'T': [],
            'C': [],
            'G': []
        }
        increment = self.p
        for i in range(0, self.k):  # goes through each column in the sequences
            colI = [col[i] for col in motifs]
            countDict['A'].append(colI.count('A') + increment)
            countDict['T'].append(colI.count('T') + increment)
            countDict['G'].append(colI.count('G') + increment)
            countDict['C'].append(colI.count('C') + increment)

        return countDict

    def findProfile(self, motifs):
        """
        Calculates the distribution profile for a set of motifs.
        Calls on counts() to generate the counts, then uses these counts to calculate the matrix profile.

        :param motifs:
        :return: a dictionary storing the probability distribution for each base at each position

        NOTE: much like the counts dictionary, the dictionary returned by findProfile() simulates a matrix.
        Example:

                Example:
                pos_1 pos_2 ... pos_k
            A:   3/8   2/8  ...  2/8
            T:   1/8   2/8  ...  2/8
            C:   2/8   2/8  ...  2/8
            G:   2/8   2/8  ...  3/8
        """
        counts = self.count(motifs)  # call on self.counts() with the same motifs
        profileDict = {
            'A': [],
            'T': [],
            'C': [],
            'G': []
        }
        total = sum([seq[0] for seq in counts.values()])
        for key in counts.keys():
            for i in range(0, len(counts[key])):
                profileDict[key].append(counts[key][i] / total)

        return profileDict

    def motifs(self, profile):
        """
        Calculates the collection of k-mers formed by the profile-most probable k-mers

        :param profile: profile distribution for a given set of motifs (calculated in profile())
        :return: list of most-probable motifs
        """
        # print('Motifs fcn has been called')
        newMotifs = [None] * len(self.Dna)  # initialize list with one element for each string in self.Dna
        for i in range(0, len(self.Dna)):
            seq = self.Dna[i]
            bestProb = 0  # stores the best probability for each motif
            for j in range(0, len(seq) - self.k + 1):
                storedMotif = seq[j:j + self.k]  # look at this motif
                prob = 1  # reset the probability for each motif
                for pos in range(0, self.k):  # for each position in this k-mer
                    base = storedMotif[pos]  # look at one base at a time
                    prob = prob * profile[base][pos]  # freq of this base at position k

                if prob > bestProb:  # look for the highest probability
                    newMotifs[i] = storedMotif  # store this motif
                    bestProb = prob  # reset bestProb with prob

        return newMotifs

    def entropyScore(self, motif):
        """
        Calculates the entropy scores for base frequencies across a given set of motifs. Uses the entropy equation
        found in the textbook.

        Calls on findProfile() to determine the probability distribution for a given set of motifs.

        :param: a list of motifs
        :return: the entropy score for those motifs
        """
        profiles = self.findProfile(motif)  # Calculate the frequency profile for the given motifs
        motifScore = 0
        for base in profiles.keys():  # for all bases
            for count in profiles[base]:  # for all positions
                if count != 0:
                    motifScore += (count * math.log2(count))
        motifScore = motifScore * -1
        return motifScore

    def randomizedMotifSearch(self):
        """
        Calls on a variety of functions to determine the best scoring list of motifs.
        This method is based on the psuedocode from the textbook and functions as follows:

        Called on by findConsensus(), and calls on randomMotifs(), findProfile(), and entropyScore(). This is so that
        only findConsensus() needs to be called in main().

        :return: a list of the best scoring motifs
        """

        for iteration in range(0, self.i):
            # Randomly find motifs
            motifs = self.randomMotifs()
            bestMotifs = motifs
            count = 0
            while True:
                count += 1
                # Find the profiles
                profile = self.findProfile(motifs)

                # Calculate the new motifs
                newMotifs = self.motifs(profile)  # reset motifs based on profile

                # Calculate the entropies based on these motifs
                newMotifsScore = self.entropyScore(newMotifs)  # calculate entropy for new motifs
                bestMotifsScore = self.entropyScore(bestMotifs)  # calculate entropy for old motifs

                # Compare the scores
                if newMotifsScore < bestMotifsScore:
                    bestMotifs = newMotifs  # reset bestMotifs to be the new motifs found

                    # Compare this score to self.bestScore
                    if newMotifsScore < self.bestScore:  # if these motifs and their score are the best through all runs
                        self.bestScore = newMotifsScore
                        self.bestMotifs = bestMotifs

                else:  # This is the best score and best motifs for this iteration
                    break
            print(count)
        return self.bestMotifs

    def findConsensus(self):
        """
        Finds the consensus motif based on the best scoring motifs returned by randomizedMotifSearch()
        Functions by going through each position i for all motifs and adds the base that appears the most to the
        consensus sequence.

        :return: consensus sequence
        """
        self.randomizedMotifSearch()
        consensus = ''

        for i in range(0, self.k):
            posBase = [base[i] for base in self.bestMotifs]
            consensus += max(posBase, key=posBase.count)
        return consensus


def main(inCL=None):
    """
    Main takes in the command line from the user and assigns it to the CommandLine() class.
    It then uses FastAreader() to read in the fasta files, whose sequences are then used to find the consensus motif
    across all sequences within the class Consensus().

    :return: prints the consensus sequence and the final entropy score
    """
    seqList = []
    if inCL is None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(inCL)

    myReader = FastAreader()
    for head, seq in myReader.readFasta():
        seqList.append(seq)

    myConsensus = Consensus(i=myCommandLine.args.i, k=myCommandLine.args.k, p=myCommandLine.args.p, Dna=seqList)
    print('Consensus sequence: ' + myConsensus.findConsensus())
    print('Final Score: ' + str(myConsensus.bestScore))


if __name__ == "__main__":
    main()
