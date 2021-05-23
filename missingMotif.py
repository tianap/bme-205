#!/usr/bin/env python3

########################################################################
# File:missingMotif.py
#  executable: missingMotif.py
# Purpose: ranks motifs gathered from sequences found in a given fasta file based on how statistically underrepresented
#          the specific motif is
#   stderr: errors and status
#   stdout:
#
# Author: Tiana Pereira
# History:      10/5/2020 Created
#
########################################################################

########################################################################
# CommandLine
########################################################################
import numpy as np
from fastaReader import FastAreader


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

        There are three arguments that are all optional, but are passed in the class Motifs():
            --minMotif : (int) choices = range(3,8), default = 3
            --maxMotif : (int) choices = range(3,8), default = 8
        """

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='python missingMotif.py --minMotif int --maxMotif int --cutoff int < input.fa > output.out'
        )

        self.parser.add_argument('--minMotif', type=int, choices=range(3, 9), action='store',
                                 help='minimum motif size to evaluate (int)', default=3)
        self.parser.add_argument('--maxMotif', type=int, choices=range(3, 9), action='store',
                                 help='maximum motif size to evaluate (int)', default=8)
        self.parser.add_argument('--cutoff', type=int, action='store', default=0,
                                 help='Z-score cutoff for motifs (int)')

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


########################################################################


# Methods not included in any classes
#
#
########################################################################
def reverseComplement(seq):
    """
    Generates the reverse complement of a given DNA sequence. Assumes the sequence passed contains valid bases only
    (A, C, G, T).

    :param seq: A DNA sequence
    :return: the reverse complement of seq
    """
    rc = ''
    for base in seq[::-1]:
        if base == 'A':
            rc += 'T'
        elif base == 'T':
            rc += 'A'
        elif base == 'C':
            rc += 'G'
        elif base == 'G':
            rc += 'C'
    return rc


def isDNA(motif):
    """
    Determines whether the parameter 'motif' contains valid bases by determining the length of the union of the motif
    and canonical bases.
    If the length of this union is 4, then the bases contained in the motif are strictly A, T, C, or G. Returns True.
    If the length is > 4, then motif contains non-canonical bases. Returns False.

    :param motif: a sequence that may or may not contain canonical bases
    :return: False if motif contains non-canonical bases, True otherwise.
    """

    bases = {'A', 'T', 'C', 'G'}
    if len(set(motif).union(bases)) > 4:
        return False
    else:
        return True

class Motifs:
    """
    This is the main class used for finding motifs of various sizes across an entire genome, and then calculating
    their counts, probability, expected values, and Z-scores.

    There are 7 methods contained within the class, and they are as follows:
        def __init__
        def findMotifs
        def probability
        def expectedValue
        def calculateZScore
        def calculateMotifScores
        def printSortedMotifDict

    There are also 8

    """

    def __init__(self, genome, min=3, max=8, cutoff=0):
        '''
        Motifs constructor.

        Takes in a list of genome sequences and outputs all k-mers from the min to max size whose Z-scores are less
        than or equal to the cutoff score.

        Attributes:
            genome: (parameter) a list of genomic sequences from which k-mers will be found.
            minSize: (parameter) the minimum int length k-mer to look for, default=3.
            maxSize: (parameter) the maximum int length k-mer to look for, default=8.
            cutoff: (parameter) the maximum value for Z-scores to be reported, default=0.
            calculatedMotifs: (set) stores all motifs that have been calculated for record-keeping.
            motifCounts: (dict) holds dictionaries that include the k-mers and their counts across the length of all
                genomic seqences.
            motifScores: (list) stores dictionaries that hold all k-mers of size k
                Ex:
                    motifScores = [dictionary for k=minSize, dictionary for k=minSize+1, ..., dictionary for k=maxSize]
                Each dictionary has the following format:
                    {motif: [reverse complement, count, expected value, z-score]}
            genomeSize: (int) the length of all sequences contained in self.genome - maxSize
                NOTE:
                    Assume that the value used to calculate probabilities and Z-scores are the length of all genomic
                    sequences - maxSize k-mer, that way the resulting scores do not over-reported.

        '''
        self.genome = genome  # list of sequences
        self.minSize = min
        self.maxSize = max
        self.cutoff = cutoff
        self.calculatedMotifs = set()
        self.motifCounts = dict()  # holds dictionaries that includes the kmers and their counts, according to their sequence
        self.motifScores = []  # List that holds dictionaries, in which each dictionary holds kmers of size k
        self.genomeSize = 0
        for seq in self.genome:
            self.genomeSize += (len(seq) - self.maxSize)

    def findMotifs(self):
        '''
        This function goes through each sequence included in self.genome and finds all the k-mers and their counts.

        Calls on the following static methods:
            isDNA() : if False, the motif and its reverse complement is discarded
            reverseComplement() : returns the reverse complement of a given motif

        NOTE:
            Functions by looping through each of the genomic sequences passed to this instantiation of the Motif class
            and finding all motifs of size k = minSize and k = maxSize.
            If the motif is new and hasn't been encountered yet, checks if the motif is DNA. If so, adds the motif as
            a key to self.motifCounts and its value to [1]. The reverse complement is also generated and its value is
            set to the same value as self.motifCounts[motif].
            In this way, if either the motif or the reverse complement is encountered, the counts of both are
            incremented, since we want to consider the counts of any motif and its reverse complement to be the same.

        :return: updates self.motifCounts
        '''
        for i in range(0, len(self.genome)):  # Loop through the each sequence
            seq = self.genome[i]
            for k in range(1, self.maxSize + 1):  # For each motif size
                for j in range(0, len(seq) - k + 1):  # Find all motifs of size k
                    motif = seq[j:j + k]
                    if motif in self.motifCounts.keys():  #
                        self.motifCounts[motif][0] += 1
                    else:
                        if isDNA(motif):
                            # if len(motif) == 2
                            self.motifCounts[motif] = [1]  # include the first count
                            rc = reverseComplement(motif)
                            self.motifCounts[rc] = self.motifCounts[motif]

    def probability(self, e):
        """
        Calculates the probability that a specific kmer appears, Pr(K).

        This calculation is based off of the markov approximation, in which Pr(K) is defined as E(K) / N, in which N
        is the length of a genome.

        :param e: (float) an expected value
        :return: (p) representing Pr(K)
        """
        p = e / self.genomeSize
        return p

    def expectedValue(self, motif):
        """
        Calculates the expected value of a given motif, based off of a markov model. Does so by considering the counts
        of the left side of the motif and the right side of the motif, divided by the counts of the middle part of the
        motif.

        :param motif: motif (string) whose expected value we're calculating <
        :return: if the left motif, right motif, and mid motif are all included in self.motifCounts, returns a float
            representing the expected value. Otherwise, returns 0.

        NOTE:
            Assumptions:
            1. if the parts of the motif are included in self.motifCounts, then their bases have already been assessed
            as canonical
            2. returns 0 if the above condition is not met. This causes the p-value to be calculated as 0 in
            self.probability(), which then results in the standard deviation being 0 in self.calculateZScore. This
            method checks to see if the standard deviation is equal to 0 in order to avoid a divide-by-zero error.
        """

        leftMotif = motif[:-1]
        rightMotif = motif[1:]
        midMotif = motif[1:-1]

        if leftMotif in self.motifCounts.keys() and rightMotif in self.motifCounts.keys() and midMotif in self.motifCounts.keys():
            num = self.motifCounts[leftMotif][0] * self.motifCounts[rightMotif][0]
            return num / self.motifCounts[midMotif][0]
        else:
            return 0

    def calculateZScore(self, s, p, n):
        """
        Calculates the Z score for a given motif, given the specific count, its probability, and the size of the genome.
        Does so by first calculating the mean and the standard deviation (sd). If sd is equal to 0,
        automatically returns 1, since any Z scores greater than the cutoff (default=0), will be discarded.

        :param s: (int) specific count to be converted
        :param p: (float) probability of success
        :param n: (int) genomic size
        :return: (float) the calculated Z score
        """
        mean = n * p
        sd = np.sqrt(n * p * (1 - p))
        if sd != 0:
            zScore = (s - mean) / sd
            return zScore
        else:
            return 1  # won't be included since it's positive

    def calculateMotifScores(self):
        """
        Calls on functions expectedValue(), probability(), and calculateZScore() to generate the Z scores for motifs
        and stores these motifs in self.motifScores iff Z-scores are <= cutoff.

        The function begins by looping through motifs of size k, beginning with k=maxSize and ending with k=minSize.
        (printSortedMotifDict() will all of these motifs in order of largest to smallest length k.)

        Then, a dictionary is created and will store the k-mers of size k whose Z-scores are <= the cutoff. This dict
        is then appended to self.motifScores (once the function completes, self.motifScores will contain a dictionary
        for each size k-mer. The format of this dictionary is as follows:

            3-mer{
                motif1 = [reverseComplement1, count, expectedValue, Z-score]
                ...
                motifN = [reverseComplementN, count, expectedValue, Z-score]
            }

        The function continues by looping through motifs throughout the sequence in self.genome.

        When a motif is found, the function checks whether this motif and its reverse complement has already been
        included in the set self.calculatedMotifs. If not, then the function calls on the previous described functions
        to ultimately the Z-score.

        If the Z-score is <= cutoff, it is then stored in the dictionary that holds the other motifs of the same length
        k. This function will also sort the motif and its reverse complement by alphabetical order. Whichever comes
        first wiill become the key in the dictionary, and the first element of its value (a list) is its reverse
        complement.

        The motif and its reverse complement is then added to self.calculatedMotifs to keep track of which motifs have
        been calculated.


        :return: updates the list self.motifScores with dictionaries, where each dictionary contains motifs of the same
        size.

        """

        # Calls findMotifs() to generate all kmers
        self.findMotifs()

        # Loop through k in descending order (max -> min)
        for k in range(self.maxSize, self.minSize - 1, -1):
            kDict = dict()  # Intialize a dictionary that holds kmers of size k
            self.motifScores.append(kDict)  # Add this dictionary to a list containing all kmers
            for i in range(len(self.genome)):  # Loop through the list of sequences
                seq = self.genome[i]
                for j in range(0, len(seq) - k + 1):  # Loop through the motifs of size k
                    motif = seq[j:j + k]
                    rc = reverseComplement(motif)
                    if motif in self.motifCounts.keys():  # Check if this motif has been recorded for this seq
                        if motif not in self.calculatedMotifs and rc not in self.calculatedMotifs:
                            kmers = [motif, rc]  # List of strings motif and its rc, to be sorted alphabetically
                            eValue = self.expectedValue(motif)  # Calculate the E-value according to the sequence
                            p = eValue / self.genomeSize
                            zScore = self.calculateZScore(self.motifCounts[motif][0], p, self.genomeSize)
                            if zScore <= self.cutoff:
                                kmers.sort()
                                # Add these scores to kDict
                                kDict[kmers[0]] = [kmers[1], self.motifCounts[motif][0], eValue, zScore]
                            self.calculatedMotifs.update([motif, rc])

    def printSortedMotifDict(self):
        """
        Prints the motifs in order of descending length of motifs, and then by ascending value of Z-score.

        The function begins by looping through dictionaries stored in self.motifScores. Since each of the dictionaries
        storing the motifs of the same lengths were added in order from max to min, the function will print the motifs
        in descending length of motifs.

        These dictionaries are then sorted by their Z-scores, and each motif and its information is printed.

        """
        # Iterate through dictionary sorted by length of kmers, and by increasing z-score
        print('N = ' + str(self.genomeSize))
        for dict in self.motifScores:  # Loop through each kmer dictionary in self.motifScores
            for i in sorted(dict.items(), key=lambda x: x[1][3]):
                seq = i[0]
                rSeq = i[1][0]
                count = i[1][1]
                E = i[1][2]
                Z = i[1][-1]
                print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}'.format(seq, rSeq, count, E, Z))


########################################################################


# Main
# Here is the main program
#
#
########################################################################

def main(inCL=None):
    """
    Reads in STDIN and calls on class Motifs() using the arguments passed in STDIN.

    NOTE:
        main() relies on FastAreader(), which is imported at the top.
        For the purposes of this assignment, FastAReader.py is a separate file that will be submitted along with this
        file.

    :param inCL: STDIN
    :return:
    """

    seqList = []
    if inCL is None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(inCL)

    myReader = FastAreader()
    for head, seq in myReader.readFasta():
        seqList.append(seq)

    myMotifs = Motifs(seqList, min=myCommandLine.args.minMotif, max=myCommandLine.args.maxMotif,
                      cutoff=myCommandLine.args.cutoff)
    myMotifs.calculateMotifScores()
    myMotifs.printSortedMotifDict()


if __name__ == "__main__":
    main()
