'''
Created on May 9, 2011
@author: vmandal
'''
# Code originally from 
# http://www.codesofmylife.com/2011/05/13/smith-waterman-algorithm-for-local-alignment-in-python/
from sys import *
import itertools
import numpy as np
import collections

'''Read the similarity matrix from BLOWSUM50 and return the similarity map of all the combination of the bases'''
def readBLOSUM50(fileName):
    '''The first row of blowsum50 matrix'''
    t = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    #similarityMatrix = np.loadtxt(fileName, delimiter='\t')
    similarityMatrix = [x.strip().split() for x in open(fileName).readlines()]
    similarityMatrix = [[int(col) for col in row[1:]] for row in similarityMatrix[1:]]
    similarityMatrixMap = dict()
    for i in range(len(t)):
        base1 = t[i]
        for j in range(len(t)):
            base2 = t[j]
            similarityMatrixMap[base1 + base2] = similarityMatrix[i][j]
    return similarityMatrixMap

def readDNA():
    similarityMatrixMap = {}
    for a, b in itertools.product('ACGTN', 'ACGTN'):
        if a == b:
            similarityMatrixMap[a+b] = 10
        else:
            similarityMatrixMap[a+b] = 0
    return similarityMatrixMap

Alignment = collections.namedtuple('Alignemnt', 'sequence aligned start end')

def computeFMatrix(seq1, seq2, gap, similarityMatrixMap):
    '''This function creates alignment score matrix
    seq1 : reference sequence
    seq2 : other sequence
    gap : gap penalty '''

    rows = len(seq1) + 1
    cols = len(seq2) + 1
    fMatrix = np.zeros((rows, cols), int)
    pointers = np.zeros((rows, cols), int)

    maxScore = 0
    iOfMax = 0
    jOfMax = 0

    fMatrix[:,0] = 0
    pointers[:,0] = 0
        
    fMatrix[0,:] = 0
    pointers[0,:] = 0

    for i in range(1, rows):
        for j in range(1, cols):
            mtch = fMatrix[i - 1, j - 1] + \
                   similarityMatrixMap[seq1[i-1], seq2[j-1]]
            delete = fMatrix[i-1, j] + gap
            insert = fMatrix[i, j-1] + gap
            fMatrix[i, j] = max(0, mtch, delete, insert)

            if(fMatrix[i, j] == 0):
                pointers[i, j] = -1

            elif(fMatrix[i, j] == delete):
                pointers[i, j] = 1

            elif(fMatrix[i, j] == insert):
                pointers[i, j] = 2

            elif(fMatrix[i, j] == mtch):
                pointers[i, j] = 3

            if fMatrix[i, j] > maxScore :
                iOfMax = i
                jOfMax = j
                maxScore = fMatrix[i, j]

    startOfAlign1 = 1
    startOfAlign2 = 1

    (aligned1, aligned2, startOfAlign1, startOfAlign2) = trackBack(pointers, seq1, seq2, gap, similarityMatrixMap, iOfMax, jOfMax)
    return ((aligned1, startOfAlign1, iOfMax), 
            (aligned2, startOfAlign2, jOfMax))

def formatSWAlignment(a1, a2):
    #    seq1, seq2, aligned1, aligned2, startOfAlign1, startOfAlign2, iOfMax, jOfMax):
    '''Some formatting for displaying alignment'''
    
    numOfSpacesToAdd1 = {True: 0, False: a2.start - a1.start}[a1.start >= a2.start]
    numOfSpacesToAdd2 = {True: 0, False: a1.start - a2.start}[a2.start >= a1.start]
    aligned1 = ' ' * numOfSpacesToAdd1 + a1.sequence[:a1.start].decode('ascii') + a1.aligned + a1.sequence[a1.end:].decode('ascii')
    aligned2 = ' ' * numOfSpacesToAdd2 + a2.sequence[:a2.start].asstring().encode('ascii') + a2.aligned + a2.sequence[a2.end:].asstring().encode('ascii')
    alignPointer = []
    for cSeq1,cSeq2 in zip(aligned1, aligned2):
        alignPointer.append({True: ':', False: ' '}[cSeq1 == cSeq2])

    return (aligned1, aligned2, ''.join(alignPointer))

def scoreSWAlignment(aligned1, aligned2):
    assert len(aligned1.aligned) == len(aligned2.aligned)
    score = 0
    for cSeq1,cSeq2 in zip(aligned1.aligned, aligned2.aligned):
        if cSeq1 == cSeq2:
            score += 1
            
    return score


def trackBack(pointers, seq1, seq2, gap, similarityMap, i, j):
    '''Tracks back to create the aligned sequence pair'''
    alignedSeq1 = []
    alignedSeq2 = []
    
    while pointers[i, j] != -1 and i > 0 and j > 0:
        if pointers[i, j] == 1:
            alignedSeq1.append(seq1[i - 1])
            alignedSeq2.append(ord('-'))
            i = i - 1
        elif pointers[i, j] == 2:
            alignedSeq1.append(ord('-'))
            alignedSeq2.append(seq2[j - 1])
            j = j - 1
        elif pointers[i, j] == 3:
            alignedSeq1.append(seq1[i - 1])
            alignedSeq2.append(seq2[j - 1])
            i = i - 1
            j = j - 1
        else:
            raise ValueError("Lost")

    alignedSeq1.reverse()
    alignedSeq2.reverse()
    return (bytes(alignedSeq1), bytes(alignedSeq2), i, j)

if __name__ == "__main__":

    similarityMatrixMap = readBLOSUM50("../blosum50.txt")
    #similarityMatrixMap = readDNA()
    
    seq1 = b'HEAGAWGHEE'
    seq2 = b'PAWHEAE'
    
    alignment1, alignment2 = computeFMatrix(seq1, seq2, -6, similarityMatrixMap)
    print(alignment1)
    print(alignment2)
    (alignedSeq1, alignedSeq2, alignPointer) = formatSWAlignment(alignment1, alignment2)
    print(alignedSeq1)
    print(alignPointer)
    print(alignedSeq2)
