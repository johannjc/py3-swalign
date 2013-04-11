from .cswalign import local_align
from .cswalign import read_matrix


def formatSWAlignment(seq1, seq2, align1, align2):
    '''Some formatting for displaying alignment

    Parameters:
      seq1 - first sequence
      seq2 - second sequence
      aligned1 - first alignment result from local_align
      aligned2 - second alignment result from local_align

      local_align returns (aligned_chars, start, stop)
    '''

    numOfSpacesToAdd1 = max(0, align2[1]-align1[1])
    numOfSpacesToAdd2 = max(0, align1[1]-align2[1])

    aligned1 = b' ' * numOfSpacesToAdd1 + \
              seq1[:align1[1]] + \
              align1[0] + \
              seq1[align1[2]:]
    aligned2 = b' ' * numOfSpacesToAdd2 + \
              seq2[:align2[1]] + \
              align2[0] + \
              seq2[align2[2]:]
    alignPointer = []
    for cSeq1,cSeq2 in zip(aligned1, aligned2):
       alignPointer.append({True: ':', False: ' '}[cSeq1 == cSeq2])

    return (aligned1, ''.join(alignPointer), aligned2)

def scoreSWAlignment(aligned1, aligned2):
    """count how many bases match
    """
    assert len(aligned1.aligned) == len(aligned2.aligned)
    score = 0
    for cSeq1,cSeq2 in zip(aligned1.aligned, aligned2.aligned):
        if cSeq1 == cSeq2:
            score += 1

    return score

