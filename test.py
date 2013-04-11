#from swalign import cswalign
from swalign import swalign
from swalign import cswalign
from swalign import formatSWAlignment, scoreSWAlignment

import numpy as np
import math
import string

import pstats, cProfile

import timeit

seq1 = b'SSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENL'

seq2 = b'SCAVPSTDDYAGKYGLQLDFQQNGTAKSVTCTYSPELNKLFCQLAKTCPLLVRVESPPPRGSILRATAVYKKSEHVAEVVKRCPHHERSVEPGEDAAPPSHLMRVEGNLQAYYMEDVNSGRHSVCVPYEGPQVGTECTTVLYNYMCNSSCMGGMNRRPILTIITLETPQGLLLGRRCFEVRVCACPGRDRRTEEDNY'

result1 = b'SSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSD-SDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEEN'
result2 = b'SCAVPSTDDYAGKYGLQLDFQQNGTAKSVTCTYSPELNKLFCQLAKTCPLLVRVESPPPRGSILRATAVYKKSEHVAEVVKRCPHHERSVEPGEDAAPPSHLMRVEGNLQAYYMEDVNSGRHSVCVPYEGPQVGTECTTVLYNYMCNSSCMGGMNRRPILTIITLETPQGLLLGRRCFEVRVCACPGRDRRTEEDN'

CSIM = cswalign.read_matrix('blosum50.txt')
CDNA = cswalign.read_matrix('dna.txt')

def test():
    pa1, pa2 = swalign.local_align(seq1, seq2, -10, CSIM)
    ca1, ca2 = cswalign.local_align(seq1, seq2, -10, CSIM)

    # is it right?
    assert ca1[0] == result1
    assert ca2[0] == result2
    # compare aligned
    assert pa1[0] == ca1[0]
    assert pa2[0] == ca2[0]
    # compare start
    assert pa1[1] == ca1[1]
    assert pa2[1] == ca2[1]
    # compare end
    assert pa1[2] == ca1[2]
    assert pa2[2] == ca2[2]

def test_formatting():
    sseq1 = b'HEAGAWGEE'
    sseq2 = b'PAWHEAE'

    ca1, ca2 = cswalign.local_align(sseq1, sseq2, -6, CSIM)
    formatted = formatSWAlignment(sseq1, sseq2, ca1, ca2)

    assert formatted[0] == sseq1
    assert formatted[1] == '    :: : '
    assert formatted[2] == b'   ' + sseq2
    
    print('--- prot ---')
    print(formatted[0].decode('ascii'))
    print(formatted[1])
    print(formatted[2].decode('ascii'))
    
def test_dna():
    seq1 = b'GGTATACC'
    seq2 = b'TATANC'
    ca1, ca2 = cswalign.local_align(seq1, seq2, -6, CDNA)
    pa1, pa2 = swalign.local_align(seq1, seq2, -6, CDNA)

    assert ca1[0] == b'TATACC'
    assert ca2[0] == b'TATANC'
    assert pa1[0] == ca1[0]
    assert pa2[0] == pa2[0]

    assert ca1[1] == 2
    assert ca2[1] == 0

    assert ca1[2] == 8
    assert ca2[2] == 6

    formatted = formatSWAlignment(seq1, seq2, ca1, ca2)
    print('--- dna ---')
    print(formatted[0].decode('ascii'))
    print(formatted[1])
    print(formatted[2].decode('ascii'))
    
def time():
    setup = 'from __main__ import  swalign, cswalign, seq1, seq2, CSIM'
    cmds = { 'py': 'swalign.local_align(seq1, seq2, -10, CSIM)',
             'cy': 'cswalign.local_align(seq1, seq2, -10, CSIM)',
             }
    seconds = 5
    for name in cmds:
        t0 = timeit.timeit(cmds[name], setup,  number=1)
        count = int(math.ceil(seconds / (max(t0, .0000001) * 3)))
        print(name, t0, count)
        trial = timeit.repeat(cmds[name], setup, repeat=3, number=count)
        normalized_trial = [ t / count for t in trial ]
        print (name, normalized_trial)

        cProfile.runctx(cmds[name], globals(), locals(), name+'.prof')

    for name in cmds:
        print('-----', name, '-----')
        s = pstats.Stats(name + '.prof')
        s.strip_dirs().sort_stats('time').print_stats()

if __name__ == '__main__':
    test()
    test_formatting()
    test_dna()
    time()
