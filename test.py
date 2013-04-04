from swalign import cswalign
from swalign import swalign
import numpy as np

seq1 = 'HEAGAWGEE'
seq2 = 'PAWHEAE'

seq1 = 'SSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENL'

seq2 = 'SCAVPSTDDYAGKYGLQLDFQQNGTAKSVTCTYSPELNKLFCQLAKTCPLLVRVESPPPRGSILRATAVYKKSEHVAEVVKRCPHHERSVEPGEDAAPPSHLMRVEGNLQAYYMEDVNSGRHSVCVPYEGPQVGTECTTVLYNYMCNSSCMGGMNRRPILTIITLETPQGLLLGRRCFEVRVCACPGRDRRTEEDNY'

def load_matrix():
    sim = swalign.readBLOSUM50('blosum50.txt')
    csim = cswalign.read_matrix('blosum50.txt')
    for pair in sim:
        pair_score = sim[pair]
        cscore = csim[ord(pair[0]), ord(pair[1])]
        assert pair_score == cscore
    return sim, csim
    
def test():
    sim, csim = load_matrix()
    #blosum = cswalign.read_matrix('blosum50.txt')
    pa1, pa2 = swalign.computeFMatrix(seq1, seq2, -6, sim)
    print (pa1, pa2)

    ca1, ca2 = cswalign.local_align(seq1, seq2, -6, csim)
    print (ca1, ca2)

if __name__ == '__main__':
    test()
    
