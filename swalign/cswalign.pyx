import numpy as np
cimport numpy as np


cimport cython
import sys
import os.path

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char *s)

cdef extern from "Python.h":
    ctypedef void PyObject
    PyObject *PyBytes_FromStringAndSize(char *, Py_ssize_t)
    int _PyBytes_Resize(PyObject **, Py_ssize_t)
    char * PyBytes_AS_STRING(PyObject *)

ctypedef np.int_t DTYPE_INT
ctypedef np.uint_t DTYPE_UINT
ctypedef np.int8_t DTYPE_BOOL

cdef size_t UP = 1, LEFT = 2, DIAG = 3, NONE = 4

cdef inline int imax2(int a, int b):
    if a >= b: return a
    return b

def read_matrix(path, dict cache={}):
    """
    so here, we read a matrix in the NCBI format and put
    it into a numpy array. so the score for a 'C' changing
    to an 'A' is stored in the matrix as:
        mat[ord('C'), ord('A')] = score
    as such, it's a direct array lookup from each pair in the alignment
    to a score. this makes if very fast. the cost is in terms of space.
    though it's usually less than 100*100.
    """
    if path in cache: return cache[path]
    cdef np.ndarray[DTYPE_INT, ndim=2] a
    cdef size_t ai = 0, i
    cdef int v, mat_size

    fh = open(path)
    headers = None
    while headers is None:
        line = fh.readline().strip()
        if line[0] == '#': continue
        headers = [ord(x) for x in line.split(' ') if x]
    mat_size = max(headers) + 1

    a = np.zeros((mat_size, mat_size), dtype=np.int)

    line = fh.readline()
    while line:
        line_vals = [int(x) for x in line[:-1].split(' ')[1:] if x]
        for ohidx, val in zip(headers, line_vals):
            a[headers[ai], ohidx] = val
        ai += 1
        line = fh.readline()
        
    cache[path] = a
    return a

@cython.boundscheck(False)
@cython.nonecheck(False)
def local_align(object _seqj, object _seqi, gap=-1, object matrix=None):
    """Perform a local sequence alignment (smith-waterson)
    """
    if matrix is None:
        raise ValueError("Matrix is actually required")
    
    return build_alignment_matrix(_seqj, _seqi, gap, matrix) 
    #return ((ai, starti, iOfMax), (aj, startj, jOfMax))
    
cdef build_alignment_matrix(object _seqi, object _seqj, gap=-1, object matrix=None):
    cdef bint flip = 0
    
    cdef char *seqi = _seqi
    cdef char *seqj = _seqj

    cdef size_t max_i = strlen(seqi) + 1
    cdef size_t max_j = strlen(seqj) + 1

    cdef np.ndarray[DTYPE_INT, ndim=2] fMatrix = \
      np.empty((max_i+1, max_j+1), dtype=np.int)
    cdef np.ndarray[DTYPE_INT, ndim=2] pointers = \
      np.empty((max_i+1, max_j+1), dtype=np.int)
    cdef np.ndarray[DTYPE_INT, ndim=2] similarity

    #similarity = read_matrix(matrix)
    similarity = matrix

    cdef int match, delete, insert, maxScore=0
    cdef int iOfMax=0, jOfMax=0
    cdef int i, j
    
    fMatrix[:,0] = 0
    pointers[:,0] = 0
    
    fMatrix[0,:] = 0
    pointers[0,:] = 0
    
    for i in range(1, max_i):
        for j in range(1, max_j):
            match = fMatrix[i-1][j-1]+similarity[seqi[i-1], seqj[j-1]]
            delete = fMatrix[i-1][j] + gap
            insert = fMatrix[i][j-1] + gap
            fMatrix[i][j] = max(0, match, delete, insert)

            if fMatrix[i][j] == 0:
                pointers[i][j] = -1
            elif fMatrix[i][j] == delete:
                pointers[i][j] = 1
            elif fMatrix[i][j] == insert:
                pointers[i][j] = 2
            elif fMatrix[i][j] == match:
                pointers[i][j] = 3

            if fMatrix[i][j] > maxScore:
                iOfMax = i
                jOfMax = j
                maxScore = fMatrix[i][j]

    ai, aj, starti, startj = trackBack(pointers, seqi, seqj, gap, similarity, iOfMax, jOfMax)
    return ((''.join(ai), starti, iOfMax), (''.join(aj), startj, jOfMax))
                

cdef trackBack(np.ndarray[DTYPE_INT, ndim=2] pointers, 
               char *seqi, 
               char *seqj, 
               int gap, 
               np.ndarray[DTYPE_INT, ndim=2] similarity, 
               int i, 
               int j):
    alignedi = []
    alignedj = []

    while pointers[i][j] != -1 and i > 0 and j > 0:
        if pointers[i][j] == 1:
            alignedi.insert(0, chr(seqi[i-1]))
            alignedj.insert(0, '-')
            i -= 1
        elif pointers[i][j] == 2:
            alignedi.insert(0, '-')
            alignedj.insert(0, chr(seqj[j-1]))
            j -= 1
        elif pointers[i][j] == 3:
            alignedi.insert(0, chr(seqi[i-1]))
            alignedj.insert(0, chr(seqj[j-1]))
            i -= 1
            j -= 1
        else:
            raise ValueError("Invalid pointer value")

    return (alignedi, alignedj, i, j)
            
