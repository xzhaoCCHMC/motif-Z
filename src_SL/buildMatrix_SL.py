# Module to build up frequency tables

import utility
import numpy as np
import time
#from collections import Counter

# Constants
AA_LIST = list('ARNDCEQGHILKMFPSTWYV')

# == Class Summary
#   Make superclass matrix for foreground, background, binomial matrices  
#   example: [[              ],
#             [              ],
#             [              ],
#             [              ]]             

# Backport Python 2.6 Counter class

from operator import itemgetter
from heapq import nlargest
from itertools import repeat, ifilter

class Counter(dict):
    '''Dict subclass for counting hashable objects.  Sometimes called a bag
    or multiset.  Elements are stored as dictionary keys and their counts
    are stored as dictionary values.

    >>> Counter('zyzygy')
    Counter({'y': 3, 'z': 2, 'g': 1})

    '''

    def __init__(self, iterable=None, **kwds):
        '''Create a new, empty Counter object.  And if given, count elements
        from an input iterable.  Or, initialize the count from another mapping
        of elements to their counts.

        >>> c = Counter()                           # a new, empty counter
        >>> c = Counter('gallahad')                 # a new counter from an iterable
        >>> c = Counter({'a': 4, 'b': 2})           # a new counter from a mapping
        >>> c = Counter(a=4, b=2)                   # a new counter from keyword args

        '''        
        self.update(iterable, **kwds)

    def __missing__(self, key):
        return 0

    def most_common(self, n=None):
        '''List the n most common elements and their counts from the most
        common to the least.  If n is None, then list all element counts.

        >>> Counter('abracadabra').most_common(3)
        [('a', 5), ('r', 2), ('b', 2)]

        '''        
        if n is None:
            return sorted(self.iteritems(), key=itemgetter(1), reverse=True)
        return nlargest(n, self.iteritems(), key=itemgetter(1))

    def elements(self):
        '''Iterator over elements repeating each as many times as its count.

        >>> c = Counter('ABCABC')
        >>> sorted(c.elements())
        ['A', 'A', 'B', 'B', 'C', 'C']

        If an element's count has been set to zero or is a negative number,
        elements() will ignore it.

        '''
        for elem, count in self.iteritems():
            for _ in repeat(None, count):
                yield elem

    # Override dict methods where the meaning changes for Counter objects.

    @classmethod
    def fromkeys(cls, iterable, v=None):
        raise NotImplementedError(
            'Counter.fromkeys() is undefined.  Use Counter(iterable) instead.')

    def update(self, iterable=None, **kwds):
        '''Like dict.update() but add counts instead of replacing them.

        Source can be an iterable, a dictionary, or another Counter instance.

        >>> c = Counter('which')
        >>> c.update('witch')           # add elements from another iterable
        >>> d = Counter('watch')
        >>> c.update(d)                 # add elements from another counter
        >>> c['h']                      # four 'h' in which, witch, and watch
        4

        '''        
        if iterable is not None:
            if hasattr(iterable, 'iteritems'):
                if self:
                    self_get = self.get
                    for elem, count in iterable.iteritems():
                        self[elem] = self_get(elem, 0) + count
                else:
                    dict.update(self, iterable) # fast path when counter is empty
            else:
                self_get = self.get
                for elem in iterable:
                    self[elem] = self_get(elem, 0) + 1
        if kwds:
            self.update(kwds)

    def copy(self):
        'Like dict.copy() but returns a Counter instance instead of a dict.'
        return Counter(self)

    def __delitem__(self, elem):
        'Like dict.__delitem__() but does not raise KeyError for missing values.'
        if elem in self:
            dict.__delitem__(self, elem)

    def __repr__(self):
        if not self:
            return '%s()' % self.__class__.__name__
        items = ', '.join(map('%r: %r'.__mod__, self.most_common()))
        return '%s({%s})' % (self.__class__.__name__, items)

    # Multiset-style mathematical operations discussed in:
    #       Knuth TAOCP Volume II section 4.6.3 exercise 19
    #       and at http://en.wikipedia.org/wiki/Multiset
    #
    # Outputs guaranteed to only include positive counts.
    #
    # To strip negative and zero counts, add-in an empty counter:
    #       c += Counter()

    def __add__(self, other):
        '''Add counts from two counters.

        >>> Counter('abbb') + Counter('bcc')
        Counter({'b': 4, 'c': 2, 'a': 1})


        '''
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] + other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __sub__(self, other):
        ''' Subtract count, but keep only results with positive counts.

        >>> Counter('abbbc') - Counter('bccd')
        Counter({'b': 2, 'a': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] - other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __or__(self, other):
        '''Union is the maximum of value in either of the input counters.

        >>> Counter('abbb') | Counter('bcc')
        Counter({'b': 3, 'c': 2, 'a': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        _max = max
        result = Counter()
        for elem in set(self) | set(other):
            newcount = _max(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result

    def __and__(self, other):
        ''' Intersection is the minimum of corresponding counts.

        >>> Counter('abbb') & Counter('bcc')
        Counter({'b': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        _min = min
        result = Counter()
        if len(self) < len(other):
            self, other = other, self
        for elem in ifilter(self.__contains__, other):
            newcount = _min(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result


class Matrix:
    '''
    Generate the basic matrix format for the speicific cases
    '''
    
    # initialize the matrix in the data structure of list of list
    # -- iPairs is the parameter to record all identified rp_pairs in a list
    def __init__(self, PEP_LENGTH, seqset):    
        self.PEP_LENGTH = PEP_LENGTH
        self.seqset = seqset
        self.matrix = np.zeros((len(AA_LIST), self.PEP_LENGTH))
        #Change the sequence set into a numpy string array
        self.pepList_str = list(seqset)
        self.total_pep = len(seqset)
        if len(self.seqset) == 0:
            print "Warning: Empty dataset to build matrix!"
        self.charMatrix = self.charMatrix()
        
        
    # return the matrix to caller
    def getMatrix(self):
        return self.pepArray_str
    
    # generate forward aligned matrix
    def charMatrix(self):
        #Change the numpy string array to a Char array with each peptide string as a row
        pepList_char = [] 
        for str_pep in self.pepList_str:
            pepList_char.append(list(str_pep))    
        pepArray_char = np.array(pepList_char)
        return pepArray_char


    # generate count matrix
    def aaCountMatrix(self):
        
        start = time.time()
        
        count_matrix_dict = dict()
        for aa in AA_LIST:
            count_matrix_dict[aa] = list()
            
        for col in range(self.PEP_LENGTH):
            count_curr_col = Counter(self.charMatrix[:, col])
            if count_curr_col.has_key(' '):
                del count_curr_col[' ']
            if count_curr_col.has_key('X'):
                del count_curr_col['X']
            for aa in AA_LIST:
                if count_curr_col.has_key(aa):
                    count_aa = count_curr_col[aa]
                else:
                    count_aa = 0
                count_matrix_dict[aa].append(count_aa)           
            
        end = time.time()
        #print "build count matrix used ", (end-start)/60, "min"
        return count_matrix_dict

 
    # generate count matrix
    def aaFreqMatrix(self):
        
        start = time.time()
        
        freq_matrix_dict = dict()
        for aa in AA_LIST:
            freq_matrix_dict[aa] = list()
            
        for col in range(self.PEP_LENGTH):
            count_curr_col = Counter(self.charMatrix[:, col])
            if count_curr_col.has_key(' '):
                del count_curr_col[' ']
            if count_curr_col.has_key('X'):
                del count_curr_col['X']
            
            for aa in AA_LIST:
                if count_curr_col.has_key(aa):
                    freq = count_curr_col[aa]/float(self.total_pep)
                    freq_matrix_dict[aa].append(freq)
                else:
                    freq_matrix_dict[aa].append(0)           
            
        end = time.time()
        #print "build frequency matrix used ", (end-start)/60, "min"
        return freq_matrix_dict
        
    # Define the binomial probability matrix
    def logbpMatrix (self, aaCountMatrix, bgMatrixWithPC, log_factorial_map):
        
        # initialize the binomial matrix as a dictionary 
        logbp_matrix = dict()
        for aa in AA_LIST:
            logbp_matrix[aa] = list()
            
        for aa in AA_LIST:
            for pos in range(self.PEP_LENGTH):
                k = aaCountMatrix[aa][pos]
                n = self.total_pep
                p = bgMatrixWithPC[aa][pos]
                # if foreground completely deplete one rp pair, set the binomial probability to 0 (its log form)
                if k == 0 or p == 0:
                    logbp_matrix[aa].append(0)
                # Fixed rp pairs from previous round would not count
                elif k == n and p == 1:
                    logbp_matrix[aa].append(0)
                else:
                    logbp_matrix[aa].append(utility.binomial_log_prob(n, k, p, log_factorial_map))
                    
        return logbp_matrix              
        
# == Class Summary
#   Make superclass matrix for foreground, background, binomial matrices  
#   example: [[              ],
#             [              ],
#             [              ],
#             [              ]]  

class bgMatrixWithPseudoCount(Matrix):
    # initialize the matrix with a data structure of dict of list
    def __init__(self, PEP_LENGTH, seqset, pseudo_count_matrix):
        Matrix.__init__(self, PEP_LENGTH, seqset)
        self.pc_matrix = pseudo_count_matrix
        self.bg_matrix = Matrix.aaFreqMatrix(self)     
        for aa in AA_LIST:
            for pos in range(self.PEP_LENGTH):
                if self.bg_matrix[aa][pos] == 0:
                #print("Warning: Zero probability in background, pseudo counts added")
#                print("The residue/position pair: " + aa + "/" + str(pos) + "from 0 change to "),
#                print(self.pc_matrix[aa][pos])
                # Add pseudo count to the background matrix 
                    self.bg_matrix[aa][pos] = self.pc_matrix[aa][pos]
                 
    def getMatrixWithPC(self):
        return self.bg_matrix
