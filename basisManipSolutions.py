__author__ = 'Sean Wilner'
'''
basisManipSolutions.py is a collection of funcion calls
'''
import fractions
from collections import  defaultdict
from diophantine import *

# Extend lcm and gcd to cover using lists of integers instead of just pairs of integers
lcm = lambda a: reduce(lambda x,y: abs(x*y)/fractions.gcd(x,y) if x*y > 0 else max(abs(x), abs(y)), a)
gcd = lambda a: reduce(lambda x,y: abs(fractions.gcd(x,y)), a)

# Perform vector addition and scalar-vector multiplication
#   (could also re-write code to use numpy arrays which natively support this)
def list_sum(a,b):
    return map(lambda (x,y): x+y, zip(a,b))
def list_mul(a,c):
    return map(lambda (a):a*c, a)

def getMinimal(a):
    temp_a = map(soft_round, list_mul(a, 1.0/gcd(a)))
    scale_val = 1
    temp_a.reverse()
    for val in temp_a:
        if val != 0:
            scale_val = val
            break
    temp_a = map(soft_round, list_mul(a, 1.0/scale_val))
    temp_a.reverse()
    return temp_a

def findManips(a):
    '''

    :param a:
    :return:
    '''
    # Create a set of all the indices which can be
    allowable = set(range(len(a[0])))
    order = {}
    to_fix = defaultdict(list)
    to_clean = defaultdict(list)
    while True:
        # Create transpose of 'a'
        c = zip(*a)

        for b_indx, c_ in enumerate(c):
            for a_indx in xrange(len(c_)):
                # None of this works if
                if c_[a_indx] != 0:
                    # If a_indx hasn't been covered yet, it needs to be
                    if a_indx in allowable:
                        if not order.has_key(b_indx):
                            order[b_indx] = a_indx
                            allowable.remove(a_indx)
                            continue
                        to_fix[b_indx].append(a_indx)
                    # If a_indx HAS been covered already, then keep track of the fact that b_indx needs clean_up along dimension a_indx
                    else:
                        to_clean[b_indx].append(a_indx)

            else:
                continue
            break
        else:
            break
    return order, to_fix, to_clean


def soft_round(x):
    if abs(x - round(x)) < 0.0001:
        return int(round(x))
    else: return x


def useXtoFixY(X,Y, index):
    scale =  - float(Y[index])/X[index]
    return list_sum(Y, list_mul(X, scale))

def fixSols(a, order, to_fix):
    while len(to_fix)>0:
        for b_indx, a_indices in to_fix.iteritems():
            for a_index in a_indices:
                a[a_index] = useXtoFixY(a[order[b_indx]], a[a_index], b_indx)
                # print "made a change"
            order, to_fix, to_clean = findManips(a)
            break

    return a

def cleanSols(a, order, to_clean):
    while len(to_clean)>0:
        for b_indx, a_indices in to_clean.iteritems():
            if b_indx > max(order):
                continue
            for a_index in a_indices:
                # a[a_index] = useXtoFixY(a[order[b_indx]], a[a_index], b_indx)
                a[a_index] = useXtoFixY(a[order[b_indx]], a[a_index], b_indx)
                # print "made a change"
            order, to_fix, to_clean = findManips(a)
            break
        else:
            break

    return a


def getManipBasis(a,b):
    assert isinstance(a, Matrix)
    assert isinstance(b, Matrix)
    a = [a._mat[i:i+a.shape[1]] for i in xrange(0, len(a._mat), a.shape[1])]
    b = b._mat
    # b = [b_[0] for b_ in b]
    a = [a_ for a_ in a if a_.reverse()==None]
    b.reverse()
    # for a_ in a:
    #     print a_
    # print b

    # c = zip(*a)
    # for c_ in c:
    #     print c_

    order, to_fix, to_clean = findManips(a)
    # print order
    # print to_fix
    # print to_clean

    a = fixSols(a, order, to_fix)



    order, to_fix, to_clean = findManips(a)
    a = cleanSols(a, order, to_clean)

    # print "-------------------------"
    # print order
    # print to_fix
    # print to_clean
    for b_index, a_index in order.iteritems():
        if a[a_index][b_index] < 0:
            a[a_index] = list_mul(a[a_index], -1)
    a = [[soft_round(a__) for a__ in a_] for a_ in a if a_.reverse()==None]
    a = [getMinimal(a_) for a_ in a]
    b.reverse()
    # c = zip(*a)
    # for c_ in c:
    #     print c_
    # # for a_ in a:
    # # #     print a_
    # dio = Diophantine()


    return a, b

if __name__ == "__main__":
    a = Matrix([[-1, 1, 0, 1, -1, 0, 0, 1, -1, 0, -1, 1, 0, 0, 0, 0, 0], [0, 0, 1, -1, 0, 0, -1, 0, 1, 0, 0, 1, -1, 0, 0, 0, 0], [-1, 1, 0, 0, 1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0], [0, -1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1], [-1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1, 0], [-1, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, -1, 1, 0, 0], [0, -1, 0, 1, 1, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, 0, 0], [0, -1, 0, 1, 1, 0, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0], [-1, 1, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0], [0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0], [-1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [-1, 1, 1, 0, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    b = Matrix([[2], [1], [2], [2], [2], [1], [1], [1], [1], [1], [0], [0], [1], [0], [-1], [0], [-1]])

    # a = Matrix([[-1, 1, 0, 1, -1, 0, 0, 1, -1, 0, -1, 1, 0, 0, 0, 0, 0], [0, 0, 1, -1, 0, 0, -1, 0, 1, 0, 0, 1, -1, 0, 0, 0, 0], [-1, 1, 0, 0, 1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0], [0, -1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1], [-1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1, 0], [-1, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, -1, 1, 0, 0], [0, -1, 0, 1, 1, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, 0, 0], [0, -1, 0, 1, 1, 0, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0], [-1, 1, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0], [0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0], [-1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [-1, 1, 1, 0, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    # b = Matrix([[-1], [1], [0], [1], [1], [2], [2], [1], [2], [1], [1], [1], [1], [1], [-1], [0], [-1]])

    # a = Matrix([[-1, 1, 0, 1, -1, 0, 0, 1, -1, 0, -1, 1, 0, 0, 0, 0, 0], [0, 0, 1, -1, 0, 0, -1, 0, 1, 0, 0, 1, -1, 0, 0, 0, 0], [-1, 1, 0, 0, 1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0], [0, -1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1], [-1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1, 0], [-1, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, -1, 1, 0, 0], [0, 0, 0, 0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0], [0, -1, 0, 1, 1, 0, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0], [-1, 1, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0], [0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0], [-1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [-1, 1, 1, 0, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    # b = Matrix([[2], [2], [2], [2], [2], [3], [3], [3], [2], [2], [2], [2], [0], [1], [1], [-1], [-2]])
    # a = Matrix([[-1, 1, 0, 1, -1, 0, 0, 1, -1, 0, -1, 1, 0, 0, 0], [0, 0, 1, -1, 0, 0, -1, 0, 1, 0, 0, 1, -1, 0, 0], [-1, 1, 0, 0, 1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0], [-1, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, -1, 1], [0, -1, 0, 1, 1, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0], [0, -1, 0, 1, 1, 0, 0, 0, -1, -1, 0, 1, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0], [-1, 1, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 1, 0, 0], [0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0], [0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0], [-1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0], [-1, 1, 1, 0, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0]])
    # b = Matrix([[120], [63], [18], [-20], [-48], [-67], [-77], [-78], [-70], [-52], [-26], [9], [56], [109], [172]])

    a,b = getManipBasis(a,b)
    print a
    print b
