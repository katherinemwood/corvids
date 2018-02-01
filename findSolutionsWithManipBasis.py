__author__ = 'Sean Wilner'
from basisManipSolutions import *
import multiprocessing as mp

def forced_neg_removal(a,b):

    force_indices = {}
    c = map(list, zip(*a))
    for indx, row in enumerate(c):
        # check if there's exactly one non-zero, and if the value of the non-zero is '1'
        if sum([1 for r in row if r!=0 ]) == 1 == sum(row):
            force_indices[indx] = row.index(1)
    for indx in  force_indices:
        b = list_sum(b, list_mul(a[force_indices[indx]], -1*b[indx]))
    return b

def dead_solution(a, b):

    for indx, val in enumerate(b):
        if val >= 0:
            continue
        if any([a_[indx]>0 for a_ in a]):
            continue
        return True
    return False

def recurse_over_solution_path(a, b, covered = set()):
    sols = []
    if len(set([indx for indx,b_ in enumerate(b) if b_<0])) == 0:
        sols.append(b)

    for a_indx, a_ in enumerate(a):
        if a_indx in covered:
            continue
        b_temp = list_sum(a_, b)
        # temp_neg =set([indx for indx,b_ in enumerate(b_temp) if b_<0])
        # if len(temp_neg) == 0:
        #     sols.append(b_temp)
        if not dead_solution(a, b):
            sols.extend(recurse_over_solution_path(a, b_temp, covered=set(covered)))
        covered.add(a_indx)

    return sols


def recurse_find_first(a, b, covered=set()):

    if len(set([indx for indx,b_ in enumerate(b) if b_<0])) == 0:
        fl_sol = [soft_round(float(v)) for v in b]
        if all([isinstance(x, int) for x in fl_sol]):
            return b

    for a_indx, a_ in enumerate(a):
        if a_indx in covered:
            continue
        b_temp = list_sum(a_, b)
        temp_neg =set([indx for indx,b_ in enumerate(b_temp) if b_<0])
        if len(temp_neg) == 0:
            fl_sol = [soft_round(float(v)) for v in b_temp]
            if all([isinstance(x, int) for x in fl_sol]):
                return b_temp
        if not dead_solution(a, b_temp):
            sol = recurse_find_first(a, b_temp, covered=set(covered))
            if sol:
                fl_sol = [soft_round(float(v)) for v in sol]
                if all([isinstance(x, int) for x in fl_sol]):
                    return sol
        covered.add(a_indx)
    return None

def multiprocess_recurse_func(X):
    a, b, depth, covered = X
    rets = []
    for a_indx, a_ in enumerate(a):
        if a_indx in covered:
            continue
        b_temp = list_sum(a_, b)
        temp_neg =set([indx for indx,b_ in enumerate(b_temp) if b_<0])
        if len(temp_neg) == 0:
            rets.append(b_temp)
        if not dead_solution(a, b):
            if depth == 0:
                rets.append((a, b_temp, 5, set(covered)))
            else:
                rets.extend(multiprocess_recurse_func((a , b_temp, depth-1, set(covered))))
        covered.add(a_indx)
    return  rets

def multiprocess_recurse_find_first(a, b, covered=set()):

    if isinstance(b, list):
            for b_inst in b:
                if len(set([indx for indx,b_ in enumerate(b_inst) if b_<0])) == 0:
                    fl_sol = [soft_round(float(v)) for v in b_inst]
                    if all([isinstance(x, int) for x in fl_sol]):
                        return b_inst
    else:
        if len(set([indx for indx,b_ in enumerate(b) if b_<0])) == 0:
            fl_sol = [soft_round(float(v)) for v in b]
            if all([isinstance(x, int) for x in fl_sol]):
                return b
    pool = mp.Pool()
    runnables = []
    if isinstance(a, list) and isinstance(b, list):
        for a_inst, b_inst in zip(a,b):
            runnables.append((a_inst, b_inst, 5, set()))
    else:
        runnables = [(a, b, 5, set())]
    # print "Runnables: ", runnables
    while len(runnables) > 0:
        return_lists = pool.map(multiprocess_recurse_func, runnables)
        runnables = []
        for returns in return_lists:
            for ret in returns:
                if not isinstance(ret[-1], set):
                    fl_sol = [soft_round(float(v)) for v in ret]
                    if not all([isinstance(x, int) for x in fl_sol]):
                        continue
                    return ret
                else:
                    runnables.append(ret)
    return None

def multiprocess_recurse_over_solution_path(a, b, covered = set()):

    solution_list = []
    if len(set([indx for indx,b_ in enumerate(b) if b_<0])) == 0:
        solution_list.append(b)
    pool = mp.Pool()
    runnables = [(a, b, 5, set())]
    while len(runnables) > 0:
        return_lists = pool.map(multiprocess_recurse_func, runnables)
        runnables = []
        for returns in return_lists:
            for ret in returns:
                if not isinstance(ret[-1], set):
                    solution_list.append(ret)
                else:
                    runnables.append(ret)
    return solution_list


if __name__ == "__main__":

    import time

    # #Constructed Normal
    # a = Matrix([[-1, 2, -1, 1, -2, 1, 0], [0, -1, 2, 0, -2, 1, 0], [-1, 2, 0, -2, 1, 0, 0], [-1, 1, 1, 0, -1, -1, 1]])
    # b = Matrix([[3], [3], [6], [4], [3], [1], [0]])

    # # Constructed Uniform
    # a = Matrix([[-1, 2, -1, 1, -2, 1, 0], [0, -1, 2, 0, -2, 1, 0], [-1, 2, 0, -2, 1, 0, 0], [-1, 1, 1, 0, -1, -1, 1]])
    # b = Matrix([[3], [3], [3], [3], [3], [3], [3]])

    # # Constructed Bimodal
    # a = Matrix([[-1, 2, -1, 1, -2, 1, 0], [0, -1, 2, 0, -2, 1, 0], [-1, 2, 0, -2, 1, 0, 0], [-1, 1, 1, 0, -1, -1, 1]])
    # b = Matrix([[5], [3], [2], [4], [2], [2], [2]])

    # # Constructed Skew
    # a = Matrix([[-1, 2, -1, 1, -2, 1, 0], [0, -1, 2, 0, -2, 1, 0], [-1, 2, 0, -2, 1, 0, 0], [-1, 1, 1, 0, -1, -1, 1]])
    # b = Matrix([[9], [6], [3], [3], [-1], [0], [0]])

    a = Matrix([[-1,  2, -1,  1, -2,  1, 0],[ 0, -1,  2,  0, -2,  1, 0],[-1,  2,  0, -2,  1,  0, 0],[-1,  1,  1,  0, -1, -1, 1]])
    b = Matrix([[0, 1, 2, 3, 1, 2, -1]])

    # a = Matrix([[-1, 1, 0, 1, -1, 0, 0, 1, -1, 0, -1, 1, 0, 0, 0, 0, 0], [0, 0, 1, -1, 0, 0, -1, 0, 1, 0, 0, 1, -1, 0, 0, 0, 0], [-1, 1, 0, 0, 1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0], [0, -1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1], [-1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1, 0], [-1, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, -1, 1, 0, 0], [0, -1, 0, 1, 1, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, 0, 0], [0, -1, 0, 1, 1, 0, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0], [-1, 1, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0], [0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0], [-1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [-1, 1, 1, 0, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    # b = Matrix([[2], [1], [2], [2], [2], [1], [1], [1], [1], [1], [0], [0], [1], [0], [-1], [0], [-1]])

    # a = Matrix([[-1, 1, 0, 1, -1, 0, 0, 1, -1, 0, -1, 1, 0, 0, 0], [0, 0, 1, -1, 0, 0, -1, 0, 1, 0, 0, 1, -1, 0, 0], [-1, 1, 0, 0, 1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0], [-1, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, -1, 1], [0, -1, 0, 1, 1, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0], [0, -1, 0, 1, 1, 0, 0, 0, -1, -1, 0, 1, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0], [-1, 1, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 1, 0, 0], [0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0], [0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0], [-1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0], [-1, 1, 1, 0, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0]])
    # b = Matrix([[-1], [0], [0], [2], [1], [2], [2], [2], [2], [1], [1], [1], [0], [0], [-1]])


    # a = Matrix([[-1, 1, 0, 1, -1, 0, 0, 1, -1, 0, -1, 1, 0, 0, 0], [0, 0, 1, -1, 0, 0, -1, 0, 1, 0, 0, 1, -1, 0, 0], [-1, 1, 0, 0, 1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0], [-1, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, -1, 1], [0, -1, 0, 1, 1, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0], [0, -1, 0, 1, 1, 0, 0, 0, -1, -1, 0, 1, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0], [-1, 1, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 1, 0, 0], [0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0], [0, 0, -1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0], [-1, 1, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0], [-1, 1, 1, 0, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0]])
    # b = Matrix([[0], [0], [0], [1], [1], [1], [1], [3], [2], [2], [1], [1], [0], [0], [0]])

    a,b = getManipBasis(a,b)
    b = forced_neg_removal(a, b)
    start = time.time()
    sols = multiprocess_recurse_over_solution_path(a, b)
    print sols
    print time.time() - start
    start = time.time()
    sols = recurse_over_solution_path(a, b)
    print sols
    print time.time() - start
    start = time.time()
    sols = multiprocess_recurse_find_first([a], [b], covered=set())
    print sols
    print time.time() - start
    start = time.time()
    sols = recurse_find_first(a, b, covered=set())
    print sols
    print time.time() - start

################## Older code that uses a worse form of multiprocessing ############################
#
# def multi_thread_recurse_func(params):
#     a_indx, a_, a, b, covered = params
#     covered.update(xrange(a_indx))
#     sols = []
#     b_temp = list_sum(a_, b)
#     temp_neg =set([indx for indx,b_ in enumerate(b_temp) if b_<0])
#     if len(temp_neg) == 0:
#         # print b_temp
#         sols.append(b_temp)
#     if not dead_solution(a, b):#temp_neg.issubset(neg_indices):
#         # print covered
#         sols.extend(recurse_over_solution_path(a, b_temp, covered=set(covered), multiprocessing=False))
#     return sols
#
# def recurse_over_solution_path(a, b, covered = set(), multiprocessing=None):
#     sols = []
#     # print b
#     # neg_indices = set([indx for indx,b_ in enumerate(b) if b_<0])
#     if multiprocessing:
#         sols.extend(list(itertools.chain.from_iterable(
#             multiprocessing.map(multi_thread_recurse_func,
#                      [(a_indx, a_, a, b, covered) for (a_indx, a_) in enumerate(a)
#                                                  if a_indx not in covered]
#                      )
#         )
#         )
#         )
#     else:
#         for a_indx, a_ in enumerate(a):
#             if a_indx in covered:
#                 continue
#             b_temp = list_sum(a_, b)
#             temp_neg =set([indx for indx,b_ in enumerate(b_temp) if b_<0])
#             if len(temp_neg) == 0:
#                 # print b_temp
#                 sols.append(b_temp)
#             if not dead_solution(a, b):#if temp_neg.issubset(neg_indices):
#                 # print covered
#                 sols.extend(recurse_over_solution_path(a, b_temp, covered=set(covered), multiprocessing=False))
#             covered.add(a_indx)
#
#     return sols
#
# def find_dependencies(a, b):
#     if isinstance(a, list):
#         a = Matrix(a)
#     if isinstance(b, list):
#         b = Matrix(b)
#     assert isinstance(a, Matrix) and isinstance(b,Matrix)
#     col_row_deps = []
#     for col_index in xrange(a.cols):
#         column = list(a.col(col_index))
#         col_row_deps.append([])
#         for i in xrange(len(column)):
#             col_row_deps[-1].append(column[i])
#     return col_row_deps