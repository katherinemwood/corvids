__author__ = 'Sean Wilner'
from findSolutionsWithManipBasis import *
from sympy import Matrix
from decimal import  Decimal
import math, random, itertools, functools, collections

from mpl_toolkits.mplot3d import Axes3D
from sys import platform as sys_pf
if sys_pf == 'darwin':
     import matplotlib
     matplotlib.use("TkAgg")
else:
    import matplotlib
import matplotlib.pyplot as plt
import numpy as np


dio = Diophantine()

def multiprocessGetManipBases(basis_and_base_vec):
    basis, base_vec = basis_and_base_vec
    manip_basis, base_vec = getManipBasis(basis, base_vec)
    # print base_vec
    # print manip_basis

    manip_base_vec = forced_neg_removal(manip_basis, base_vec)
    return (manip_basis, manip_base_vec)

def multiprocessGetSolutionSpace(min_score, max_score, num_samples, mean_and_variance,
                                 check_val=None, poss_vals = None, debug=True):
    '''
    A function call to mimic the one on the RecreateData object while remaining multiprocess compatible
        (can't pickle an object method within the object namespace)

    This function handles getting a valid initial solution (possibly with negative values) and a corresponding vector
        space of valid transformations to that solution.

    Typical users should never need to call this function as it is called from within the recreateData() method on
        the RecreateData object.
    :param min_score: minimum value in the range of possible values for the dataset
    :param max_score: maximum value in the range of possible values for the dataset
    :param num_samples: total number of samples to find solutions for
    :param mean_and_variance: a tuple containing the mean and variance of the dataset e.g. (mean,variance)
    :param check_val: values we want to assume the solution space MUST contain presented as either a single integer
                            or as a list of integers, or as a dictionary whose keys is the integers assumed to be in
                            the dataset with corresponding values of the number of times the integer is assumed to
                            appear in the dataset.  That is, this value lets the user check for specific values in
                            solutions.
    :param poss_vals: an iterable (eg list) containing all potential values to consider when constructing viable
                            datasets.  That is, this value lets the user remove given values from consideration for
                            solution spaces (by presenting all other potential values).
    :param debug: Boolean to indicate if print statements are allowed.
    :return: Either:
                            A 5-tuple containing an initial (potentially negative) solution, the basis for a vector space of
                            transformations, the two matrices used to calculate it, and the mean_and_variance tuple
                            passed as an argument if such a solution exists

             Or:
                            None if no solution exists
    '''
    if not poss_vals:
        poss_vals = xrange(min_score, max_score+1)

    mean, variance = mean_and_variance
    param_tuple = (mean, variance)

    mean *=num_samples
    mean = int(round(mean))
    variance *=(num_samples-1)
    variance *=num_samples**2
    variance = int(Decimal(str(variance)))

    A_list = []
    for i in poss_vals:
        coef = []
        coef.append(1) # participants
        coef.append( (i)) # scaled mean -- total_sum
        coef.append(((num_samples * i)-mean)**2) # variance
        A_list.append(coef)
    A = Matrix(A_list).T
    if check_val:
        if isinstance(check_val,int):
            variance -= (num_samples*check_val - mean)**2
            mean -= check_val
            num_samples -= 1
        elif isinstance(check_val, list):
            for val in check_val:
                variance -= (num_samples*val - mean)**2
                mean -= val
                num_samples -= 1
        elif isinstance(check_val, dict):
            for val, num in check_val.iteritems():
                variance -= num*(num_samples*val - mean)**2
                mean -= val*num
                num_samples -= num
        else:
            raise TypeError
    b = Matrix([num_samples, mean, variance])
    basis = Matrix(dio.getBasis(A, b))

    try:
        base_vec = basis[-1,:]
        basis = basis[:-1,:]
        if debug:
            print "found potential at: " + str(mean_and_variance)
        return base_vec, basis, A, b, param_tuple
    except IndexError:
        return None

class RecreateData:
    '''
    An object which conatins all the relevant information about a given set of summary statistics and allows methods
        to discover all potential solutions
    '''


    def __init__(self, min_score, max_score, num_samples, mean, variance, debug=True, mean_precision=0.0, variance_precision=0.0):

        self.simpleData = defaultdict(list)
        self.debug = debug
        self.absolute_min = min_score
        self.min_score = min_score
        self.absolute_max = max_score
        self.max_score = max_score
        self.num_samples = num_samples
        self.mean = mean
        self.variance = variance
        self.un_mut_num_samples = num_samples
        self.un_mut_mean = mean
        self.un_mut_variance = variance
        self.sols = None
        # self.sols_as_samples
        self.mean_precision = mean_precision
        self.variance_precision = variance_precision
        self.extended_poss_vals = None

    def adjust_sol(self, solution, num_adjustments):
        # print "hoping for ", num_adjustments
        if num_adjustments == 0:
            return solution
        new_sol = deepcopy(solution)
        new_sol.sort()
        num_adjusted = 0
        while num_adjusted < num_adjustments:
            correction = num_adjustments - num_adjusted
            items = list(set(new_sol))
            items.sort()
            if correction == 1:
                candidates = [item for item, count in collections.Counter(new_sol).items() if count > 1 and (item + 1) in self.poss_vals and (item - 1) in self.poss_vals]
                if len(candidates)==0:
                    # print "None1 was encountered ", num_adjustments
                    return None
                potential_diff = candidates[0]
                new_sol.remove(potential_diff)
                new_sol.remove(potential_diff)
                new_sol.append(potential_diff + 1)
                new_sol.append(potential_diff - 1)
                new_sol.sort()
                return new_sol
            potential_diffs = []
            temp_items = []
            if len(new_sol)>len(items):
                for val in items:
                    if new_sol.count(val) == 1:
                        continue
                    temp_items.append(val)
            items.extend(temp_items)
            items.sort()
            for val_1, val_2 in itertools.combinations(items, 2):
                #since we sorted items, val_1 is the smaller one
                for size in xrange(1, min((max(self.poss_vals) - val_2), val_1 - min(self.poss_vals))):
                    if val_1 - size not in self.poss_vals or val_2 + size not in self.poss_vals:
                        continue
                    adjustment = (math.fabs(val_1 - val_2)*size + size**2)
                    if adjustment <= correction:
                        potential_diffs.append((adjustment, val_1, val_2, size))
            # if len(potential_diffs) == 0:
            #     candidates = [item for item, count in collections.Counter(new_sol).items() if count > 1 and (item + 1) in self.poss_vals and (item - 1) in self.poss_vals]
            #     if len(candidates)==0:
            #         return None
            #     candidate = candidates[0]
            #     new_sol.remove(candidate)
            #     new_sol.remove(candidate)
            #     new_sol.append(candidate + 1)
            #     new_sol.append(candidate - 1)
            #     new_sol.sort()
            #     num_adjusted += 1
            #     continue
            if len(potential_diffs) == 0:
                return
            potential_diffs.sort()
            potential_diff = potential_diffs[-1] #Greed algorithm here to choose the biggest adjustment possible
            num_adjusted += potential_diff[0]
            new_sol.remove(potential_diff[1])
            new_sol.remove(potential_diff[2])
            new_sol.append(potential_diff[2] + potential_diff[3])
            new_sol.append(potential_diff[1] - potential_diff[3])
            new_sol.sort()
        return new_sol







    def validMeansVariances(self, findFirst):
        solutions = defaultdict(list)
        means_list = []
        for i in xrange(int(math.ceil((self.mean - self.mean_precision)*self.num_samples)),
                        int(math.floor((self.mean + self.mean_precision)*self.num_samples))+1):
            means_list.append(float(i)/self.num_samples)
        mean_variances = []
        step_size = 2*self.num_samples**2
        for m in means_list:
            # construct the minimum variance (not within our range or anything, we are just going to use this for granularity purposes)
            initial_var = ((self.num_samples - m*self.num_samples%self.num_samples)*(math.floor(m)-m)**2 +
                           (m*self.num_samples%self.num_samples)*(math.ceil(m)-m)**2)/(self.num_samples-1)
            initial_mean_valid = [int(math.floor(m))] * (int((self.num_samples - m*self.num_samples%self.num_samples))) + [int(math.ceil(m))] * int((m*self.num_samples%self.num_samples))
            # initial_var = np.std([int(math.floor(m))] * (int((self.num_samples - m*self.num_samples%self.num_samples))) + [int(math.ceil(m))] * int((m*self.num_samples%self.num_samples)), ddof=1, dtype=np.float128)**2
            # print "Initial Variance: " + str(initial_var) + " with " + str((self.num_samples - m*self.num_samples%self.num_samples)) + " values at " + str(math.floor(m)) + " and "+ str((m*self.num_samples%self.num_samples)) + " values at "  + str((math.ceil(m)))
            initial_adjusted_var = initial_var *(self.num_samples - 1)*self.num_samples**2
            # print "Initial Adjusted: " + str(initial_adjusted_var)
            i = int(math.ceil((self.variance - self.variance_precision)*((self.num_samples-1)*self.num_samples**2)))
            while i <(math.floor((self.variance + self.variance_precision)*((self.num_samples-1)*self.num_samples**2)))+1:
                if (i-initial_adjusted_var)%step_size != 0:
                    i +=step_size - (i-initial_adjusted_var)%step_size
                    continue
                mean_variances.append((m, float(i)/((self.num_samples - 1)*self.num_samples**2)))
                solution = self.adjust_sol(initial_mean_valid, int((i-initial_adjusted_var)/step_size))
                if solution:
                    solutions[(m, float(i)/((self.num_samples - 1)*self.num_samples**2))].append(solution)
                i +=step_size
        if findFirst:
            if len(solutions) > 0:
                self.simpleData.update(solutions)
                if self.debug:
                    print str(sum([len(x) for x in self.simpleData.itervalues()])) + " unique solutions found simulatenously (not neccessarily complete!!)."
                    index = 0
                    for params in self.simpleData:
                        if index > 100:
                            break
                        print "At mean, variance", params, ":"
                        for simpleSol in self.simpleData[params]:
                            if index > 100:
                                break
                            index += 1
                            print simpleSol

                print "Done."
                return
        if self.debug:
            print "Total potential mean/variance pairs to consider: " + str(len(mean_variances))
        return mean_variances


    # def validMeansVariances(self):
    #     means_list = []
    #     for i in xrange(int(math.ceil((self.mean - self.mean_precision)*self.num_samples)),
    #                     int(math.floor((self.mean + self.mean_precision)*self.num_samples))+1):
    #         means_list.append(float(i)/self.num_samples)
    #     variances_list = []
    #     for i in xrange(int(math.ceil((self.variance - self.variance_precision)*((self.num_samples-1)*self.num_samples**2))),
    #                     int(math.floor((self.variance + self.variance_precision)*((self.num_samples-1)*self.num_samples**2)))+1):
    #         variances_list.append(float(i)/((self.num_samples - 1)*self.num_samples**2))
    #
    #     if self.debug:
    #         print str(len(means_list) * len(variances_list)) + " possible mean and variance combinations to consider."
    #
    #     return means_list, variances_list

    def getSolutionSpace(self, check_val=None, poss_vals=None):


        mean = self.mean
        variance = self.variance

        if not poss_vals:
            poss_vals = xrange(self.absolute_min, self.absolute_max+1)
        else:
            self.min_score = min(poss_vals)
            self.max_score = max(poss_vals)

        mean *=self.num_samples
        mean = int(round(mean))
        variance *=(self.num_samples-1)
        variance *=self.num_samples**2
        variance = int(Decimal(str(variance)))

        A_list = []
        for i in poss_vals:
            coef = []
            coef.append(1) # participants
            coef.append( (i)) # scaled mean -- total_sum
            coef.append((self.num_samples * (i)-mean)**2) # variance
            A_list.append(coef)
        A = Matrix(A_list).T
        if check_val:
            if isinstance(check_val,int):
                variance -= (self.num_samples*check_val - mean)**2
                mean -= check_val
                self.num_samples -= 1
            elif isinstance(check_val, list):
                for val in check_val:
                    variance -= (self.num_samples*val - mean)**2
                    mean -= val
                    self.num_samples -= 1
            elif isinstance(check_val, dict):
                for val, num in check_val.iteritems():
                    variance -= num*(self.num_samples*val - mean)**2
                    mean -= val*num
                    self.num_samples -= num
            else:
                raise TypeError
        b = Matrix([self.num_samples, mean, variance])
        basis = Matrix(dio.getBasis(A, b))
        self.num_samples = self.un_mut_num_samples
        self.mean = self.un_mut_mean
        self.variance = self.un_mut_variance
        try:
            base_vec = basis[-1,:]
            basis = basis[:-1,:]
            return base_vec, basis, A, b
        except IndexError:
            if self.debug:
                print "No solutions exist (postive or otherwise)"
            return None


    def _recreateData_piece_1(self, check_val=None, poss_vals=None, multiprocess=True, find_first=False):
        means_list = [self.mean]
        variances_list = [self.variance]

        if not poss_vals:
            poss_vals = range(self.absolute_min, self.absolute_max+1)
        else:
            self.min_score = min(poss_vals)
            self.max_score = max(poss_vals)
        self.poss_vals = poss_vals


        # if self.variance_precision or self.mean_precision > 0:
        #     means_list, variances_list = self.validMeansVariances()
        #
        # mean_variance_pairs = []
        # for pair in itertools.product(means_list, variances_list):
        #     mean_variance_pairs.append((pair[0], pair[1]))

        mean_variance_pairs = self.validMeansVariances(find_first)

        return mean_variance_pairs

    def _recreateData_piece_2(self, mean_variance_pairs, check_val=None, poss_vals=None, multiprocess=True, find_first=False):
        if self.debug:
            print "Checking for potential solution spaces."
        if multiprocess:
            pool = mp.Pool()
            func = functools.partial(multiprocessGetSolutionSpace,self.min_score, self.max_score, self.num_samples,
                                     check_val=check_val, poss_vals=poss_vals, debug=self.debug)
            solution_spaces = pool.map(func, mean_variance_pairs)
            pool.close()
            pool.join()
        else:
            solution_spaces = []
            for mean_variance_pair in mean_variance_pairs:
                solution_spaces.append(multiprocessGetSolutionSpace(self.min_score, self.max_score, self.num_samples,
                                                                    mean_variance_pair,
                                                                    check_val=check_val, poss_vals=poss_vals,
                                                                    debug=self.debug))
        return solution_spaces

    # def _findAll_piece_1_sng_proc(self, solution_spaces, check_val=None, poss_vals=None, multiprocess=True, find_first=False):
    #     self.sols = {}
    #     init_base_vecs = []
    #     init_bases = []
    #     param_tuples = []
    #     for solution_space in solution_spaces:
    #         if solution_space == None:
    #             continue
    #
    #         base_vec, basis, _, _, param_tuple = solution_space
    #         if len(basis) == 0:
    #             if not any([val<0 for val in base_vec._mat]):
    #                 sol = base_vec._mat
    #                 temp_sol =[int(v) for v in sol]
    #                 if check_val:
    #                     if isinstance(check_val,int):
    #                         temp_sol[poss_vals.index(check_val)]+=1
    #
    #                     elif isinstance(check_val, list):
    #                         for val in check_val:
    #                             temp_sol[poss_vals.index(val)]+=1
    #                     elif isinstance(check_val, dict):
    #                         for val, num in check_val.iteritems():
    #                             temp_sol[poss_vals.index(val)]+=num
    #                     else:
    #                         raise TypeError
    #                 self.sols[param_tuple] = [temp_sol]
    #             continue
    #
    #         manip_basis, base_vec = getManipBasis(basis, base_vec)
    #         manip_base_vec = forced_neg_removal(manip_basis, base_vec)
    #         single_set_sols = recurse_over_solution_path(manip_basis, manip_base_vec, covered=set())
    #         # The returned solutions are lists of numpy data structures, so we cast them back to integers for ease of use
    #
    #         temp_sols = []
    #         for sol in single_set_sols:
    #             temp_sols.append([int(v) for v in sol])
    #             if check_val:
    #                 if isinstance(check_val,int):
    #                     temp_sols[-1][poss_vals.index(check_val)]+=1
    #
    #                 elif isinstance(check_val, list):
    #                     for val in check_val:
    #                         temp_sols[-1][poss_vals.index(val)]+=1
    #                 elif isinstance(check_val, dict):
    #                     for val, num in check_val.iteritems():
    #                         temp_sols[-1][poss_vals.index(val)]+=num
    #                 else:
    #                     raise TypeError
    #         self.sols[param_tuple] = temp_sols
    #     if self.debug:
    #         print "Done."
    #     temp_sols =self.sols
    #     self.sols = {}
    #     for key, value in temp_sols.iteritems():
    #         if len(value)>0:
    #             self.sols[key] = value
    #     return self.sols

    def _findAll_piece_1_multi_proc(self, solution_spaces, check_val=None, poss_vals=None, multiprocess=True, find_first=False):
        self.sols = {}
        init_base_vecs = []
        init_bases = []
        param_tuples = []
        for solution_space in solution_spaces:
            if solution_space == None:
                continue
            base_vec, basis, _, _, param_tuple = solution_space
            if len(basis) == 0:
                if not any([val<0 for val in base_vec._mat]):
                    sol = base_vec._mat
                    temp_sol =[int(v) for v in sol]
                    if check_val:
                        if isinstance(check_val,int):
                            try:
                                temp_sol[poss_vals.index(check_val)]+=1
                            except ValueError:
                                temp_sol.append(1)
                                poss_vals.append(check_val)
                                # temp_sol[poss_vals.index(val)]=1

                        elif isinstance(check_val, list):
                            for val in check_val:
                                try:
                                    temp_sol[poss_vals.index(val)]+=1
                                except ValueError:
                                    temp_sol.append(1)
                                    poss_vals.append(val)
                                    # temp_sol[poss_vals.index(val)]=1
                        elif isinstance(check_val, dict):
                            for val, num in check_val.iteritems():
                                try:
                                    temp_sol[poss_vals.index(val)]+=1
                                except ValueError:
                                    poss_vals.append(val)
                                    temp_sol.append(1)
                                    # temp_sol[poss_vals.index(val)]=1
                        else:
                            raise TypeError
                    self.sols[param_tuple] = [temp_sol]

                continue
            init_base_vecs.append(base_vec)
            init_bases.append(basis)
            param_tuples.append(param_tuple)
        if self.debug:
            print "Found " + str(len(param_tuples) + len(self.sols)) + " potentially viable mean/variance pairs." #Get Katherine to UPDATE this!
            print "Manipulating Bases and Initial Vectors for Complete Search Guarantee"
        return init_bases, init_base_vecs, param_tuples

    def _findAll_piece_2_multi_proc(self, init_bases, init_base_vecs, param_tuples, check_val=None, poss_vals=None, multiprocess=True, find_first=False):
        pool = mp.Pool()
        # bases_and_inits = pool.map(multiprocessGetManipBases, zip(init_bases, init_base_vecs))
        bases_and_inits= []
        for X in zip(init_bases, init_base_vecs):
            bases_and_inits.append(multiprocessGetManipBases(X))

        for basis_and_init, param_tuple in zip(bases_and_inits, param_tuples):
            if self.debug:
                print "Checking for solutions at: " + str(param_tuple)
            manip_base_vec = basis_and_init[1]
            manip_basis = basis_and_init[0]
            # print manip_base_vec
            # print manip_basis
            single_set_sols = multiprocess_recurse_over_solution_path(manip_basis, manip_base_vec)
            temp_sols = []
            for sol in single_set_sols:
                fl_sol = [float(v) for v in sol]
                if not all([x.is_integer() for x in fl_sol]):
                    continue
                temp_sols.append([int(v) for v in fl_sol])
                if check_val:
                    if isinstance(check_val,int):
                        try:
                            temp_sols[-1][poss_vals.index(check_val)]+=1
                        except ValueError:
                            for temp_sol in temp_sols:
                                temp_sol.append(0)
                            poss_vals.append(check_val)
                            temp_sols[-1][poss_vals.index(check_val)]=1

                    elif isinstance(check_val, list):
                        for val in check_val:
                            try:
                                temp_sols[-1][poss_vals.index(val)]+=1
                            except ValueError:
                                for temp_sol in temp_sols:
                                    temp_sol.append(0)
                                poss_vals.append(val)
                                temp_sols[-1][poss_vals.index(val)]=1
                    elif isinstance(check_val, dict):
                        for val, num in check_val.iteritems():
                            try:
                                temp_sols[-1][poss_vals.index(val)]+=1
                            except ValueError:
                                poss_vals.append(val)
                                for temp_sol in temp_sols:
                                    temp_sol.append(0)
                                temp_sols[-1][poss_vals.index(val)]=1
                    else:
                        raise TypeError
            self.sols[param_tuple] = temp_sols

        if self.debug:
            print "Done."
        temp_sols =self.sols
        self.sols = {}
        for key, value in temp_sols.iteritems():
            if len(value)>0:
                self.sols[key] = value
        self.extended_poss_vals = poss_vals
        return self.sols


    def _findFirst_piece_1(self, solution_spaces, check_val=None, poss_vals=None, multiprocess=True, find_first=False):
        base_vecs = []
        bases = []

        if multiprocess:
            init_base_vecs = []
            init_bases = []
            for solution_space in solution_spaces:
                if solution_space == None:
                    continue
                base_vec, basis, _, _, param_tuple = solution_space
                if len(basis) == 0:
                    if not any([val<0 for val in base_vec._mat]):
                        sol = base_vec._mat
                        temp_sol =[int(v) for v in sol]
                        if check_val:
                            if isinstance(check_val,int):
                                try:
                                    temp_sol[poss_vals.index(check_val)]+=1
                                except ValueError:
                                    temp_sol.append(1)
                                    poss_vals.append(check_val)
                                    # temp_sol[poss_vals.index(val)]=1

                            elif isinstance(check_val, list):
                                for val in check_val:
                                    try:
                                        temp_sol[poss_vals.index(val)]+=1
                                    except ValueError:
                                        temp_sol.append(1)
                                        poss_vals.append(val)
                                        # temp_sol[poss_vals.index(val)]=1
                            elif isinstance(check_val, dict):
                                for val, num in check_val.iteritems():
                                    try:
                                        temp_sol[poss_vals.index(val)]+=1
                                    except ValueError:
                                        poss_vals.append(val)
                                        temp_sol.append(1)
                                        # temp_sol[poss_vals.index(val)]=1
                            else:
                                raise TypeError
                        self.sols[param_tuple] = [temp_sol]
                        self.extended_poss_vals = poss_vals
                        return self.sols
                init_base_vecs.append(base_vec)
                init_bases.append(basis)
            pool = mp.Pool()
            bases_and_inits = pool.map(multiprocessGetManipBases, zip(init_bases, init_base_vecs))
            for basis_and_init in bases_and_inits:
                base_vecs.append(basis_and_init[1])
                bases.append(basis_and_init[0])
            sol = multiprocess_recurse_find_first(bases, base_vecs, covered=set())
            # print sol
        else:
            for solution_space in solution_spaces:
                if solution_space == None:
                    continue
                base_vec, basis, _, _, param_tuple = solution_space
                if len(basis) == 0:
                    if not any([val<0 for val in base_vec._mat]):
                        sol = base_vec._mat
                        temp_sol =[int(v) for v in sol]
                        if check_val:
                            if isinstance(check_val,int):
                                try:
                                    temp_sol[poss_vals.index(check_val)]+=1
                                except ValueError:
                                    temp_sol.append(1)
                                    poss_vals.append(check_val)
                                    # temp_sol[poss_vals.index(val)]=1

                            elif isinstance(check_val, list):
                                for val in check_val:
                                    try:
                                        temp_sol[poss_vals.index(val)]+=1
                                    except ValueError:
                                        temp_sol.append(1)
                                        poss_vals.append(val)
                                        # temp_sol[poss_vals.index(val)]=1
                            elif isinstance(check_val, dict):
                                for val, num in check_val.iteritems():
                                    try:
                                        temp_sol[poss_vals.index(val)]+=1
                                    except ValueError:
                                        poss_vals.append(val)
                                        temp_sol.append(1)
                                        # temp_sol[poss_vals.index(val)]=1
                            else:
                                raise TypeError
                        self.sols[param_tuple] = [temp_sol]
                        self.extended_poss_vals = poss_vals
                        return self.sols
                manip_basis, base_vec = getManipBasis(basis, base_vec)
                manip_base_vec = forced_neg_removal(manip_basis, base_vec)
                base_vecs.append(manip_base_vec)
                bases.append(manip_basis)

            for basis, base_vec in zip(bases, base_vecs):
                sol = recurse_find_first(basis, base_vec)
                if sol:
                    break
        if not sol:
            self.extended_poss_vals = poss_vals
            return None
        # print "Solution: ", type(sol[0]), sol
        sol = [int(v) for v in sol]
        # print "Solution: ", type(sol[0]), sol

        if check_val:
            if isinstance(check_val,int):
                try:
                    sol[poss_vals.index(check_val)]+=1
                except ValueError:
                    sol.append(1)
                    poss_vals.append(check_val)
                    # temp_sol[poss_vals.index(val)]=1

            elif isinstance(check_val, list):
                for val in check_val:
                    try:
                        sol[poss_vals.index(val)]+=1
                    except ValueError:
                        sol.append(1)
                        poss_vals.append(val)
                        # temp_sol[poss_vals.index(val)]=1
            elif isinstance(check_val, dict):
                for val, num in check_val.iteritems():
                    try:
                        sol[poss_vals.index(val)]+=1
                    except ValueError:
                        poss_vals.append(val)
                        sol.append(1)
                        # temp_sol[poss_vals.index(val)]=1
            else:
                raise TypeError

        self.sols = {'_':[sol]}
        if self.debug:
            print "Done."
        self.extended_poss_vals = poss_vals
        return self.sols

    def recreateData(self, check_val=None, poss_vals=None, multiprocess=True, find_first=False):

        means_list = [self.mean]
        variances_list = [self.variance]

        if not poss_vals:
            poss_vals = range(self.absolute_min, self.absolute_max+1)
        else:
            self.min_score = min(poss_vals)
            self.max_score = max(poss_vals)
        self.poss_vals = poss_vals

        if self.variance_precision or self.mean_precision > 0:
            means_list, variances_list = self.validMeansVariances()

        mean_variance_pairs = []
        for pair in itertools.product(means_list, variances_list):
            mean_variance_pairs.append((pair[0], pair[1]))

        if self.debug:
            print str(len(mean_variance_pairs)) + " total mean and variance pairs to check."

        if multiprocess:
            pool = mp.Pool()
            func = functools.partial(multiprocessGetSolutionSpace, self.min_score, self.max_score, self.num_samples,
                                     check_val=check_val, poss_vals=poss_vals, debug=self.debug)
            solution_spaces = pool.map(func, mean_variance_pairs)
            pool.close()
            pool.join()
        else:
            solution_spaces = []
            for mean_variance_pair in mean_variance_pairs:
                solution_spaces.append(multiprocessGetSolutionSpace(self.min_score, self.max_score, self.num_samples,
                                                                    mean_variance_pair,
                                                                    check_val=check_val, poss_vals=poss_vals,
                                                                    debug=self.debug))
        # try:
        #     base_vec, basis, _, _, _ = self.getSolutionSpace(check_val=check_val, poss_vals=poss_vals)
        #     if self.debug:
        #         print "At least one solution exists (possibly invalid due to negative values)"
        # except TypeError as E:
        #     return None

        if find_first:
            base_vecs = []
            bases = []

            if multiprocess:
                init_base_vecs = []
                init_bases = []
                for solution_space in solution_spaces:
                    if solution_space == None:
                        continue
                    base_vec, basis, _, _, param_tuple = solution_space
                    if len(basis) == 0:
                        if not any([val<0 for val in base_vec._mat]):
                            sol = base_vec._mat
                            temp_sol =[int(v) for v in sol]
                            if check_val:
                                if isinstance(check_val,int):
                                    try:
                                        temp_sol[poss_vals.index(check_val)]+=1
                                    except ValueError:
                                        temp_sol.append(1)
                                        poss_vals.append(check_val)
                                        # temp_sol[poss_vals.index(val)]=1

                                elif isinstance(check_val, list):
                                    for val in check_val:
                                        try:
                                            temp_sol[poss_vals.index(val)]+=1
                                        except ValueError:
                                            temp_sol.append(1)
                                            poss_vals.append(val)
                                            # temp_sol[poss_vals.index(val)]=1
                                elif isinstance(check_val, dict):
                                    for val, num in check_val.iteritems():
                                        try:
                                            temp_sol[poss_vals.index(val)]+=1
                                        except ValueError:
                                            poss_vals.append(val)
                                            temp_sol.append(1)
                                            # temp_sol[poss_vals.index(val)]=1
                                else:
                                    raise TypeError

                            self.sols[param_tuple] = [temp_sol]
                            self.extended_poss_vals = poss_vals
                            return self.sols

                    init_base_vecs.append(base_vec)
                    init_bases.append(basis)
                pool = mp.Pool()
                bases_and_inits = pool.map(multiprocessGetManipBases, zip(init_bases, init_base_vecs))
                for basis_and_init in bases_and_inits:
                    base_vecs.append(basis_and_init[1])
                    bases.append(basis_and_init[0])
                sol = multiprocess_recurse_find_first(bases, base_vecs, covered=set())
            else:
                for solution_space in solution_spaces:
                    if solution_space == None:
                        continue
                    base_vec, basis, _, _, param_tuple = solution_space
                    if len(basis) == 0:
                        if not any([val<0 for val in base_vec._mat]):
                            sol = base_vec._mat
                            temp_sol =[int(v) for v in sol]
                            if check_val:
                                if isinstance(check_val,int):
                                    try:
                                        temp_sol[poss_vals.index(check_val)]+=1
                                    except ValueError:
                                        temp_sol.append(1)
                                        poss_vals.append(check_val)
                                        # temp_sol[poss_vals.index(val)]=1

                                elif isinstance(check_val, list):
                                    for val in check_val:
                                        try:
                                            temp_sol[poss_vals.index(val)]+=1
                                        except ValueError:
                                            temp_sol.append(1)
                                            poss_vals.append(val)
                                            # temp_sol[poss_vals.index(val)]=1
                                elif isinstance(check_val, dict):
                                    for val, num in check_val.iteritems():
                                        try:
                                            temp_sol[poss_vals.index(val)]+=1
                                        except ValueError:
                                            poss_vals.append(val)
                                            temp_sol.append(1)
                                            # temp_sol[poss_vals.index(val)]=1
                                else:
                                    raise TypeError
                            self.sols[param_tuple] = [temp_sol]
                            self.extended_poss_vals = poss_vals
                            return self.sols
                    manip_basis, base_vec = getManipBasis(basis, base_vec)
                    manip_base_vec = forced_neg_removal(manip_basis, base_vec)
                    base_vecs.append(manip_base_vec)
                    bases.append(manip_basis)

                for basis, base_vec in zip(bases, base_vecs):
                    sol = recurse_find_first(basis, base_vec)
                    if sol:
                        break
            if not sol:
                self.extended_poss_vals = poss_vals
                return None
            # print "Solution: ", type(sol[0]), sol
            sol = [int(v) for v in sol]
            # print "Solution: ", type(sol[0]), sol

            if check_val:
                if isinstance(check_val,int):
                    try:
                        sol[poss_vals.index(check_val)]+=1
                    except ValueError:
                        sol.append(1)
                        poss_vals.append(check_val)
                        # temp_sol[poss_vals.index(val)]=1

                elif isinstance(check_val, list):
                    for val in check_val:
                        try:
                            sol[poss_vals.index(val)]+=1
                        except ValueError:
                            sol.append(1)
                            poss_vals.append(val)
                            # temp_sol[poss_vals.index(val)]=1
                elif isinstance(check_val, dict):
                    for val, num in check_val.iteritems():
                        try:
                            sol[poss_vals.index(val)]+=1
                        except ValueError:
                            poss_vals.append(val)
                            sol.append(1)
                            # temp_sol[poss_vals.index(val)]=1
                else:
                    raise TypeError

            self.sols = {'_':[sol]}
            if self.debug:
                print "Done."
            self.extended_poss_vals = poss_vals
            return self.sols

        self.sols = {}
        init_base_vecs = []
        init_bases = []
        param_tuples = []
        if multiprocess:
            for solution_space in solution_spaces:
                if solution_space == None:
                    continue
                base_vec, basis, _, _, param_tuple = solution_space
                if len(basis) == 0:
                    if not any([val<0 for val in base_vec._mat]):
                        sol = base_vec._mat
                        temp_sol =[int(v) for v in sol]
                        if check_val:
                            if isinstance(check_val,int):
                                try:
                                    temp_sol[poss_vals.index(check_val)]+=1
                                except ValueError:
                                    temp_sol.append(1)
                                    poss_vals.append(check_val)
                                    # temp_sol[poss_vals.index(val)]=1

                            elif isinstance(check_val, list):
                                for val in check_val:
                                    try:
                                        temp_sol[poss_vals.index(val)]+=1
                                    except ValueError:
                                        temp_sol.append(1)
                                        poss_vals.append(val)
                                        # temp_sol[poss_vals.index(val)]=1
                            elif isinstance(check_val, dict):
                                for val, num in check_val.iteritems():
                                    try:
                                        temp_sol[poss_vals.index(val)]+=1
                                    except ValueError:
                                        poss_vals.append(val)
                                        temp_sol.append(1)
                                        # temp_sol[poss_vals.index(val)]=1
                            else:
                                raise TypeError
                        self.sols[param_tuple] = [temp_sol]

                    continue
                init_base_vecs.append(base_vec)
                init_bases.append(basis)
                param_tuples.append(param_tuple)
            pool = mp.Pool()
            if self.debug:
                print "Found " + str(len(param_tuples)) + " potentially viable mean/variance pairs found."
                print "Manipulating Bases and Initial Vectors for Complete Search Guarantee"
            bases_and_inits = pool.map(multiprocessGetManipBases, zip(init_bases, init_base_vecs))
            for basis_and_init, param_tuple in zip(bases_and_inits, param_tuples):
                if self.debug:
                    print "Checking for solutions at: " + str(param_tuple)
                manip_base_vec = basis_and_init[1]
                manip_basis = basis_and_init[0]
                single_set_sols = multiprocess_recurse_over_solution_path(manip_basis, manip_base_vec)
                temp_sols = []
                for sol in single_set_sols:
                    fl_sol = [float(v) for v in sol]
                    if not all([x.is_integer() for x in fl_sol]):
                        continue
                    temp_sols.append([int(v) for v in fl_sol])
                    if check_val:
                        if isinstance(check_val,int):
                            try:
                                temp_sols[-1][poss_vals.index(check_val)]+=1
                            except ValueError:
                                for temp_sol in temp_sols:
                                    temp_sol.append(0)
                                poss_vals.append(check_val)
                                temp_sols[-1][poss_vals.index(check_val)]=1

                        elif isinstance(check_val, list):
                            for val in check_val:
                                try:
                                    temp_sols[-1][poss_vals.index(val)]+=1
                                except ValueError:
                                    for temp_sol in temp_sols:
                                        temp_sol.append(0)
                                    poss_vals.append(val)
                                    temp_sols[-1][poss_vals.index(val)]=1
                        elif isinstance(check_val, dict):
                            for val, num in check_val.iteritems():
                                try:
                                    temp_sols[-1][poss_vals.index(val)]+=1
                                except ValueError:
                                    poss_vals.append(val)
                                    for temp_sol in temp_sols:
                                        temp_sol.append(0)
                                    temp_sols[-1][poss_vals.index(val)]=1
                        else:
                            raise TypeError
                self.sols[param_tuple] = temp_sols
        else:
            for solution_space in solution_spaces:
                if solution_space == None:
                    continue

                base_vec, basis, _, _, param_tuple = solution_space
                if len(basis) == 0:
                    if not any([val<0 for val in base_vec._mat]):
                        sol = base_vec._mat
                        temp_sol =[int(v) for v in sol]
                        if check_val:
                            if isinstance(check_val,int):
                                try:
                                    temp_sol[poss_vals.index(check_val)]+=1
                                except ValueError:
                                    temp_sol.append(1)
                                    poss_vals.append(check_val)
                                    # temp_sol[poss_vals.index(val)]=1

                            elif isinstance(check_val, list):
                                for val in check_val:
                                    try:
                                        temp_sol[poss_vals.index(val)]+=1
                                    except ValueError:
                                        temp_sol.append(1)
                                        poss_vals.append(val)
                                        # temp_sol[poss_vals.index(val)]=1
                            elif isinstance(check_val, dict):
                                for val, num in check_val.iteritems():
                                    try:
                                        temp_sol[poss_vals.index(val)]+=1
                                    except ValueError:
                                        poss_vals.append(val)
                                        temp_sol.append(1)
                                        # temp_sol[poss_vals.index(val)]=1
                            else:
                                raise TypeError
                        self.sols[param_tuple] = [temp_sol]
                    continue

                manip_basis, base_vec = getManipBasis(basis, base_vec)
                manip_base_vec = forced_neg_removal(manip_basis, base_vec)
                single_set_sols = recurse_over_solution_path(manip_basis, manip_base_vec, covered=set())
                # The returned solutions are lists of numpy data structures, so we cast them back to integers for ease of use

                temp_sols = []
                for sol in single_set_sols:
                    fl_sol = [float(v) for v in sol]
                    if not all([x.is_integer() for x in fl_sol]):
                        continue
                    temp_sols.append([int(v) for v in fl_sol])
                    if check_val:
                        if isinstance(check_val,int):
                            try:
                                temp_sols[-1][poss_vals.index(check_val)]+=1
                            except ValueError:
                                for temp_sol in temp_sols:
                                    temp_sol.append(0)
                                poss_vals.append(check_val)
                                temp_sols[-1][poss_vals.index(check_val)]=1

                        elif isinstance(check_val, list):
                            for val in check_val:
                                try:
                                    temp_sols[-1][poss_vals.index(val)]+=1
                                except ValueError:
                                    for temp_sol in temp_sols:
                                        temp_sol.append(0)
                                    poss_vals.append(val)
                                    temp_sols[-1][poss_vals.index(val)]=1
                        elif isinstance(check_val, dict):
                            for val, num in check_val.iteritems():
                                try:
                                    temp_sols[-1][poss_vals.index(val)]+=1
                                except ValueError:
                                    poss_vals.append(val)
                                    for temp_sol in temp_sols:
                                        temp_sol.append(0)
                                    temp_sols[-1][poss_vals.index(val)]=1
                        else:
                            raise TypeError
                self.sols[param_tuple] = temp_sols
        if self.debug:
            print "Done."
        temp_sols =self.sols
        self.sols = {}
        for key, value in temp_sols.iteritems():
            if len(value)>0:
                self.sols[key] = value
        self.extended_poss_vals = poss_vals
        return self.sols



    def analyzeSkew(self):

        import scipy.stats.stats
        import numpy as np

        if not self.sols:
            if self.debug:
                print "No solutions to run analysis over.  NB: recreateData() must be run before analyzeSkew()"
            raise ValueError

        skews = []
        for sol_set in self.sols.values():
            for dist in sol_set:
                sol = []
                for i, x in enumerate(dist):
                    sol += [i + 1]*x
                skews.append(scipy.stats.stats.skew(sol))

        return skews, np.mean(skews), np.std(skews, dtype=np.float64, ddof=1.0)**2


    def graphData(self, max_samples=40):

        # if self.simpleData and not self.sols:
        #     self.sols = {}
        #     for key, value in self.simpleData.iteritems():
        #         sol =


        if not self.sols:
            if self.debug:
                print "No solutions to run analysis over.  NB: graphData() must be run before analyzeSkew()"
            raise ValueError




        sols = []
        all_sols = []
        for sol_set in self.sols.values():
            all_sols.extend(sol_set)
        if max_samples == -1:
            sols = all_sols
        else:
            sols = random.sample(all_sols, min(max_samples, len(all_sols)))
        if len(sols) == 1:
            sols.append([0]*len(sols[0]))

        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')

        xpos = []
        ypos = []
        dz = []
        for index1, sol in enumerate(sols):
            for index2,val in enumerate(sol):
                ypos.append(self.poss_vals[index2] - .5)
                xpos.append(index1 + .5) #added + .5
                dz.append(val)

        num_elements = len(xpos)
        zpos = [0]*num_elements
        dx = np.ones(num_elements)
        dy = np.ones(num_elements)

        plt.ylim(min(ypos), max(ypos) + 1)

        ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color='#00ceaa')
        ax1.set_xlabel('Solution Number')
        ax1.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1))
        ax1.set_ylabel('Response Value')
        ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1))
        ax1.set_zlabel('Frequency')
        ax1.zaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1))
        plt.show()

    def getDataSimple(self):
        if not self.sols:
            if self.debug:
                print "No solutions to run analysis over.  NB: graphData() must be run before analyzeSkew()"
            raise ValueError
        for param, sol_list in self.sols.iteritems():
            for sol in sol_list:
                simple_sol = []
                poss_vals = self.poss_vals
                if self.extended_poss_vals:
                    poss_vals = self.extended_poss_vals
                for value, num_instances in zip(poss_vals, sol):
                    simple_sol += [value]*num_instances
                self.simpleData[param].append(simple_sol)
        return self.simpleData


if __name__ == "__main__":

    import sys

    # sys.stdout.write('Minimum Value: ')
    # min_score = int(raw_input())
    # print
    #
    # sys.stdout.write('Maximum Value: ')
    # max_score = int(raw_input())
    # print
    #
    # sys.stdout.write('Number of Samples: ')
    # num_samples = int(raw_input())
    # print
    #
    # sys.stdout.write('Mean Value: ')
    # mean = float(eval(raw_input()))
    # print
    #
    # sys.stdout.write('Variance Value: ')
    # variance = float(eval(raw_input()))
    # print

    RD = RecreateData(1,7,10,3,4)#(min_score, max_score, num_samples, mean, variance, mean_precision=0.0, variance_precision=0.0)
    RD.recreateData(multiprocess=True, find_first=False)
    print RD.getDataSimple()

    # for sol_set in RD.sols.values():
    #     for sol in sol_set:
    #         print sol

    print RD.analyzeSkew()

    RD.graphData()