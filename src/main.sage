#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


import datetime
import time
import numpy
import matplotlib.pyplot as plt

from colorama import Fore, Back, Style


class Coppersmith(object):

    def __init__(self, bit_length, delta, m, jsonoutput, debug):
        # generate a complete crt-rsa keyset with the given parameters
        keys = generate_keys(
                bit_length=bit_length,
                delta=delta,
                m=m
                )

        # sort parameters to class variables
        self.__m = m
        self.__delta = delta
        self.__bit_length = bit_length
        self.__N = keys['N']
        self.__e = keys['e']
        self.__p = keys['p']
        self.__q = keys['q']
        self.__kp = keys['kp']
        self.__kq = keys['kq']
        self.__dp = keys['dp']
        self.__dq = keys['dq']
        self.__keys = keys

        # calculate bounds for xp1, xp2, xq1, xq2, yp, yq
        self.__X = ceil(keys['e'] * keys['N']^(delta - 0.5))
        self.__Y = 2*ceil(keys['N']^(0.5))

        self.__matrix = None
        self.__reduced_matrix = None
        self.__col_indice = None
        self.__inverted_col_indice = None
        self.__polynomials_tuple = None
        self.__row_index = None
        self.__detB = None
        self.__detL = None
        self.__tau = None
        self.__polynom_vector = []
        self.__results = None
        self.__recoverable = False

        self.__duration = 0
        self.__duration_lll = 0
        self.__duration_groebner = 0

        self.__jsonoutput = jsonoutput
        self.__debug = debug

    def print_stats(self):
        # print some stats
        if not self.__jsonoutput:
            pprint('generated crt-rsa parameters')
            pprint('N = {}'.format(self.__N))
            pprint('e = {}'.format(self.__e))
            pprint('dp = {}'.format(self.__dp))
            pprint('dq = {}'.format(self.__dq))
            pprint('p = {}'.format(self.__p))
            pprint('q = {}'.format(self.__q))
            pprint('kp = {}'.format(self.__kp))
            pprint('kq = {}'.format(self.__kq))
            pprint('|N| = {} Bit'.format(self.__bit_length))
            pprint('m = {}'.format(self.__m))

    def generate_matrix(self):
        # generate the lattice for our parameters
        self.__matrix, \
        self.__col_indice, \
        self.__polynomials_tuple, \
        self.__row_index, \
        self.__detB, \
        self.__tau = generate_lattice(
                self.__N,
                self.__e,
                self.__X, self.__Y,
                self.__m,
                self.__tau,
                self.__debug,
                self.__jsonoutput
                )

        # invert the colum index relation
        self.__inverted_col_indice \
                = dict((v,k) for k,v in self.__col_indice.iteritems())


    def pre_process_matrix(self):
        # sort the matrix triangular
        self.__matrix, \
        self.__inverted_col_indice, \
        self.__row_index, \
        self.__polynomials_tuple = matrix_sort_triangle(
                self.__matrix,
                self.__inverted_col_indice,
                self.__row_index,
                self.__polynomials_tuple,
                self.__jsonoutput
                )

        self.__matrix = array_to_matrix(self.__matrix)

        # substitute N from the diagonal
        self.__matrix = substitute_N(
                self.__matrix,
                self.__N,
                self.__e,
                self.__m,
                self.__X,
                self.__Y
                )

        # create g(xX, yY)
        self.__matrix = set_upper_bound(
                self.__matrix,
                self.__X,
                self.__Y,
                self.__inverted_col_indice,
                self.__e,
                self.__m
                )

        self.__detL = self.__matrix.determinant()

    def root_check(self):
        substituted_polynomials = []

        # define the polynomial ring
        R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

        # iterate through the rows of the reduced coefficient matrix
        for row in self.__matrix.rows():
            # initialize the polunomial for each row
            p = 0

            # reset the column index
            col_index = 0

            # iterate through all values of the row
            for i in range(len(row)):
                # get the monome grade from the column index
                monom_grade = self.__inverted_col_indice[i]

                # get the coefficient
                coefficient = row[i]

                # set the correct grade to the variables of each monomial
                x_p_1 = xp1^int(monom_grade[0])
                x_p_2 = xp2^int(monom_grade[1])
                x_q_1 = xq1^int(monom_grade[2])
                x_q_2 = xq2^int(monom_grade[3])
                y_p = yp^int(monom_grade[4])
                y_q = yq^int(monom_grade[5])

                exp_xp1 = int(monom_grade[0])
                exp_xp2 = int(monom_grade[1])
                exp_xq1 = int(monom_grade[2])
                exp_xq2 = int(monom_grade[3])
                exp_yp = int(monom_grade[4])
                exp_yq = int(monom_grade[5])

                x_bound = self.__X^(exp_xp1 + exp_xp2 + exp_xq1 + exp_xq2)
                y_bound = self.__Y^(exp_yp + exp_yq)
                
                x_bound_inverse = inverse_mod(x_bound, self.__e**self.__m)
                y_bound_inverse = inverse_mod(y_bound, self.__e**self.__m)

                coefficient = (coefficient * x_bound_inverse * y_bound_inverse) % self.__e**self.__m

                # add the monomial myltyplied with its coefficient to the polunomial
                p += coefficient * x_p_1 * x_p_2 * x_q_1 * x_q_2 * y_p * y_q

            # append the polynomials to the reduced polynomials array
            substituted_polynomials.append(p)

        # check if the roots match
        if root_check(substituted_polynomials, self.__keys, self.__debug, self.__m):
            if not self.__jsonoutput:
                pprint("root check              [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
        else:
            if not self.__jsonoutput:
                pprint("root check              [" + Fore.RED + " failed " + Fore.RESET + "]") 

    def reduction(self):
        # start time for duration
        start_time = time.mktime(datetime.datetime.now().timetuple())

        self.__reduced_matrix = self.__matrix.LLL()

        # start time for duration
        end_time = time.mktime(datetime.datetime.now().timetuple())
        
        # set duration
        self.__duration_lll = end_time - start_time

        # print some logs
        if not self.__jsonoutput:
            pprint("LLL-reduction           [" + Fore.GREEN + " passed " + Fore.RESET + "]") 

    def build_polynomial_vector(self):
        # initialize an array for the polynomials
        reduced_polynomials = []

        # define the polynomial ring
        R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

        # iterate through the rows of the reduced coefficient matrix
        for row in self.__reduced_matrix.rows():
            # initialize the polunomial for each row
            p = 0

            # initialize variables for howgrave graham theorem
            howgrave_sum = 0

            monome_count = 0

            # iterate through all values of the row
            for i in range(len(row)):
                # get the monome grade from the column index
                monom_grade = self.__inverted_col_indice[i]

                # get the coefficient
                coefficient = row[i]

                # increase the monom count by 1 if we find a monom != 0
                if coefficient != 0:
                    monome_count += 1

                # set the correct grade to the variables of each monomial
                x_p_1 = (xp1)^int(monom_grade[0])
                x_p_2 = (xp2)^int(monom_grade[1])
                x_q_1 = (xq1)^int(monom_grade[2])
                x_q_2 = (xq2)^int(monom_grade[3])
                y_p = (yp)^int(monom_grade[4])
                y_q = (yq)^int(monom_grade[5])

                # add the absolute value of the coefficient squared for howgrave
                # grahams theorem howgrave_sum += abs(coefficient)^2
                howgrave_sum += abs(coefficient)^2

                exp_xp1 = int(monom_grade[0])
                exp_xp2 = int(monom_grade[1])
                exp_xq1 = int(monom_grade[2])
                exp_xq2 = int(monom_grade[3])
                exp_yp = int(monom_grade[4])
                exp_yq = int(monom_grade[5])

                x_bound = self.__X^(exp_xp1 + exp_xp2 + exp_xq1 + exp_xq2)
                y_bound = self.__Y^(exp_yp + exp_yq)

                coefficient = ((coefficient) / (x_bound * y_bound))

                # add the monomial multiplied with its coefficient to the polunomial
                p += (coefficient * x_p_1 * x_p_2 * x_q_1 * x_q_2 * y_p * y_q)

            # append the polynomials to the reduced polynomials array
            reduced_polynomials.append(p)

            # get the square root for the euclidian norm of the vector
            howgrave = floor(sqrt(howgrave_sum))

            # howgrave-grahams lemma in order to get polynoms, which are
            # small enough
            if howgrave < ((self.__e**self.__m) / sqrt(monome_count)):
                # append the polynomial to the polunomial vector for calculating the
                # groebner basis
                self.__polynom_vector.append(p)

    def reduced_root_check(self):
        # check the roots of the reduced polynomials
        if reduced_root_check(self.__polynom_vector, self.__keys, self.__debug, self.__m):
            if not self.__jsonoutput:
                pprint("reduced root check      [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
                pprint('create ideal with {} polynoms'.format(len(self.__polynom_vector)))
        else:
            if not self.__jsonoutput:
                pprint("reduced root check      [" + Fore.RED + " failed " + Fore.RESET + "]") 
                pprint('create ideal with {} polynoms'.format(len(self.__polynom_vector)))

    def get_roots(self):
        # create an ideal from the polynom vector
        I = Ideal(self.__polynom_vector)

        # start time for duration
        start_time = time.mktime(datetime.datetime.now().timetuple())

        # calculate the groebner basis with verbose on or off
        if not self.__jsonoutput:
            pprint('calculate groebner basis')
            B = I.groebner_basis(algorithm='libsingular:groebner', prot=True)
        else:
            B = I.groebner_basis(algorithm='libsingular:groebner', prot=False)

        # start time for duration
        end_time = time.mktime(datetime.datetime.now().timetuple())
        
        # set duration
        self.__duration_groebner = end_time - start_time

        # print the basis
        if not self.__jsonoutput:
            pprint('basis length: {}'.format(len(B)))

        equations = []
        x1, x2, x3, x4, y1, y2 = var('x1 x2 x3 x4 y1 y2')

        #if len(B) >= 6:
        for row in B:

            eq1 = 0

            for monom in row.dict():
                #print("{} - {}".format(monom, row.dict()[monom]))
                mon = row.dict()[monom] * x1^monom[0] * x2^monom[1] * y1^monom[4] * y2^monom[5] * x3^monom[2] * x4^monom[3]
                eq1 += mon

            eq = eq1 == 0

            equations.append(eq)

        # solve the linear equations
        self.__results = solve(equations, x1, x2, x3, x4, y1, y2)

    def factor_check(self):
        # check if p or q is in the solutions
        try:
            for res in [ r.rhs() for r in self.__results ]:
                if res == self.__p or res == self.__q:
                    self.__recoverable = True
        except:
            for result in self.__results:
                for res in [ r.rhs() for r in result ]:
                    if res == self.__p or res == self.__q:
                        self.__recoverable = True

        # check if its recoverable
        if self.__recoverable and not self.__jsonoutput:
            pprint('got p and q!')

    def print_json(self):
        # print json output
        if self.__jsonoutput:
            output = {}
            output['m'] = self.__m
            output['tau'] = self.__tau
            output['delta'] = self.__delta
            output['matrix_dimension'] = self.__reduced_matrix.ncols()
            output['bit_length'] = self.__bit_length
            output['recoverable'] = self.__recoverable
            output['duration'] = self.__duration
            output['duration_lll'] = self.__duration_lll
            output['duration_groebner'] = self.__duration_groebner
            print(output)

    def print_heatmap(self, matrix, filename):
        # get matrix with entrys log to base 2 for bitsize
        heatmap_matrix = matrix_to_sizematrix(matrix)

        # plot heatmap of the matrix to file
        print_heatmap(heatmap_matrix, filename)

    def print_timings(self):
        # print time stats
        if not self.__jsonoutput:
            pprint('runtime LLL       {} sec'.format(self.__duration_lll))
            pprint('runtime groebner  {} sec'.format(self.__duration_groebner))
            pprint('runtime total     {} sec'.format(self.__duration))

    def get_matrix(self):
        return self.__matrix

    def get_reduced_matrix(self):
        return self.__reduced_matrix

    def set_duration(self, duration):
        self.__duration = duration


def instantiate_attack(delta, m, bit_length, tau, debug, printmatrix, jsonoutput):
    c = Coppersmith(
            bit_length,
            delta,
            m,
            jsonoutput,
            debug
            )

    c.print_stats()

    # start time for duration
    start_time = time.mktime(datetime.datetime.now().timetuple())

    # bulid up the correct matrix
    c.generate_matrix()

    # print heatmap to file
    c.print_heatmap(
            c.get_matrix(),
            'matrix.png'
            )

    c.pre_process_matrix()
    c.root_check()

    # print heatmap to file
    c.print_heatmap(
            c.get_matrix(),
            'matrix_before_reduction.png'
            )

    # reduce it
    c.reduction()

    # print heatmap to file
    c.print_heatmap(
            c.get_reduced_matrix(),
            'matrix_after_reduction.png'
            )

    # get the polynomials back from the matrix
    c.build_polynomial_vector()
    c.reduced_root_check()

    # groebner
    c.get_roots()

    # end time for duration and set it
    end_time = time.mktime(datetime.datetime.now().timetuple())
    c.set_duration(end_time - start_time)

    # check if we have the correct factors
    c.factor_check()

    # print json and some stats
    c.print_json()
    c.print_timings()


# execute the main function, if the script is called on the commandline
if __name__ == '__main__':
    # import sage modules
    load('./argparser.sage')
    load('./coppersmith.sage')
    load('./keygen.sage')
    load('./print.sage')
    load('./polynomials.sage')
    load('./shiftpolynomials.sage')
    load('./sort.sage')
    load('./substitute.sage')
    load('./test.sage')

    # parse arguments
    delta, m, bit_length, tau, debug, printmatrix, jsonoutput = parse_args()

    # call the main function
    instantiate_attack(
            delta,
            m,
            bit_length,
            tau,
            debug,
            printmatrix,
            jsonoutput
            )
