#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


import sys


def substitute_y(polynomial, N, e, m):
    """
    Substitute yp*yq to N
    """
    # get the polynomial as dict with multigrade => coefficient relation
    polynomial_dict = polynomial.dict()

    # initialize the output polynomial
    output_polynomial = 0

    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lp')

    # iterate through the multigrades
    for multigrade in polynomial_dict:
        # save the coefficient to the monome
        monome = polynomial_dict[multigrade]

        # save the exponents of yq and yp for the substitution
        exp_xp1 = multigrade[0]
        exp_xp2 = multigrade[1]
        exp_xq1 = multigrade[2]
        exp_xq2 = multigrade[3]
        exp_yp = multigrade[4]
        exp_yq = multigrade[5]

        #print('[sub|y] -> {} {} {} {} {} {}'.format(exp_xp1, exp_xp2, exp_yq, exp_yp, exp_xq1, exp_xq2))

        # check if both exponents are greater than 0
        if exp_yq != 0 and exp_yp != 0:
            #print('[sub|y] -> forwarding to substitution')
            # if we have more yqs
            if exp_yp <= exp_yq:
                #print('[sub|y] -> substitute yp')
                # substitute all yp times yq*qp=N
                exp_N = exp_yp
                monome *= N^exp_N

                # set the new multigrade
                exp_yq = exp_yq - exp_yp
                exp_yp = 0
                #print('[sub|y] N^{} * yq^{} * yp^{}'.format(exp_N, exp_yq, exp_yp))
            elif exp_yp > exp_yq:
                #print('[sub|y] -> substitute yq')
                # substitute all yq times yq*qp=N
                exp_N = exp_yq
                monome *= N^exp_N

                # set the new multigrade
                exp_yp = exp_yp - exp_yq
                exp_yq = 0
                #print('[sub|y] N^{} * yq^{} * yp^{}'.format(exp_N, exp_yq, exp_yp))
        #print('[sub|y] -> {} {} {} {} {} {} after substitution'.format(exp_xp1, exp_xp2, exp_yq, exp_yp, exp_xq1, exp_xq2))

        # set the exponents related to the multigrades to our momone
        monome *= xp1^exp_xp1
        monome *= xp2^exp_xp2
        monome *= yq^exp_yq
        monome *= yp^exp_yp
        monome *= xq1^exp_xq1
        monome *= xq2^exp_xq2

        # add the monome to the polynomial
        output_polynomial += monome

    # return it
    return output_polynomial


def substitute_x(polynomial, e, m):
    """
    Substitute xp1 = xq1 - 1
               xp2 = xq2 + 1
    """
    if polynomial == 0:
        return 0

    # get the polynomial as dict with multigrade => coefficient relation
    polynomial_dict = polynomial.dict()

    # initialize the output polynomial
    output_polynomial = 0

    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lp')

    # iterate through every monomes multigrade
    for multigrade in polynomial_dict:
        # save the coefficient to the monome
        monome = polynomial_dict[multigrade]

        # sort exponents from multigrade
        exp_xp1 = multigrade[0]
        exp_xp2 = multigrade[1]
        exp_xq1 = multigrade[2]
        exp_xq2 = multigrade[3]
        exp_yp = multigrade[4]
        exp_yq = multigrade[5]

        # initial grades for xp and xq
        xp = xp1^exp_xp1 * xp2^exp_xp2
        xq = xq1^exp_xq1 * xq2^exp_xq2

        # check if the monome contains p terms only
        if exp_yq == 0 and exp_yp > 0:
            # if xq1 or xq2 exist in this monome
            if exp_xq1 != 0 or exp_xq2 != 0:
                # replace xq in this monome by xp
                xp = (xp1 + 1)^exp_xq1 * (xp2 - 1)^exp_xq2 * xp1^exp_xp1 * xp2^exp_xp2

                # set xq to 1
                xq = 1
        # check if the monome contains q terms only
        elif exp_yq > 0 and exp_yp == 0:
            # if xp1 or xp2 exist in this monome
            if exp_xp1 != 0 or exp_xp2 != 0:
                # replace xp in this monome by xq
                xq = (xq1 - 1)^exp_xp1 * (xq2 + 1)^exp_xp2 * xq1^exp_xq1 * xq2^exp_xq2

                # set xp to 1
                xp = 1
        else:
            # if xp1 or xp2 exist in this monome
            if exp_xp1 != 0 or exp_xp2 != 0:
                # replace xp in this monome by xq
                xq = (xq1 - 1)^exp_xp1 * (xq2 + 1)^exp_xp2 * xq1^exp_xq1 * xq2^exp_xq2
        
                # set xp to 1
                xp = 1
        
        # build the correct grade
        monome *= xp
        monome *= yq^exp_yq
        monome *= yp^exp_yp
        monome *= xq

        # add monome to the output polynomial
        output_polynomial += monome

    # return it
    return output_polynomial


def substitute_xp(polynomial):
    """
    Substitute xq1 = xp1 + 1
               xq2 = xp2 - 1
    """
    # get the polynomial as dict with multigrade => coefficient relation
    polynomial_dict = polynomial.dict()

    # initialize the output polynomial
    output_polynomial = 0

    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lp')

    # iterate through every monomes multigrade
    for multigrade in polynomial_dict:
        # save the coefficient to the monome
        monome = polynomial_dict[multigrade]

        # sort exponents from multigrade
        exp_xp1 = multigrade[0]
        exp_xp2 = multigrade[1]
        exp_xq1 = multigrade[2]
        exp_xq2 = multigrade[3]
        exp_yp = multigrade[4]
        exp_yq = multigrade[5]


        # replace xp in this monome by xq
        xq = (xq1 - 1)^exp_xp1 * (xq2 + 1)^exp_xp2 * xq1^exp_xq1 * xq2^exp_xq2

        # set xp to 1
        xp = 1
        #else:
        
        # build the correct grade
        monome *= xp
        monome *= yq^exp_yq
        monome *= yp^exp_yp
        monome *= xq

        # add monome to the output polynomial
        output_polynomial += monome

    # return it
    return output_polynomial


def substitute_N(matrix, N, e, m, X, Y):
    # calculate e^m for the inverse
    em = e^m

    # get the inverse of N in e^m
    N_inv = inverse_mod(N, em)
    N1_inv = inverse_mod(N - 1, em)

    # iterate through the matrix
    for i in range(0, matrix.nrows()):
        # while the matrix diagonals contains multiple of N
        while (matrix[i, i] % N) == 0:
            # iterate through the columns
            for j in range(i + 1):
                if (matrix[i, j] % N) == 0 and matrix[i, j] != 0:
                    matrix[i, j] = matrix[i, j] / N
                else:
                    matrix[i, j] = matrix[i, j] * N_inv


        # while the matrix diagonals contains multiple of N
        while (matrix[i, i] % (N - 1)) == 0:
            # iterate through the columns
            for j in range(i + 1):
                if (matrix[i, j] % (N - 1)) == 0 and matrix[i, j] != 0:
                    matrix[i, j] = matrix[i, j] / (N - 1)
                else:
                    matrix[i, j] = matrix[i, j] * N1_inv
    # return the matrix
    return matrix
