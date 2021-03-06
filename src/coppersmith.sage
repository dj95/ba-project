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

from colorama import Fore, Back, Style


def determinant_paper(X, Y, e, m, tau):
    """
    Calculate the determinant according to the formula
    in the paper on page 25.
    """
    # calculate determinant according to the paper
    n  = 0
    sx = 0
    sy = 0
    se = 0

    # initialize the generator for the last index set
    I_x = index_set_x(m, tau)

    # iterate through the last first set
    for (i1, i2, j1, j2, u) in I_x:
        n += 1                              # dimension
                                            # determinant exponents below
        sx += (i1 + i2 + j1 + j2 + 2*u)     # exponent for X
        if (i1 + i2) % 2 == 0:              #
            sy += int(floor((i1 + i2) / 2)) # exponent for Y with even (i1 + i2)
        else:                               #
            sy += int(ceil((i1 + i2) / 2))  # exponent for Y with odd (i1 + i2)
        se += (m - (i1 + i2 + u))           # exponent for e

    # initialize the generator for the last index set
    I_y_p = index_set_y_p(m, tau)

    # iterate through the second index set
    for (i1, i2, j1) in I_y_p:
        n += 1                              # dimension
                                            # determinant exponents below
        sx += (i1 + i2)                     # exponent for X
        sy += (ceil((i1 + i2) / 2) + j1)    # exponent for Y
        se += (m - (i1 + i2))               # exponent for e

    # initialize the generator for the last index set
    I_y_q = index_set_y_q(m, tau)

    # iterate through the last index set
    for (i1, i2, j2) in I_y_q:
        n += 1                              # dimension
                                            # determinant exponents below
        sx += (i1 + i2)                     # exponent for X
        sy += (floor((i1 + i2) / 2) + j2)   # exponent for Y
        se += (m - (i1 + i2))               # exponent for e

    # calculate the bound according to the papers formula
    detB = (X**sx) * (Y**sy) * (e**se)

    # return the determinant
    return detB


def calculate_theoretical_delta(m, tau, N, e):
    # initialize variables
    n = sx = sy = se = 0

    # calculate alpha in e = N**alpha
    alpha = log(e, N)

    # initialize the generator for the last index set
    I_x = index_set_x(m, tau)

    # iterate through the last first set
    for (i1, i2, j1, j2, u) in I_x:
        n += 1                              # dimension
                                            # determinant exponents below
        sx += (i1 + i2 + j1 + j2 + 2*u)     # exponent for X
        if (i1 + i2) % 2 == 0:              #
            sy += int(floor((i1 + i2) / 2)) # exponent for Y with even (i1 + i2)
        else:                               #
            sy += int(ceil((i1 + i2) / 2))  # exponent for Y with odd (i1 + i2)
        se += (m - (i1 + i2 + u))           # exponent for e

    # initialize the generator for the last index set
    I_y_p = index_set_y_p(m, tau)

    # iterate through the second index set
    for (i1, i2, j1) in I_y_p:
        n += 1                              # dimension
                                            # determinant exponents below
        sx += (i1 + i2)                     # exponent for X
        sy += (ceil((i1 + i2) / 2) + j1)    # exponent for Y
        se += (m - (i1 + i2))               # exponent for e

    # initialize the generator for the last index set
    I_y_q = index_set_y_q(m, tau)

    # iterate through the last index set
    for (i1, i2, j2) in I_y_q:
        n += 1                              # dimension
                                            # determinant exponents below
        sx += (i1 + i2)                     # exponent for X
        sy += (floor((i1 + i2) / 2) + j2)   # exponent for Y
        se += (m - (i1 + i2))               # exponent for e

    # calulate explicit delta
    delta = ((((1/2) - alpha)*sx) + (n*m) - ((1/2)*sy) - alpha*se) / sx
    asymptotic_delta = ((((1/2) - 1)*sx) + (n*m) - ((1/2)*sy) - 1*se) / sx

    return delta, asymptotic_delta


def optimize_tau(e, m, X, Y, N, jsonoutput):
    """
    Brute force the optimal tau for the lattice attack.
    Calculate the determinant and check if its smaller
    than e**nm. If not, exit the program.
    """
    # initialize variables
    n = sx = sy = se = 0
    saved_delta = saved_tau = 0
    tau = 0.5

    # calculate alpha in e = N**alpha
    alpha = log(e, N)

    # iterate 0.5 <= tau <= 1
    while tau <= 1:
        # initialize the generator for the last index set
        I_x = index_set_x(m, tau)

        # iterate through the last first set
        for (i1, i2, j1, j2, u) in I_x:
            n += 1                              # dimension
                                                # determinant exponents below
            sx += (i1 + i2 + j1 + j2 + 2*u)     # exponent for X
            if (i1 + i2) % 2 == 0:              #
                sy += int(floor((i1 + i2) / 2)) # exponent for Y with even (i1 + i2)
            else:                               #
                sy += int(ceil((i1 + i2) / 2))  # exponent for Y with odd (i1 + i2)
            se += (m - (i1 + i2 + u))           # exponent for e

        # initialize the generator for the last index set
        I_y_p = index_set_y_p(m, tau)

        # iterate through the second index set
        for (i1, i2, j1) in I_y_p:
            n += 1                              # dimension
                                                # determinant exponents below
            sx += (i1 + i2)                     # exponent for X
            sy += (ceil((i1 + i2) / 2) + j1)    # exponent for Y
            se += (m - (i1 + i2))               # exponent for e

        # initialize the generator for the last index set
        I_y_q = index_set_y_q(m, tau)

        # iterate through the last index set
        for (i1, i2, j2) in I_y_q:
            n += 1                              # dimension
                                                # determinant exponents below
            sx += (i1 + i2)                     # exponent for X
            sy += (floor((i1 + i2) / 2) + j2)   # exponent for Y
            se += (m - (i1 + i2))               # exponent for e

        # calulate explicit delta
        delta = ((((1/2) - alpha)*sx) + (n*m) - ((1/2)*sy) - alpha*se) / sx

        # check if the delta is greater than our save one
        if delta > saved_delta:
            saved_delta = delta
            saved_tau = tau

        # increase tau
        tau += 0.001

    if not jsonoutput:
        pprint('tau = {} (optimized)'.format(tau))
        pprint("upper theoretical bound: {}".format(saved_delta))

    # calculate the bound according to the papers formula
    detB = (X**sx) * (Y**sy) * (e**se)

    # if the determinant doesnt fulfill howgrave graham...
    if not (detB < e**(n*m)) and not jsonoutput:
        # ...exit with an error
        pprint("[ " + Fore.YELLOW + "WARNING" + Fore.RESET + " ] determinant greater than e**nm")

    # otherwise return the optimal tau
    return tau


def generate_lattice(N, e, X, Y, m=8, tau=0.75, debug=False, jsonoutput=False):
    """
    Generate the lattice for the coppersmith variant of
    the given attack.
    """
    # load required functions from other sage files
    load('./equations.sage')
    load('./index_set.sage')
    load('./shiftpolynomials.sage')
    load('./substitute.sage')
    load('./utils.sage')

    # get ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

    # get an optimized tau
    tau = optimize_tau(e, m, X, Y, N, jsonoutput)

    # initial values
    coeffs = {}
    count = 0
    polynomials = []

    # initialize the generator for the last index set
    I_x = index_set_x(m, tau)

    # iterate through the last first set
    for (i1, i2, j1, j2, u) in I_x:
        #print('  -> (i1, i2, j1, j2, u) = ({}, {}, {}, {}, {})'.format(i1, i2, j1, j2, u))

        # calculate the polynomial for the index set
        eq_p = h_eq(i1, i2, j1, j2, u, N, e ,m, X, Y)

        # polynomial from substitution
        p = g(i1, i2, j1, j2, u, N, e, m)

        # check if substitution matches the equations
        if eq_p != p and not jsonoutput:
            pprint('[ ERROR ] polynomial is not equal to equation')

        # avoid rows with 0 only
        if p == 0:
            print(p)
            continue

        polynomials.append((p, 1, (i1, i2, j1, j2, u)))

        coeffs[count] = {}

        # save the multigrade as dict in base-m
        for monom in p.dict():
            coeffs[count][tupel_to_string(monom)] = p.dict()[monom]

            # check if monom grades have the correct bound
            if monom[4] > max(ceil((i1 + i2) / 2), 0):
                print("Error")
            if monom[5] > floor((i1 + i2) / 2):
                print("Error")

        # increase the counter
        count += 1

    # initialize the generator for the last index set
    I_y_p = index_set_y_p(m, tau)

    # iterate through the second index set
    for (i1, i2, j1) in I_y_p:
        #print('  -> (i1, i2, j1) = ({}, {}, {})'.format(i1, i2, j1))
        
        # calculate the polynomial for the index set
        eq_p = g_p(i1, i2, j1, N, e, m, X, Y)
        
        # polynomial from substitution
        p = gp(i1, i2, j1, N, e, m)

        # check if substitution matches the equations
        if eq_p != p and not jsonoutput:
            pprint('[ ERROR ] polynomial is not equal to equation')

        # avoid rows with 0 only
        if p == 0:
            print(p)
            continue

        polynomials.append((p, 2, (i1, i2, j1)))

        coeffs[count] = {}

        # save the multigrade as dict in base-m
        for monom in p.dict():
            coeffs[count][tupel_to_string(monom)] = p.dict()[monom]

            # check if monom grades have the correct bound
            if monom[4] > max(ceil((i1 + i2) / 2) + j1, 0):
                print("Error3")
            if monom[5] > max(((i1 + i2) / 2) - j1, 0):
                print("Error4")

        # increase the counter
        count += 1

    # initialize the generator for the last index set
    I_y_q = index_set_y_q(m, tau)

    # iterate through the last index set
    for (i1, i2, j2) in I_y_q:
        #print('  -> (i1, i2, j2) = ({}, {}, {})'.format(i1, i2, j2))

        # calculate the polynomial for the index set
        eq_p = g_q(i1, i2, j2, N, e ,m, X, Y)
        
        # polynomial from substitution
        p = gq(i1, i2, j2, N, e, m)

        # check if substitution matches the equations
        if eq_p != p and not jsonoutput:
            pprint('[ ERROR ] polynomial is not equal to equation')

        # avoid rows with 0 only
        if p == 0:
            print(p)
            continue

        polynomials.append((p, 3, (i1, i2, j2)))

        coeffs[count] = {}

        # save the multigrade as dict in base-m
        for monom in p.dict():
            # save the multigrade
            coeffs[count][tupel_to_string(monom)] = p.dict()[monom]

            # check if monom grades have the correct bound
            if monom[4] > max(ceil((i1 + i2) / 2) - j2, 0):
                print("Error5")
            if monom[5] > max(((i1 + i2) / 2) + j2, 0):
                print("Error6")

        # increase the counter
        count += 1

    if not jsonoutput:
        pprint('got {} index sets'.format(count))

    col_indice = {}
    col_index = 0
    # get the polynom count for every polynom with at least one monom
    # this prevents rows with 0 only
    for polynom in sorted(list(coeffs.keys())):
        for monom in coeffs[polynom]:
            if monom not in col_indice:
                col_indice[monom] = col_index
                col_index += 1

    # initialize the matrix
    matrix = Matrix(len(coeffs), col_index)

    row_index = []

    # iterate through the polynoms we save before
    for row in range(len(coeffs)):
        row_index.append(polynomials[row][1])
        # check if we have monoms
        if len(coeffs[polynom]) > 0:
            # iterate through every monom of the polynom
            for monom in coeffs[row]:
                # get the index from the monom grade
                col = col_indice[monom]

                # set the depending cell of the matrix to the coefficient
                matrix[row, col] = coeffs[row][monom]

    # calculate the determinant according to the paper
    detB = determinant_paper(X, Y, e, m, tau)

    # return the lattice and its multigrade-column-relation
    return matrix, col_indice, polynomials, row_index, detB, tau
