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

from numpy import floor, ceil
from sympy import Symbol, Poly, QQ, pprint


def generate_lattice(N, e, m=8, tau=0.75, debug=False):
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

    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lex')

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
        p = h_eq(i1, i2, j1, j2, u, N, e ,m)

        #print(p)
        #p = h_u_check(i1, i2, j1, j2, u, N, e ,m)
        #p = g(i1, i2, j1, j2, u, N, e, m)

        new_red_p = substitute_N(p, N)
        while True:
            red_p = substitute_N(new_red_p, N)

            if new_red_p != red_p:
                new_red_p = red_p
            else:
                break
        p = new_red_p

        # avoid rows with 0 only
        if p == 0:
            print(p)
            continue

        polynomials.append((p, 1, (i1, i2, j1, j2, u)))

        coeffs[count] = {}

        # save the multigrade as dict in base-m
        for monom in p.dict():
            coeffs[count][tupel_to_string(monom)] = p.dict()[monom]

        # increase the counter
        count += 1

    # initialize the generator for the last index set
    I_y_p = index_set_y_p(m, tau)

    # iterate through the second index set
    for (i1, i2, j1) in I_y_p:
        #print('  -> (i1, i2, j1) = ({}, {}, {})'.format(i1, i2, j1))
        
        # calculate the polynomial for the index set
        p = g_p(i1, i2, j1, N, e, m)
        #print(p)
        
        #p = g_p_check(i1, i2, j1, N, e, m, True)

        #p = gp(i1, i2, j1, N, e, m)
        #print(p)
        #print('')

        new_red_p = substitute_N(p, N)
        while True:
            red_p = substitute_N(new_red_p, N)

            if new_red_p != red_p:
                new_red_p = red_p
            else:
                break
        p = new_red_p

        # avoid rows with 0 only
        if p == 0:
            print(p)
            continue

        polynomials.append((p, 2, (i1, i2, j1)))

        coeffs[count] = {}

        # save the multigrade as dict in base-m
        for monom in p.dict():
            coeffs[count][tupel_to_string(monom)] = p.dict()[monom]

        # increase the counter
        count += 1

    # initialize the generator for the last index set
    I_y_q = index_set_y_q(m, tau)

    # iterate through the last index set
    for (i1, i2, j2) in I_y_q:
        #print('  -> (i1, i2, j2) = ({}, {}, {})'.format(i1, i2, j2))

        # calculate the polynomial for the index set
        p = g_q(i1, i2, j2, N, e ,m)
        #print(p)
        
        #p = g_q_check(i1, i2, j2, N, e, m)
        #p = gq(i1, i2, j2, N, e, m)
        #print(p)

        new_red_p = substitute_N(p, N)
        while True:
            red_p = substitute_N(new_red_p, N)

            if new_red_p != red_p:
                new_red_p = red_p
            else:
                break
        p = new_red_p

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

        # increase the counter
        count += 1

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

    row = 0
    # iterate through the polynoms we save before
    for polynom in coeffs:
        # check if we have monoms
        if len(coeffs[polynom]) > 0:
            # iterate through every monom of the polynom
            for monom in coeffs[polynom]:
                # get the index from the monom grade
                col = col_indice[monom]

                #TODO: check shift polynomials dependencies
                
                # set the depending cell of the matrix to the coefficient
                matrix[row, col] = long(coeffs[polynom][monom])
            row += 1

    # eliminate cols with 0
    cols = matrix.ncols()
    rows = matrix.nrows()
    delete_index = []
    index = 0

    # iterate through columns
    for col in matrix.columns():
        # iterate through every row of the column
        for i in range(rows):
            # if we have no 0, skip to the next column
            if col[i] != 0:
                break
            # if we have the last rows value in one column and its 0
            if col[i] == 0 and i == (rows - 1):
                # delete the column because it contains 0 only
                delete_index.append(index)
        index += 1

    if len(delete_index) == 0:
        pprint("column check            [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
    else:
        pprint("column check            [" + Fore.RED + " failed " + Fore.RESET + "]") 

    # print some stats
    matrix = matrix.delete_columns(delete_index)

    # return the lattice and its multigrade-column-relation
    return matrix, col_indice, polynomials
