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


def generate_lattice(N, e, m=8, tau=0.75):
    """
    Generate the lattice for the coppersmith variant of
    the given attack.
    """
    print('==> generate lattice')

    # load required functions from other sage files
    load('./equations.sage')
    load('./index_set.sage')
    load('./utils.sage')

    R.<xp1, xp2, yq, yp, xq1, xq2> = PolynomialRing(ZZ)

    # initial values
    coeffs = {}
    count = 0

    print('==> Processing index set x')

    # initialize the generator for the last index set
    I_x = index_set_x(m, tau)

    # iterate through the last first set
    for (i1, i2, j1, j2, u) in I_x:
        #print('  -> (i1, i2, j1, j2, u) = ({}, {}, {}, {}, {})'.format(i1, i2, j1, j2, u))

        # calculate the polynomial for the index set
        p = h(i1, i2, j1, j2, u, N, e ,m)

        # avoid rows with 0 only
        if p == 0:
            print(p)
            continue

        coeffs[count] = {}

        # save the multigrade as dict in base-m
        for monom in p.dict():
            coeffs[count][tupel_to_string(monom)] = p.dict()[monom]

        # increase the counter
        count += 1

    print('==> Processing index set y p')

    # initialize the generator for the last index set
    I_y_p = index_set_y_p(m, tau)

    # iterate through the second index set
    for (i1, i2, j1) in I_y_p:
        #print('  -> (i1, i2, j1) = ({}, {}, {})'.format(i1, i2, j1))
        
        # calculate the polynomial for the index set
        p = g_p(i1, i2, j1, N, e, m)

        # avoid rows with 0 only
        if p == 0:
            print(p)
            continue

        coeffs[count] = {}

        # save the multigrade as dict in base-m
        for monom in p.dict():
            # sort the multigrade that we have the grades vor xp1 xp2 yp1 yp2 xq1 xq2
            monom_string = tupel_to_string(monom)

            coeffs[count][monom_string] = p.dict()[monom]

        # increase the counter
        count += 1

    print('==> Processing index set y q\n')

    # initialize the generator for the last index set
    I_y_q = index_set_y_q(m, tau)

    # iterate through the last index set
    for (i1, i2, j2) in I_y_q:
        #print('  -> (i1, i2, j2) = ({}, {}, {})'.format(i1, i2, j2))
       
        # calculate the polynomial for the index set
        p = g_q(i1, i2, j2, N, e ,m)

        # avoid rows with 0 only
        if p == 0:
            print(p)
            continue

        coeffs[count] = {}

        # save the multigrade as dict in base-m
        for monom in p.dict():
            # sort the multigrade that we have the grades vor xp1 xp2 yp1 yp2 xq1 xq2
            monom_string = tupel_to_string(monom)

            # save the mulitgrade
            coeffs[count][monom_string] = p.dict()[monom]

        # increase the counter
        count += 1

    print('==> Got {} index sets'.format(count))

    c = 0
    col_indice = {}
    col_index = 0
    # get the polynom count for every polynom with at least one monom
    # this prevents rows with 0 only
    for polynom in sorted(list(coeffs.keys())):
        if len(coeffs[polynom]) > 0:
            c += 1
        for monom in coeffs[polynom]:
            if monom not in col_indice:
                col_indice[monom] = col_index
                col_index += 1

    # initialize the matrix
    matrix = Matrix(c, col_index)

    row = 0
    # iterate through the polynoms we save before
    for polynom in coeffs:
        # check if we have monoms
        if len(coeffs[polynom]) > 0:
            # iterate through every monom of the polynom
            for monom in coeffs[polynom]:
                # get the index from the monom grade
                col = col_indice[monom]
                
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
        # debug
        sys.stdout.write('\r==> Checking column {}/{}'.format(index + 1, cols))
        sys.stdout.flush()
        
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

    # print some stats
    print('\n==> deleting {} columns with 0'.format(len(delete_index)))
    matrix = matrix.delete_columns(delete_index)
    print('==> finished deleting colums')

    # return the lattice and its multigrade-column-relation
    return matrix, col_indice
