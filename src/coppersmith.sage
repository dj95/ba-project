#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2017 - Daniel Jankowski

import sys 

from numpy import floor, ceil
from sympy import Symbol, Poly, QQ, pprint
from sympy.polys.polymatrix import PolyMatrix

def coppersmith():
    """
    Apply the coppersmith-algorithm
    """

    pass


def get_polynoms(N, e):
    # define symbols for equations
    R.<xp1, xp2, xq1, xq2, yq, yp> = PolynomialRing(ZZ, order='lex')

    # create the polynoms
    h1 = (N - 1) * xp1 * xp2 + xp1 + N * xp2
    h2 = (N - 1) * xq1 * xq2 + xq1 + N * xq2
    fp1 = N + xp1 * (N - yp)
    fp2 = 1 + xp2 * (yp - 1)
    fq1 = 1 + xq1 * (yq - 1)
    fq2 = N + xq2 * (N - yq)

    return h1, h2, fp1, fq1, fp2, fq2


def index_set_x_stages(m, t, x):
    """
    Generate the first index set I_x in a generator with
    the first, second, third or fourth set.
    """
    if x == 1:
        # initial values
        i1, i2, j1, j2 = 0, 0, 0, 0
        u = -1

        # let the generator work until we hit the boundary
        while i1 <= floor(m / 2):
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while i2 <= floor(m / 2):
                # if j2 is in the boundary, return
                while u <= (min(int(floor(m / 2)) - i1, int(floor(m / 2)) - i2) - 1):
                    u = u + 1
                    yield i1, i2, j1, j2, u

                # calculate the values for the next call
                i2, u = i2 + 1, -1
            # calculate the values for the next call
            i1, i2, u = i1 + 1, 0, -1
    elif x == 2:
        # initial values
        i1, i2, j1, j2 = 0, 1, 1, 0
        u = -1

        # let the generator work until we hit the boundary
        while i1 <= floor(m / 2) - 1:
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while i2 <= floor(m / 2):
                # if j2 is in the boundary, return
                while u <= (min(int(floor(m / 2)) - i1 - 1, int(floor(m / 2)) - i2) - 1):
                    u = u + 1
                    yield i1, i2, j1, j2, u

                # calculate the values for the next call
                i2, u = i2 + 1, -1
            # calculate the values for the next call
            i1, i2, u = i1 + 1, 0, -1
    elif x == 3:
        # initial values
        i1, i2, j1, j2 = 0, 0, 1, 0
        u = -1

        # let the generator work until we hit the boundary
        while i1 <= floor(m / 2) - 1:
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while j1 <= (floor(m / 2) - i1):
                # if j2 is in the boundary, return
                while u <= (int(int(floor(m / 2)) - i1 - j1) - 1):
                    u = u + 1
                    yield i1, i2, j1, j2, u

                # calculate the values for the next call
                j1, u = j1 + 1, -1
            # calculate the values for the next call
            i1, j1, u = i1 + 1, 0, -1
    elif x == 4:
        # initial values
        i1, i2, j1, j2 = 0, 0, 0, 1
        u = -1

        # let the generator work until we hit the boundary
        while i2 <= floor(m / 2) - 1:
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while j2 <= (floor(m / 2) - i2):
                # if j2 is in the boundary, return
                while u <= (int(int(floor(m / 2)) - i2 - j2) - 1):
                    u = u + 1
                    yield i1, i2, j1, j2, u

                # calculate the values for the next call
                j2, u = j2 + 1, -1
            # calculate the values for the next call
            i2, j2, u = i2 + 1, 0, -1


def index_set_x(m, t):
    for set in index_set_x_stages(m, t, 1):
        yield set
    for set in index_set_x_stages(m, t, 2):
        yield set
    for set in index_set_x_stages(m, t, 3):
        yield set
    for set in index_set_x_stages(m, t, 4):
        yield set


def index_set_y_p(m, t):
    """
    Generate the second index set I_y,p in a generator.
    """
    # initial values
    i1, i2 = 0, 0
    j1 = 0

    # let the generator work until we hit the boundary
    while i1 <= floor(m / 2):
        # iterate through all i2 to the boundary and increase i1 if we hit it
        while i2 <= floor(m / 2):
            # if j2 is in the boundary, return
            while j1 <= (int(floor(t * (i1 + i2)) - ceil((i1 + i2) / 2)) - 1):
                j1 = j1 + 1
                yield i1, i2, j1

            # calculate the values for the next call
            i2, j1 = i2 + 1, 0
        # calculate the values for the next call
        i1, i2, j1 = i1 + 1, 0, 0


def index_set_y_q(m, t):
    """
    Generate the third index set I_y,q in a generator.
    """
    # initial values
    i1, i2 = 0, 0
    j2 = 0

    # let the generator work until we hit the boundary
    while i1 <= floor(m / 2):
        # iterate through all i2 to the boundary and increase i1 if we hit it
        while i2 <= floor(m / 2):
            # if j2 is in the boundary, return
            while j2 <= (int(floor(t * (i1 + i2)) - floor((i1 + i2) / 2)) - 1):
                j2 = j2 + 1
                yield i1, i2, j2

            # calculate the values for the next call
            i2, j2 = i2 + 1, 0
        # calculate the values for the next call
        i1, i2, j2 = i1 + 1, 0, 0


def tupel_to_string(tupel):
    output = ''

    for value in tupel:
        output = '{}{}'.format(output, value)

    return output


def generate_lattice(N, e, m):
    """
    Generate the lattice for the coppersmith variant of
    the given attack.
    """
    print('==> generate lattice')

    R.<xp1, xp2, yq, yp, xq1, xq2> = PolynomialRing(ZZ)

    # get the polynomials
    h1, h2, fp1, fq1, fp2, fq2 = get_polynoms(N, e)

    # initial values
    coeffs = {}
    count = 0

    print('==> Processing index set x')

    # initialize the generator for the last index set
    I_x = index_set_x(8, 0.75)

    # iterate through the last first set
    for (i1, i2, j1, j2, u) in I_x:
        #print('  -> (i1, i2, j1, j2, u) = ({}, {}, {}, {}, {})'.format(i1, i2, j1, j2, u))

        # insert the indices into the polynom
        p = xp1**j1 * xp2**j2 * yq**(int(floor((i1 + i2) / 2))) * fp1^i1 * fp2^i2 * h2^u * e^(m - (i1 + i2))

        coeffs[count] = {}

        # save the multigrade as dict in base-m
        for monom in p.dict():
            coeffs[count][tupel_to_string(monom)] = p.dict()[monom]

        # increase the counter
        count += 1

    print('==> Processing index set y p')

    # initialize the generator for the last index set
    I_y_p = index_set_y_p(8, 0.75)

    # iterate through the second index set
    for (i1, i2, j1) in I_y_p:
        #print('  -> (i1, i2, j1) = ({}, {}, {})'.format(i1, i2, j1))
        
        # calculate the exponent for y_q
        exponent = int(floor((i1 + i2) / 2) + j1)
        
        # insert the indices into the polynom
        p = yp**exponent * fp1**i1 * fp2**i2 * e**(m - (i1 + i2))

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
    I_y_q = index_set_y_q(8, 0.75)

    # iterate through the last index set
    for (i1, i2, j2) in I_y_q:
        #print('  -> (i1, i2, j2) = ({}, {}, {})'.format(i1, i2, j2))

        # calculate the exponent for y_q
        exponent = int(floor((i1 + i2) / 2) + j2)

        # insert the indices into the polynom
        p = yq**exponent * fq1**i1 * fq2**i2 * e**(m - (i1 + i2))

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
                
                #TODO: rework polynoms with sage instead of sympy
                # set the depending cell of the matrix to the coefficient
                matrix[row, col] = long(coeffs[polynom][monom])
            row += 1

    # eleminate cols and rows with 0
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

    print('\n==> deleting {} columns with 0'.format(len(delete_index)))
    matrix = matrix.delete_columns(delete_index)
    print('==> finished deleting colums')

    #print(matrix)
    return matrix, col_indice
