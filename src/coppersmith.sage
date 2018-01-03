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


def index_set_x_stages(m, t, x):
    """
    Generate the first index set I_x in a generator with
    the first, second, third or fourth subset.
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
    """
    Generate all index sets of the different subsets of index set I_x
    with one generator expression.
    """
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
    """
    Concatenate all values of a tupel to one string.
    """
    output = ''

    # iterate through all values
    for value in tupel:
        # concatenate the value to the output string
        output = '{}{}'.format(output, value)

    # return it
    return output


def h_odd(i1, i2, j1, j2, u, N, e ,m):
    # define the integer ring and its variables
    R.<xp1, xp2, yq, yp, xq1, xq2> = PolynomialRing(ZZ)

    # initialize the polynom
    p = 0

    # for loops simulate the sums
    for ip in range(int(floor((i1 + i2) / 2)) + 1, i1 + i2 + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            for l in range(0, u + 1):
                for kp1 in range(0, i1 - ip1 + l + 1):
                    for kp2 in range(0, i2 - ip + ip1 + u - l + 1):
                        # initialize the monom
                        monom = 1

                        # set the correct sign
                        monom *= (-1)^(2*ip1 + i2 - ip + u - l + kp2)

                        # multiply all of the binomial coefficients
                        monom *= binomial(i2, ip - ip1)
                        monom *= binomial(i1, ip1)
                        monom *= binomial(u, l)
                        monom *= binomial(i1 - ip1 + l, kp1)
                        monom *= binomial(i2 - ip + ip1 + u - l, kp2)

                        # multiply N and e to the coefficient
                        monom *= N^(i1 - ip1 + int(floor((i1 + i2) / 2)) + l)
                        monom *= e^(m - (i1 + i2 + u))

                        # multiply with the correct grades for xp1, xp2, yp
                        monom *= xp1^(i1 + j1 + u - kp1)
                        monom *= xp2^(i2 + j2 + u - kp2)
                        monom *= yp^(ip - int(floor((i1 + i2) / 2)))

                        # add the monom to the polynom
                        p += monom
    # return it
    return p


def h_even(i1, i2, j1, j2, u, N, e ,m):
    # define the integer ring and its variables
    R.<xp1, xp2, yq, yp, xq1, xq2> = PolynomialRing(ZZ)

    # initialize the polynom
    p = 0

    # for loops simulate the sums
    for ip in range(0, int(floor((i1 + i2) / 2)) + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            for l in range(0, u + 1):
                for kq1 in range(0, ip1 + j1 + u - l + 1):
                    for kq2 in range(0, ip - ip1 + j2 + l + 1):
                        # initialize the monom
                        monom = 1

                        # set the correct sign
                        monom *= (-1)^(2*ip1 + i2 - ip + u - l + kq1)

                        # multiply all of the binomial coefficients
                        monom *= binomial(i2, ip - ip1)
                        monom *= binomial(i1, ip1)
                        monom *= binomial(u, l)
                        monom *= binomial(ip1 + j1 + u - l, kq2)
                        monom *= binomial(ip + j2 + l - ip1, kq2)

                        # multiply N and e to the coefficient
                        monom *= N^(i1 - ip1 + ip + l)
                        monom *= e^(m - (i1 + i2 + u))

                        # multiply with the correct grades for xp1, xp2, yp
                        monom *= xq1^(i1 + j1 + u - kq1)
                        monom *= xq2^(i2 + j2 + u - kq2)
                        monom *= yq^(int(floor((i1 + i2) / 2)) - ip)

                        # add the monom to the polynom
                        p += monom
    # return it
    return p


def g_p(i1, i2, j1, N, e, m):
    """
    Calculate the second shift polunomial with all substitutions
    like in the proof on page 31/32.
    """
    # define the integer ring and its variables
    R.<xp1, xp2, yq, yp, xq1, xq2> = PolynomialRing(ZZ)

    # initialize the polynom
    p = 0
    
    # for loops simulate the sums
    for ip in range(int(floor((i1 + i2) / 2) - j1 + 1), i1 + i2 + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            for kp1 in range(0, i1 - ip1 + 1):
                for kp2 in range(0, i2 - ip + ip1 + 1):
                    # initialize the monom
                    monom = 1

                    # set the correct sign
                    monom *= (-1)^(2*ip1 + i2 - ip + kp2)

                    # multiply all of the binomial coefficients
                    monom *= binomial(i2, ip - ip1)
                    monom *= binomial(i1, ip1)
                    monom *= binomial(i1 - ip1, kp1)
                    monom *= binomial(i2 - ip + ip1, kp2)

                    # multiply N and e to the coefficient
                    monom *= N^(i1 - ip1 + int(floor((i1 + i2) / 2)))
                    monom *= e^(m - (i1 + i2))
        
                    # multiply with the correct grades for xp1, xp2, yp
                    monom *= xp1^(i1 - kp1)
                    monom *= xp2^(i2 - kp2)
                    monom *= yp^(ip - int(floor((i1 + i2) / 2)) + j1)

                    # add the monom to the polynom
                    p += monom

    # for loops simulate the sums
    for ip in range(0, int(floor((i1 + i2) / 2)) - j1 + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            for kq1 in range(0, ip1 + 1):
                for kq2 in range(0, ip - ip1 + 1):
                    # initialize the monom
                    monom = 1

                    # set the correct sign
                    monom *= (-1)^(2*ip1 + i2 - ip + kq1)

                    # multiply all of the binomial coefficients
                    monom *= binomial(i2, ip - ip1)
                    monom *= binomial(i1, ip1)
                    monom *= binomial(ip1, kq1)
                    monom *= binomial(ip - ip1, kq2)

                    # multiply N and e to the coefficient
                    monom *= N^(i1 - ip1 + ip)
                    monom *= e^(m - (i1 + i2))
        
                    # multiply with the correct grades for xp1, xp2, yp
                    monom *= xq1^(i1 - kq1)
                    monom *= xq2^(i2 - kq2)
                    monom *= yq^(int(floor((i1 + i2) / 2)) - j1 - ip)

                    # add the monom to the polynom
                    p += monom
    # return it
    return p


def g_q(i1, i2, j2, N, e, m):
    """
    Calculate the third shift polunomial with all substitutions
    like in the proof on page 31/32., just with -j1 substituted
    by + j2 for the third polynomial instead of the second one.
    """
    # define the integer ring and its variables
    R.<xp1, xp2, yq, yp, xq1, xq2> = PolynomialRing(ZZ)

    # initialize the polynom
    p = 0
    
    # for loops simulate the sums
    for ip in range(int(floor((i1 + i2) / 2) + j2 + 1), i1 + i2 + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            for kp1 in range(0, i1 - ip1 + 1):
                for kp2 in range(0, i2 - ip + ip1 + 1):
                    # initialize the monom
                    monom = 1

                    # set the correct sign
                    monom *= (-1)^(2*ip1 + i2 - ip + kp2)

                    # multiply all of the binomial coefficients
                    monom *= binomial(i2, ip - ip1)
                    monom *= binomial(i1, ip1)
                    monom *= binomial(i1 - ip1, kp1)
                    monom *= binomial(i2 - ip + ip1, kp2)

                    # multiply N and e to the coefficient
                    monom *= N^(i1 - ip1 + int(floor((i1 + i2) / 2)))
                    monom *= e^(m - (i1 + i2))
        
                    # multiply with the correct grades for xp1, xp2, yp
                    monom *= xp1^(i1 - kp1)
                    monom *= xp2^(i2 - kp2)
                    monom *= yp^(ip - int(floor((i1 + i2) / 2)) - j2)

                    # add the monom to the polynom
                    p += monom

    # for loops simulate the sums
    for ip in range(0, int(floor((i1 + i2) / 2)) + j2 + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            for kq1 in range(0, ip1 + 1):
                for kq2 in range(0, ip - ip1 + 1):
                    # initialize the monom
                    monom = 1

                    # set the correct sign
                    monom *= (-1)^(2*ip1 + i2 - ip + kq1)

                    # multiply all of the binomial coefficients
                    monom *= binomial(i2, ip - ip1)
                    monom *= binomial(i1, ip1)
                    monom *= binomial(ip1, kq1)
                    monom *= binomial(ip - ip1, kq2)

                    # multiply N and e to the coefficient
                    monom *= N^(i1 - ip1 + ip)
                    monom *= e^(m - (i1 + i2))
        
                    # multiply with the correct grades for xp1, xp2, yp
                    monom *= xq1^(i1 - kq1)
                    monom *= xq2^(i2 - kq2)
                    monom *= yq^(int(floor((i1 + i2) / 2)) + j2 - ip)

                    # add the monom to the polynom
                    p += monom
    # return it
    return p


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

        # calculate the polynomial for the index set
        p = h_odd(i1, i2, j1, j2, u, N, e ,m)
        p += h_even(i1, i2, j1, j2, u, N, e ,m)

        # avoid rows with 0 only
        if p == 0:
            continue

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
        
        # calculate the polynomial for the index set
        p = g_p(i1, i2, j1, N, e, m)

        # avoid rows with 0 only
        if p == 0:
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
    I_y_q = index_set_y_q(8, 0.75)

    # iterate through the last index set
    for (i1, i2, j2) in I_y_q:
        #print('  -> (i1, i2, j2) = ({}, {}, {})'.format(i1, i2, j2))
       
        # calculate the polynomial for the index set
        p = g_q(i1, i2, j2, N, e ,m)

        # avoid rows with 0 only
        if p == 0:
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
