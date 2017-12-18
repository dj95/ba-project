#!/bin/env sage
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
from sympy import Symbol, Poly, ZZ, QQ, pprint
from sympy.polys.polymatrix import PolyMatrix

def coppersmith():
    """
    Apply the coppersmith-algorithm
    """

    pass


def get_polynoms(N, e):
    # define symbols for equations
    x_p_1, x_q_1 = Symbol('x_p_1'), Symbol('x_q_1')
    x_p_2, x_q_2 = Symbol('x_p_2'), Symbol('x_q_2')
    y_p, y_q = Symbol('y_p'), Symbol('y_q')

    # create the polynoms
    h_1 = Poly(
            (N - 1) * x_p_1 * x_p_2 + x_p_1 + N * x_p_2,
            modulus=e
            )
    h_2 = Poly(
            (N - 1) * x_q_1 * x_q_2 + N * x_q_1 + x_q_2,
            modulus=e
            )
    
    f_p_1 = Poly(
            N + x_p_1 * (N - y_p),
            modulus=e
            )
    f_q_1 = Poly(
            1 + x_q_1 * (y_q -1),
            modulus=e
            )
    f_p_2 = Poly(
            1 + x_p_2 * (y_p -1),
            modulus=e
            )
    f_q_2= Poly(
            N + x_q_2 * (N - y_q),
            modulus=e
            )

    return h_1, h_2, f_p_1, f_q_1, f_p_2, f_q_2


def index_set_x_stages(m, t, x):
    """
    Generate the first index set I_x in a generator with
    the first, second, third or fourth set.
    """
    if x == 1:
        # initial values
        i1, i2, j1, j2 = 0, 0, 0, 0
        u = 0

        # let the generator work until we hit the boundary
        while i1 < floor(m / 2):
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while i2 < floor(m / 2):
                # if j2 is in the boundary, return
                if min(int(floor(m / 2)) - i1, int(floor(m / 2)) - i2) > 0:
                    yield i1, i2, j1, j2, min(int(floor(m / 2)) - i1, int(floor(m / 2)) - i2)

                # calculate the values for the next call
                i2 = i2 + 1
            # calculate the values for the next call
            i1, i2 = i1 + 1, 0
    elif x == 2:
        # initial values
        i1, i2, j1, j2 = 0, 1, 1, 0
        u = 0

        # let the generator work until we hit the boundary
        while i1 < floor(m / 2) - 1:
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while i2 < floor(m / 2):
                # if j2 is in the boundary, return
                if min(int(floor(m / 2)) - i1 - 1, int(floor(m / 2)) - i2) > 0:
                    yield i1, i2, j1, j2, min(int(floor(m / 2)) - i1 - 1, int(floor(m / 2)) - i2)

                # calculate the values for the next call
                i2 = i2 + 1
            # calculate the values for the next call
            i1, i2 = i1 + 1, 0
    elif x == 3:
        # initial values
        i1, i2, j1, j2 = 0, 0, 1, 0
        u = 0

        # let the generator work until we hit the boundary
        while i1 < floor(m / 2) - 1:
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while j1 < (floor(m / 2) - i1):
                # if j2 is in the boundary, return
                if int(int(floor(m / 2)) - i1 - j1) > 0:
                    yield i1, i2, j1, j2, int(int(floor(m / 2)) - i1 - j1)

                # calculate the values for the next call
                j1 = j1 + 1
            # calculate the values for the next call
            i1, j1 = i1 + 1, 0
    elif x == 4:
        # initial values
        i1, i2, j1, j2 = 0, 0, 0, 1
        u = 0

        # let the generator work until we hit the boundary
        while i2 < floor(m / 2) - 1:
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while j2 < (floor(m / 2) - i2):
                # if j2 is in the boundary, return
                if int(int(floor(m / 2)) - i2 - j2) > 0:
                    yield i1, i2, j1, j2, int(int(floor(m / 2)) - i2 - j2)

                # calculate the values for the next call
                j2 = j2 + 1
            # calculate the values for the next call
            i2, j2 = i2 + 1, 0


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
    i1, i2, j2 = 0, 0, 1

    # let the generator work until we hit the boundary
    while i1 < floor(m / 2):
        # iterate through all i2 to the boundary and increase i1 if we hit it
        while i2 < floor(m / 2):
            # if j2 is in the boundary, return
            if int(floor(t * (i1 + i2)) - ceil((i1 + i2) / 2)) > 0:
                yield i1, i2, int(floor(t * (i1 + i2)) - ceil((i1 + i2) / 2))

            # calculate the values for the next call
            i2 = i2 + 1
        # calculate the values for the next call
        i1, i2 = i1 + 1, 0


def index_set_y_q(m, t):
    """
    Generate the third index set I_y,q in a generator.
    """
    # initial values
    i1, i2, j2 = 0, 0, 1

    # let the generator work until we hit the boundary
    while i1 < floor(m / 2):
        # iterate through all i2 to the boundary and increase i1 if we hit it
        while i2 < floor(m / 2):
            # if j2 is in the boundary, return
            if int(floor(t * (i1 + i2)) - floor((i1 + i2) / 2)) > 0:
                yield i1, i2, int(floor(t * (i1 + i2)) - floor((i1 + i2) / 2))

            # calculate the values for the next call
            i2 = i2 + 1
        # calculate the values for the next call
        i1, i2 = i1 + 1, 0


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
    # set up symbols
    x_p_1, x_q_1 = Symbol('x_p_1'), Symbol('x_q_1')
    x_p_2, x_q_2 = Symbol('x_p_2'), Symbol('x_q_2')
    y_p, y_q = Symbol('y_p'), Symbol('y_q')

    # get the polynomials
    h_1, h_2, f_p_1, f_q_1, f_p_2, f_q_2 = get_polynoms(N, e)

    # get the polynoms degree for the matrix
    degree = h_1.degree()

    # initialize the matrix with 0 values
    #lattice = Matrix(degree, degree)

    print('==> Processing index set x')

    coeffs = {}
    count = 0

    I_x = index_set_x(8, 0.75)
    for (i1, i2, j1, j2, u) in I_x:
        #print('  -> (i1, i2, j1, j2, u) = ({}, {}, {}, {}, {})'.format(i1, i2, j1, j2, u))
        
        p = Poly(x_p_1**j1 * x_p_2**j2 * y_q**(int(floor((i1 + i2) / 2))), x_p_1, x_p_2, y_q) * f_p_1.mul_ground(i1) * f_p_2.mul_ground(i2) * h_2.mul_ground(u) * e**(m - (i1 + i2))

        coeffs[count] = {}

        for monom in p.as_dict():
            #print('{} - {}'.format(monom, p.as_dict()[monom]))
            coeffs[count][tupel_to_string(monom)] = p.as_dict()[monom]

        count += 1

    print('==> Processing index set y p')

    #TODO: test indices, check if indices are calculated correctly
    I_y_p = index_set_y_p(8, 0.75)
    for (i1, i2, j1) in I_y_p:
        #print('  -> (i1, i2, j1) = ({}, {}, {})'.format(i1, i2, j1))
        exponent = int(floor((i1 + i2) / 2) + j1)
        
        #TODO: correct function?
        p = Poly(y_p**exponent, y_p) * f_p_1.mul_ground(i1) * f_p_2.mul_ground(i2) * e**(m - (i1 + i2))

        coeffs[count] = {}

        for monom in p.as_dict():
            #print('{}000 - {}'.format(tupel_to_string(monom), p.as_dict()[monom]))
            coeffs[count][tupel_to_string(monom) + '000'] = p.as_dict()[monom]

        count += 1

    print('==> Processing index set y q\n')

    I_y_q = index_set_y_q(8, 0.75)
    for (i1, i2, j2) in I_y_q:
        #print('  -> (i1, i2, j2) = ({}, {}, {})'.format(i1, i2, j2))
        exponent = int(floor((i1 + i2) / 2) + j2)

        p = Poly(y_q**exponent, y_q) * f_q_1.mul_ground(i1) * f_q_2.mul_ground(i2) * e**(m - (i1 + i2))

        coeffs[count] = {}

        for monom in p.as_dict():
            #print('000{} - {}'.format(tupel_to_string(monom), p.as_dict()[monom]))
            coeffs[count]['000' + tupel_to_string(monom)] = p.as_dict()[monom]

        count += 1

    c = 0
    for polynom in sorted(list(coeffs.keys())):
        if len(coeffs[polynom]) > 0:
            #print('{} - {}'.format(polynom, coeffs[polynom]))
            c += 1

    test_matrix = Matrix((int(m)**6), c)

    y = 0
    for polynom in coeffs:
        if len(coeffs[polynom]) > 0:
            for monom in coeffs[polynom]:
                x = int(monom, base=8)
                
                #TODO: rework polynoms with sage instead of sympy
                test_matrix[x, y] = long(coeffs[polynom][monom])
            y += 1

    #print(test_matrix)
    return test_matrix
