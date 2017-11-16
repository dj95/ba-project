#!/bin/env python3
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2017 - Daniel Jankowski

import sys 


from matrix import Matrix
from numpy import *
from sympy import Symbol, Poly, ZZ, QQ, pprint
from sympy.polys.polymatrix import PolyMatrix


def coppersmith():
    """
    Apply the coppersmith-algorithm
    """

    pass


def index_set_y_p(m, t):
    """
    Generate the second index set I_y,p in a generator.
    """
    # initial values
    i1, i2, j1 = 0, 0, 1

    # let the generator work until we hit the boundary
    while i1 < floor(m / 2):
        # if j1 is in theboundary, return it, too
        if j1 <= int(t * (floor(m / 2) * 2)) - int(ceil((floor(m / 2) * 2) / 2)):
            yield i1, i2, j1
        else: # else return 0 for j1
            yield i1, i2, 0

        # calculate the values for the next call
        i1, i2, j1 = i1 + 1, i2 + 1, j1 + 1


def index_set_y_q(m, t):
    """
    Generate the third index set I_y,q in a generator.
    """
    # initial values
    i1, i2, j2 = 0, 0, 1

    # let the generator work until we hit the boundary
    while i1 < floor(m / 2):
        # if j2 is in theboundary, return it, too
        if j2 <= int(t * (floor(m / 2) * 2)) - int(floor((floor(m / 2) * 2) / 2)):
            yield i1, i2, j2
        else: # elxe return 0 for j2
            yield i1, i2, 0

        # calculate the values for the next call
        i1, i2, j2 = i1 + 1, i2 + 1, j2 + 1


def generate_lattice(N, e, m):
    """
    Generate the lattice for the coppersmith variant of
    the given attack.
    """
    print('==> generate lattice')

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
            1 + x_p_1 * (y_p -1),
            modulus=e
            )
    f_q_2= Poly(
            N + x_q_2 * (N - y_q),
            modulus=e
            )

    # get the polynoms degree for the matrix
    degree = h_1.degree()

    # initialize the matrix with 0 values
    lattice = Matrix(degree, degree)

    #TODO: test indices, check if indices are calculated correctly
    I_y_p = index_set_y_p(8, 0.75)
    for (i1, i2, j1) in I_y_p:
        print('  -> (i1, i2, j1) = ({}, {}, {})'.format(i1, i2, j1))
        exponent = int(floor((i1 + i2) / 2) - j1)
        #TODO: correct function?
        p = Poly(y_q) * f_p_1.mul_ground(i1) * f_p_2.mul_ground(i2) * e**(m - (i1 + i2))
        pprint(p) 
        print()

    I_y_q = index_set_y_q(8, 0.75)
    for (i1, i2, j2) in I_y_q:
        print('  -> (i1, i2, j1) = ({}, {}, {})'.format(i1, i2, j2))
        exponent = int(floor((i1 + i2) / 2) + j2)
        p = Poly(y_q**exponent, y_q) * f_p_1.mul_ground(i1) * f_p_2.mul_ground(i2) * e**(m - (i1 + i2))
        pprint(p)
        pprint(p.as_dict()) 
        print()
