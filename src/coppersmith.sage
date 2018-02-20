#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


def optimize_tau(e, m, X, Y, N):
    # calculate determinant according to the paper
    n  = 0
    sx = 0
    sy = 0
    se = 0

    alpha = log(e, N)
    print(alpha)
    tau = 0.5
    saved_delta = 0
    saved_tau = 0

    while tau <= 1:
        # initialize the generator for the last index set
        I_x = index_set_x(m, tau)

        # iterate through the last first set
        for (i1, i2, j1, j2, u) in I_x:
            n += 1
            sx += (i1 + i2 + j1 + j2 + 2*u) 
            if (i1 + i2) % 2 == 0:
                sy += int(floor((i1 + i2) / 2))
            else:
                sy += int(ceil((i1 + i2) / 2))
            se += (m - (i1 + i2 + u))

        # initialize the generator for the last index set
        I_y_p = index_set_y_p(m, tau)

        # iterate through the second index set
        for (i1, i2, j1) in I_y_p:
            n += 1
            sx += (i1 + i2)
            sy += (ceil((i1 + i2) / 2) + j1)
            se += (m - (i1 + i2))

        # initialize the generator for the last index set
        I_y_q = index_set_y_q(m, tau)

        # iterate through the last index set
        for (i1, i2, j2) in I_y_q:
            n += 1
            sx += (i1 + i2)
            sy += (floor((i1 + i2) / 2) + j2)
            se += (m - (i1 + i2))

        # (alpha + delta - 1/2)sx + 1/2sy + alphase = nm
        delta = ((((1/2) - alpha)*sx) + (n*m) - ((1/2)*sy) - alpha*se) / sx
        print(delta)
        if delta > saved_delta:
            saved_delta = delta
            saved_tau = tau
        tau += 0.001

    detB = (X^sx) * (Y^sy) * (e^se)
    print(detB < e^(n*m))
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

    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lp')

    tau = optimize_tau(e, m, X, Y, N)

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

    # print some stats
    #matrix = matrix.delete_columns(delete_index)

    # calculate determinant according to the paper
    sx = 0
    sy = 0
    se = 0

    # initialize the generator for the last index set
    I_x = index_set_x(m, tau)

    # iterate through the last first set
    for (i1, i2, j1, j2, u) in I_x:
        sx += (i1 + i2 + j1 + j2 + 2*u) 
        if (i1 + i2) % 2 == 0:
            sy += int(floor((i1 + i2) / 2))
        else:
            sy += int(ceil((i1 + i2) / 2))
        se += (m - (i1 + i2 + u))

    # initialize the generator for the last index set
    I_y_p = index_set_y_p(m, tau)

    # iterate through the second index set
    for (i1, i2, j1) in I_y_p:
        sx += (i1 + i2)
        sy += (ceil((i1 + i2) / 2) + j1)
        se += (m - (i1 + i2))

    # initialize the generator for the last index set
    I_y_q = index_set_y_q(m, tau)

    # iterate through the last index set
    for (i1, i2, j2) in I_y_q:
        sx += (i1 + i2)
        sy += (floor((i1 + i2) / 2) + j2)
        se += (m - (i1 + i2))

    detB = (X^sx) * (Y^sy) * (e^se)
    print(detB < e^(m*matrix.ncols()))

    # return the lattice and its multigrade-column-relation
    return matrix, col_indice, polynomials, row_index, detB, sx, sy, se
