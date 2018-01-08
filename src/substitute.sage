#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


def substitute_y(polynomial, N):
    """
    Substitute yp*yq to N
    """
    # get the polynomial as dict with multigrade => coefficient relation
    polynomial_dict = polynomial.dict()

    # initialize the output polynomial
    output_polynomial = 0

    # define the polynomial ring
    R.<xp1, xp2, yq, yp, xq1, xq2> = PolynomialRing(ZZ)

    # iterate through the multigrades
    for multigrade in polynomial_dict:
        # save the coefficient to the monome
        monome = polynomial_dict[multigrade]

        # save the exponents of yq and yp for the substitution
        exp_yq = multigrade[2]
        exp_yp = multigrade[3]         

        # check if both exponents are greater than 0
        if exp_yq != 0 and exp_yp != 0:
            # if we have more yqs
            if exp_yp <= exp_yq:
                # substitute all yp times yq*qp=N
                monome *= exp_yp * N

                # set the new multigrade
                exp_yq = exp_yq - exp_yp
                exp_yp = 0
            elif exp_yp > exp_yq:
                # substitute all yq times yq*qp=N
                monome *= exp_yq * N

                # set the new multigrade
                exp_yp = exp_yp - exp_yq
                exp_yq = 0

        # set the exponents related to the multigrades to our momone
        monome *= xp1^multigrade[0]
        monome *= xp2^multigrade[1]
        monome *= yq^exp_yq
        monome *= yp^exp_yp
        monome *= xq1^multigrade[4]
        monome *= xq2^multigrade[5]

        # add the monome to the polynomial
        output_polynomial += monome

    # return it
    return output_polynomial


def substitute_x(polynomial):
    """
    Substitute xp1 = xq1 - 1
               xp2 = xq2 + 1
    """
    # get the polynomial as dict with multigrade => coefficient relation
    polynomial_dict = polynomial.dict()

    # initialize the output polynomial
    output_polynomial = 0

    # define the polynomial ring
    R.<xp1, xp2, yq, yp, xq1, xq2> = PolynomialRing(ZZ)

    # iterate through every monomes multigrade
    for multigrade in polynomial_dict:
        # save the coefficient to the monome
        monome = polynomial_dict[multigrade]

        # sort exponents from multigrade
        exp_xp1 = multigrade[0]
        exp_xp2 = multigrade[1]
        exp_yq = multigrade[2]
        exp_yp = multigrade[3]
        exp_xq1 = multigrade[4]
        exp_xq2 = multigrade[5]

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

        # build the correct grade
        monome *= xp
        monome *= yq^exp_yq
        monome *= yp^exp_yp
        monome *= xq

        # add monome to the output polynomial
        output_polynomial += monome

    # return it
    return output_polynomial