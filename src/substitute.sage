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
