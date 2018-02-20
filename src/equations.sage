#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


def h_eq(i1, i2, j1, j2, u, N, e ,m, X, Y):
    """
    The first shift polunomial with the substitutions.
    """
    # define the integer ring and its variables
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

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
                        monom *= (xp1)^(i1 + j1 + u - kp1)
                        monom *= (xp2)^(i2 + j2 + u - kp2)
                        monom *= (yp)^(ip - int(floor((i1 + i2) / 2)))

                        # add the monom to the polynom
                        p += monom

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
                        monom *= binomial(ip1 + j1 + u - l, kq1)
                        monom *= binomial(ip + j2 + l - ip1, kq2)

                        # multiply N and e to the coefficient
                        monom *= N^(i1 - ip1 + ip + l)
                        monom *= e^(m - (i1 + i2 + u))

                        # multiply with the correct grades for xp1, xp2, yp
                        monom *= (xq1)^(i1 + j1 + u - kq1)
                        monom *= (xq2)^(i2 + j2 + u - kq2)
                        monom *= (yq)^(int(floor((i1 + i2) / 2)) - ip)

                        # add the monom to the polynom
                        p += monom
    # return it
    return p


def h_u_check(i1, i2, j1, j2, u, N, e, m):
    # define the integer ring and its variables
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

    # initialize the polynom
    p = 0

    for ip in range(int(floor((i1 + i2) / 2)) + 1, i1 + i2 + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            for l in range(0, u + 1):
                # initialize the monom
                monom = 1

                monom *= (-1)^(2*ip1 + i2 - ip + u - l)

                monom *= binomial(i2, ip - ip1)
                monom *= binomial(i1, ip1)
                monom *= binomial(u, l)

                monom *= N^(i1 - ip1 + int(floor((i1 + i2) / 2)) + l)
                monom *= e^(m - (i1 + i2 + u))

                monom *= xp1^(ip1 + j1 + u - l)
                monom *= xp2^(ip - ip1 + j2 + l)
                monom *= (xp1 + 1)^(i1 - ip1 + l)
                monom *= (xp2 - 1)^(i2 - ip + ip1 + u - l)
                monom *= yp^(ip - int(floor((i1 + i2) / 2)))

                p += monom

    for ip in range(0, int(floor((i1 + i2) / 2)) + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            for l in range(0, u + 1):
                # initialize the monom
                monom = 1
    
                monom *= (-1)^(2*ip1 + i2 - ip + u - l)
    
                monom *= binomial(i2, ip - ip1)
                monom *= binomial(i1, ip1)
                monom *= binomial(u, l)
    
                monom *= N^(i1 - ip1 + ip + l)
                monom *= e^(m - (i1 + i2 + u))
    
                monom *= (xq1 - 1)^(ip1 + j1 + u - l)
                monom *= (xq2 + 1)^(ip - ip1 + j2 + l)
                monom *= xq1^(i1 - ip1 + l)
                monom *= xq2^(i2 - ip + ip1 + u - l)
                monom *= yq^(int(floor((i1 + i2) / 2)) - ip)
    
                p += monom

    return p


def g_p_check_before_y_sub(i1, i2, j1, N, e, m):
    # define the integer ring and its variables
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

    # initialize the polynom
    p = 0

    for ip in range(0, i1 + i2 + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            monom = 1

            monom *= (-1)^(2*ip1 + i2 - ip)

            monom *= binomial(i2, ip - ip1)
            monom *= binomial(i1, ip1)

            monom *= N^(i1 - ip1)
            monom *= e^(m - (i1 + i2))

            monom *= xp1^(ip1)
            monom *= xp2^(ip - ip1)
            monom *= xq1^(i1 - ip1)
            monom *= xq2^(i2 - ip + ip1)
            monom *= yp^(ip)
            monom *= yq^(int(floor((i1 + i2) / 2)) - j1)

            p += monom
    return p


def g_p_check(i1, i2, j1, N, e, m, substitute=True):
    # define the integer ring and its variables
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

    # initialize the polynom
    p = 0

    for ip in range(int(floor((i1 + i2) / 2) - j1 + 1), i1 + i2 + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            monom = 1

            monom *= (-1)^(2*ip1 + i2 - ip)

            monom *= binomial(i2, ip - ip1)
            monom *= binomial(i1, ip1)

            #TODO: fehler
            #monom *= N^(i1 - ip1 + int(floor((i1 + i2) / 2)))

            # korrekter Exponent fuer N
            monom *= N^(i1 - ip1 + int(floor((i1 + i2) / 2)) - j1)
            monom *= e^(m - (i1 + i2))

            monom *= xp1^(ip1)
            monom *= xp2^(ip - ip1)

            if substitute:
                monom *= (xp1 + 1)^(i1 - ip1)
                monom *= (xp2 - 1)^(i2 - ip + ip1)
            else:
                monom *= xq1^(i1 - ip1)
                monom *= xq2^(i2 - ip + ip1)

            monom *= yp^(ip - int(floor((i1 + i2) / 2)) + j1)

            p += monom

    for ip in range(0, int(floor((i1 + i2) / 2)) - j1 + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            monom = 1

            monom *= (-1)^(2*ip1 + i2 - ip)

            monom *= binomial(i2, ip - ip1)
            monom *= binomial(i1, ip1)

            monom *= N^(i1 - ip1 + ip)
            monom *= e^(m - (i1 + i2))

            if substitute:
                monom *= (xq1 - 1)^(ip1)
                monom *= (xq2 + 1)^(ip - ip1)
            else:
                monom *= xp1^(ip1)
                monom *= xp2^(ip - ip1)

            monom *= xq1^(i1 - ip1)
            monom *= xq2^(i2 - ip + ip1)

            monom *= yq^(int(floor((i1 + i2) / 2)) - ip - j1)

            p += monom
    return p


def g_p(i1, i2, j1, N, e, m, X, Y): 
    """
    Calculate the second shift polunomial with all substitutions
    like in the proof on page 31/32.
    """
    # define the integer ring and its variables
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

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
                    monom *= N^(i1 - ip1 + int(floor((i1 + i2) / 2)) - j1)
                    monom *= e^(m - (i1 + i2))
        
                    # multiply with the correct grades for xp1, xp2, yp
                    monom *= (xp1)^(i1 - kp1)
                    monom *= (xp2)^(i2 - kp2)
                    monom *= (yp)^(ip - int(floor((i1 + i2) / 2)) + j1)

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
                    monom *= (xq1)^(i1 - kq1)
                    monom *= (xq2)^(i2 - kq2)
                    monom *= (yq)^(int(floor((i1 + i2) / 2)) - j1 - ip)

                    # add the monom to the polynom
                    p += monom
    # return it
    return p


def g_q_check(i1, i2, j2, N, e, m):
    # define the integer ring and its variables
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')

    # initialize the polynom
    p = 0

    for ip in range(int(floor((i1 + i2) / 2) + j2 + 1), i1 + i2 + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            monom = 1

            monom *= (-1)^(2*ip1 + i2 - ip)

            monom *= binomial(i2, ip - ip1)
            monom *= binomial(i1, ip1)

            #TODO: fehler
            #monom *= N^(i1 - ip1 + int(floor((i1 + i2) / 2)))

            # korrekter Exponent fuer N
            monom *= N^(i1 - ip1 + int(floor((i1 + i2) / 2)) + j2)
            monom *= e^(m - (i1 + i2))

            monom *= xp1^(ip1)
            monom *= xp2^(ip - ip1)
            monom *= xq1^(i1 - ip1)
            monom *= xq2^(i2 - ip + ip1)
            monom *= yp^(ip - int(floor((i1 + i2) / 2)) - j2)

            p += monom

    for ip in range(0, int(floor((i1 + i2) / 2)) + j2 + 1):
        for ip1 in range(max(0, ip - i2), min(i1, ip) + 1):
            monom = 1

            monom *= (-1)^(2*ip1 + i2 - ip)

            monom *= binomial(i2, ip - ip1)
            monom *= binomial(i1, ip1)

            monom *= N^(i1 - ip1 + ip)
            monom *= e^(m - (i1 + i2))

            monom *= xp1^(ip1)
            monom *= xp2^(ip - ip1)
            monom *= xq1^(i1 - ip1)
            monom *= xq2^(i2 - ip + ip1)
            monom *= yq^(int(floor((i1 + i2) / 2)) - ip + j2)

            p += monom
    return p


def g_q(i1, i2, j2, N, e, m, X, Y):
    """
    Calculate the third shift polunomial with all substitutions
    like in the proof on page 31/32., just with -j1 substituted
    by + j2 for the third polynomial instead of the second one.
    """
    # define the integer ring and its variables
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='deglex')


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
                    #monom *= N^(i1 - ip1 + int(floor((i1 + i2) / 2)))
                    monom *= N^(i1 - ip1 + int(floor((i1 + i2) / 2)) + j2)
                    monom *= e^(m - (i1 + i2))
        
                    # multiply with the correct grades for xp1, xp2, yp
                    monom *= (xp1)^(i1 - kp1)
                    monom *= (xp2)^(i2 - kp2)
                    monom *= (yp)^(ip - int(floor((i1 + i2) / 2)) - j2)

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
                    monom *= (xq1)^(i1 - kq1)
                    monom *= (xq2)^(i2 - kq2)
                    monom *= (yq)^(int(floor((i1 + i2) / 2)) + j2 - ip)

                    # add the monom to the polynom
                    p += monom
    # return it
    return p
