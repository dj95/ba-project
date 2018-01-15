#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


def h(N):
    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lex')

    #TODO: validate by hand
    polynomial = (N*xp2*xq1 - xp1*xq2)

    return polynomial


def fp1(N):
    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lex')

    polynomial = N + xp1*(N - yp)

    return polynomial


def fp2(N):
    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lex')

    polynomial = 1 + xp2*(yp - 1)

    return polynomial


def fq1(N):
    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lex')

    polynomial = 1 + xq1*(yq - 1)

    return polynomial


def fq2(N):
    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lex')

    polynomial = N + xq2*(N - yq)

    return polynomial


def fp1_sub_xp1(N):
    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lex')

    polynomial = (N*xq1 - xp1*yp)

    return polynomial


def fp2_sub_xp2(N):
    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lex')

    polynomial = (-xq2 + xp2*yp)

    return polynomial
