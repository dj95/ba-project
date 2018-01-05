#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2017 - Daniel Jankowski


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


def sqm(a, e, m):
    """
    Calculate a**e mod m with square and multiply.
    """
    bitlen = len(bin(e)[2:])
    c = a
    bitlen -= 1

    while bitlen > 0:
        bitlen -= 1
        c = (c * c) % m
        mul = ((e >> bitlen) & 0x1)
        if mul:
            c = (c * a) % m
    return c
