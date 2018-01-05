#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


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
    # get the bitlength of the exponent
    bitlen = len(bin(e)[2:])
    c = a
    bitlen -= 1

    # iterate through every bit of the exponent
    while bitlen > 0:
        bitlen -= 1
        # square the result
        c = (c * c) % m

        # get the actual bit of the exponent
        mul = ((e >> bitlen) & 0x1)

        # if the bit is 1...
        if mul:
            # ...multiply the result with the number
            c = (c * a) % m
    
    # return the result
    return c


def matrix_to_ones(matrix, N):
    """
    Create a matrix that indicates with 1, -1 and 0 where positive
    and where negative values are.
    """
    # get stats from the old matrix
    cols = matrix.ncols()
    rows = matrix.nrows()

    # initialize new matrix
    output_matrix = []

    counter = 0
    # iterate through every row of the matrixx
    for row in matrix.rows():
        # append an empty row to the output matrix
        output_matrix.append([])

        # iterate through the values of each row
        for i in range(len(row)):
            # if the value is 0...
            if row[i] == 0:
                # ...append a 0
                output_matrix[counter].append(0)
            elif row[i] < 0: # if its less than 0...
                # ...append -1
                output_matrix[counter].append(-1)
            elif row[i] > 0: # if its greater than 0...
                # ...append 1
                output_matrix[counter].append(1)
        counter += 1

    # return it
    return output_matrix
