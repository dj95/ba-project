#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


from colorama import Fore, Back, Style


def pprint(text):
    """
    Pretty print the text with colors
    """
    print(Fore.BLUE + '==>' + Fore.RESET + ' ' + text)


def root_check(polynomials, keys, debug):
    """
    Check if every polynomial in polynomials has the roots at
    the correct places in mod e
    """
    # initialize helper variables
    polynomial_count, correct_count = len(polynomials), 0

    # iterate through polynomials
    for p in polynomials:
        # calculate the y value with key parameters
        y = p(
            xp1=keys['kq'] - 1,
            xp2=keys['kp'],
            xq1=keys['kq'],
            xq2=keys['kq'] - 1,
            yp=keys['p'],
            yq=keys['q']
            ) % keys['e']

        # if we found a root, raise the counter
        if y == 0:
            correct_count += 1
        else:
            if debug:
                print(p)

    # if all polynomials share the roots
    if correct_count == polynomial_count:
        # return true
        return True

    # log output
    pprint('{} wrong roots'.format(polynomial_count - correct_count))

    # else return false
    return False


def matrix_sort_stairs(matrix):
    """
    Sort the matrix rows that they form steps or a triangle.
    """
    # initialize some variables for work
    work_matrix = []
    row_index = {}

    counter = -1
    # convert matrix object to array matrix
    for row in matrix:
        # append a new row for every row
        work_matrix.append([])
        counter += 1

        # append the values to each row
        for value in row:
            work_matrix[counter].append(value)

    # append an index to each row
    for row in work_matrix:
        index = 0

        # iterate through the columns
        for i in range(len(row)):
            # set the index to the column count of the value in order to save
            # the highest column count of an existing monom
            if row[i] != 0:
                index = i

        # append the index to each row
        row.append(index)

    sorted_matrix = []

    # sort the matrix with the index
    # iterate through all possible indices
    for i in range(len(work_matrix[0])):
        # iterate through all columns of the matrix
        for row in work_matrix:
            # check if the index is i
            if row[-1:][0] == i:
                # remove the index and append the row
                row = row[:-1]
                sorted_matrix.append(row)

    # return it
    return sorted_matrix


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
