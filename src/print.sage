#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


def print_matrix(matrix):
    """
    Print the matrix in TeX-format to a file.
    """
    # start the matrix in tex
    tex_string = '\\begin{bmatrix}\n'

    # iterate through every row of the matrix
    for row in matrix:
        # count the values of each row
        counter = 0

        # iterate through every value
        for value in row:
            # append the value with the divider to the tex string
            tex_string = '{} {} &'.format(tex_string, value)

            # increase the counter by 1
            counter += 1

        # after each row, insert a new line and remove the last dividing character
        tex_string = tex_string[:-1] + ' \\\\\n'


    # after the values, remove the last new line and close the matrix
    tex_string = tex_string[:-3] + '\n\\end{bmatrix}'

    # open the tex file...
    with open('../values/matrix.tex', 'w') as fp:
        # ...and write the tex string to the file
        fp.write(tex_string)
        print('==> wrote matrix to file ./values/matrix.tex')

    #TODO: write matrix dimension to a seperate file in order to handle
    #      big matrices in LaTeX properly
