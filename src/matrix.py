#!/bin/env python3
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2017 - Daniel Jankowski


from fpylll import IntegerMatrix
from numpy import *


class Matrix(IntegerMatrix):
    """
    An abbrevatoin of fpylll.IntegerMatrix with additional functions
    like calculating the determinant converting to an array matrix.
    """

    def toArrayMatrix(self):
        """
        Convert the IntegerMatrix to a array matrix.

        :return: returns an array representation of the matrix

        Example:
            > A = Matrix([[1, 2], [3, 4]])
            > A.toArrayMatrix()
            [[1, 2], [3, 4]]
            > type(A.toArrayMatrix())
            <class 'list'>
        """
        # the output array, that consists of the row-arrays
        output_matrix = []

        # iterate through every row of the matrix
        for row in range(self.nrows):
            # clear the array for each new row of the output
            # array
            new_row = []

            # iterate through every value of each row
            for value in range(self.ncols):
                # append the value to the new row...
                new_row.append(self.__getitem__((row, value)))
            
            # ...and append the new row to the output matrix
            output_matrix.append(new_row)

        # return the output matrix
        return output_matrix 

    def det(self):
        """
        Calculate the determinant of the matrix with numpy.

        :return: An integer value of the determinant

        Example:
            > A = Matrix.from_matrix([[1, 2], [3, 4]])
            > A.det()
            -2.0
        """
        # convert the IntegerMatrix to a matrix consisting of arrays
        array_matrix = self.toArrayMatrix()

        # calculate the determinant with numpy
        return numpy.linalg.det(array_matrix)
