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
from sage.modules.free_module_integer import IntegerLattice


def create_row_index(matrix):
    row_index = []

    # append an index to each row
    for row in matrix:
        index = 0

        # iterate through the columns
        for i in range(len(row)):
            # set the index to the column count of the value in order to save
            # the highest column count of an existing monom
            if row[i] != 0:
                index += 1

        # append the index to each row
        row_index.append(index)

    return row_index


def create_column_index(matrix):
    column_index = []

    # iterate through every column
    for i in range(len(matrix[0])):
        index = 0

        # iterate through every value of each column
        for j in range(len(matrix)):
            # if the value is not zero, increase the index
            if matrix[j][i] != 0:
                index += 1

        # append the index to all indices
        column_index.append(index)

    # return it
    return column_index


def matrix_to_array(matrix):
    work_matrix = []

    zero_count = 0
    non_zero_count = 0

    counter = -1
    # convert matrix object to array matrix
    for row in matrix:
        # append a new row for every row
        work_matrix.append([])
        counter += 1

        # append the values to each row
        for value in row:
            if value == 0:
                zero_count += 1
            else:
                non_zero_count += 1
            work_matrix[counter].append(value)

    if non_zero_count > len(matrix[0]) and non_zero_count < zero_count:
        pprint("triangle possible       [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
    else:
        pprint("triangle possible       [" + Fore.RED + " failed " + Fore.RESET + "]") 

    return work_matrix


def check_right_values(matrix, i):
    j = len(matrix) - 1 - i
    # check if top left and bottom right values are 1
    if matrix[j][j] == 0:
        return False

    # check if top row and right colum contains 0 only
    for k in range(i, j):
        if matrix[k][j] != 0:
            return False

    # return True if it is the case
    return True


def check_top_values(matrix, i):
    # check if top left and bottom right values are 1
    if matrix[i][i] == 0:
        return False

    # check if top row and right colum contains 0 only
    for k in range(i + 1, len(matrix) - i):
        if matrix[i][k] != 0:
            return False

    # return True if it is the case
    return True


def check_outer_values(matrix, i):
    j = len(matrix) - 1 - i
    # check if top left and bottom right values are 1
    if matrix[i][i] == 0 or matrix[j][j] == 0:
        return False

    # check if top row and right colum contains 0 only
    for k in range(i, len(matrix) - 1 - i):
        if matrix[i][k + 1] != 0 or matrix[k][j] != 0:
            return False

    # return True if it is the case
    return True


def update_index(matrix, i, row_index, col_index):
    # iterate through the sublattice
    for j in range(i, len(matrix) - i):
        # if the top and bottom row != 0... 
        if matrix[i][j] != 0:
            # ...reduce the colum index by 1
            col_index[j] -= 1

        if matrix[len(matrix) - 1 - i][j] != 0:
            # ...reduce the colum index by 1
            col_index[j] -= 1

        # if the left or right column != 0...
        if matrix[j][i] != 0:
            # ...reduce the row index by 1
            row_index[j] -= 1

        if matrix[j][len(matrix) - 1 - i] != 0:
            # ...reduce the row index by 1
            row_index[j] -= 1

    # return the updated indices
    return row_index, col_index


def swap_index(matrix, i, j):
    # save the first row to swap
    save_row = matrix[i]

    # overwrite the first row with the other one
    matrix[i] = matrix[j]

    # overwrite the other row with the saved one
    matrix[j] = save_row

    # return the matrix
    return matrix


def swap_rows(matrix, i, j):
    # save the first row to swap
    save_row = matrix[i]

    # overwrite the first row with the other one
    matrix[i] = matrix[j]

    # overwrite the other row with the saved one
    matrix[j] = save_row

    # return the matrix
    return matrix


def swap_cols(matrix, i, j):
    # iterate through every row
    for k in range(len(matrix)):
        # save the ith col value in kth row
        save_value = matrix[k][i]

        # overwrite the saved value with the one to swap
        matrix[k][i] = matrix[k][j]

        # overwrite the swapped value with the saved one
        matrix[k][j] = save_value

    # return the matrix
    return matrix


def matrix_sort_triangle(matrix, monom_index, polynom_index):
    """
    Sort the matrix rows that they form steps or a triangle.
    """
    # convert the matrix into an array matrix
    work_matrix = matrix_to_array(matrix)

    # get an index for rows and cols, that counts values != 0
    row_index = create_row_index(work_matrix)
    col_index = create_column_index(work_matrix)

    # iterate through the half of the matrix row lengths
    for i in range(int(floor(len(work_matrix) / 2))):
        # get the ith value from the right side
        j = len(work_matrix) - 1 - i

        # check if top-left and bottom-right corners != 0 and top+right row equals 0
        if check_outer_values(work_matrix, i):
            # update the row and col index for the next iteration
            row_index, col_index = update_index(work_matrix, i, row_index, col_index)
            continue

        # iterate throygh the row index
        for k in range(i, j + 1):
            # it there is only 1 value in the row...
            if row_index[k] == 1:
                # swap the row to the top
                work_matrix = swap_rows(work_matrix, k, i)

                # swap the row index...
                row_index = swap_index(row_index, k, i)

                # ...and the polynom_index, too
                polynom_index = swap_index(polynom_index, k, i)

                # if the top left value == 0
                if not check_top_values(work_matrix, i):
                    # iterate through the columes
                    for l in range(i + 1, j + 1):
                        # find the column != 0 in this row
                        if work_matrix[i][l] != 0:
                            # and swap it to the left
                            work_matrix = swap_cols(work_matrix, l, i)

                            # update the colum index...
                            col_index = swap_index(col_index, l, i)

                            # ...and the monom index as well
                            monom_index = swap_index(monom_index, l, i)

                            # escape the for loop since we got our result
                            break

                # escape the for loop since we got our result
                break

        # iterate throygh the col index
        for k in range(i + 1, j + 1):
            # it there is only 1 value in the row...
            if col_index[k] == 1:
                # swap the row to the right
                work_matrix = swap_cols(work_matrix, k, j)

                # swap the col index
                col_index = swap_index(col_index, k, j)

                # and swap the monom index
                monom_index = swap_index(monom_index, k, j)

                # if the bottom right corner is != 0
                if not check_right_values(work_matrix, i):
                    # iterate through every row
                    for l in range(i, j):
                        # if we found the value != 0
                        if work_matrix[l][j] != 0:
                            # swap it to the bottom
                            work_matrix = swap_rows(work_matrix, l, j)

                            # swap the row index
                            row_index = swap_index(row_index, l, j)

                            # and swap the polynom index as well
                            polynom_index = swap_index(polynom_index, l, j)

                            # escape the for loop since we got our result
                            break

                # escape the for loop since we got our result
                break

        # update the row and col index for the next iteration
        row_index, col_index = update_index(work_matrix, i, row_index, col_index)

    # return it
    return Matrix(work_matrix), monom_index, polynom_index
