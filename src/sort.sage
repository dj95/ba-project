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
    for j in range(i, len(matrix) - i):

        if matrix[i][j] != 0:
            col_index[j] -= 1

        if matrix[len(matrix) - 1 - i][j] != 0:
            col_index[j] -= 1

    for j in range(i, len(matrix) - i):

        if matrix[j][i] != 0:
            row_index[j] -= 1

        if matrix[j][len(matrix) - 1 - i] != 0:
            row_index[j] -= 1

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
    work_matrix = matrix_to_array(matrix)

    row_index = create_row_index(work_matrix)
    col_index = create_column_index(work_matrix)

    sorted_matrix = []

    for i in range(int(floor(len(work_matrix) / 2))):
        j = len(work_matrix) - 1 - i
        if check_outer_values(work_matrix, i):
            row_index, col_index = update_index(work_matrix, i, row_index, col_index)
            continue

        for k in range(i, j + 1):
            if row_index[k] == 1:
                swap_rows(work_matrix, k, i)
                row_index = swap_index(row_index, k, i)
                polynom_index = swap_index(polynom_index, k, i)

                if not check_top_values(work_matrix, i):
                    for l in range(i + 1, j + 1):
                        if work_matrix[i][l] != 0:
                            swap_cols(work_matrix, l, i)
                            col_index = swap_index(col_index, l, i)
                            monom_index = swap_index(monom_index, l, i)
                            break
                
                break

        for k in range(i + 1, j + 1):
            if col_index[k] == 1:
                swap_cols(work_matrix, k, j)
                col_index = swap_index(col_index, k, j)
                monom_index = swap_index(monom_index, k, j)

                if not check_right_values(work_matrix, i):
                    for l in range(i, j):
                        if work_matrix[l][j] != 0:
                            swap_rows(work_matrix, l, j)
                            row_index = swap_index(row_index, l, j)
                            polynom_index = swap_index(polynom_index, l, j)
                            break

                break


        row_index, col_index = update_index(work_matrix, i, row_index, col_index)

    # return it
    return Matrix(work_matrix), monom_index, polynom_index
