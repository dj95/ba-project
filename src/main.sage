#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2017 - Daniel Jankowski


def main():
    """
    The main function of this application and the entrypoint for it.
    """
    # print a logo
    print('#####################################')
    print('# Attack an small CRT-RSA exponents #')
    print('#####################################')
    print('')

    # import sage modules
    load('./coppersmith.sage')
    load('./test.sage')
    load('./keygen.sage')

    # test the lll algorithm
    test_lll_reduction()

    # generate a complete crt-rsa keyset with the given parameters
    keys = generate_keys(
            bit_length=1024,
            delta=0.010
            )

    # generate the lattice for our parameters
    matrix, col_indice = generate_lattice(
            keys['N'],
            keys['e'],
            8
            )

    #NOTE: check the dimension: in the paper its 180 for m=8
    print('\n==> Got a {}x{} matrix'.format(matrix.nrows(), matrix.ncols()))
    print('==> Reducing matrix with LLL-algorithm')

    # invert the colum index relation
    inverted_col_indice = {}
    for index in col_indice:
        inverted_col_indice[col_indice[index]] = index
    inverted_col_indice

    # reduce it
    reduced_matrix = matrix.LLL()

    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yq, yp> = PolynomialRing(ZZ)

    # initialize an array for the polynomials
    polynom_vector = []

    # iterate through the rows of the reduced coefficient matrix
    for row in reduced_matrix:
        # initialize the polunomial for each row
        p = 1

        # reset the column index
        col_index = 0

        # iterate through all values of the row
        for coefficient in row:
            # get the monome grade from the column index
            monom_grade = inverted_col_indice[col_index]

            # increase the column index for the next column
            col_index += 1

            # set the correct grade to the variables of each monomial
            x_p_1 = xp1^int(monom_grade[0])
            x_p_2 = xp2^int(monom_grade[1])
            x_q_1 = xq1^int(monom_grade[2])
            x_q_2 = xq2^int(monom_grade[3])
            y_p = yp^int(monom_grade[4])
            y_q = yq^int(monom_grade[5])

            # add the monomial myltyplied with its coefficient to the polunomial
            p += coefficient * x_p_1 * x_p_2 * x_q_1 * x_q_2 * y_p * y_q

        # append the polunomial to the polunomial vector for calculating the 
        # groebner basis
        polynom_vector.append(p)

    # create an ideal out of the polunomial vector
    I = Ideal(polynom_vector)

    # calculate the groebner basis
    B = I.groebner_basis()

    # print the basis
    print(B)


# execute the main function, if the script is called on the commandline
if __name__ == '__main__':
    main()
