#!/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2017 - Daniel Jankowski


def main():
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

    keys = generate_keys(1024)

    # generate the lattice for our parameters
    matrix, col_indice = generate_lattice(
            keys['N'],
            keys['e'],
            8
            )

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
    R.<xp1, xp2, yq, yp, xq1, xq2> = PolynomialRing(ZZ)

    polynom_vector = []

    for row in reduced_matrix:
        p = 1
        col_index = 0

        for coefficient in row:
            monom_grade = inverted_col_indice[col_index]
            col_index += 1

            x_p_1 = xp1^int(monom_grade[0])
            x_p_2 = xp2^int(monom_grade[1])
            x_q_1 = xq1^int(monom_grade[2])
            x_q_2 = xq2^int(monom_grade[3])
            y_p = yp^int(monom_grade[4])
            y_q = yq^int(monom_grade[5])

            p += coefficient * x_p_1 * x_p_2 * x_q_1 * x_q_2 * y_p * y_q
        polynom_vector.append(p)

    I = Ideal(polynom_vector)
    B = I.groebner_basis()
    print(B)


    #TODO: use groebner
    pass


if __name__ == '__main__':
    main()
