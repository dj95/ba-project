#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski

def main():
    """
    The main function of this application and the entrypoint for it.
    """
    # import sage modules
    load('./argparser.sage')
    load('./coppersmith.sage')
    load('./keygen.sage')
    load('./print.sage')
    load('./polynomials.sage')
    load('./shiftpolynomials.sage')
    load('./substitute.sage')
    load('./test.sage')

    # parse arguments
    delta, m, bit_length, tau = parse_args()

    # print a logo
    print('#####################################')
    print('# Attack an small CRT-RSA exponents #')
    print('#####################################')
    print('')

    # test the lll algorithm
    #test_lll_reduction()

    # generate a complete crt-rsa keyset with the given parameters
    keys = generate_keys(
            bit_length=bit_length,
            delta=delta
            )

    print('==> e*dp = 1 + {}(p - 1)'.format(keys['kp']))
    print('==> e*dq = 1 + {}(q - 1)'.format(keys['kq']))

    # generate the lattice for our parameters
    matrix, col_indice, polynomials = generate_lattice(
            keys['N'],
            keys['e'],
            m,
            tau
            )

    #NOTE: check the dimension: in the paper its 180 for m=8
    print('\n==> Got a {}x{} matrix'.format(matrix.nrows(), matrix.ncols()))

    # define the polynomial ring
    R.<xp1, xp2, yq, yp, xq1, xq2> = PolynomialRing(ZZ, order='lex')

    # checking polynomials
    polynomial_count, correct_count = len(polynomials), 0
    for p, set_count, index_set in polynomials:
        y = p(
            xp1=keys['kq'] - 1,
            xp2=keys['kp'],
            xq1=keys['kq'],
            xq2=keys['kq'] - 1,
            yp=keys['p'],
            yq=keys['q']
            ) % keys['e']
        correct_count += 1
        #print('\r==> [{}/{}] polynomials have the correct zero'.format(correct_count, polynomial_count)),
        #print('  -> {} {} {} : {}'.format(correct_count, set_count, index_set, y))
    print('==> {}/{} polynoms have the correct zero\n'.format(correct_count, polynomial_count))


    # invert the colum index relation
    inverted_col_indice = {}
    for index in col_indice:
        inverted_col_indice[col_indice[index]] = index
    inverted_col_indice

    #NOTE: debugging purpose
    ones_matrix = matrix_to_ones(matrix, N)
    sorted_matrix = matrix_sort_stairs(ones_matrix)
    print_matrix(sorted_matrix)

    # reduce it
    print('==> Reducing matrix with LLL-algorithm')
    reduced_matrix = matrix.LLL()
    print('==> Got the reduced matrix')

    # initialize an array for the polynomials
    polynom_vector = []

    # iterate through the rows of the reduced coefficient matrix
    for row in reduced_matrix:
        # initialize the polunomial for each row
        p = 1

        # reset the column index
        col_index = 0

        # initialize variables for howgrave graham theorem
        howgrave_sum = 0

        monome_count = 0

        # iterate through all values of the row
        for coefficient in row:
            # get the monome grade from the column index
            monom_grade = inverted_col_indice[col_index]

            # increase the column index for the next column
            col_index += 1

            if coefficient != 0:
                monome_count += 1

            # set the correct grade to the variables of each monomial
            x_p_1 = xp1^int(monom_grade[0])
            x_p_2 = xp2^int(monom_grade[1])
            y_q = yq^int(monom_grade[2])
            y_p = yp^int(monom_grade[3])
            x_q_1 = xq1^int(monom_grade[4])
            x_q_2 = xq2^int(monom_grade[5])

            # add the monomial myltyplied with its coefficient to the polunomial
            p += coefficient * x_p_1 * x_p_2 * x_q_1 * x_q_2 * y_p * y_q

            # add the absolute value of the coefficient squared for howgrave
            # grahams theorem
            howgrave_sum += abs(coefficient)^2

        #print(p)
        #print('==> Checking {} < {}'.format(sqrt(howgrave_sum), ((e^m) / sqrt(reduced_matrix.ncols()))))

        # howgrave-grahams lemma in order to get polynoms, which are
        # small enough
        if sqrt(howgrave_sum) < ((keys['e']^m) / sqrt(monome_count)):
            # append the polynomial to the polunomial vector for calculating the
            # groebner basis
            polynom_vector.append(p)

    #polynom_vector = polynom_vector[:len(polynom_vector)/10]

    p_v = []
    for i in range(6):
        p_v.append(polynom_vector[i])


    #print(polynom_vector)

    print('==> Create ideal from {} polynoms'.format(len(polynom_vector)))

    # create an ideal out of the polunomial vector
    I = ideal(p_v)

    print('==> Calculate groebner basis')

    # calculate the groebner basis
    B = I.groebner_basis(algorithm='libsingular:slimgb', prot=True)

    # print the basis
    print(B)


# execute the main function, if the script is called on the commandline
if __name__ == '__main__':
    main()
