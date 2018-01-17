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
    delta, m, bit_length, tau, debug, test, nogroebner, noreduction = parse_args()

    # test the lll algorithm
    if test:
        test_lll_reduction()

    # generate a complete crt-rsa keyset with the given parameters
    keys = generate_keys(
            bit_length=bit_length,
            delta=delta
            )

    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq> = PolynomialRing(ZZ, order='lex')

    # print some stats
    pprint('generated crt-rsa parameters')
    pprint('N = {}'.format(keys['N']))
    pprint('e = {}'.format(keys['e']))
    pprint('dp = {}'.format(keys['dp']))
    pprint('dq = {}'.format(keys['dq']))
    pprint('p = {}'.format(keys['p']))
    pprint('q = {}'.format(keys['q']))
    pprint('|N| = {} Bit'.format(bit_length))

    pprint('starting coppersmith with')
    pprint('m = {}    tau = {}'.format(m, tau))

    # generate the lattice for our parameters
    matrix, col_indice, polynomials_tuple = generate_lattice(
            keys['N'],
            keys['e'],
            m,
            tau,
            debug
            )

    #NOTE: check the dimension: in the paper its 180 for m=8
    if matrix.nrows() == matrix.ncols():
        pprint("squared matrix          [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
    else:
        pprint("squared matrix          [" + Fore.RED + " failed " + Fore.RESET + "]") 
        pprint('[' + Fore.RED + ' ERROR ' + Fore.RESET + '] got a {}x{} matrix'.format(matrix.nrows(), matrix.ncols()))

    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lex')

    # get polynomials from polynomials tuple
    polynomials = [value[0] for value in polynomials_tuple]

    # check if the roots match
    if root_check(polynomials, keys, debug):
        pprint("root check              [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
    else:
        pprint("root check              [" + Fore.RED + " failed " + Fore.RESET + "]") 

    # invert the colum index relation
    inverted_col_indice = {}
    for index in col_indice:
        inverted_col_indice[col_indice[index]] = index
    inverted_col_indice

    # reduce it
    if not noreduction:
        try:
            reduced_matrix = matrix.LLL()
        except:
            pprint("LLL-reduction           [" + Fore.RED + " failed " + Fore.RESET + "]") 
            return

        pprint("LLL-reduction           [" + Fore.GREEN + " passed " + Fore.RESET + "]") 

        if reduced_matrix.is_LLL_reduced():
            pprint("is LLL reduced          [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
        else:
            pprint("is LLL reduced          [" + Fore.RED + " failed " + Fore.RESET + "]") 
    else:
        pprint("LLL-reduction           [" + Fore.YELLOW + "  skip  " + Fore.RESET + "]") 
        reduced_matrix = matrix

    #NOTE: debugging purpose
    #TODO: debug flag
    #ones_matrix = matrix_to_ones(matrix, N)
    #sorted_matrix = matrix_sort_stairs(ones_matrix)
    #print_matrix(sorted_matrix)

    # initialize an array for the polynomials
    polynom_vector, reduced_polynomials = [], []

    counter = -1

    # iterate through the rows of the reduced coefficient matrix
    for row in reduced_matrix.rows():
        # initialize the polunomial for each row
        p = 0

        counter += 1

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
            x_q_1 = xq1^int(monom_grade[2])
            x_q_2 = xq2^int(monom_grade[3])
            y_p = yp^int(monom_grade[4])
            y_q = yq^int(monom_grade[5])

            # add the monomial myltyplied with its coefficient to the polunomial
            p += coefficient * x_p_1 * x_p_2 * x_q_1 * x_q_2 * y_p * y_q

            # add the absolute value of the coefficient squared for howgrave
            # grahams theorem howgrave_sum += abs(coefficient)^2
            howgrave_sum += abs(coefficient)^2

        # append the polynomials to the reduced polynomials array
        reduced_polynomials.append(p)

        # howgrave-grahams lemma in order to get polynoms, which are
        # small enough
        if sqrt(howgrave_sum) < ((keys['e']^m) / sqrt(reduced_matrix.ncols())):
            # append the polynomial to the polunomial vector for calculating the
            # groebner basis
            polynom_vector.append(p)

    # check the roots of the reduced polynomials
    if root_check(reduced_polynomials, keys, debug):
        pprint("reduced root check      [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
    else:
        pprint("reduced root check      [" + Fore.RED + " failed " + Fore.RESET + "]") 

    # create an ideal out of the polunomial vector
    pprint('create ideal with {} polynoms'.format(len(polynom_vector)))
    I = ideal(polynom_vector)

    # calculate the groebner basis
    if not nogroebner:
        pprint('calculate groebner basis')
        B = I.groebner_basis(algorithm='libsingular:slimgb', prot=True)

        # print the basis
        print(B)


# execute the main function, if the script is called on the commandline
if __name__ == '__main__':
    main()
