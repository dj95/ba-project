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
    load('./sort.sage')
    load('./substitute.sage')
    load('./test.sage')

    # parse arguments
    delta, m, bit_length, tau, debug, test, nogroebner, noreduction, forcetriangle, printmatrix, jsonoutput = parse_args()

    # test the lll algorithm
    if test:
        test_lll_reduction()

    # generate a complete crt-rsa keyset with the given parameters
    keys = generate_keys(
            bit_length=bit_length,
            delta=delta
            )

    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq> = PolynomialRing(ZZ, order='lp')

    # calculate upper bounds for roots
    alpha = log(keys['e'], keys['N'])

    if keys['p'] < keys['q']:
        print(keys['q'] < 2 * keys['p'])
    else:
        print(keys['p'] < 2 * keys['q'])

    print(alpha + delta - 0.25)
    X = ceil(keys['N']^(alpha + delta - 0.5))
    Y = 2*ceil(1 + keys['N']^(0.5))

    #X = long(max(keys['kp'] + 1, keys['kq'] + 1))
    Y = long(max(keys['p'] + 1, keys['q'] + 1))

    print(keys['dp'] < keys['N']^delta)
    print(keys['dq'] < keys['N']^delta)

    # print some stats
    if not jsonoutput:
        pprint('generated crt-rsa parameters')
        pprint('N = {}'.format(keys['N']))
        pprint('e = {}'.format(keys['e']))
        pprint('dp = {}'.format(keys['dp']))
        pprint('dq = {}'.format(keys['dq']))
        pprint('p = {}'.format(keys['p']))
        pprint('q = {}'.format(keys['q']))
        pprint('kp = {}'.format(keys['kp']))
        pprint('kq = {}'.format(keys['kq']))
        pprint('|N| = {} Bit'.format(bit_length))
        pprint('kp < X {}'.format(keys['kp'] < X))
        pprint('kq < X {}'.format(keys['kq'] < X))
        pprint('p < Y {}'.format(keys['p'] < Y))
        pprint('q < Y {}'.format(keys['q'] < Y))

        pprint('{}'.format(
            keys['e']*keys['dp'] == 1 + keys['kp'] * (keys['p'] - 1)
                    ))
        pprint('{}'.format(
                keys['e']*keys['dq'] == 1 + keys['kq'] * (keys['q'] - 1)
                    ))

        pprint('{}'.format(
            (X^4 * Y^2) < keys['e']^m
                    ))

        pprint('starting coppersmith with')
        pprint('m = {}    tau = {}'.format(m, tau))

    pprint('alpha >= 7/16: {}'.format(alpha >= 7/16))
    pprint('delta < 1/2 - sqrt(alpha/7): {}'.format(delta < (1/2) - sqrt(alpha / 7)))

    # generate the lattice for our parameters
    matrix, col_indice, polynomials_tuple, row_index = generate_lattice(
            keys['N'],
            keys['e'],
            X, Y,
            m,
            tau,
            debug,
            jsonoutput
            )

    #NOTE: check the dimension: in the paper its 180 for m=8
    if matrix.nrows() == matrix.ncols():
        if not jsonoutput:
            pprint("squared matrix          [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
    else:
        if not jsonoutput:
            pprint("squared matrix          [" + Fore.RED + " failed " + Fore.RESET + "]") 
            pprint('[' + Fore.RED + ' ERROR ' + Fore.RESET + '] got a {}x{} matrix'.format(matrix.nrows(), matrix.ncols()))

    # define the polynomial ring
    R.<xp1, xp2, xq1, xq2, yp, yq > = PolynomialRing(ZZ, order='lp')

    # get polynomials from polynomials tuple
    polynomials = [value[0] for value in polynomials_tuple]

    # invert the colum index relation
    inverted_col_indice = dict((v,k) for k,v in col_indice.iteritems())

    # sort the matrix triangular
    matrix, inverted_col_indice, row_index = matrix_sort_triangle(matrix, inverted_col_indice, row_index)

    matrix = set_upper_bound(matrix, X, Y, inverted_col_indice, keys['e'], m)

    # substitute N from the diagonal
    #TODO: hier steck der Fehler!!!
    matrix = substitute_N(matrix, keys['N'], keys['e'], m, X, Y)

    substituted_polynomials = []

    # iterate through the rows of the reduced coefficient matrix
    for row in matrix.rows():
        # initialize the polunomial for each row
        p = 0

        # reset the column index
        col_index = 0

        # iterate through all values of the row
        for i in range(len(row)):
            # get the monome grade from the column index
            monom_grade = inverted_col_indice[i]

            # get the coefficient
            coefficient = row[i]

            # set the correct grade to the variables of each monomial
            x_p_1 = xp1^int(monom_grade[0])
            x_p_2 = xp2^int(monom_grade[1])
            x_q_1 = xq1^int(monom_grade[2])
            x_q_2 = xq2^int(monom_grade[3])
            y_p = yp^int(monom_grade[4])
            y_q = yq^int(monom_grade[5])

            exp_xp1 = int(monom_grade[0])
            exp_xp2 = int(monom_grade[1])
            exp_xq1 = int(monom_grade[2])
            exp_xq2 = int(monom_grade[3])
            exp_yp = int(monom_grade[4])
            exp_yq = int(monom_grade[5])

            x_bound = X^(exp_xp1 + exp_xp2 + exp_xq1 + exp_xq2)
            y_bound = Y^(exp_yp + exp_yq)
            
            x_bound_inverse = inverse_mod(x_bound, keys['e']^m)
            y_bound_inverse = inverse_mod(y_bound, keys['e']^m)

            #if coefficient % X == 0:
            #    print('X')
            #if coefficient % Y == 0:
            #    print('Y')
            #if coefficient % keys['N'] == 0:
            #    print('N')
            #if coefficient % keys['e']^m == 0:
            #    print('em')

            test = coefficient
            eX = 0
            eY = 0
            while test % X == 0 and test != 0:
                test /= X
                eX += 1

            while test % Y == 0 and test != 0:
                test /= Y
                eY += 1

            #if coefficient != 0:
            #    print(coefficient)
            #    print(coefficient % keys['e']^m)
            #    print(monom_grade)
            #    print(eX)
            #    print(eY)

            coefficient = (coefficient * x_bound_inverse * y_bound_inverse) % keys['e']^m

            # add the monomial myltyplied with its coefficient to the polunomial
            p += coefficient * x_p_1 * x_p_2 * x_q_1 * x_q_2 * y_p * y_q

        # append the polynomials to the reduced polynomials array
        substituted_polynomials.append(p)

    # check if the roots match
    if root_check(substituted_polynomials, keys, debug, m):
        if not jsonoutput:
            pprint("root check              [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
    else:
        if not jsonoutput:
            pprint("root check              [" + Fore.RED + " failed " + Fore.RESET + "]") 

    # reduce it
    if not noreduction:
        try:
            reduced_matrix = matrix.LLL()
        except Exception as e:
            if not jsonoutput:
                pprint("LLL-reduction           [" + Fore.RED + " failed " + Fore.RESET + "]") 
                print(e)
            return

        if not jsonoutput:
            pprint("LLL-reduction           [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
    else:
        if not jsonoutput:
            pprint("LLL-reduction           [" + Fore.YELLOW + "  skip  " + Fore.RESET + "]") 
        reduced_matrix = matrix

    # check if determinant is lower than e^nm
    if (reduced_matrix.determinant() < keys['e']^(reduced_matrix.ncols() * m)):
        pprint("det(B) < e^nm           [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
    else:
        pprint("det(B) < e^nm           [" + Fore.RED + " failed " + Fore.RESET + "]") 

    # it arg --print is true, print the matrix to tex file
    if printmatrix:
        # substitute values != 0 by 1 in order to make the matrix readable
        ones_matrix = matrix_to_ones(reduced_matrix, keys['N'])

        # print the matrix to file
        print_matrix(ones_matrix, inverted_col_indice, row_index)

    # initialize an array for the polynomials
    polynom_vector, reduced_polynomials = [], []

    detL = reduced_matrix.determinant()
    print(detL)
    detL = sqrt(detL, reduced_matrix.ncols())
    detL = long(ceil(detL))
    print(detL)

    pprint(' Check ||v|| < det(L)^(1/n)')

    # iterate through the rows of the reduced coefficient matrix
    for row in reduced_matrix.rows():
        # initialize the polunomial for each row
        p = 0

        # initialize variables for howgrave graham theorem
        howgrave_sum = 0

        monome_count = 0

        # iterate through all values of the row
        for i in range(len(row)):
            # get the monome grade from the column index
            monom_grade = inverted_col_indice[i]

            # get the coefficient
            coefficient = row[i]

            # increase the monom count by 1 if we find a monom != 0
            if coefficient != 0:
                monome_count += 1

            # set the correct grade to the variables of each monomial
            x_p_1 = (xp1 / X)^int(monom_grade[0])
            x_p_2 = (xp2 / X)^int(monom_grade[1])
            x_q_1 = (xq1 / X)^int(monom_grade[2])
            x_q_2 = (xq2 / X)^int(monom_grade[3])
            y_p = (yp / Y)^int(monom_grade[4])
            y_q = (yq / Y)^int(monom_grade[5])

            # add the absolute value of the coefficient squared for howgrave
            # grahams theorem howgrave_sum += abs(coefficient)^2
            howgrave_sum += abs(coefficient)^2

            # add the monomial multiplied with its coefficient to the polunomial
            p += (coefficient * x_p_1 * x_p_2 * x_q_1 * x_q_2 * y_p * y_q)

        # append the polynomials to the reduced polynomials array
        reduced_polynomials.append(p)

        howgrave = floor(sqrt(howgrave_sum))

        if howgrave < detL:
            print('  -> True')
        else:
            print('  -> False')

        # howgrave-grahams lemma in order to get polynoms, which are
        # small enough
        if howgrave < keys['e']^m:
            # append the polynomial to the polunomial vector for calculating the
            # groebner basis
            polynom_vector.append(p)

    # check the roots of the reduced polynomials
    if reduced_root_check(reduced_polynomials, keys, debug, m):
        if not jsonoutput:
            pprint("reduced root check      [" + Fore.GREEN + " passed " + Fore.RESET + "]") 
    else:
        if not jsonoutput:
            pprint("reduced root check      [" + Fore.RED + " failed " + Fore.RESET + "]") 

    if not jsonoutput:
    #    pprint("det(M) = {}".format(reduced_matrix.det()))
    #    pprint("howgrave = {}".format(reduced_matrix.det() < keys['e']^(reduced_matrix.ncols() * m)))

        # create an ideal out of the polunomial vector
        pprint('create ideal with {} polynoms'.format(len(polynom_vector)))

    #new_pv = []
    #for p in polynom_vector:
    #    new_pv.append(p(xp1 = (keys['kq'] - 1)))
    #polynom_vector = new_pv

    #print(polynom_vector)


    # calculate the groebner basis
    if not nogroebner:
        # create an ideal from the polynom vector
        I = Ideal(polynom_vector)
        if not jsonoutput:

            pprint('calculate groebner basis')
            B = I.groebner_basis(algorithm='libsingular:slimgb', prot=True)
        else:
            B = I.groebner_basis(algorithm='libsingular:slimgb', prot=False)

        # print the basis
        #print(B)
        if not jsonoutput:
            pprint('basis length: {}'.format(len(B)))
        # print json output
        else:
            output = {}
            output['m'] = m
            output['tau'] = tau
            output['delta'] = delta
            output['matrix_dimension'] = reduced_matrix.ncols()
            output['no_errors'] = True
            output['bit_length'] = bit_length
            print(output)

        equations = []
        x1, x2, x3, x4, y1, y2 = var('x1 x2 x3 x4 y1 y2')

        #if len(B) >= 6:
        for row in B:
            if not jsonoutput:
                print(row)


            #print(row.dict())
            
            #if len(row.dict()) == 1:
            #    print(row % keys['e'])

            eq1 = 0

            for monom in row.dict():
                print("{} - {}".format(monom, row.dict()[monom]))
                m = row.dict()[monom] * x1^monom[0] * x2^monom[1] * y1^monom[4] * y2^monom[5] * x3^monom[2] * x4^monom[3]
                eq1 += m

            eq = eq1 == 0

            equations.append(eq)

        for row in equations:
            print(row)


        print(solve(equations, x1, x2, y1, y2))


# execute the main function, if the script is called on the commandline
if __name__ == '__main__':
    main()
