#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


def matrix_to_sizematrix(matrix):
    # initialize numpy array with float64 for the heatmap
    output_matrix = numpy.array(matrix, dtype=numpy.float64)

    # iterate through every row
    for i in range(len(output_matrix)):
        # iterate through every value of each row
        for j in range(len(output_matrix[i])):
            # if the entry is greater than 0
            if output_matrix[i][j] > 0:
                # set the value to log base 2
                output_matrix[i][j] = math.log(output_matrix[i][j], 2)
            # if the entry is lower than 0
            elif output_matrix[i][j] < 0:
                # set the value to - log base 2
                output_matrix[i][j] = - math.log(abs(output_matrix[i][j]), 2)

    # return the matrix
    return output_matrix


def print_heatmap(heatmap_matrix, filename):
    # initialize the figure to plot in to
    fig = plt.figure()

    # get the coordinate system
    ax = fig.add_subplot(111)

    # set labels for x and y
    ax.set_xlabel('Monomials in $\{x_{p_1}, x_{p_2}, x_{q_1}, x_{q_2}, y_p, y_q \}$')
    ax.set_ylabel('Polynomials (Sum of Monomials)')

    # create the heatmap from the matrix 
    im1 = ax.matshow(heatmap_matrix, cmap='coolwarm')

    # set title for the axis
    ax.set_title(
            "Heatmap of monomial's size log(2) in " \
            + "$%3d\\times%3d$ coefficient matrix after sort." \
            %(len(heatmap_matrix), len(heatmap_matrix)),
            fontsize="18"
            )

    # print the colorbar for the heatmap into the figure
    fig.colorbar(im1, ax=ax)

    # the the figure to a png
    fig.savefig(filename, bbox_inches='tight')


def print_matrix(matrix, inverted_col_indice, row_index):
    """
    Print the matrix in TeX-format to a file.
    """
    # start the matrix in tex
    tex_string = '\\begin{bmatrix}\n'

    # iterate through every value
    #for i in range(len(matrix[0])):
        # append the value with the divider to the tex string
    #    tex_string = '{} {} &'.format(tex_string, inverted_col_indice[i])

        # increase the counter by 1
        #counter += 1

    # after each row, insert a new line and remove the last dividing character
    #tex_string = tex_string[:-1] + ' \\\\\n'
    
    r_counter = 0

    # iterate through every row of the matrix
    for row in matrix:
        # count the values of each row
        counter = 0

        # iterate through every value
        for value in row:
            if counter == r_counter:
                # append the value with the divider to the tex string
                if value > 0:
                    tex_string = '{} {} &'.format(tex_string, '+')
                elif value < 0:
                    tex_string = '{} {} &'.format(tex_string, '-')
                else:
                    tex_string = '{} {} &'.format(tex_string, '.')
            else:
                # append the value with the divider to the tex string
                if value > 0:
                    tex_string = '{} {} &'.format(tex_string, '+')
                elif value < 0:
                    tex_string = '{} {} &'.format(tex_string, '-')
                else:
                    tex_string = '{} {} &'.format(tex_string, '.')

            # increase the counter by 1
            counter += 1

        # append the value with the divider to the tex string
        tex_string = '{} {} &'.format(tex_string, row_index[r_counter])

        r_counter += 1

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
