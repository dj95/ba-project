#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


def index_set_x_stages(m, t, x):
    """
    Generate the first index set I_x in a generator with
    the first, second, third or fourth subset.
    """
    if x == 1:
        # initial values
        i1, i2, j1, j2 = 0, 0, 0, 0
        u = -1

        # let the generator work until we hit the boundary
        while i1 <= floor(m / 2):
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while i2 <= floor(m / 2):
                # if j2 is in the boundary, return
                while u <= (min(int(floor(m / 2)) - i1, int(floor(m / 2)) - i2) - 1):
                    u = u + 1
                    yield i1, i2, j1, j2, u

                # calculate the values for the next call
                i2, u = i2 + 1, -1
            # calculate the values for the next call
            i1, i2, u = i1 + 1, 0, -1
    elif x == 2:
        # initial values
        i1, i2, j1, j2 = 0, 1, 1, 0
        u = -1

        # let the generator work until we hit the boundary
        while i1 <= floor(m / 2) - 1:
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while i2 <= floor(m / 2):
                # if j2 is in the boundary, return
                while u <= (min(int(floor(m / 2)) - i1 - 1, int(floor(m / 2)) - i2) - 1):
                    u = u + 1
                    yield i1, i2, j1, j2, u

                # calculate the values for the next call
                i2, u = i2 + 1, -1
            # calculate the values for the next call
            i1, i2, u = i1 + 1, 1, -1
    elif x == 3:
        # initial values
        i1, i2, j1, j2 = 0, 0, 1, 0
        u = -1

        # let the generator work until we hit the boundary
        while i1 <= floor(m / 2) - 1:
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while j1 <= (floor(m / 2) - i1):
                # if j2 is in the boundary, return
                while u <= (int(int(floor(m / 2)) - i1 - j1) - 1):
                    u = u + 1
                    yield i1, i2, j1, j2, u

                # calculate the values for the next call
                j1, u = j1 + 1, -1
            # calculate the values for the next call
            i1, j1, u = i1 + 1, 1, -1
    elif x == 4:
        # initial values
        i1, i2, j1, j2 = 0, 0, 0, 1
        u = -1

        # let the generator work until we hit the boundary
        while i2 <= floor(m / 2) - 1:
            # iterate through all i2 to the boundary and increase i1 if we hit it
            while j2 <= (floor(m / 2) - i2):
                # if j2 is in the boundary, return
                while u <= (int(int(floor(m / 2)) - i2 - j2) - 1):
                    u = u + 1
                    yield i1, i2, j1, j2, u

                # calculate the values for the next call
                j2, u = j2 + 1, -1
            # calculate the values for the next call
            i2, j2, u = i2 + 1, 1, -1


def index_set_x(m, t):
    """
    Generate all index sets of the different subsets of index set I_x
    with one generator expression.
    """
    for set in index_set_x_stages(m, t, 1):
        yield set
    for set in index_set_x_stages(m, t, 2):
        yield set
    for set in index_set_x_stages(m, t, 3):
        yield set
    for set in index_set_x_stages(m, t, 4):
        yield set


def index_set_y_p(m, t):
    """
    Generate the second index set I_y,p in a generator.
    """
    # initial values
    i1, i2 = 0, 0
    j1 = 0

    # let the generator work until we hit the boundary
    while i1 <= floor(m / 2):
        # iterate through all i2 to the boundary and increase i1 if we hit it
        while i2 <= floor(m / 2):
            # if j2 is in the boundary, return
            while j1 <= (int(floor(t * (i1 + i2)) - ceil((i1 + i2) / 2)) - 1):
                j1 = j1 + 1
                yield i1, i2, j1

            # calculate the values for the next call
            i2, j1 = i2 + 1, 0
        # calculate the values for the next call
        i1, i2, j1 = i1 + 1, 0, 0


def index_set_y_q(m, t):
    """
    Generate the third index set I_y,q in a generator.
    """
    # initial values
    i1, i2 = 0, 0
    j2 = 0

    # let the generator work until we hit the boundary
    while i1 <= floor(m / 2):
        # iterate through all i2 to the boundary and increase i1 if we hit it
        while i2 <= floor(m / 2):
            # if j2 is in the boundary, return
            while j2 <= (int(floor(t * (i1 + i2)) - floor((i1 + i2) / 2)) - 1):
                j2 = j2 + 1
                yield i1, i2, j2

            # calculate the values for the next call
            i2, j2 = i2 + 1, 0
        # calculate the values for the next call
        i1, i2, j2 = i1 + 1, 0, 0

