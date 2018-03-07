#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski


import numpy

from Crypto.Util import number
from functools import reduce
from math import log


def generate_keys(bit_length=1024, delta=0.020, m=4):
    """
    Generate p, q, N, dp, dq, d, e for CRT-RSA with
    small dp and dq for the given bit size.
    """
    # load required functions from other sage files
    load('./utils.sage')

    # generate primes for rsa
    p = number.getPrime(int(bit_length // 2))
    q = number.getPrime(int(bit_length // 2))

    # calculate the N
    N = p * q

    # get our boundary for the attack so we are able to compute
    # small dp's and dq's
    #TODO: find a way to calculate this with big N
    boundary = N ** delta

    # get random smal dp and dq with the calculated boundary
    dp = ZZ.random_element(2, int(boundary) - 1)
    dq = ZZ.random_element(2, int(boundary) - 1)


    # get the related e to our d
    try:
        d = crt(dp, dq, p - 1, q - 1)

        e = inverse_mod(d, (p - 1) * (q - 1))

        # make sure N and (N - 1) can be inverted in e^m
        if gcd(e**m, N) != 1 or gcd(e**m, N - 1) != 1:
            keys = generate_keys(bit_length, delta, m)

        # calculate the kq and kp  
        kp = ((e*dp) - 1) /  (p - 1)
        kq = ((e*dq) - 1) /  (q - 1)

        keys = {
                'p': p,
                'q': q,
                'N': p * q,
                'e': e,
                'd': d,
                'dp': dp,
                'dq': dq,
                'kp': kp,
                'kq': kq
                }

        if kp == 0 or kq == 0:
            keys = generate_keys(bit_length, delta, m)

    except ZeroDivisionError:
        # if no inverse exists, start all over again
        keys = generate_keys(bit_length, delta, m)
    except ValueError:
        # if crt does not work
        keys = generate_keys(bit_length, delta, m)

    return keys


def check_parameters(N, e, d):
    """
    Check the given RSA-parameters N, e, d. This function
    generates a random messge, encrypts and decrypts it in order
    to verify the RSA-parameters. If the decrypted and generated
    message are equal, this function return True.
    """
    # generate a random message
    message = number.getRandomInteger(number.size(N) - int(2))

    # encrypt the message with e
    encrypted_message = sqm(message, e, N)

    # encrypt the message with d
    decrypted_message = sqm(encrypted_message, d, N)

    # check if encrypted and decrypted message are equal
    if message == decrypted_message:
        return True
    else:
        return False


def main():
    """
    The main function of the module in order to get
    CRT-RSA-parameters as text values and in order to
    test this script.
    """
    # print a small logo
    print(' ______ ______ _______       ______ _______ _______        _______ _______ _______ ')
    print('|      |   __ \_     _|_____|   __ \     __|   _   |______|     __|    ___|    |  |')
    print('|   ---|      < |   ||______|      <__     |       |______|    |  |    ___|       |')
    print('|______|___|__| |___|       |___|__|_______|___|___|      |_______|_______|__|____|')
    print('\n')

    # generate a parameter set for crt rsa
    keys = generate_keys()

    # calculate the bit sizes of the generated values
    bit_length_dp = number.size(keys['dp'])
    bit_length_dq = number.size(keys['dq'])
    bit_length_e = number.size(keys['e'])
    bit_length_d = number.size(keys['d'])
    bit_length_N = number.size(keys['N'])

    # print the sizes
    print('Bitlength N = {}'.format(bit_length_N))
    print('Bitlength e = {}'.format(bit_length_e))
    print('Bitlength d = {}'.format(bit_length_d))
    print('Bitlength dp = {}'.format(bit_length_dp))
    print('Bitlength dq = {}'.format(bit_length_dq))

    # print the parameters
    for key in keys:
        print('')
        print('{} - {}'.format(key, keys[key]))

    # check if the parameters work for rsa
    check = check_parameters(keys['N'], keys['e'], keys['d'])
    print('\n==> Parameters working: {}'.format(check))
    pass


#if __name__ == '__main__':
#    main()
