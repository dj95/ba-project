#!/usr/bin/env sage
#
# keygen for crt rsa
#
# (c) 2017 Daniel Jankowski


from Crypto.Util import number
from functools import reduce
from math import log


def sqm(a, e, m):
    """
    Calculate a**e mod m with square and multiply.
    """
    bitlen = len(bin(e)[2:])
    c = a
    bitlen -= 1

    while bitlen > 0:
        bitlen -= 1
        c = (c * c) % m
        mul = ((e >> bitlen) & 0x1)
        if mul:
            c = (c * a) % m
    return c


def generate_keys(bit_length=1024):
    """
    Generate p, q, N, dp, dq, d, e for CRT-RSA with
    small dp and dq for the given bit size.
    """
    # generate primes for rsa
    p = number.getStrongPrime(bit_length // 2, 6)
    q = number.getStrongPrime(bit_length // 2, 6)

    # calculate the N
    N = p * q

    # get our boundary for the attack so we are able to compute
    # small dp's and dq's
    boundary = N ** 0.090

    # ger random smal dp and dq with the calculated boundary
    dp = number.getRandomInteger(number.size(int(boundary)))
    dq = number.getRandomInteger(number.size(int(boundary)))

    d = crt(dp, dq, p, q)

    # get the related e to our d
    try:
        e = inverse_mod(d, (p - 1) * (q - 1))

        keys = {
                'p': p,
                'q': q,
                'N': p * q,
                'e': e,
                'd': d,
                'dp': dp,
                'dq': dq
                }
    except ZeroDivisionError:
        keys = generate_keys(bit_length)

    return keys


def check_parameters(N, e, d):
    """
    Check the given RSA-parameters N, e, d. This function
    generates a random messge, encrypts and decrypts it in order
    to verify the RSA-parameters. If the decrypted and generated
    message are equal, this function return True.
    """
    message = number.getRandomInteger(number.size(N) - int(2))

    encrypted_message = sqm(message, e, N)

    decrypted_message = sqm(encrypted_message, d, N)

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
