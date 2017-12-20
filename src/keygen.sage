#!/bin/env sage
#
# keygen for crt rsa
#
# (c) 2017 Daniel Jankowski


from Crypto.Util import number
from functools import reduce
from math import log


def eea(a, b):
    """
    The extended euclidian algorithm with a and b as input.
    """
    if (a==0): return (b, 0, 1)
    else: g, y, x = eea(b % a, a); return (g, x-(b // a)*y,y)


def modinv(a, m):
    """
    Returns the modular mulitplicative inverse of a in m.
    Return 'nicht existent' if the inverse doesn't exist.
    """
    gcd, x, y = eea(a, m)
    if (gcd!=1): return "nicht existent"
    else: return x % m


def chinese_remainder_theorem(x, moduln):
    """
    Solve modular equations:

        z = x[0] mod moduln[0]
        z = x[1] mod moduln[2]
        ...
        z = x[n] mod moduln[n]

    and return z.
    """
    # initial values
    phi, z = 1, 0

    # get the common modulus
    prod = reduce(lambda a, b: a*b, moduln)

    # calculate the crt
    for i in range(len(x)):
        phi = prod // moduln[i]
        z += x[i]  * modinv(phi, moduln[i]) * phi

    # return our value we searched for
    return z % prod


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

    # calculate our d with the random dp and dq
    d = chinese_remainder_theorem(
            [dp, dq],
            [p, q]
            )

    # get the related e to our d
    e = modinv(d, (p - 1) * (q - 1))
    
    # check if the modular inverse exist and we have all our parameters
    if e != 'nicht existent':
        return {
                'p': p,
                'q': q,
                'N': p * q,
                'e': e,
                'd': d,
                'dp': dp,
                'dq': dq
                }
    else:
        # if e is non existant because of the modular inverse
        # generate new values
        return generate_keys(bit_length)


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
