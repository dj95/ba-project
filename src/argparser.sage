#!/usr/bin/env sage
#
# ba-proj
#
# This project is the implementation of the small
# CRT-RSA attack of Takayusa, Lu and Peng utilizing
# the coppersmith and L3 algorithm.
#
# (c) 2018 - Daniel Jankowski



import argparse


def parse_args():
    """
    Parse the command line arguments for this project.
    """
    # initialize the argument parser
    parser = argparse.ArgumentParser()

    # add the different options
    parser.add_argument(
            '-s',
            help='bit size of N [default=1024]',
            default=1024,
            type=int
            )
    parser.add_argument(
            '-d',
            help='delta [default=0.010]',
            default=0.010,
            type=float
            )
    parser.add_argument(
            '-m',
            help='m [default=8]',
            default=8,
            type=int
            )
    parser.add_argument(
            '-t',
            help='tau [default=0.75]',
            default=0.75,
            type=float
            )
    parser.add_argument(
            '--debug',
            help='Print debug messages',
            dest='debug',
            action='store_const',
            const=True,
            default=False
            )
    parser.add_argument(
            '--test',
            help='Test the LLL algorithm',
            dest='test',
            action='store_const',
            const=True,
            default=False
            )
    parser.add_argument(
            '--nogroebner',
            help='Skip the groebner basis',
            dest='nogroebner',
            action='store_const',
            const=True,
            default=False
            )
    parser.add_argument(
            '--noreduction',
            help='Skip the lll reduction',
            dest='noreduction',
            action='store_const',
            const=True,
            default=False
            )
    
    #TODO: add testing flag

    # parse the arguments
    args = parser.parse_args()

    # and return them
    return args.d, args.m, args.s, args.t, args.debug, args.test, args.nogroebner, args.noreduction
