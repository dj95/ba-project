# ba-project

This repository contains the project for my bachelor thesis.


### Summary

This is an implementation of the modified coppersmith attack
on CRT-RSA with small decoding-exponents by Takayasu, Lu and Peng ([Link to their paper](https://eprint.iacr.org/2017/092.pdf)).


This implemtation of the attack if for educational purposes only!
Please bear in mind that this project is still unfinished!


### Requirements

- sagemath 8.1
- Python 2.7
- Python 2.7 modules:
  - colorama
  - crypto


### Usage

  ```
    usage: main.sage.py [-h] [-s S] [-d D] [-m M] [-t T] [--debug] [--test]                                                                                                                                                         
                        [--nogroebner] [--noreduction] [--forcetriangle] [--print]
                        [--json]

    optional arguments:
      -h, --help       show this help message and exit
      -s S             bit size of N [default=1024]
      -d D             delta [default=0.010]
      -m M             m [default=8]
      -t T             tau [default=0.75]
      --debug          Print debug messages
      --test           Test the LLL algorithm
      --nogroebner     Skip the groebner basis
      --noreduction    Skip the LLL reduction
      --forcetriangle  Force a triangular matrix
      --print          Print the matrix to tech
      --json           JSON output
  ```

