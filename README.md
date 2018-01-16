# ba-project

This repository contiains the project for my bachelor thesis.


### Requirements

- sagemath 8.1
- Python 2.7
- Python 2.7 modules:
  - colorama
  - crypto


### Usage

  ```
  usage: main.sage.py [-h] [-s S] [-d D] [-m M] [-t T] [--debug] [--test]
                    [--nogroebner] [--noreduction]

  optional arguments:
    -h, --help     show this help message and exit
    -s S           bit size of N [default=1024]
    -d D           delta [default=0.010]
    -m M           m [default=8]
    -t T           tau [default=0.75]
    --debug        Print debug messages
    --test         Test the LLL algorithm
    --nogroebner   Skip the groebner basis
    --noreduction  Skip the LLL reduction

  ```

