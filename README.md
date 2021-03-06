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
  - time
  - datetime


### Further tasks

- [x] code capsulation for loops
- [ ] better tests and plots
- [ ] more object orientated programming
- [ ] parameter check
- [ ] d in [dp,dq]: disp. hamming weight
- [ ] convert to pure python
- [ ] sanity checks
- [x] N^delta for big N
- [ ] dp == dq
- [x] solver timing + single timings
- [ ] lattice as csv(reduced/unreduced)
- [x] lattice as heatmap
- [ ] optimize matrix
- [ ] optimize polynomials
- [ ] optimize bound
- [ ] code clean up and more readability


### Usage

```
usage: ./main.sage [-h] [-s S] [-d D] [-m M] [-t T] [--debug] [--print]
                   [--json] [--deltagen]

optional arguments:
  -h, --help  show this help message and exit
  -s S        bit size of N [default=1024]
  -d D        delta [default=0.010]
  -m M        m [default=8]
  -t T        tau [default=0.75]
  --debug     Print debug messages
  --print     Print the matrix to tech
  --json      JSON output
  --deltagen  Generate table with theoretical deltas
```


### Usage with docker

Run `docker-compose up`.


If you want to modify the parameters, check out the
`./docker-compose.yml`.
