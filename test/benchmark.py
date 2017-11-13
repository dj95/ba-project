from fpylll import IntegerMatrix, LLL
from multiprocessing import Pool

d, workers, tasks = 30, 4, 128

def run_it(p, f, A, prefix=""):
    """Print status during parallel execution."""
    import sys
    r = []
    for i, retval in enumerate(p.imap_unordered(f, A, 1)):
        r.append(retval)
        sys.stderr.write('\r{0} done: {1:.2%}'.format(prefix, float(i)/len(A)))
        sys.stderr.flush()
    sys.stderr.write('\r{0} done {1:.2%}\n'.format(prefix, float(i+1)/len(A)))
    return r

A = [IntegerMatrix.random(d, "uniform", bits=30) for _ in range(tasks)]
A = run_it(Pool(workers), LLL.reduction, A)
