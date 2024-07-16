from decimal import Decimal as deci
import numpy as np
from scipy.special import comb

delta = 10

ns = (512, 3488, 4608, 6688, 6960, 8192)
ts = ( 20,   64,   96,  128,  119,  128)
ms = (  9,   12,   13,   13,   13,   13)
ds = ( 18,   22,   22,   22,   22,   22)

print("\t" + "\t".join([f"{n}" for n in ns]))

accuracies = [0.963, 0.962, 0.961, 0.960,
              0.959, 0.958, 0.957, 0.956, 0.955, 0.954, 0.953, 0.952, 0.951, 0.950,
              0.949, 0.948, 0.947, 0.946, 0.945, 0.944, 0.943, 0.942, 0.941, 0.940,
              0.939, 0.938, 0.937, 0.936, 0.935, 0.934, 0.933, 0.932, 0.931, 0.930,
              0.929, 0.928, 0.927, 0.926, 0.925]

for accuracy in accuracies:
    print(accuracy, end="\t")
    for n, t, m, d in zip(ns, ts, ms, ds):
        if m == 9:
            template_accuracy = 0.9714
        elif m == 12:
            template_accuracy = 0.9684
        else:
            template_accuracy = 0.9996
        p = template_accuracy * (accuracy ** d)
        cs   = [comb(n, i, exact=True)     for i in range(m*t+delta, n+1)]
        pos  = [deci(p) ** deci(i)         for i in range(m*t+delta, n+1)]
        npos = [deci(1 - p) ** deci(n - i) for i in range(m*t+delta, n+1)]
        success_rate = 100 * sum([c * po * npo for c, po, npo in zip(cs, pos, npos)])
        print(f"{success_rate:0.3f}", end="\t")
    print()

