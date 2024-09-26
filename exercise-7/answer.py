#!/usr/bin/env python3

import scipy
from numpy import cos, sin, exp


a = 1.3
b = 2.2
alpha = 0.0
beta = 5.0 / 6.0


def f(x):
    return (
        4.0 * cos(0.5 * x) * exp(-5.0 * x / 4.0)
        + 2.0 * sin(4.5 * x) * exp(x / 8.0)
        + 2.0
    )


ans = scipy.integrate.quad(f, a, b)
ansp = scipy.integrate.quad(f, a, b, weight="alg", wvar=(alpha, beta))

print(f"p=  1, ans={ans[0]:20.16f}, err={ans[1]:20.16}")
print(f"p=alg, ans={ansp[0]:20.16f}, err={ansp[1]:20.16}")
