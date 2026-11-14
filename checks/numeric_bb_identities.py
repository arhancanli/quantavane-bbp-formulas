"""Numerical check (independent of the certificates) of the Borwein-Broadhurst identities (28), (29), (30), (37):
  sum_k c_k Cl_2(k theta_D) = Z_D = |D| sum_{k=1}^{|D|/2} chi_D(k) Cl_2(2 pi k/|D|),
  theta_D = arg w_D, w_-11 = (-5+sqrt-11)/6, w_-15 = (-1+sqrt-15)/4, w_-20 = (-2+sqrt-5)/3, w_-24 = (-1+2sqrt-6)/5,
i.e. theta = 2 arctan sqrt11, 2 arctan sqrt(5/3), 2 arctan sqrt5, 2 arctan sqrt(3/2). Red control: one coefficient + 1.
Usage: python -I numeric_bb_identities.py digits"""
import sys
import mpmath as mp
dig = int(sys.argv[1]) if len(sys.argv) > 1 else 220
mp.mp.dps = dig + 20


def kron(d, n):
    from sympy import jacobi_symbol, factorint
    r = 1
    for p, e in factorint(n).items():
        if p == 2:
            v = 0 if d % 2 == 0 else (1 if d % 8 in (1, 7) else -1)
        else:
            v = 0 if d % p == 0 else jacobi_symbol(d % p, p)
        r *= v ** e
    return r


cases = {-11: (2 * mp.atan(mp.sqrt(11)), {1: 15, 2: -10, 5: 1}),
         -15: (2 * mp.atan(mp.sqrt(mp.mpf(5) / 3)), {1: 24, 2: -12, 3: -8, 4: 6}),
         -20: (2 * mp.atan(mp.sqrt(5)), {1: 36, 2: -30, 3: 4, 4: 3}),
         -24: (2 * mp.atan(mp.sqrt(mp.mpf(3) / 2)), {1: 60, 2: -18, 3: -4, 4: 3})}
for D, (th, cf) in cases.items():
    n = -D
    Z = n * mp.fsum(kron(D, k) * mp.clsin(2, 2 * mp.pi * k / n) for k in range(1, n // 2 + 1))
    lhs = mp.fsum(c * mp.clsin(2, k * th) for k, c in cf.items())
    red = lhs + mp.clsin(2, th) - Z
    print(f"D = {D}: |lhs - Z_D| = {mp.nstr(abs(lhs - Z), 3)} at {dig} digits; red control (c_1 + 1): {mp.nstr(red, 8)}")
