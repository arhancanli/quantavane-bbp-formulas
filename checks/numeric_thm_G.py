"""Numerical check of Theorem 6.3 (written for v4; independent of the search and certificate programs).
  D(z_1) - D(z_7) + D(z_9) + D(z_17) = D(i/2)/2 + G/5,           z_a = exp(i pi a/20)/sqrt2
  Im Li2(z_1) - Im Li2(z_7) + Im Li2(z_9) + Im Li2(z_17) - Ti2(1/2)/2 = G/5
  G = 5 sum_{q>=0} (-1/1024)^q sum_{j=1}^{20} g_j/(20q+j)^2
Usage: python -I numeric_thm_G.py digits"""
import sys
import mpmath as mp

dig = int(sys.argv[1]) if len(sys.argv) > 1 else 1100
mp.mp.dps = dig + 20


def D(z):
    return mp.im(mp.polylog(2, z)) + mp.arg(1 - z) * mp.log(abs(z))


z = {a: mp.expjpi(mp.mpf(a) / 20) / mp.sqrt(2) for a in (1, 7, 9, 17)}
G = mp.catalan
h = mp.mpc(0, mp.mpf(1) / 2)
lhs_D = D(z[1]) - D(z[7]) + D(z[9]) + D(z[17])
r1 = lhs_D - D(h) / 2 - G / 5
ImL = lambda w: mp.im(mp.polylog(2, w))
Ti2_half = mp.nsum(lambda k: (-1) ** k * mp.mpf(2) ** (-(2 * k + 1)) / (2 * k + 1) ** 2, [0, mp.inf])
r2 = ImL(z[1]) - ImL(z[7]) + ImL(z[9]) + ImL(z[17]) - Ti2_half / 2 - G / 5
r2b = ImL(h) - Ti2_half
g = [mp.mpf(x) for x in ('1/2', '-3/2', '1/4', '0', '1/2', '3/8', '-1/16', '0', '1/32', '1/16', '1/64', '0',
                         '-1/128', '3/128', '1/64', '0', '1/512', '-3/512', '1/1024', '0')] if False else None
from fractions import Fraction as Fr
gf = [Fr(1, 2), Fr(-3, 2), Fr(1, 4), 0, Fr(1, 2), Fr(3, 8), Fr(-1, 16), 0, Fr(1, 32), Fr(1, 16), Fr(1, 64), 0,
      Fr(-1, 128), Fr(3, 128), Fr(1, 64), 0, Fr(1, 512), Fr(-3, 512), Fr(1, 1024), 0]
gm = [mp.mpf(Fr(x).numerator) / Fr(x).denominator for x in gf]
nq = int(dig / 3.0) + 10
S = mp.fsum(mp.mpf(-1) ** q / mp.mpf(1024) ** q * mp.fsum(gm[j - 1] / (20 * q + j) ** 2 for j in range(1, 21))
            for q in range(nq))
r3 = 5 * S - G
red = lhs_D - D(h) / 2 - G / 4
fmt = lambda x: mp.nstr(abs(x), 3) if x != 0 else '0'
print(f"digits {dig}:")
print(f"  D form residual        {fmt(r1)}")
print(f"  Im Li2 form residual   {fmt(r2)}   (Im Li2(i/2) - Ti2(1/2) = {fmt(r2b)})")
print(f"  series form residual   {fmt(r3)}   ({nq} terms)")
print(f"  red control (G/5 -> G/4): {mp.nstr(red, 6)}")
