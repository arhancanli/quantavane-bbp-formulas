"""Exact finite steps of the proof of Theorem 6.3 (the chi_-4 direction on the orbit of zeta_40/sqrt2), in Q(zeta_40).

(a) h = (1 - z_1)(1 - z_9)(1 - z_17) conj(1 - z_7) conj(1 - i/2) is real (h = conj h exactly); its sign is printed.
(b) t_m = sin(pi m/20) - sin(7 pi m/20) + sin(9 pi m/20) + sin(17 pi m/20) - 2 sin(pi m/4) [m even] satisfies, exactly,
    t_m 2^(-m/2) = g_m and t_(m+20) 2^(-(m+20)/2) = -2^(-10) g_m for m = 1..20, with the g of Theorem 6.3.
Pure Python (fractions); run: python3 -I prove_thm_G.py
"""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from fractions import Fraction as Fr
from cyclo import CycloField

K = CycloField(40)
Z = K.zeta
ONE = K.one()
I = Z(10)
SQ2 = K.add(Z(5), Z(-5))                       # zeta_8 + zeta_8^-1 = sqrt 2
assert K.is_zero(K.sub(K.mul(SQ2, SQ2), K.elt([2])))
ISQ2 = K.inv(SQ2)


def conj(x):
    r = K.elt([0])
    for i, c in enumerate(x):
        if c:
            r = K.add(r, K.scal(c, Z(-i)))
    return r


def z(a):
    return K.mul(Z(a), ISQ2)


half_i = K.scal(Fr(1, 2), I)
h = K.mul(K.mul(K.mul(K.sub(ONE, z(1)), K.sub(ONE, z(9))), K.sub(ONE, z(17))),
          K.mul(conj(K.sub(ONE, z(7))), conj(K.sub(ONE, half_i))))
real = K.is_zero(K.sub(h, conj(h)))
print("(a) h real exactly:", real, "; h = 5/8 exactly:", K.is_zero(K.sub(h, K.elt([Fr(5, 8)]))))

INV2I = K.inv(K.scal(2, I))


def sin_m(a, m):            # sin(pi a m / 20) = (zeta^(a m) - zeta^(-a m)) / (2 i)
    return K.mul(K.sub(Z(a * m), Z(-a * m)), INV2I)


def t(m):
    v = K.add(K.sub(sin_m(1, m), sin_m(7, m)), K.add(sin_m(9, m), sin_m(17, m)))
    if m % 2 == 0:
        v = K.sub(v, K.scal(2, sin_m(5, m)))       # sin(pi m/4) = sin(pi * 5m / 20)
    return v


def pow_isq2(m):
    r = ONE
    for _ in range(m):
        r = K.mul(r, ISQ2)
    return r


g = [Fr(1, 2), Fr(-3, 2), Fr(1, 4), 0, Fr(1, 2), Fr(3, 8), Fr(-1, 16), 0, Fr(1, 32), Fr(1, 16), Fr(1, 64), 0,
     Fr(-1, 128), Fr(3, 128), Fr(1, 64), 0, Fr(1, 512), Fr(-3, 512), Fr(1, 1024), 0]
ok = True
for m in range(1, 21):
    lhs1 = K.mul(t(m), pow_isq2(m))
    lhs2 = K.mul(t(m + 20), pow_isq2(m + 20))
    ok &= K.is_zero(K.sub(lhs1, K.elt([g[m - 1]])))
    ok &= K.is_zero(K.sub(lhs2, K.elt([-Fr(g[m - 1]) / 1024])))
print("(b) t_m 2^(-m/2) = g_m and t_(m+20) 2^(-(m+20)/2) = -g_m/1024 for m = 1..20 (exact):", ok)
print("THEOREM 6.3 STEPS VERIFIED EXACTLY" if (real and ok and K.is_zero(K.sub(h, K.elt([Fr(5, 8)])))) else "FAILED")
