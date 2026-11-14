# Exact facts used to fix the branches of arg in the proofs of Theorems 1.2 and 4.2 (version of 12 Nov 2026).
# z_a = zeta^a/sqrt2 with zeta = exp(pi i/20).  All checks are exact in Q(zeta_40) (cyclo.py); floating point is used
# only to print the approximate values of (b) and (c) for orientation.
#   (a) g = (1 - z_1)(1 - z_7)(1 - z_9) conj(1 - z_17) satisfies g * i = sqrt5/4;
#   (b) (1 - z_7)(1 - z_17) zeta^6 is real;   (c) (1 - z_1) conj(1 - z_17) zeta is real.
# With arg(1 - z_a) in (-pi/2, 0) (because Im z_a > 0 and |z_a| < 1), these fix the argument sums to -pi/2, -3pi/10
# and -pi/20 (see the proofs in the paper).
import sys
from fractions import Fraction as F
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from cyclo import CycloField
K = CycloField(40)
z = K.zeta
ONE = K.one()
I = z(10)
SQ2 = K.add(z(5), z(-5))
SQ5 = K.add(K.scal(F(2), K.add(z(8), z(-8))), ONE)       # sqrt5 = 1 + 4 cos(2 pi/5)
assert K.is_zero(K.sub(K.mul(SQ5, SQ5), K.scal(F(5), ONE)))
assert K.is_zero(K.sub(K.mul(SQ2, SQ2), K.scal(F(2), ONE)))
def conj(x):
    r = K.elt([0])
    for i, c in enumerate(x):
        if c:
            r = K.add(r, K.scal(c, z(-i)))
    return r
def is_real(x):
    return K.is_zero(K.sub(x, conj(x)))
za = lambda a: K.div(z(a), SQ2)
om = lambda a: K.sub(ONE, za(a))
g = K.mul(K.mul(K.mul(om(1), om(7)), om(9)), conj(om(17)))
ok_a = K.is_zero(K.sub(K.mul(g, I), K.scal(F(1, 4), SQ5)))
b = K.mul(K.mul(om(7), om(17)), z(6))
c = K.mul(K.mul(om(1), conj(om(17))), z(1))
ok_b, ok_c = is_real(b), is_real(c)
print('(a) g * i = sqrt5/4 exactly:', ok_a)
print('(b) (1 - z_7)(1 - z_17) zeta^6 is real:', ok_b, ' (numerically', round(K.num(b).real, 6), ')')
print('(c) (1 - z_1) conj(1 - z_17) zeta is real:', ok_c, ' (numerically', round(K.num(c).real, 6), ')')
print('ALL BRANCH FACTS VERIFIED EXACTLY' if ok_a and ok_b and ok_c else 'FAILED')
