"""Exact checks of the elementary steps (arguments, Fourier coefficients, BBP coefficients) in Q(zeta_40).
z_a = zeta_40^a / sqrt2.  All identities below are checked in exact rational arithmetic (power basis of Q(zeta_40));
numeric values are used only to fix the branch of an argument (a sign test far from 0).
Usage: python3 -I verify_steps.py"""
import sys, cmath, math
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from fractions import Fraction as Fr
from cyclo import CycloField
K = CycloField(40); T = K.zeta; ONE = K.one(); ZERO = K.elt([0])
def conj(x):
    r = ZERO
    for i, c in enumerate(x):
        if c: r = K.add(r, K.scal(c, T(-i)))
    return r
SQ2 = K.add(T(5), T(-5))                                  # zeta_8 + zeta_8^-1
I = T(10)
SQ5 = K.sub(K.add(T(8), T(32)), K.add(T(16), T(24)))      # Gauss sum: z5 - z5^2 - z5^3 + z5^4
assert K.is_zero(K.sub(K.mul(SQ2, SQ2), K.scal(2, ONE))) and K.is_zero(K.sub(K.mul(SQ5, SQ5), K.scal(5, ONE)))
INV_SQ2 = K.scal(Fr(1, 2), SQ2)
z = lambda a: K.mul(T(a), INV_SQ2)
ok = True
def check(name, cond):
    global ok
    print(("OK   " if cond else "FAIL ") + name); ok &= bool(cond)

# Step 2: argument sums.  W real-positive after rotation  <=>  W*rot - conj(W*rot) = 0 and numeric value > 0.
def argcheck(name, factors, rot_exp):
    W = ONE
    for a, s in factors:
        f = K.sub(ONE, z(a))
        W = K.mul(W, f if s > 0 else conj(f))
    R = K.mul(W, T(rot_exp))
    real = K.is_zero(K.sub(R, conj(R)))
    v = K.num(R)
    # principal args of 1 - z lie in (-pi/2, pi/2) because Re(1 - z) > 0 for |z| < 1
    S = sum(s * cmath.phase(1 - K.num(z(a))) for a, s in factors)
    target = -2 * math.pi * rot_exp / 40
    check(f"{name}: (prod) * zeta_40^{rot_exp} is real (exact) and positive ({v.real:.3f}); sum of args = {S/math.pi:.12f} pi, claimed {target/math.pi:.4f} pi",
          real and v.real > 0 and abs(S - target) < 1e-9)
argcheck("A : arg(1-z1)+arg(1-z7)+arg(1-z9)-arg(1-z17) = -pi/2", [(1, 1), (7, 1), (9, 1), (17, -1)], 10)
argcheck("B1: arg(1-z7)+arg(1-z17) = -3pi/10", [(7, 1), (17, 1)], 6)
argcheck("B2: arg(1-z1)-arg(1-z17) = -pi/20", [(1, 1), (17, -1)], 1)
argcheck("B3: arg(1-z9)-arg(1-z17) = -3pi/20", [(9, 1), (17, -1)], 3)

# Step 3: sum_{b in {1,3,7,9}} sin(2 pi b n/20) = sqrt5 chi_-20(n) for n = 1..20 (period 20)  [exact]
chi = {1: 1, 3: 1, 7: 1, 9: 1, 11: -1, 13: -1, 17: -1, 19: -1}
good = True
for n in range(1, 21):
    lhs = ZERO
    for b in (1, 3, 7, 9):
        lhs = K.add(lhs, K.sub(T(2 * b * n), T(-2 * b * n)))      # 2i sin(2 pi b n /20)
    rhs = K.scal(2 * chi.get(n % 20, 0), K.mul(I, SQ5))
    good &= K.is_zero(K.sub(lhs, rhs))
check("Step 3: sum_{b=1,3,7,9} sin(pi b n/10) = sqrt5 * chi_-20(n), n = 1..20", good)

# Step 4: s_k = sin(pi k/20)+sin(7 pi k/20)+sin(9 pi k/20)-sin(17 pi k/20) = sqrt5 * c_k * 2^(k/2 - ceil(k/2)) ... i.e.
#   s_k = sqrt5 c_k (k even),  s_k = (sqrt10/2) c_k (k odd),  c periodic mod 40 with c_{k+20} = -c_k
c20 = [1, 1, -1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, -1, 1, 1, 0]
c = c20 + [-x for x in c20]
good = True
for k in range(1, 41):
    lhs = ZERO
    for a, s in ((1, 1), (7, 1), (9, 1), (17, -1)):
        lhs = K.add(lhs, K.scal(s, K.sub(T(a * k), T(-a * k))))     # 2i s_k
    fac = SQ5 if k % 2 == 0 else K.scal(Fr(1, 2), K.mul(SQ2, SQ5))
    rhs = K.scal(2 * c[k - 1], K.mul(I, fac))
    good &= K.is_zero(K.sub(lhs, rhs))
check("Step 4: s_k = sqrt5 c_k (k even), (sqrt10/2) c_k (k odd), k = 1..40", good)
print("ALL STEPS VERIFIED EXACTLY" if ok else "SOME STEP FAILED")
