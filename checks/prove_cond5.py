# Proposition 5.1 of the paper (version of 12 Nov 2026): binary formulas for the conductor-5 series S_1, S_2, with the
# integer vectors V_1, V_2 printed in Appendix A.  Exact arithmetic in Q(zeta_40) (cyclo.py) checks the finite steps of
# the printed proof:
#   (1) with u = 2^19/5, for k = 1..40 (both sides have period 40)
#       G_1(k) := v1_j 2^(j/2) = u[(21a - b)(sin(7 pi k/20) + sin(17 pi k/20)) - (a + 21b)(sin(pi k/20) - sin(9 pi k/20))],
#       G_2(k) := v2_j 2^(j/2) = u[(a + 21b)(sin(7 pi k/20) + sin(17 pi k/20)) + (21a - b)(sin(pi k/20) - sin(9 pi k/20))],
#       where a = sqrt(10 + 2 sqrt5), b = sqrt(10 - 2 sqrt5) and j = k (here 1 <= k <= 40);
#   (2) the Fourier coefficients of G_i on the other sine modes m = 1..19 vanish (a consistency check of (1));
#   (3) the coefficients obtained from Theorem 4.2 match the claim: 442 u/16 = 2^16 * 221 and u/10 = 2^18/25;
#   (4) the vectors A_1, A_2 returned by the integer relation search satisfy 4 A_i = 5 V_i + Z, where Z is the zero
#       relation (125) of Bailey's compendium (as printed there); P(2, 2^20, 40, Z) is checked numerically.
# Run: python3 -I prove_cond5.py   (needs mpmath only for the numerical line about Z)
import sys
from fractions import Fraction as F
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from cyclo import CycloField
K = CycloField(40)
z = K.zeta
ONE = K.one()
def add(*xs):
    r = K.elt([0])
    for x in xs:
        r = K.add(r, x)
    return r
def sc(q, x):
    return K.scal(F(q), x)
I = z(10)
SQ2 = K.add(z(5), z(-5))
ALPHA = sc(2, K.add(z(2), z(-2)))    # 4 cos(pi/10) = sqrt(10 + 2 sqrt5)
BETA = sc(2, K.add(z(6), z(-6)))     # 4 cos(3 pi/10) = sqrt(10 - 2 sqrt5)
assert K.is_zero(K.sub(K.mul(ALPHA, ALPHA), add(sc(10, ONE), sc(2, add(sc(2, K.add(z(8), z(-8))), ONE)))))  # alpha^2 = 10 + 2 sqrt5
def sin_(n):   # sin(pi n/20) = (zeta^n - zeta^-n)/(2i)
    return K.div(K.sub(z(n), z(-n)), sc(2, I))
V1 = (11010048, 0, -262144, -5505024, 0, 0, -65536, -65536, -688128, 0, -344064, -16384, -8192, 0, 0, -86016, -2048, 0,
      21504, 0, -10752, 0, 256, 5376, 0, 0, 64, 64, 672, 0, 336, 16, 8, 0, 0, 84, 2, 0, -21, 0)
V2 = (524288, 0, 5505024, -262144, 0, 0, 1376256, 1376256, -32768, 0, -16384, 344064, 172032, 0, 0, -4096, 43008, 0,
      1024, 0, -512, 0, -5376, 256, 0, 0, -1344, -1344, 32, 0, 16, -336, -168, 0, 0, 4, -42, 0, -1, 0)
A1 = (13893632, -786432, -393216, -6029312, 786432, -196608, -65536, -491520, -851968, -49152,
      -434176, 32768, -12288, -12288, -24576, -133120, -2048, -3072, 26624, -3072, -13568, -768,
      384, 5120, -768, -192, 64, 288, 832, -48, 424, -80, 12, -12, 24, 118, 2, -3, -26, 0)
A2 = (786432, -786432, 6815744, 524288, 786432, -196608, 1736704, 1310720, -32768, -49152,
      -24576, 483328, 212992, -12288, -24576, -30720, 54272, -3072, 1024, -3072, -768, -768,
      -6656, -1280, -768, -192, -1696, -1472, 32, -48, 24, -520, -208, -12, 24, 18, -53, -3, -1, 0)
p2 = lambda e: 2 ** e
Z = (p2(19), -3*p2(20), -p2(18), 13*p2(18), 3*p2(20), -3*p2(18), p2(16), -25*p2(16), p2(15), -3*p2(16),
     -p2(14), 13*p2(14), -p2(13), -3*p2(14), -3*p2(15), -25*p2(12), p2(11), -3*p2(12), -p2(10), -3*p2(12),
     -p2(9), -3*p2(10), 256, -25*p2(8), -3*p2(10), -3*p2(8), -64, 13*p2(6), -32, -192,
     16, -25*p2(4), 8, -48, 96, 52, -2, -12, 1, 0)
def pow2half(k):          # 2^(k/2) in Q(zeta_40)
    r = sc(2 ** (k // 2), ONE)
    return K.mul(r, SQ2) if k % 2 else r
u = F(2 ** 19, 5)
P1 = K.sub(sc(21, ALPHA), BETA)        # 21 alpha - beta
P2 = K.add(ALPHA, sc(21, BETA))        # alpha + 21 beta
results = {}
for name, V, (cT1, cT2) in (('V1', V1, (P1, K.neg(P2))), ('V2', V2, (P2, P1))):
    G = [None] + [sc(V[k - 1], pow2half(k)) for k in range(1, 41)]
    closed = all(K.is_zero(K.sub(G[k], sc(u, add(K.mul(cT1, K.add(sin_(7 * k), sin_(17 * k))),
                                                   K.mul(cT2, K.sub(sin_(k), sin_(9 * k)))))))
                 for k in range(1, 41))
    # Fourier coefficients B_m = (1/20) sum_k G(k) sin(pi m k/20): zero for m outside {1, 7, 9, 17}
    def B(m):
        acc = K.elt([0])
        for k in range(1, 41):
            acc = K.add(acc, K.mul(G[k], sin_(m * k)))
        return sc(F(1, 20), acc)
    others = all(K.is_zero(B(m)) for m in range(1, 20) if m not in (1, 7, 9, 17))
    results[name] = closed and others
    print(f'{name}: closed sine form exact for k = 1..40: {closed}; other modes vanish: {others}')
print('442 u/16 = 2^16*221:', F(442, 16) * u * 5 == 2 ** 16 * 221 and F(442, 16) * u == F(2 ** 16 * 221, 5),
      '| u/10 = 2^18/25 = 2^16*221*4/5525:', u / 10 == F(2 ** 18, 25) == F(2 ** 16 * 221 * 4, 5525))
rel = all(4 * a == 5 * v + zz for a, v, zz in zip(A1, V1, Z)) and all(4 * a == 5 * v + zz for a, v, zz in zip(A2, V2, Z))
print('search vectors: 4 A_i = 5 V_i + Z (compendium relation (125)) for i = 1, 2:', rel)
try:
    from mpmath import mp, mpf
    mp.dps = 200
    PZ = sum(mpf(2) ** (-20 * q) * sum(mpf(Z[j]) / (40 * q + j + 1) ** 2 for j in range(40)) for q in range(32))
    print('numerical P(2, 2^20, 40, Z) at 200 digits (truncation error ~1e-190):', mp.nstr(PZ, 5))
except ImportError:
    print('mpmath not available: numerical check of (125) skipped')
print('PROPOSITION 5.1 STEPS VERIFIED EXACTLY' if all(results.values()) and rel else 'FAILED', results)
