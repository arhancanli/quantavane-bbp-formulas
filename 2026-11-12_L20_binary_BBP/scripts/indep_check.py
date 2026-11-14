"""Independent recomputation + red controls for
   L(2,chi_-20) - (sqrt5/15) pi log2 = (4/3) sum_{k>=1} c_k / (2^ceil(k/2) k^2),  c periodic mod 40, c_{k+20} = -c_k.
BBP side: pure Python big-integer fixed point, summed in REVERSE order (no mpmath).
Constant side: mpmath (Hurwitz zeta for L, mp.pi, mp.log).  Usage: python -I indep_check.py digits"""
import sys, time
import mpmath as mp
P = int(sys.argv[1])
c20 = [1, 1, -1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, -1, 1, 1, 0]
def c(k):
    r = (k - 1) % 40
    return c20[r] if r < 20 else -c20[r - 20]
def bbp_fixed(coef, P):
    """sum_k coef(k) / (2^ceil(k/2) k^2) as an integer scaled by 10^(P+10), reverse order."""
    scale = 10 ** (P + 10)
    K = int((P + 12) * 3.3219280948873626 * 2) + 4
    tot = 0
    for k in range(K, 0, -1):
        ck = coef(k)
        if ck:
            tot += ck * (scale // ((1 << ((k + 1) // 2)) * k * k))
    return tot
t0 = time.time()
S = bbp_fixed(c, P)
print(f"BBP fixed-point sum ({time.time()-t0:.1f}s)")
mp.mp.dps = P + 30
chi = {1: 1, 3: 1, 7: 1, 9: 1, 11: -1, 13: -1, 17: -1, 19: -1}
L = mp.fsum(v * mp.zeta(2, mp.mpf(r) / 20) for r, v in chi.items()) / 400
C = L - mp.sqrt(5) / 15 * mp.pi * mp.log(2)
lhs = mp.mpf(4) / 3 * mp.mpf(S) / mp.mpf(10) ** (P + 10)
d = C - lhs
print(f"main identity at {P} digits: log10|diff| = {float(mp.log10(abs(d))) if d else '-inf'}")
# degree-1 sibling (Adegoke-Layeni 2013, compendium formula 22): sum c_k/(2^ceil(k/2) k) = pi sqrt5 / 10
mp.mp.dps = 120
s1 = mp.fsum(c(k) / (mp.mpf(2) ** ((k + 1) // 2) * k) for k in range(1, 900))
print(f"degree-1 sibling: sum c_k/(2^ceil(k/2) k) - pi sqrt5/10 = {mp.nstr(s1 - mp.pi * mp.sqrt(5) / 10, 5)}")
# red controls (at 300 digits)
mp.mp.dps = 330
Cs = C
def val(coef, P=300):
    return mp.mpf(4) / 3 * mp.mpf(bbp_fixed(coef, P)) / mp.mpf(10) ** (P + 10)
base = val(c)
print("red controls (each must FAIL, i.e. |diff| far above 1e-300):")
for j in (1, 3, 19):
    cc = lambda k, j=j: c(k) + (1 if (k - 1) % 40 == j - 1 else 0)
    print(f"  c_{j} -> c_{j}+1 (period 40): |diff| = {mp.nstr(abs(Cs - val(cc)), 5)}")
print(f"  companion sqrt5/15 -> sqrt5/16: |diff| = {mp.nstr(abs(L - mp.sqrt(5)/16*mp.pi*mp.log(2) - base), 5)}")
print(f"  drop companion: |diff| = {mp.nstr(abs(L - base), 5)}")
print(f"  prefactor 4/3 -> 5/4: |diff| = {mp.nstr(abs(Cs - base * mp.mpf(15) / 16), 5)}")
print(f"  unperturbed at 300 digits: |diff| = {mp.nstr(abs(Cs - base), 5)}")
