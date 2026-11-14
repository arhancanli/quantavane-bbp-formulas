import sys, mpmath as mp
from math import gcd
from sympy import factorint
mp.mp.dps = 50
M = int(sys.argv[1]); b = int(sys.argv[2]); pts = [tuple(int(x) for x in p.split(':')) for p in sys.argv[3].split(',')]
# sqrt(b) as element: use Galois action via character of Q(sqrt b): compute sign by evaluating Gauss sum numerically
def sb_conj(t):
    # sigma_t(sqrt b) where sqrt b = sum over its Gauss-sum expression; use numeric: sqrt b in Q(zeta_M) via 2cos for b=2 (M%8==0)
    if b == 2: return mp.sqrt(2) * (1 if t % 8 in (1, 7) else -1)
    if b == 3: return mp.sqrt(3) * (1 if t % 12 in (1, 11) else -1)
    if b == 5: return mp.sqrt(5) * (1 if t % 5 in (1, 4) else -1)
    if b == 6: return mp.sqrt(6) * (1 if t % 24 in (1, 5, 19, 23) else -1)
    if b == 7: return mp.sqrt(7) * (1 if t % 28 in (1, 3, 9, 19, 25, 27) else -1)
    raise
for a, k in pts:
    pr = mp.mpf(1)
    for t in range(1, M):
        if gcd(t, M) != 1: continue
        pr *= abs(sb_conj(t) ** k - mp.expjpi(mp.mpf(2 * a * t) / M))
    v = int(mp.nint(pr))
    z = mp.expjpi(mp.mpf(2*a)/M)/mp.sqrt(b)**k
    print(f"a={a} k={k}: N(b^(k/2) - zeta_{M}^{a}) = {v} = {factorint(v)}  arg(1-z)/pi = {mp.nstr(mp.arg(1-z)/mp.pi, 12)}")
