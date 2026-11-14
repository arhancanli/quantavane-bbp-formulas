"""Exact verification of a five-term-relation certificate for a Bloch-Wigner identity in Q(zeta_M).

Input JSON (from cert_search.py): M, cert = [{coef, x, y}], recipes (exact construction of each element), target =
[[a, k, c], ...] meaning sum of c * D(zeta_M^a * (1/sqrt2)^k).  For each certificate entry the five-term relation
  D(x) + D(y) + D((1-x)/(1-xy)) + D(1-xy) + D((1-y)/(1-xy)) = 0
holds (Bloch-Wigner; valid for all complex x, y off the degenerate set, which is checked).  Every argument is reduced
EXACTLY (rational coordinates in the power basis of Q(zeta_M)) to a class representative using
  D(w) = D(1-1/w) = D(1/(1-w)) = -D(1/w) = -D(1-w) = -D(w/(w-1)),  D(conj w) = -D(w),  D(real) = 0,
and the weighted sum of all relations is compared, class by class and in exact arithmetic, with the target.
Usage: python3 -I verify_cert_exact.py cert.json
"""
import sys, json
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from fractions import Fraction as Fr
from cyclo import CycloField

data = json.load(open(sys.argv[1]))
M = data['M']
K = CycloField(M)
T = K.zeta
ONE = K.one()
ZERO = K.elt([0])
HALF_SQ2 = K.scal(Fr(1, 2), K.add(T(M // 8), T(-M // 8)))   # 1/sqrt2 = (zeta_8 + zeta_8^-1)/2

def conj(x):
    res = ZERO
    for i, c in enumerate(x):
        if c:
            res = K.add(res, K.scal(c, T(-i)))
    return res

def bbp_point(a, k):
    x = T(a)
    for _ in range(k):
        x = K.mul(x, HALF_SQ2)
    return x

recipes = data['recipes']
memo = {}
def build(i):
    if i in memo:
        return memo[i]
    r = recipes[i]
    if r[0] == 'bbp':
        x = bbp_point(r[1], r[2])
    elif r[0] == 'sym':
        w = build(r[1]); si = r[2]
        if si == 0: x = w
        elif si == 1: x = K.sub(ONE, K.inv(w))
        elif si == 2: x = K.inv(K.sub(ONE, w))
        elif si == 3: x = K.inv(w)
        elif si == 4: x = K.sub(ONE, w)
        elif si == 5: x = K.div(w, K.sub(w, ONE))
        if r[3]:
            x = conj(x)
    else:
        raise ValueError(r)
    memo[i] = x
    return x

def numkey(x):
    v = K.num(x)
    return (round(v.real, 9) + 0.0, round(v.imag, 9) + 0.0)

def canon(w):
    if K.is_zero(K.sub(w, conj(w))):
        return None, 0
    inv = K.inv(w); omw = K.sub(ONE, w)
    imgs = [(w, 1), (K.sub(ONE, inv), 1), (K.inv(omw), 1), (inv, -1), (omw, -1), (K.div(w, K.sub(w, ONE)), -1)]
    best = None
    for v, s in imgs:
        for vv, ss in ((v, s), (conj(v), -s)):
            k = numkey(vv)
            if best is None or k < best[0]:
                best = (k, vv, ss)
    return best[1], best[2]

total = {}
def addterm(w, c):
    cv, s = canon(w)
    if cv is not None:
        total[cv] = total.get(cv, Fr(0)) + c * s

for item in data['cert']:
    c = Fr(item['coef'])
    x, y = build(item['x']), build(item['y'])
    xy = K.mul(x, y)
    for bad in (x, y, K.sub(ONE, x), K.sub(ONE, y), K.sub(ONE, xy)):
        assert not K.is_zero(bad), "degenerate five-term relation"
    d = K.sub(ONE, xy)
    for a in (x, y, K.div(K.sub(ONE, x), d), d, K.div(K.sub(ONE, y), d)):
        addterm(a, c)

claim = {}
for a, k, c in data['target']:
    cv, s = canon(bbp_point(a, k))
    if cv is not None:
        claim[cv] = claim.get(cv, Fr(0)) + Fr(c) * s
keys = set(total) | set(claim)
bad = {k: total.get(k, Fr(0)) - claim.get(k, Fr(0)) for k in keys}
bad = {k: v for k, v in bad.items() if v != 0}
print(f"field Q(zeta_{M}), relations {len(data['cert'])}, classes involved {len(keys)}")
if not bad:
    print("CERTIFICATE VERIFIED EXACTLY: sum_i c_i (five-term relation_i) = target combination")
else:
    print("MISMATCH in", len(bad), "classes")
