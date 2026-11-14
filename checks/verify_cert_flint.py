"""Second exact verifier for the certificates of Theorems 6.3 and 6.4 (any Q(zeta_M)), written separately from
verify_cert_exact.py.

Differences from the first verifier: field arithmetic in Q(zeta_40) = Q[x]/(Phi_40) with FLINT polynomials (fmpq_poly,
inverses by the extended Euclidean algorithm), and orbit representatives chosen by an exact total order (the
lexicographically smallest coefficient vector among the 12 images), with no floating point anywhere. It rebuilds every
point from the recipes in the certificate file, forms the five arguments of each relation, checks non-degeneracy,
reduces every symbol, reduces the target the same way, and compares.
Usage: python -I verify_cert_flint.py cert.json
"""
import sys, json
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq, fmpz_poly

data = json.load(open(sys.argv[1]))
M = data['M']
PHI = fmpq_poly(fmpz_poly.cyclotomic(M))
DEG = PHI.degree()
X = fmpq_poly([0, 1])
ONE = fmpq_poly([1])


def red(p):
    return p % PHI


def zeta(k):
    k %= M
    return red(X ** k)


def inv(p):
    g, s, t = p.xgcd(PHI)          # s p + t PHI = g (a nonzero constant)
    assert g.degree() == 0
    return red(s * fmpq_poly([1 / g[0]]))


def conj(p):
    # x -> x^-1 = x^39
    r = fmpq_poly([0])
    for i in range(p.degree() + 1):
        c = p[i]
        if c != 0:
            r += fmpq_poly([c]) * zeta(-i)
    return red(r)


def key(p):
    return tuple(Fr(int(p[i].p), int(p[i].q)) for i in range(DEG))


if M % 8 == 0:
    SQRT2 = red(zeta(M // 8) + zeta(-M // 8))
    assert red(SQRT2 * SQRT2) == fmpq_poly([2])
    ISQ2 = inv(SQRT2)
else:
    ISQ2 = None

recipes = data['recipes']
memo = {}


def build(i):
    if i in memo:
        return memo[i]
    r = recipes[i]
    if r[0] == 'bbp':
        x = zeta(r[1])
        for _ in range(r[2]):
            x = red(x * ISQ2)
    elif r[0] == 'sym':
        w = build(r[1]); si = r[2]
        if si == 0:
            x = w
        elif si == 1:
            x = red(ONE - inv(w))
        elif si == 2:
            x = inv(red(ONE - w))
        elif si == 3:
            x = inv(w)
        elif si == 4:
            x = red(ONE - w)
        else:
            x = red(w * inv(red(w - ONE)))
        if r[3]:
            x = conj(x)
    elif r[0] == 'pow':
        x = fmpq_poly([0])
        for j, c in r[1].items():
            cf = Fr(c)
            x += fmpq_poly([fmpq(cf.numerator, cf.denominator)]) * zeta(int(j))
        x = red(x)
    elif r[0] == 'rel':
        a, b, j = build(r[1]), build(r[2]), r[3]
        dd = red(ONE - a * b)
        x = {2: red((ONE - a) * inv(dd)), 3: dd, 4: red((ONE - b) * inv(dd))}[j]
    else:
        raise ValueError(r)
    memo[i] = x
    return x


def canon(w):
    if w == conj(w):
        return None, 0
    iw = inv(w); omw = red(ONE - w)
    imgs = [(w, 1), (red(ONE - iw), 1), (inv(omw), 1), (iw, -1), (omw, -1), (red(w * inv(red(w - ONE))), -1)]
    best = None
    for v, sg in imgs:
        for vv, ss in ((v, sg), (conj(v), -sg)):
            k = key(vv)
            if best is None or k < best[0]:
                best = (k, ss)
    return best


total = {}
for item in data['cert']:
    c = Fr(item['coef'])
    x, y = build(item['x']), build(item['y'])
    xy = red(x * y)
    for bad in (x, y, red(ONE - x), red(ONE - y), red(ONE - xy)):
        assert bad != 0, "degenerate relation"
    dd = red(ONE - xy)
    idd = inv(dd)
    for a in (x, y, red((ONE - x) * idd), dd, red((ONE - y) * idd)):
        k, sg = canon(a)
        if k is not None:
            total[k] = total.get(k, Fr(0)) + c * sg
claim = {}
for item in data['target']:
    if item[0] == 'pow':
        w = fmpq_poly([0])
        for j, cc in item[1].items():
            cf = Fr(cc)
            w += fmpq_poly([fmpq(cf.numerator, cf.denominator)]) * zeta(int(j))
        w = red(w)
        kk, sg = canon(w)
        if kk is not None:
            claim[kk] = claim.get(kk, Fr(0)) + Fr(item[2]) * sg
        continue
    a, k, c = item
    w = zeta(a)
    for _ in range(k):
        w = red(w * ISQ2)
    kk, sg = canon(w)
    if kk is not None:
        claim[kk] = claim.get(kk, Fr(0)) + Fr(c) * sg
keys = set(total) | set(claim)
bad = [k for k in keys if total.get(k, 0) != claim.get(k, 0)]
print(f"Q(zeta_{M}) via FLINT polynomials: {len(data['cert'])} relations, {len(keys)} classes")
print("CERTIFICATE VERIFIED EXACTLY (second verifier)" if not bad else f"MISMATCH in {len(bad)} classes")
