# Independent exact verifier for five-term-relation certificates in Q(zeta_8) = Q(i, sqrt2), written for
# papers/bbp-formulas-2026 (12 Nov 2026). It shares no code with the lane's scripts: elements are tuples
# (a, b, c, d) of Fractions meaning (a + b*sqrt2) + i*(c + d*sqrt2); orbit representatives are chosen by an exact total
# order, and every equality is decided exactly. Floating point is used only to look up which exact element a stored
# certificate entry denotes (the lane stores points by rounded numerical keys; an assertion checks the match).
# Input: the unit list (exponent vectors over zeta8, 1+sqrt2, 1-zeta8, 1+sqrt(-2), 1-sqrt(-2), with a float key) and
# a certificate (coef, float keys of x and y), as stored by the lane.
# Usage: python3 -I verify_q8_independent.py units.json cert.json claim   (claim = L8 or BB27)
import sys, json
from fractions import Fraction as F
H = F(1, 2)
def mul(x, y):
    a, b, c, d = x; e, f, g, h = y
    rr = (a * e + 2 * b * f - (c * g + 2 * d * h), a * f + b * e - (c * h + d * g))
    ii = (a * g + 2 * b * h + c * e + 2 * d * f, a * h + b * g + c * f + d * e)
    return (rr[0], rr[1], ii[0], ii[1])
def conj(x): return (x[0], x[1], -x[2], -x[3])
def inv(x):
    a, b, c, d = x
    p = a * a + 2 * b * b + c * c + 2 * d * d          # |x|^2 = p + q sqrt2
    q = 2 * a * b + 2 * c * d
    den = p * p - 2 * q * q
    ip, iq = p / den, -q / den                         # 1/(p + q sqrt2)
    xc = conj(x)
    return mul(xc, (ip, iq, F(0), F(0)))
def sub(x, y): return tuple(s - t for s, t in zip(x, y))
ONE = (F(1), F(0), F(0), F(0)); ZERO = (F(0),) * 4
def div(x, y): return mul(x, inv(y))
def power(x, e):
    r = ONE; base = x if e >= 0 else inv(x)
    for _ in range(abs(e)): r = mul(r, base)
    return r
Z8 = (F(0), H, F(0), H)
GENS = [(F(1), F(1), F(0), F(0)), (F(1), -H, F(0), -H), (F(1), F(0), F(0), F(1)), (F(1), F(0), F(0), F(-1))]
def from_exps(e):
    x = power(Z8, e[0])
    for g, p in zip(GENS, e[1:]): x = mul(x, power(g, p))
    return x
def num(x):
    s2 = 2 ** 0.5
    return complex(float(x[0]) + float(x[1]) * s2, float(x[2]) + float(x[3]) * s2)
def images(w):
    iw = inv(w); omw = sub(ONE, w)
    base = [(w, 1), (sub(ONE, iw), 1), (inv(omw), 1), (iw, -1), (omw, -1), (div(w, sub(w, ONE)), -1)]
    out = []
    for v, s in base: out += [(v, s), (conj(v), -s)]
    return out
def canon(w):
    if w[2] == 0 and w[3] == 0: return None, 0       # real: D = 0
    return min(images(w), key=lambda t: t[0])           # exact lexicographic order on the Fraction tuple
units = json.load(open(sys.argv[1])); cert = json.load(open(sys.argv[2])); claim_name = sys.argv[3]
bykey = {(round(re, 9) + 0.0, round(im, 9) + 0.0): e for e, (re, im) in units}
total = {}
def add(w, c):
    r, s = canon(w)
    if r is not None: total[r] = total.get(r, F(0)) + c * s
for item in cert:
    c = F(item['coef'])
    x = from_exps(bykey[(round(item['x'][0], 9) + 0.0, round(item['x'][1], 9) + 0.0)])
    y = from_exps(bykey[(round(item['y'][0], 9) + 0.0, round(item['y'][1], 9) + 0.0)])
    assert abs(num(x) - complex(*item['x'])) < 1e-9 and abs(num(y) - complex(*item['y'])) < 1e-9
    xy = mul(x, y)
    for t in (x, y, sub(ONE, x), sub(ONE, y), sub(ONE, xy)):
        assert t != ZERO, 'degenerate relation'
    d = sub(ONE, xy)
    for a in (x, y, div(sub(ONE, x), d), d, div(sub(ONE, y), d)): add(a, c)
if claim_name == 'L8':     # 6[i/sqrt2] + [i/(2 sqrt2)] - 4[zeta8] - 4[zeta8^3]
    claim_pts = [((F(0), F(0), F(0), H), F(6)), ((F(0), F(0), F(0), F(1, 4)), F(1)), (Z8, F(-4)), (power(Z8, 3), F(-4))]
else:                      # (27/2)[w] - (9/2)[w^2] + (1/2)[w^3] - 8[zeta8] - 8[zeta8^3], w = (1+sqrt-2)/(1-sqrt-2)
    w = div(GENS[2], GENS[3])
    claim_pts = [(w, F(27, 2)), (power(w, 2), F(-9, 2)), (power(w, 3), F(1, 2)), (Z8, F(-8)), (power(Z8, 3), F(-8))]
claim = {}
for p, c in claim_pts:
    r, s = canon(p)
    if r is not None: claim[r] = claim.get(r, F(0)) + c * s
keys = set(total) | set(claim)
bad = {k: total.get(k, F(0)) - claim.get(k, F(0)) for k in keys}
bad = {k: v for k, v in bad.items() if v != 0}
print(f'{claim_name}: relations {len(cert)}, orbit classes {len(keys)}:', 'VERIFIED EXACTLY' if not bad else f'MISMATCH in {len(bad)} classes')
