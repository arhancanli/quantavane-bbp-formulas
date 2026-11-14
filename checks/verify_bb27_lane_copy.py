"""Exact verification of a five-term-relation certificate in Q(zeta_8) for the Borwein-Broadhurst identity (27) (copy of the L8 lane verify_cert.py with the claim replaced; papers/bbp-formulas-2026, 12 Nov 2026).

Claim:  6 D(sqrt(-2)/2) + D(sqrt(-2)/4) - 4 D(zeta_8) - 4 D(zeta_8^3) = 0,
where D is the Bloch-Wigner dilogarithm at the identity embedding (zeta_8 = exp(2 pi i/8)).

Certificate: rational coefficients c_i and pairs (x_i, y_i) of elements of Q(zeta_8).  For each pair the
five-term relation  D(x)+D(y)+D((1-x)/(1-xy))+D(1-xy)+D((1-y)/(1-xy)) = 0  holds (Bloch-Wigner).  Each
argument w is rewritten with the exact symmetries D(w) = D(1-1/w) = D(1/(1-w)) = -D(1/w) = -D(1-w)
= -D(w/(w-1)), D(conj w) = -D(w), and D(w) = 0 for real w.  We check in exact arithmetic that
sum_i c_i * (relation_i) equals the claim, class by class.
Usage: python -I verify_cert.py exc.json cert.json
"""
import sys, json, cmath, math
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from fractions import Fraction as Fr
from cyclo import CycloField

K = CycloField(8)
T = K.zeta
SQ2 = K.sub(T(1), T(3))
SQM2 = K.add(T(1), T(3))
GEN = {'eps': K.add(K.one(), SQ2), 'p2': K.sub(K.one(), T(1)), 'p3': K.add(K.one(), SQM2), 'q3': K.sub(K.one(), SQM2)}
ONE = K.one()


def power(ex, e):
    res = K.one()
    base = ex if e >= 0 else K.inv(ex)
    for _ in range(abs(e)):
        res = K.mul(res, base)
    return res


def elt_from_exps(e):
    k, a, b, c, d = e
    x = T(k)
    for g, p in zip(('eps', 'p2', 'p3', 'q3'), (a, b, c, d)):
        x = K.mul(x, power(GEN[g], p))
    return x


def conj(x):
    # complex conjugation = automorphism zeta -> zeta^-1 ; apply to power basis
    res = K.elt([0])
    for i, c in enumerate(x):
        if c:
            res = K.add(res, K.scal(c, T(-i)))
    return res


def is_real(x):
    return K.is_zero(K.sub(x, conj(x)))


def numkey(x):
    v = K.num(x)
    return (round(v.real, 9) + 0.0, round(v.imag, 9) + 0.0)


def images(w):
    """(image, sign) for the 12 symmetries."""
    one = ONE
    out = []
    inv = K.inv(w)
    omw = K.sub(one, w)
    base = [(w, 1), (K.sub(one, inv), 1), (K.inv(omw), 1), (inv, -1), (omw, -1), (K.div(w, K.sub(w, one)), -1)]
    for v, s in base:
        out.append((v, s))
        out.append((conj(v), -s))
    return out


def canon(w):
    if is_real(w):
        return None, 0
    best = None
    for v, s in images(w):
        k = numkey(v)
        if best is None or k < best[0]:
            best = (k, v, s)
    return best[1], best[2]


def main():
    exc = json.load(open(sys.argv[1]))
    cert = json.load(open(sys.argv[2]))
    bykey = {}
    for e, (re, im) in exc:
        bykey[(round(re, 9) + 0.0, round(im, 9) + 0.0)] = e
    total = {}
    def addterm(w, coef):
        cv, s = canon(w)
        if cv is None:
            return
        total[cv] = total.get(cv, Fr(0)) + coef * s
    for item in cert:
        c = Fr(item['coef'])
        kx = (round(item['x'][0], 9) + 0.0, round(item['x'][1], 9) + 0.0)
        ky = (round(item['y'][0], 9) + 0.0, round(item['y'][1], 9) + 0.0)
        x = elt_from_exps(bykey[kx])
        y = elt_from_exps(bykey[ky])
        # sanity: exact element evaluates to the recorded number
        assert abs(K.num(x) - complex(*item['x'])) < 1e-9 and abs(K.num(y) - complex(*item['y'])) < 1e-9
        xy = K.mul(x, y)
        omxy = K.sub(ONE, xy)
        args = [x, y, K.div(K.sub(ONE, x), omxy), omxy, K.div(K.sub(ONE, y), omxy)]
        for a in args:
            addterm(a, c)
    claim = {}
    # claim (27) of Borwein-Broadhurst, moved to one side: (27/2)[w] - (9/2)[w^2] + (1/2)[w^3] - 8[z8] - 8[z8^3],
    # w = (1 + sqrt(-2))/(1 - sqrt(-2)) = exp(2i arctan sqrt2)
    w1 = K.div(GEN['p3'], GEN['q3'])
    w2 = K.mul(w1, w1)
    w3 = K.mul(w2, w1)
    for w, coef in ((w1, Fr(27, 2)), (w2, Fr(-9, 2)), (w3, Fr(1, 2)), (T(1), Fr(-8)), (T(3), Fr(-8))):
        cv, s = canon(w)
        claim[cv] = claim.get(cv, Fr(0)) + coef * s
    keys = set(total) | set(claim)
    diff = {k: total.get(k, Fr(0)) - claim.get(k, Fr(0)) for k in keys}
    nonzero = {k: v for k, v in diff.items() if v != 0}
    print("classes involved:", len(keys))
    if not nonzero:
        print("CERTIFICATE VERIFIED EXACTLY: sum_i c_i (five-term relation_i) = claim")
    else:
        print("MISMATCH in", len(nonzero), "classes")
        for k, v in list(nonzero.items())[:10]:
            print("  ", K.num(k), v)


if __name__ == '__main__':
    main()
