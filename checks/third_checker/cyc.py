# Minimal exact arithmetic in Q(zeta_N), written independently by the referee.
from fractions import Fraction as F


def cyclotomic_coeffs(N):
    # Phi_N via prod_{d|N} (x^d - 1)^{mu(N/d)} using integer polynomial mult/div
    def mobius(n):
        res, p, m = 1, 2, n
        while p * p <= m:
            if m % p == 0:
                m //= p
                if m % p == 0:
                    return 0
                res = -res
            p += 1
        if m > 1:
            res = -res
        return res

    def pmul(a, b):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                r[i + j] += x * y
        return r

    def pdiv(a, b):  # exact division, b monic-ish (leading coeff 1)
        a = a[:]
        q = [0] * (len(a) - len(b) + 1)
        for i in range(len(a) - len(b), -1, -1):
            c = a[i + len(b) - 1] // b[-1]
            q[i] = c
            for j, y in enumerate(b):
                a[i + j] -= c * y
        assert all(x == 0 for x in a), "non-exact division"
        return q

    num, den = [1], [1]
    for d in range(1, N + 1):
        if N % d == 0:
            mu = mobius(N // d)
            poly = [-1] + [0] * (d - 1) + [1]
            if mu == 1:
                num = pmul(num, poly)
            elif mu == -1:
                den = pmul(den, poly)
    return pdiv(num, den)


class Field:
    def __init__(self, N):
        self.N = N
        self.phi = cyclotomic_coeffs(N)  # low -> high, monic
        self.d = len(self.phi) - 1
        self.pow_cache = {}

    def reduce(self, poly):
        poly = list(poly)
        d = self.d
        for m in range(len(poly) - 1, d - 1, -1):
            c = poly[m]
            if c:
                for j in range(d + 1):
                    poly[m - d + j] -= c * self.phi[j]
        poly = poly[:d] + [F(0)] * max(0, d - len(poly))
        return tuple(F(x) for x in poly[:d])

    def zpow(self, k):
        k %= self.N
        if k not in self.pow_cache:
            p = [F(0)] * (k + 1)
            p[k] = F(1)
            self.pow_cache[k] = self.reduce(p)
        return Elem(self, self.pow_cache[k])

    def const(self, r):
        v = [F(0)] * self.d
        v[0] = F(r)
        return Elem(self, tuple(v))


class Elem:
    __slots__ = ("K", "v")

    def __init__(self, K, v):
        self.K = K
        self.v = v

    def _coerce(self, o):
        if isinstance(o, Elem):
            return o
        return self.K.const(o)

    def __add__(self, o):
        o = self._coerce(o)
        return Elem(self.K, tuple(a + b for a, b in zip(self.v, o.v)))

    __radd__ = __add__

    def __neg__(self):
        return Elem(self.K, tuple(-a for a in self.v))

    def __sub__(self, o):
        return self + (-self._coerce(o))

    def __rsub__(self, o):
        return self._coerce(o) - self

    def __mul__(self, o):
        o = self._coerce(o)
        d = self.K.d
        r = [F(0)] * (2 * d - 1)
        for i, a in enumerate(self.v):
            if a:
                for j, b in enumerate(o.v):
                    if b:
                        r[i + j] += a * b
        return Elem(self.K, self.K.reduce(r))

    __rmul__ = __mul__

    def is_zero(self):
        return all(a == 0 for a in self.v)

    def inv(self):
        K = self.K
        d = K.d
        assert not self.is_zero(), "inverse of zero"
        # columns: self * zeta^j
        cols = [(self * K.zpow(j)).v for j in range(d)]
        M = [[cols[j][i] for j in range(d)] + [F(1) if i == 0 else F(0)] for i in range(d)]
        # Gaussian elimination
        for col in range(d):
            piv = next(r for r in range(col, d) if M[r][col] != 0)
            M[col], M[piv] = M[piv], M[col]
            pv = M[col][col]
            M[col] = [x / pv for x in M[col]]
            for r in range(d):
                if r != col and M[r][col] != 0:
                    f = M[r][col]
                    M[r] = [x - f * y for x, y in zip(M[r], M[col])]
        return Elem(K, tuple(M[i][d] for i in range(d)))

    def __truediv__(self, o):
        return self * self._coerce(o).inv()

    def __rtruediv__(self, o):
        return self._coerce(o) * self.inv()

    def conj(self):
        K = self.K
        r = K.const(0)
        for k, a in enumerate(self.v):
            if a:
                r = r + K.zpow(-k) * a
        return r

    def __eq__(self, o):
        return isinstance(o, Elem) and self.v == o.v

    def __hash__(self):
        return hash(self.v)

    def key(self):
        return self.v

    def to_complex(self):
        import cmath
        N = self.K.N
        return sum(float(a) * cmath.exp(2j * cmath.pi * k / N) for k, a in enumerate(self.v))


def orbit(z):
    """The 12 images of z under <z->1-z, z->1/z, z->conj z>, with the sign eps."""
    one = 1
    iz = z.inv()
    omz = one - z
    iomz = omz.inv()
    els = [(z, 1), (one - iz, 1), (iomz, 1), (iz, -1), (omz, -1), (-(z * iomz), -1)]
    els += [(e.conj(), -s) for e, s in els]
    return els


def reduce_symbol(z, cache):
    """Return (rep, sign) with D(z) = sign*D(rep), or None if D(z)=0 is forced."""
    if z in cache:
        return cache[z]
    if z == z.conj():
        cache[z] = None
        return None
    els = orbit(z)
    signs = {}
    for e, s in els:
        signs.setdefault(e, set()).add(s)
    if any(len(s) > 1 for s in signs.values()):
        cache[z] = None
        return None
    rep = min(signs.keys(), key=lambda e: e.key())
    s = next(iter(signs[rep]))
    cache[z] = (rep, s)
    return cache[z]


def five_args(x, y):
    one = 1
    xy = x * y
    for t, name in ((x, "x"), (y, "y"), (one - x, "1-x"), (one - y, "1-y"), (one - xy, "1-xy")):
        assert not t.is_zero(), f"degenerate relation: {name}=0"
    inv = (one - xy).inv()
    return [x, y, (one - x) * inv, one - xy, (one - y) * inv]


def reduce_combo(terms, cache):
    """terms: list of (coeff, element). Return dict rep -> coeff (nonzero only)."""
    out = {}
    for c, z in terms:
        r = reduce_symbol(z, cache)
        if r is None:
            continue
        rep, s = r
        out[rep] = out.get(rep, F(0)) + F(c) * s
    return {k: v for k, v in out.items() if v != 0}


def check_certificate(rels, target, label):
    cache = {}
    terms = []
    for c, x, y in rels:
        for a in five_args(x, y):
            terms.append((c, a))
    got = reduce_combo(terms, cache)
    want = reduce_combo(target, cache)
    keys = set(got) | set(want)
    diff = {k: got.get(k, F(0)) - want.get(k, F(0)) for k in keys}
    diff = {k: v for k, v in diff.items() if v != 0}
    print(f"[{label}] relations={len(rels)} reduced symbols={len(got)} target symbols={len(want)} ->",
          "EXACT MATCH" if not diff else f"MISMATCH in {len(diff)} symbols")
    return got, want, diff
