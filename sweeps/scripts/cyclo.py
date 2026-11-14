"""Exact arithmetic in the cyclotomic field Q(zeta_N) (elements = rational coefficient vectors in the
power basis 1, t, ..., t^(phi(N)-1) with t = exp(2 pi i / N)), plus numeric evaluation.
Small and slow, but exact: used to certify that the arguments of functional equations are what we claim."""
from fractions import Fraction as Fr
import cmath, math


def cyclotomic_poly(N):
    """Coefficients (low -> high) of Phi_N as integers."""
    # Phi_N = prod_{d | N} (x^d - 1)^{mu(N/d)}
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
            if x:
                for j, y in enumerate(b):
                    r[i + j] += x * y
        return r
    def pdiv(a, b):  # exact division of integer polys
        a = a[:]
        q = [0] * (len(a) - len(b) + 1)
        for i in range(len(q) - 1, -1, -1):
            c = a[i + len(b) - 1] // b[-1]
            q[i] = c
            for j, y in enumerate(b):
                a[i + j] -= c * y
        return q
    num, den = [1], [1]
    for d in range(1, N + 1):
        if N % d == 0:
            mu = mobius(N // d)
            f = [-1] + [0] * (d - 1) + [1]
            if mu == 1:
                num = pmul(num, f)
            elif mu == -1:
                den = pmul(den, f)
    return pdiv(num, den)


class CycloField:
    def __init__(self, N):
        self.N = N
        self.phi = cyclotomic_poly(N)
        self.deg = len(self.phi) - 1
        self.t = cmath.exp(2j * math.pi / N)

    def reduce(self, coeffs):
        c = [Fr(x) for x in coeffs]
        d = self.deg
        for i in range(len(c) - 1, d - 1, -1):
            if c[i] != 0:
                f = c[i]
                for j in range(d + 1):
                    c[i - d + j] -= f * self.phi[j]
        c = c[:d] + [Fr(0)] * max(0, d - len(c))
        return tuple(c[:d])

    def elt(self, coeffs):
        return self.reduce(list(coeffs))

    def one(self):
        return self.elt([1])

    def zeta(self, k):
        k %= self.N
        return self.elt([0] * k + [1])

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple(x - y for x, y in zip(a, b))

    def neg(self, a):
        return tuple(-x for x in a)

    def mul(self, a, b):
        r = [Fr(0)] * (2 * self.deg - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y:
                        r[i + j] += x * y
        return self.reduce(r)

    def scal(self, q, a):
        return tuple(Fr(q) * x for x in a)

    def inv(self, a):
        # solve M x = e1 where M = multiplication-by-a matrix
        d = self.deg
        cols = []
        e = [Fr(0)] * d
        for j in range(d):
            ej = [Fr(0)] * d
            ej[j] = Fr(1)
            cols.append(self.mul(a, tuple(ej)))
        M = [[cols[j][i] for j in range(d)] + [Fr(1) if i == 0 else Fr(0)] for i in range(d)]
        # Gaussian elimination
        for c in range(d):
            p = next(r for r in range(c, d) if M[r][c] != 0)
            M[c], M[p] = M[p], M[c]
            pv = M[c][c]
            M[c] = [x / pv for x in M[c]]
            for r in range(d):
                if r != c and M[r][c] != 0:
                    f = M[r][c]
                    M[r] = [x - f * y for x, y in zip(M[r], M[c])]
        return tuple(M[i][d] for i in range(d))

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    def num(self, a):
        return sum(complex(float(x)) * self.t ** i for i, x in enumerate(a))

    def is_zero(self, a):
        return all(x == 0 for x in a)
