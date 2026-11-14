"""Exact recomputation of the relative norms nu used by survey2.py (no floating point).

For each orbit and discriminant listed in the survey JSON, compute in exact arithmetic in Q(zeta_M)
  nu = b^E * prod_{s in G, chi_d(s) = +1, one s per orbit point} (1 - s(z)),   z = zeta_M^j b^(-k/2),
and check that nu = (x + y g)/2 with the survey's integers x, y, where g = sum_t chi_d(t) zeta_|d|^t is the Gauss sum
(g^2 = d, g = i sqrt|d|).  Also re-derives the obstruction vector omega from x, y with an independent routine
(sympy factorization).  Usage: python -I check_nu_exact.py survey.json
"""
import sys, json
from math import gcd
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from cyclo import CycloField
from sympy import factorint, legendre_symbol

data = json.load(open(sys.argv[1]))
b = data['b']
DISC_B = {2: 8, 3: 12, 5: 5}[b]


def kron(d, n):
    res = 1
    for p, e in factorint(n).items():
        if p == 2:
            v = 0 if d % 2 == 0 else (1 if d % 8 in (1, 7) else -1)
        else:
            v = 0 if d % p == 0 else legendre_symbol(d % p, p)
        res *= v ** e
    return res


fields = {}
def field(M):
    if M not in fields:
        K = CycloField(M)
        T = K.zeta
        cond = {2: 8, 3: 12, 5: 5}[b]
        sb = None
        if M % cond == 0:
            if b == 2:
                sb = K.add(T(M // 8), T(-M // 8))                     # sqrt 2 = zeta_8 + zeta_8^-1
            elif b == 3:
                sb = K.add(T(M // 12), T(-M // 12))                   # sqrt 3 = zeta_12 + zeta_12^-1
            else:
                sb = K.add(K.one(), K.scal(2, K.add(T(M // 5), T(-M // 5))))   # sqrt 5 = 1 + 2(zeta_5 + zeta_5^-1)
            assert K.is_zero(K.sub(K.mul(sb, sb), K.elt([b])))
        fields[M] = (K, sb)
    return fields[M]


def act(t, j, k, M):
    e = kron(DISC_B, t) if k % 2 else 1
    return (j * t + (M // 2 if e == -1 else 0)) % M


def root_mod(d, p, E):
    """s with s^2 = d mod p^E: Hensel lifting one power at a time (p odd: Newton; p = 2: bit by bit, s = 1 mod 4)."""
    if p == 2:
        s = 1
        for e in range(3, E):
            if (s * s - d) % 2 ** (e + 1):
                s += 2 ** (e - 1)
        assert (s * s - d) % 2 ** E == 0
        return s
    s = next(r for r in range(p) if (r * r - d) % p == 0) if p < 10 ** 6 else None
    if s is None:
        from sympy import sqrt_mod
        s = min(sqrt_mod(d % p, p, all_roots=True))
    for e in range(1, E):
        pe1 = p ** (e + 1)
        s = (s - (s * s - d) * pow(2 * s, -1, pe1)) % pe1
    assert (s * s - d) % p ** E == 0
    return s


def omega(d, x, y):
    n = (x * x - d * y * y) // 4
    out = {}
    for p, e in factorint(n).items():
        if kron(d, p) != 1:
            continue
        E = e + 6
        s = root_mod(d, p, E)
        mod = p ** (E - 1)
        sh = 1 if p == 2 else 0
        def v(u):
            u %= mod
            c = 0
            while u % p == 0 and c < E:
                u //= p; c += 1
            return c - sh
        vp_, vm_ = v(x + y * s), v(x - y * s)
        assert vp_ + vm_ == e, (p, e, vp_, vm_)
        if vp_ != vm_:
            out[str(p)] = vp_ - vm_
    return out


nchk = 0
for d, s in data['summary'].items():
    d = int(d)
    for ob in s['obstructions']:
        M, j, k = ob['M'], ob['j'], ob['k']
        K, sb = field(M)
        G = [t for t in range(1, M) if gcd(t, M) == 1]
        pts = {}
        for t in G:
            pts.setdefault(act(t, j, k, M), t)
        inv_sb = K.inv(sb) if k % 2 else None
        prod = K.one(); h = 0
        for jj, t in pts.items():
            if kron(d, t) != 1:
                continue
            h += 1
            z = K.zeta(jj)
            if k % 2:
                z = K.mul(z, inv_sb)
            z = K.scal(Fr(1, b ** (k // 2)), z)
            prod = K.mul(prod, K.sub(K.one(), z))
        E = (k * h + 1) // 2
        nu = K.scal(Fr(b) ** E, prod)
        g = K.elt([0])
        for t in range(1, -d):
            c = kron(d, t)
            if c:
                g = K.add(g, K.scal(c, K.zeta((M // -d) * t)))
        rhs = K.scal(Fr(1, 2), K.add(K.elt([ob['x']]), K.scal(ob['y'], g)))
        assert K.is_zero(K.sub(nu, rhs)), (d, M, j, k)
        assert omega(d, ob['x'], ob['y']) == ob['omega'], (d, M, j, k, omega(d, ob['x'], ob['y']), ob['omega'])
        nchk += 1
print(f"base {b}: {nchk} relative norms nu recomputed exactly in Q(zeta_M) and equal to (x + y sqrt d)/2; "
      f"obstruction vectors re-derived: all equal")
