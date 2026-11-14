"""Complete search with the twisted Bloch criterion, for base-b points z = zeta_M^j * b^(-k/2).

For every Galois orbit O of such points (z in Q(zeta_M) with M <= MFMAX, 1 <= k <= KMAX) and every fundamental
discriminant d < 0 with sqrt(d) in Q(z):
  xi_{O,d} = sum_{s in Gal(Q(z)/Q)} chi_d(s) [s z]          (chi_d-twisted orbit sum)
  nu_O     = N_{Q(z)/Q(sqrt d)}(1 - z) times a power of b (an algebraic integer of Q(sqrt d), computed exactly
             by rounding a 45-digit evaluation; recomputed in exact arithmetic by check_nu_exact.py)
  omega_O  = (v_P(nu) - v_Pbar(nu))_{p split in Q(sqrt d)}    (obstruction vector)
Criterion (Theorem 3.3 of the paper): sum_O B_O xi_{O,d} has a five-term certificate iff sum_O B_O k_O omega_O = 0.
Output per d: the kernel K_d of B -> sum B_O k_O omega_O, the subspace Dist_d spanned by the distribution relations
D(z^n) = n sum_{w^n=1} D(w z) among the points, and for a basis of a complement of Dist_d in K_d the numbers
  q = D(xi)/(sqrt|d| L(2,chi_d)),   r = (1/(2 pi)) sum_O B_O k_O A_O,  A_O = sum chi_d(s) arg(1 - s z),
so that  sum_O B_O sum_s chi_d(s) Im Li_2(s z) = q sqrt|d| L(2,chi_d) + r pi log b.
Usage: python -I survey2.py b KMAX MFMAX [json_out]
"""
import sys, math, time, json
from math import gcd
from fractions import Fraction as Fr
import mpmath as mp
from flint import fmpz, fmpq_mat, fmpq

b = int(sys.argv[1]); KMAX = int(sys.argv[2]); MFMAX = int(sys.argv[3])
JOUT = sys.argv[4] if len(sys.argv) > 4 else None
FORMAL_OFF = len(sys.argv) > 5 and sys.argv[5] == 'noformal'
mp.mp.dps = 45
DISC_B = {2: 8, 3: 12, 5: 5, 6: 24, 7: 28}[b]
COND_B = abs(DISC_B)


def factor(n):
    return {int(p): int(e) for p, e in fmpz(n).factor()} if n > 1 else {}


_kc = {}
def kron(d, n):
    """Kronecker symbol (d/n) for n > 0."""
    if (d, n) in _kc:
        return _kc[(d, n)]
    res = 1
    for p, e in factor(n).items():
        if p == 2:
            v = 0 if d % 2 == 0 else (1 if d % 8 in (1, 7) else -1)
        else:
            r = d % p
            v = 0 if r == 0 else (1 if pow(r, (p - 1) // 2, p) == 1 else -1)
        res *= v ** e
    _kc[(d, n)] = res
    return res


def is_fund_neg(d):
    if d >= 0:
        return False
    if d % 4 == 1:
        return all(e == 1 for e in factor(-d).values())
    if d % 4 == 0:
        m = d // 4
        return m % 4 in (2, 3) and all(e == 1 for e in factor(-m).values())
    return False


def units(M):
    return [t for t in range(1, M) if gcd(t, M) == 1]


def act(t, j, k, M):
    e = kron(DISC_B, t) if k % 2 else 1
    return (j * t + (M // 2 if e == -1 else 0)) % M


def point(j, k, M):
    return mp.expjpi(mp.mpf(2 * j) / M) / mp.sqrt(b) ** k


_cache = {}
def DA(j, k, M):
    key = (k, Fr(j, M))
    if key not in _cache:
        z = point(j, k, M)
        _cache[key] = (mp.im(mp.polylog(2, z)) + mp.arg(1 - z) * mp.log(abs(z)), mp.arg(1 - z))
    return _cache[key]


def Lval(d):
    n = -d
    return mp.fsum(kron(d, a) * mp.zeta(2, mp.mpf(a) / n) for a in range(1, n + 1)) / n ** 2


def vp(n, p):
    if n == 0:
        return 10 ** 9
    v = 0
    while n % p == 0:
        n //= p; v += 1
    return v


def padic_sqrt(d, p, E):
    if p == 2:
        s = 1
        for e in range(4, E + 1):
            if (s * s - d) % (2 ** e):
                s += 2 ** (e - 2)
        assert (s * s - d) % (2 ** E) == 0
        return s
    from sympy.ntheory import sqrt_mod
    r = min(sqrt_mod(d % p, p, all_roots=True))
    s, pe = r, p
    while pe < p ** E:
        pe = min(pe * pe, p ** E)
        s = (s - (s * s - d) * pow(2 * s, -1, pe)) % pe
    assert (s * s - d) % (p ** E) == 0
    return s


def omega(d, x, y):
    """nu = (x + y sqrt d)/2.  Returns ({p: v_P - v_Pbar} over split p with nonzero difference, N(nu))."""
    n = (x * x - d * y * y)
    assert n % 4 == 0 and n > 0
    n //= 4
    out = {}
    for p, e in factor(n).items():
        if kron(d, p) != 1:
            continue
        E = e + 6
        s = padic_sqrt(d, p, E)
        mod = p ** E
        sh = 1 if p == 2 else 0
        vplus = vp((x + y * s) % mod, p) - sh
        vminus = vp((x - y * s) % mod, p) - sh
        assert vplus + vminus == e, (d, p, e, vplus, vminus)
        if vplus != vminus:
            out[p] = vplus - vminus
    return out, n


def to_fmpq(c):
    return fmpq(c.numerator, c.denominator)


def from_fmpq(x):
    return Fr(int(x.p), int(x.q))


def nullspace(rows, ncols):
    """basis of {v : rows v = 0} (exact, FLINT rref)."""
    if not rows:
        return [[Fr(int(i == j)) for i in range(ncols)] for j in range(ncols)]
    A = fmpq_mat(len(rows), ncols, [to_fmpq(c) for r in rows for c in r])
    R, rk = A.rref()
    Rl = [[from_fmpq(R[i, j]) for j in range(ncols)] for i in range(rk)]
    piv = [next(j for j in range(ncols) if Rl[i][j] != 0) for i in range(rk)]
    basis = []
    for f in (c for c in range(ncols) if c not in piv):
        v = [Fr(0)] * ncols
        v[f] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -Rl[i][f]
        basis.append(v)
    return basis


def nullspace_slow(rows, ncols):
    A = [r[:] for r in rows]
    piv = []
    ri = 0
    for c in range(ncols):
        pr = next((i for i in range(ri, len(A)) if A[i][c] != 0), None)
        if pr is None:
            continue
        A[ri], A[pr] = A[pr], A[ri]
        pv = A[ri][c]
        A[ri] = [v / pv for v in A[ri]]
        for i in range(len(A)):
            if i != ri and A[i][c] != 0:
                f = A[i][c]
                A[i] = [vi - f * vr for vi, vr in zip(A[i], A[ri])]
        piv.append(c); ri += 1
    basis = []
    for f in (c for c in range(ncols) if c not in piv):
        v = [Fr(0)] * ncols
        v[f] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -A[i][f]
        basis.append(v)
    return basis


def rank(vecs):
    A = [v[:] for v in vecs]
    r = 0
    ncols = len(A[0]) if A else 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        for i in range(r + 1, len(A)):
            if A[i][c] != 0:
                f = A[i][c] / A[r][c]
                A[i] = [vi - f * vr for vi, vr in zip(A[i], A[r])]
        r += 1
    return r


class Echelon:
    """incremental row echelon form over Q; add(v) returns True if v was independent."""
    def __init__(self):
        self.rows = {}          # pivot column -> row (normalized, pivot 1)
    def reduce(self, v):
        v = v[:]
        for c in sorted(self.rows):
            if v[c] != 0:
                f = v[c]
                rw = self.rows[c]
                v = [a - f * bb for a, bb in zip(v, rw)]
        return v
    def add(self, v):
        v = self.reduce(v)
        piv = next((i for i, a in enumerate(v) if a != 0), None)
        if piv is None:
            return False
        pv = v[piv]
        v = [a / pv for a in v]
        # keep reduced: eliminate piv from existing rows
        for c in list(self.rows):
            rw = self.rows[c]
            if rw[piv] != 0:
                f = rw[piv]
                self.rows[c] = [a - f * bb for a, bb in zip(rw, v)]
        self.rows[piv] = v
        return True


def ratid(x, den=10 ** 8):
    if abs(x) < mp.mpf(10) ** -35:
        return Fr(0)
    f = Fr(mp.nstr(x, 42, strip_zeros=False).replace(' ', '')).limit_denominator(den)
    if abs(x - mp.mpf(f.numerator) / f.denominator) < mp.mpf(10) ** -35:
        return f
    return None


# ---------------- enumerate orbits ----------------
t0 = time.time()
seen = {}
orbits = []
for M in range(4, MFMAX + 1, 2):
    G = units(M)
    for k in range(1, KMAX + 1):
        if k % 2 and M % COND_B:
            continue
        for j in range(1, M):
            if 2 * j == M:
                continue                       # real point
            key = (k, Fr(j, M))
            if key in seen:
                continue
            # minimal level: the point and sqrt(b)^k must lie in Q(zeta_M); skip if a smaller even level works
            m = M // gcd(j, M)
            Mmin = m * (COND_B // gcd(m, COND_B) if k % 2 else 1)
            Mmin = Mmin * 2 // gcd(Mmin, 2)
            if Mmin != M:
                continue
            pts = {}
            for t in G:
                pts.setdefault(act(t, j, k, M), t)
            stab = [t for t in G if act(t, j, k, M) == j]
            oi = len(orbits)
            for jj in pts:
                seen[(k, Fr(jj, M))] = oi
            orbits.append(dict(k=k, M=M, j=j, pts=pts, stab=stab))
npts = len(seen)
print(f"base {b}: {len(orbits)} orbits, {npts} non-real points zeta_M^j {b}^(-k/2) with M <= {MFMAX}, "
      f"k <= {KMAX} ({time.time()-t0:.1f}s)", flush=True)


def locate(k, num, den):
    """point zeta^(num/den) b^(-k/2) -> (orbit index, M, j) or None if real / outside the set."""
    f = Fr(num, den) % 1
    if f == 0 or f == Fr(1, 2):
        return None
    oi = seen.get((k, f))
    if oi is None:
        return 'missing'
    o = orbits[oi]
    return oi, (f * o['M']).numerator


# distribution relations, as point-level relations: [z^n] - n sum_{w^n = 1} [w z]
dist_rel = []
nmiss = 0
for (k, f), oi in seen.items():
    for n in range(2, KMAX // k + 1):
        terms = []
        tgt = locate(k * n, f.numerator * n, f.denominator)
        ok = tgt != 'missing'
        if tgt not in (None, 'missing'):
            terms.append((tgt, Fr(1)))
        for w in range(n):
            src = locate(k, f.numerator * n + w * f.denominator, f.denominator * n)
            if src == 'missing':
                ok = False; break
            if src is not None:
                terms.append((src, Fr(-n)))
        if ok and terms:
            dist_rel.append(terms)
        elif not ok:
            nmiss += 1
print(f"distribution relations inside the point set: {len(dist_rel)} (skipped {nmiss} reaching outside)", flush=True)

ds = sorted(set(d for o in orbits for d in range(-o['M'], 0) if o['M'] % (-d) == 0 and is_fund_neg(d)), reverse=True)
summary = {}
for d in ds:
    t1 = time.time()
    rows = []
    for oi, o in enumerate(orbits):
        M, j, k = o['M'], o['j'], o['k']
        if M % (-d):
            continue
        if any(kron(d, t) != 1 for t in o['stab']):
            continue
        coef = {jj: kron(d, t) for jj, t in o['pts'].items()}
        h = sum(1 for c in coef.values() if c == 1)
        E = (k * h + 1) // 2
        prec = 40 + int(h * math.log10(b ** (k / 2.0) + 1) + E * math.log10(b)) + 5
        with mp.workdps(prec):
            nup = mp.fprod(1 - point(jj, k, M) for jj, c in coef.items() if c == 1)
            nui = nup * mp.mpf(b) ** E
            X = 2 * mp.re(nui); Y = 2 * mp.im(nui) / mp.sqrt(-d)
            x = int(mp.nint(X)); y = int(mp.nint(Y))
            assert abs(X - x) < mp.mpf(10) ** -30 and abs(Y - y) < mp.mpf(10) ** -30, (d, o['M'], o['j'], X, Y)
        om, n = omega(d, x, y)
        rows.append(dict(oi=oi, k=k, M=M, j=j, size=len(coef), coef=coef, x=x, y=y, n=n, om=om))
    if not rows:
        continue
    idx = {r['oi']: i for i, r in enumerate(rows)}
    primes = sorted(set(p for r in rows for p in r['om']))
    eqs = [[Fr(r['k'] * r['om'].get(p, 0)) for r in rows] for p in primes]
    K = nullspace(eqs, len(rows))
    # distribution relations projected to the chi_d-isotypic part, in the basis xi_{O,d}
    Dv = []
    for terms in dist_rel:
        v = [Fr(0)] * len(rows)
        for (oi, jj), c in terms:
            if oi not in idx:
                continue        # sqrt d not in Q(z): projection 0
            r = rows[idx[oi]]
            v[idx[oi]] += c * r['coef'][jj] / r['size']
        if any(v):
            Dv.append(v)
    n_ = len(rows)
    # sanity: every distribution relation satisfies the criterion (exact)
    if Dv and eqs:
        Wm = fmpq_mat(len(eqs), n_, [to_fmpq(c) for r in eqs for c in r])
        Dm = fmpq_mat(len(Dv), n_, [to_fmpq(c) for r in Dv for c in r])
        P = Wm * Dm.transpose()
        assert all(P[i, j] == 0 for i in range(P.nrows()) for j in range(P.ncols())), "distribution relation violates criterion"
    if Dv:
        Dm = fmpq_mat(len(Dv), n_, [to_fmpq(c) for r in Dv for c in r])
        RD, rD = Dm.rref()
        basis = [[RD[i, j] for j in range(n_)] for i in range(rD)]
    else:
        rD, basis = 0, []
    comp = []
    cur = rD
    for v in K:
        trial = basis + [[to_fmpq(c) for c in v]]
        if fmpq_mat(len(trial), n_, [c for r in trial for c in r]).rank() > cur:
            basis = trial; comp.append(v); cur += 1
    L = Lval(d); sq = mp.sqrt(-d)
    res = []
    for v in comp:
        Dk = mp.mpf(0); Ak = mp.mpf(0)
        for c, r in zip(v, rows):
            if c:
                cc = mp.mpf(c.numerator) / c.denominator
                Dk += cc * mp.fsum(s * DA(jj, r['k'], r['M'])[0] for jj, s in r['coef'].items())
                Ak += cc * r['k'] * mp.fsum(s * DA(jj, r['k'], r['M'])[1] for jj, s in r['coef'].items())
        q = Dk / (sq * L); rr = Ak / (2 * mp.pi)
        res.append(dict(v=v, q=q, r=rr, qf=ratid(q), rf=ratid(rr)))
    # complement vectors with q = r = 0 numerically: prove q = 0 exactly by showing that the Im Li_2 combination
    # vanishes formally (every coefficient of b^(-kn/2)/n^2 is 0) and r = 0 (r is rational with small denominator)
    nformal = 0; unexplained = []
    for x in res:
        if abs(x['q']) > mp.mpf(10) ** -30 or abs(x['r']) > mp.mpf(10) ** -30:
            continue
        if FORMAL_OFF:
            unexplained.append(x); continue
        ok = True
        # coefficient of b^(-e/2)/e^2 is sum_{k | e} k^2 S_k(e/k), S_k(n) = sum over rows with that k
        sup = [(mp.mpf(c.numerator) / c.denominator, rows[i]) for i, c in enumerate(x['v']) if c]
        per = 1
        for c, rw in sup:
            pk = rw['k'] * rw['M']
            per = per * pk // gcd(per, pk)
        for e in range(1, per + 1):
            tot = mp.mpf(0)
            for c, rw in sup:
                kk = rw['k']
                if e % kk:
                    continue
                n = e // kk
                tot += kk * kk * c * mp.fsum(sgn * mp.sinpi(mp.mpf(2 * jj * n) / rw['M']) for jj, sgn in rw['coef'].items())
            if abs(tot) > mp.mpf(10) ** -25:
                ok = False; break
        if ok:
            nformal += 1
        else:
            unexplained.append(x)
    print(f"   zero complement vectors: {nformal} vanish formally (so q = 0 exactly); unexplained: {len(unexplained)}",
          flush=True)
    # rank of (q, r) on the complement
    qr = [(x['q'], x['r']) for x in res]
    rk = 0
    nzq = [t for t in qr if abs(t[0]) > mp.mpf(10) ** -30 or abs(t[1]) > mp.mpf(10) ** -30]
    if nzq:
        rk = 1
        for t in nzq[1:]:
            if abs(nzq[0][0] * t[1] - t[0] * nzq[0][1]) > mp.mpf(10) ** -30:
                rk = 2
    single = [r for r in rows if not r['om']]
    print(f"d={d}: orbits {len(rows)} (single-orbit certifiable: {len(single)}), dim K = {len(K)}, "
          f"dim Dist = {rD}, dim K/Dist = {len(comp)}, rank(q,r) = {rk} ({time.time()-t1:.1f}s)", flush=True)
    for x in res:
        sup = [(rows[i]['M'], rows[i]['j'], rows[i]['k'], str(c)) for i, c in enumerate(x['v']) if c]
        print(f"   COMP {sup}: q = {mp.nstr(x['q'], 20)} [{x['qf']}], r = {mp.nstr(x['r'], 20)} [{x['rf']}]", flush=True)
    summary[d] = dict(orbits=len(rows), single=[(r['M'], r['j'], r['k']) for r in single], dimK=len(K), dimDist=rD,
                      zero_formal=nformal, zero_unexplained=len(unexplained),
                      comp=[dict(sup=[(rows[i]['M'], rows[i]['j'], rows[i]['k'], str(c)) for i, c in enumerate(x['v']) if c],
                                 q=str(x['qf']), r=str(x['rf']), qnum=mp.nstr(x['q'], 30), rnum=mp.nstr(x['r'], 30))
                            for x in res], rank=rk,
                      obstructions=[dict(M=r['M'], j=r['j'], k=r['k'], x=r['x'], y=r['y'], N=r['n'],
                                         omega={str(p): v for p, v in r['om'].items()}) for r in rows])
print(f"done ({time.time()-t0:.1f}s)", flush=True)
print("POSITIVE (some certifiable combination with q != 0):",
      [d for d, s in summary.items() if any(c['q'] not in ('0', 'None') for c in s['comp'])])
if JOUT:
    json.dump(dict(b=b, KMAX=KMAX, MFMAX=MFMAX, orbits=len(orbits), points=npts, summary=summary), open(JOUT, 'w'),
              indent=1)
