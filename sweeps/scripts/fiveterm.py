"""Search for five-term relations of the Bloch-Wigner dilogarithm D among points zeta_N^r * 2^(k/2)
(N = 24), and certify each found instance exactly in Q(zeta_24).

D(x) + D(y) + D((1-x)/(1-xy)) + D(1-xy) + D((1-y)/(1-xy)) = 0          (five-term relation)
D(z) = D(1-1/z) = D(1/(1-z)) = -D(1/z) = -D(1-z) = -D(z/(z-1)),  D(conj z) = -D(z)
Canonical unknowns: D(p), p = zeta^r 2^(-k/2), 1 <= r <= N/2 - 1, 0 <= k <= KMAX.
Output: JSON list of exact relations (dict canonical-key -> rational coefficient).
"""
import sys, json, cmath, math, itertools
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from fractions import Fraction as Fr
from cyclo import CycloField

N = 24
KMAX = 6     # |k| range for canonical points (override in main)
KX = 4       # |k| range for x, y seeds
K = CycloField(N)
S2 = K.add(K.zeta(N // 8), K.zeta(-(N // 8)))           # sqrt 2 in Q(zeta_24)
HALF = K.inv(S2)                                       # 1/sqrt2
TOL = 1e-9


def point_exact(r, k):
    """zeta^r * 2^(k/2) exactly."""
    e = K.zeta(r)
    f = S2 if k >= 0 else HALF
    for _ in range(abs(k)):
        e = K.mul(e, f)
    return e


def ident_point(w):
    """If w (complex) equals zeta^r 2^(k/2) numerically, return (r, k); else None."""
    if abs(w) < 1e-12:
        return None
    k2 = 2 * math.log2(abs(w))
    k = round(k2)
    if abs(k2 - k) > 1e-7 or abs(k) > KMAX:
        return None
    a = cmath.phase(w) * N / (2 * math.pi)
    r = round(a)
    if abs(a - r) > 1e-7:
        return None
    return (r % N, k)


ORBIT = [  # (map, sign): D(map(z)) = sign * D(z)
    (lambda z: z, 1), (lambda z: 1 - 1 / z, 1), (lambda z: 1 / (1 - z), 1),
    (lambda z: 1 / z, -1), (lambda z: 1 - z, -1), (lambda z: z / (z - 1), -1)]


def canon(w):
    """Return (key, sign) with D(w) = sign * D(canonical point), or (None, 0) if D(w) = 0 trivially,
    or None if w is not in the orbit of a point."""
    if abs(w.imag) < 1e-12 * max(1, abs(w)):
        return (None, 0)                       # real argument: D = 0
    for f, sgn in ORBIT:
        try:
            v = f(w)
        except ZeroDivisionError:
            continue
        # D(w) = sgn * D(v)?  D(f(z)) = sgn D(z) => with z = w: D(f(w)) = sgn D(w) => D(w) = sgn D(f(w))
        p = ident_point(v)
        if p is None:
            continue
        r, k = p
        s = sgn
        if k > 0:                               # D(1/p) = -D(p)
            r, k, s = (-r) % N, -k, -s
        if r > N // 2:                          # D(conj p) = -D(p)
            r, s = N - r, -s
        if r == 0 or r == N // 2:
            return (None, 0)
        return ((r, -k), s)                     # key (r, kk) means zeta^r 2^(-kk/2), kk >= 0
    return None


def main():
    global KMAX, KX
    KMAX = int(sys.argv[1]); KX = int(sys.argv[2])
    seeds = []
    for r in range(N):
        for k in range(-KX, KX + 1):
            z = cmath.exp(2j * math.pi * r / N) * 2 ** (k / 2)
            for f, _ in ORBIT:
                try:
                    v = f(z)
                except ZeroDivisionError:
                    continue
                seeds.append(v)
    # dedupe numerically
    uniq = {}
    for v in seeds:
        uniq[(round(v.real, 9), round(v.imag, 9))] = v
    seeds = list(uniq.values())
    print(f"seeds: {len(seeds)}", flush=True)
    rels = {}
    for x in seeds:
        if abs(x) < 1e-9 or abs(x - 1) < 1e-9:
            continue
        cx = canon(x)
        if cx is None:
            continue
        for y in seeds:
            if abs(y) < 1e-9 or abs(y - 1) < 1e-9 or abs(1 - x * y) < 1e-9:
                continue
            cy = canon(y)
            if cy is None:
                continue
            w3 = (1 - x) / (1 - x * y)
            w4 = 1 - x * y
            w5 = (1 - y) / (1 - x * y)
            cs = [cx, cy]
            ok = True
            for w in (w3, w4, w5):
                c = canon(w)
                if c is None:
                    ok = False
                    break
                cs.append(c)
            if not ok:
                continue
            vec = {}
            for key, s in cs:
                if key is None:
                    continue
                vec[key] = vec.get(key, 0) + s
            vec = {k: v for k, v in vec.items() if v != 0}
            if not vec:
                continue
            sig = tuple(sorted(vec.items()))
            if sig not in rels:
                rels[sig] = (x, y)
    print(f"distinct nontrivial five-term relations: {len(rels)}", flush=True)
    out = [{"x": [rels[s][0].real, rels[s][0].imag], "y": [rels[s][1].real, rels[s][1].imag],
            "rel": [[list(k), v] for k, v in s]} for s in rels]
    json.dump(out, open(sys.argv[3] if len(sys.argv) > 3 else 'fiveterm_rels.json', 'w'))


if __name__ == '__main__':
    main()
