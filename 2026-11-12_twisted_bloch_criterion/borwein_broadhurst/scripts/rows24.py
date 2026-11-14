"""Readable rows for the 14-relation certificate of BB98 (37) in Q(zeta_24): each x, y written as
(a + b sqrt2 + c sqrt3 + d sqrt6) + i (e + f sqrt2 + g sqrt3 + h sqrt6) with rationals, checked exactly."""
import json, sys, cmath, math
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from cyclo import CycloField
d = json.load(open(sys.argv[1]))
M = d['M']; assert M == 24
K = CycloField(M); T = K.zeta; ONE = K.one()
SQ2 = K.add(T(3), T(-3)); SQ3 = K.add(T(2), T(-2)); I = T(6)
SQ6 = K.mul(SQ2, SQ3)
assert K.is_zero(K.sub(K.mul(SQ2, SQ2), K.elt([2]))) and K.is_zero(K.sub(K.mul(SQ3, SQ3), K.elt([3])))
basis = [ONE, SQ2, SQ3, SQ6]
basis = basis + [K.mul(I, b) for b in basis]
recipes = d['recipes']; memo = {}
def conj(x):
    r = K.elt([0])
    for i, c in enumerate(x):
        if c: r = K.add(r, K.scal(c, T(-i)))
    return r
def build(i):
    if i in memo: return memo[i]
    r = recipes[i]
    if r[0] == 'bbp':
        x = T(r[1]); assert r[2] == 0
    elif r[0] == 'pow':
        x = K.elt([0])
        for j, c in r[1].items(): x = K.add(x, K.scal(Fr(c), T(int(j))))
    elif r[0] == 'sym':
        w = build(r[1]); si = r[2]
        x = [w, K.sub(ONE, K.inv(w)), K.inv(K.sub(ONE, w)), K.inv(w), K.sub(ONE, w), K.div(w, K.sub(w, ONE))][si]
        if r[3]: x = conj(x)
    elif r[0] == 'rel':
        a, b, j = build(r[1]), build(r[2]), r[3]; dd = K.sub(ONE, K.mul(a, b))
        x = {2: K.div(K.sub(ONE, a), dd), 3: dd, 4: K.div(K.sub(ONE, b), dd)}[j]
    memo[i] = x; return x
# coordinates: solve exactly via linear algebra over Q on power-basis coordinates
import itertools
def coords(x):
    # basis matrix columns = basis elements' power-basis coordinates (8x8), solve B c = x
    n = 8
    A = [[basis[j][i] for j in range(n)] + [x[i]] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [v / pv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return [A[i][n] for i in range(n)]
def fmt_real(cs):
    names = ['', r'\sqrt2', r'\sqrt3', r'\sqrt6']
    terms = []
    for c, nm in zip(cs, names):
        if c == 0: continue
        sgn = '-' if c < 0 else '+'
        a = abs(c)
        if nm == '':
            body = (r'\frac{%d}{%d}' % (a.numerator, a.denominator)) if a.denominator != 1 else str(a.numerator)
        else:
            if a == 1: body = nm
            elif a.denominator == 1: body = '%d%s' % (a.numerator, nm)
            else: body = r'\frac{%d}{%d}%s' % (a.numerator, a.denominator, nm) if a.numerator != 1 else r'\frac{%s}{%d}' % (nm, a.denominator)
        terms.append((sgn, body))
    if not terms: return None
    s = ''.join(('' if i == 0 and sg == '+' else sg) + b for i, (sg, b) in enumerate(terms))
    return s, len(terms)
def fmt(x):
    c = coords(x)
    # exact check
    y = K.elt([0])
    for cc, b in zip(c, basis): y = K.add(y, K.scal(cc, b))
    assert K.is_zero(K.sub(x, y))
    re, im = fmt_real(c[:4]), fmt_real(c[4:])
    parts = []
    if re: parts.append(re[0] if re[1] == 1 else r'\bigl(%s\bigr)' % re[0])
    if im:
        body = im[0]
        if im[1] == 1:
            if body == '1': s_im = 'i'
            elif body == '-1': s_im = '-i'
            else: s_im = ('-' if body.startswith('-') else '') + 'i' + body.lstrip('-')
        else:
            s_im = r'i\bigl(%s\bigr)' % body
        if parts and not s_im.startswith('-'): parts.append('+' + s_im)
        else: parts.append(s_im)
    return ''.join(parts) if parts else '0'
rows = []
for item in d['cert']:
    rows.append((Fr(item['coef']), fmt(build(item['x'])), fmt(build(item['y']))))
for c, x, y in rows:
    cs = (r'\frac{%d}{%d}' % (c.numerator, c.denominator)) if c.denominator != 1 else str(c.numerator)
    if c < 0 and c.denominator != 1: cs = r'-\frac{%d}{%d}' % (-c.numerator, c.denominator)
    print('$%s$ & $%s$ & $%s$\\\\[2pt]' % (cs, x, y))
