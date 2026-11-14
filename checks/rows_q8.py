# Generates the LaTeX rows of Tables 1 and 5 (certificates in Q(zeta_8)) from the exact elements, and checks the
# printed text by evaluating it back in exact arithmetic (revision of 12 Nov 2026).
# Elements are (a, b, c, d) = (a + b sqrt2) + i (c + d sqrt2) with rational a, b, c, d.  The stored certificates name
# their points by rounded numerical keys; each key is looked up in the unit list (exponent vectors over zeta_8,
# 1 + sqrt2, 1 - zeta_8, 1 + sqrt(-2), 1 - sqrt(-2)), the exact element is built, and an assertion checks that it has
# the stored numerical value.  Everything after the lookup is exact.
# Usage: python3 -I rows_q8.py UNITS.json CERT.json OUT.tex
import sys, json, re
from fractions import Fraction as F
H = F(1, 2)
def mul(x, y):
    a, b, c, d = x; e, f, g, h = y
    return (a*e + 2*b*f - c*g - 2*d*h, a*f + b*e - c*h - d*g, a*g + 2*b*h + c*e + 2*d*f, a*h + b*g + c*f + d*e)
def inv(x):
    a, b, c, d = x
    p = a*a + 2*b*b + c*c + 2*d*d; q = 2*a*b + 2*c*d; den = p*p - 2*q*q
    return mul((a, b, -c, -d), (p/den, -q/den, F(0), F(0)))
ONE = (F(1), F(0), F(0), F(0))
def power(x, e):
    r = ONE; base = x if e >= 0 else inv(x)
    for _ in range(abs(e)):
        r = mul(r, base)
    return r
Z8 = (F(0), H, F(0), H)
GENS = [(F(1), F(1), F(0), F(0)), (F(1), -H, F(0), -H), (F(1), F(0), F(0), F(1)), (F(1), F(0), F(0), F(-1))]
def from_exps(e):
    x = power(Z8, e[0])
    for g, p in zip(GENS, e[1:]):
        x = mul(x, power(g, p))
    return x
def num(x):
    s2 = 2 ** 0.5
    return complex(float(x[0]) + float(x[1]) * s2, float(x[2]) + float(x[3]) * s2)
# ---------- formatting in the style of Table 1 ----------
def frac(p, q):
    return f'\\frac{{{p}}}{{{q}}}'
def term(t, root):                 # |t| (times sqrt2 if root), without sign
    t = abs(t)
    if not root:
        return str(t.numerator) if t.denominator == 1 else frac(t.numerator, t.denominator)
    if t.denominator == 1:
        return '\\sqrt2' if t.numerator == 1 else f'{t.numerator}\\sqrt2'
    return frac('\\sqrt2' if t.numerator == 1 else f'{t.numerator}\\sqrt2', t.denominator)
def lin(a, b):                     # a + b sqrt2 -> (string, number of terms, has fraction)
    parts = []
    for t, root in ((a, False), (b, True)):
        if t != 0:
            parts.append(('-' if t < 0 else '+') + term(t, root))
    s = ''.join(parts)
    if s.startswith('+'):
        s = s[1:]
    return s, len(parts), any(t != 0 and t.denominator != 1 for t in (a, b))
def imag_single(t, root):          # i times a single term, with sign
    sg = '-' if t < 0 else '+'
    u = abs(t)
    if not root:
        if u == 1:
            body = 'i'
        elif u.denominator == 1:
            body = f'{u.numerator}i'
        else:
            body = frac('i' if u.numerator == 1 else f'{u.numerator}i', u.denominator)
    else:
        body = 'i' + term(u, True)
    return sg, body
def tex(x):
    a, b, c, d = x
    re_s, re_n, re_f = lin(a, b)
    if c == 0 and d == 0:
        return re_s if re_s else '0'
    im_s, im_n, im_f = lin(c, d)
    if im_n == 1:
        sg, body = imag_single(c if c != 0 else d, d != 0)
        im_part = (sg if (re_s or sg == '-') else '') + body
    else:
        im_part = '+i' + ('\\bigl(' + im_s + '\\bigr)' if im_f else '(' + im_s + ')')
    if not re_s:
        return im_part.lstrip('+')
    if re_n == 2:
        re_s = ('\\bigl(' + re_s + '\\bigr)') if re_f else ('(' + re_s + ')')
    return re_s + im_part
def coef_tex(c):
    c = F(c)
    s = '-' if c < 0 else ''
    u = abs(c)
    return s + (str(u.numerator) if u.denominator == 1 else frac(u.numerator, u.denominator))
# ---------- exact evaluation of the printed text (round trip) ----------
class E:
    def __init__(self, v): self.v = v
    def __add__(self, o): return E(tuple(p + q for p, q in zip(self.v, lift(o).v)))
    __radd__ = __add__
    def __sub__(self, o): return E(tuple(p - q for p, q in zip(self.v, lift(o).v)))
    def __rsub__(self, o): return lift(o) - self
    def __mul__(self, o): return E(mul(self.v, lift(o).v))
    __rmul__ = __mul__
    def __truediv__(self, o): return E(mul(self.v, inv(lift(o).v)))
    def __rtruediv__(self, o): return lift(o) / self
    def __neg__(self): return E(tuple(-p for p in self.v))
    def __pos__(self): return self
def lift(o):
    return o if isinstance(o, E) else E((F(o), F(0), F(0), F(0)))
SQ, II = E((F(0), F(1), F(0), F(0))), E((F(0), F(0), F(1), F(0)))
def evaluate(s):
    """Exact value of the printed text: a small recursive-descent parser (no eval) for the grammar of our output:
    sums and differences of products/quotients of integers, sqrt2 (S), i (I) and parenthesised groups."""
    s = s.replace('\\bigl(', '(').replace('\\bigr)', ')')
    while '\\frac{' in s:          # arguments of \frac contain no braces in our output
        s = re.sub(r'\\frac\{([^{}]*)\}\{([^{}]*)\}', r'((\1)/(\2))', s)
    s = s.replace('\\sqrt2', 'S')
    toks = re.findall(r'\d+|[SI()+\-/]|i', s)
    assert ''.join(toks) == s.replace(' ', ''), ('unparsed text', s)
    toks = ['I' if t == 'i' else t for t in toks]
    pos = [0]
    def peek():
        return toks[pos[0]] if pos[0] < len(toks) else None
    def take():
        pos[0] += 1
        return toks[pos[0] - 1]
    def factor():
        t = take()
        if t == '-':
            return -factor()
        if t == '+':
            return factor()
        if t == '(':
            v = expr()
            assert take() == ')'
            return v
        if t == 'S':
            return SQ
        if t == 'I':
            return II
        assert t.isdigit(), t
        return lift(int(t))
    def term():
        v = factor()
        while True:
            t = peek()
            if t == '/':
                take(); v = v / factor()
            elif t is not None and (t.isdigit() or t in ('S', 'I', '(')):
                v = v * factor()          # implicit multiplication
            else:
                return v
    def expr():
        v = term()
        while peek() in ('+', '-'):
            v = v + term() if take() == '+' else v - term()
        return v
    v = expr()
    assert pos[0] == len(toks)
    return v.v
units = json.load(open(sys.argv[1])); cert = json.load(open(sys.argv[2]))
bykey = {(round(re_, 9) + 0.0, round(im_, 9) + 0.0): e for e, (re_, im_) in units}
rows, ok = [], True
for item in cert:
    x = from_exps(bykey[(round(item['x'][0], 9) + 0.0, round(item['x'][1], 9) + 0.0)])
    y = from_exps(bykey[(round(item['y'][0], 9) + 0.0, round(item['y'][1], 9) + 0.0)])
    assert abs(num(x) - complex(*item['x'])) < 1e-9 and abs(num(y) - complex(*item['y'])) < 1e-9
    tx, ty, tc = tex(x), tex(y), coef_tex(item['coef'])
    ok &= evaluate(tx) == x and evaluate(ty) == y and F(evaluate(tc)[0]) == F(item['coef'])
    rows.append(f'${tc}$ & ${tx}$ & ${ty}$\\\\[2pt]')
open(sys.argv[3], 'w').write('\n'.join(rows) + '\n')
print(len(rows), 'rows written to', sys.argv[3], '| printed text evaluates back to the exact elements:', ok)
