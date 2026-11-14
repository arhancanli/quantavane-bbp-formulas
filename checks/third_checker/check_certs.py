# Referee's exact check of the four certificates, transcribed by hand from main.tex.
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from cyc import Field, check_certificate, five_args, reduce_symbol, orbit

# ---------------- Q(zeta_8) ----------------
K8 = Field(8)
z8 = K8.zpow
I8 = z8(2)
S2_8 = z8(1) + z8(-1)
assert S2_8 * S2_8 == K8.const(2)


def q(a, b, c, d):
    """(a + b*sqrt2) + i*(c + d*sqrt2) with rational a,b,c,d."""
    return (S2_8 * F(b) + F(a)) + I8 * (S2_8 * F(d) + F(c))


h = F(1, 2)
t = F(1, 3)
tab1 = [
    (-2, q(h, 0, F(3, 2), 1), q(F(1, 6), -t, -F(1, 6), F(1, 6))),
    (-2, q(h, 0, h, h), q(t, 0, 0, -t)),
    (2, q(0, h, -1, -h), q(0, 0, -1, 1)),
    (-2, q(h, 0, h, h), q(F(2, 3), 0, 0, -t)),
    (2, q(F(3, 2), h, h, 0), q(1, h, 0, -h)),
    (2, q(-1, 1, 2, -1), q(h, h, -h, -h)),
    (-1, q(0, 0, 0, F(1, 4)), q(1, -1, -1, 0)),
    (-2, q(-1, 1, 2, -1), q(0, F(2, 3), -t, 0)),
]
sqm2 = I8 * S2_8  # sqrt(-2)
target1 = [(6, sqm2 / 2), (1, sqm2 / 4), (-4, z8(1)), (-4, z8(3))]
got1, want1, diff1 = check_certificate(tab1, target1, "Table 1 / Thm 3.1, Q(zeta_8)")

tab5 = [
    (4, q(t, 0, 0, t), q(F(4, 3), -F(2, 3), t, -t)),
    (4, q(-1, 1, 2, -1), q(h, h, -h, -h)),
    (4, q(h, 0, 0, -h), q(1, 0, -1, 0)),
    (-4, q(1, -h, 0, -h), q(h, 0, 0, h)),
    (4, q(h, -h, -h, h), q(-h, -1, -h, -h)),
    (4, q(h, -h, -h, h), q(0, -h, 0, 0)),
    (2, q(0, 0, 0, F(1, 4)), q(1, -1, -1, 0)),
    (-2, q(-t, 0, 0, -F(2, 3)), q(F(4, 3), 0, 0, t)),
    (4, q(h, F(1, 6), F(1, 6), 0), q(-1, F(3, 2), 0, h)),
    (4, q(F(1, 4), F(1, 4), -F(1, 4), 0), q(2, -1, 0, 0)),
    (-1, q(-F(1, 9), 0, 0, -F(2, 9)), q(2, 0, 0, h)),
    (F(15, 2), q(t, 0, 0, -t), q(1, 0, 0, -1)),
    (4, q(h, -h, -h, h), q(1, h, 0, -h)),
    (4, q(-1, 1, 0, 0), q(0, h, -1, h)),
    (-4, q(F(3, 4), -F(1, 4), F(1, 4), 0), q(F(2, 3), -t, F(2, 3), t)),
    (4, q(-h, h, h, 0), q(1, h, 0, -h)),
    (-4, q(F(3, 2), -1, -h, h), q(1, h, 1, h)),
    (4, q(1, -h, -2, F(3, 2)), q(1, h, 0, h)),
    (4, q(t, 0, 0, t), q(1, h, 0, -h)),
    (-4, q(F(1, 4), F(1, 4), -F(1, 4), 0), q(h, 0, F(3, 2), -h)),
    (-4, q(-t, -F(2, 3), -t, -t), q(F(2, 3), -F(1, 6), -F(2, 3), -F(1, 6))),
    (1, q(-F(1, 9), 0, 0, -F(2, 9)), q(1, 0, 0, 1)),
    (-F(9, 2), q(-F(7, 9), 0, 0, -F(4, 9)), q(F(1, 9), 0, 0, F(2, 9))),
    (h, q(F(23, 27), 0, 0, -F(10, 27)), q(F(5, 9), 0, 0, F(1, 9))),
]
w = (1 + sqm2) / (1 - sqm2)
assert w == q(-t, 0, 0, F(2, 3)), "w = (-1+2sqrt(-2))/3"
target5 = [(F(27, 2), w), (-F(9, 2), w * w), (h, w * w * w), (-8, z8(1)), (-8, z8(3))]
got5, want5, diff5 = check_certificate(tab5, target5, "Table 5 / Thm 1.3 (BB (27)), Q(zeta_8)")

# ---------------- Q(zeta_40) ----------------
K = Field(40)
z = K.zpow
S2 = z(5) + z(-5)
assert S2 * S2 == K.const(2)


def inv(a):
    return a.inv()


def bar(a):
    return a.conj()


def za(a):
    return z(a) / S2


tab2 = [
    (-F(1, 4), 1 - z(-12), inv(1 - z(18))),
    (F(1, 2), z(30), bar(inv(1 - z(4)))),
    (-F(1, 2), S2 * z(-1), inv(1 - z(14))),
    (-F(1, 2), z(36), S2 * z(-1)),
    (F(1, 4), bar(inv(1 - z(5) / S2)), inv(1 - z(12))),
    (-F(1, 2), z(17) / S2, 1 - z(-2)),
    (-F(1, 2), z(35) / S2, S2 * z(-7)),
    (F(1, 4), z(4), z(10)),
    (F(1, 4), z(10), z(36)),
    (-F(1, 2), z(31) / S2, 1 - z(-6)),
    (F(1, 2), z(10), bar(1 - S2 * z(-17))),
    (F(1, 2), z(7) / S2, bar(1 - z(-18))),
]
target2 = [(1, za(1)), (1, za(7)), (1, za(9)), (-1, za(17))] + [(-F(3, 4), z(2 * b)) for b in (1, 3, 7, 9)]
got2, want2, diff2 = check_certificate(tab2, target2, "Table 2 / Thm 4.1, Q(zeta_40)")

tab3 = [
    (-F(1, 4), bar(1 - z(-2)), bar(inv(1 - z(2)))),
    (-F(7, 16), z(6), bar(inv(1 - z(14)))),
    (F(9, 16), bar(1 - z(-2)), inv(1 - z(16))),
    (-F(9, 16), z(28), z(38)),
    (-F(9, 8), z(34), 1 - z(-8)),
    (F(1, 4), z(17) / S2, bar(inv(1 - z(5) / S2))),
    (F(3, 16), z(4), bar(inv(1 - z(16)))),
    (-F(1, 8), 1 - z(-4), inv(1 - z(16))),
    (F(1, 4), z(17) / S2, inv(1 - z(7) / S2)),
    (-F(1, 4), z(30), (z(7) / S2) / (z(7) / S2 - 1)),
    (F(1, 4), bar(inv(1 - z(12))), inv(1 - z(18))),
    (-F(1, 4), inv(1 - z(6)), bar(inv(1 - z(7) / S2))),
    (F(1, 4), z(5) / S2, z(18)),
    (-F(1, 4), z(2), 1 - z(-12)),
    (-F(1, 4), 1 - z(-2), bar(inv(1 - z(7) / S2))),
    (-F(1, 4), z(26), bar(1 - S2 * z(-17))),
]
target3 = [(1, za(7)), (1, za(17)), (-F(21, 16), z(8)), (F(1, 16), z(16))]
got3, want3, diff3 = check_certificate(tab3, target3, "Table 3 / Thm 4.3 first, Q(zeta_40)")

tab4 = [
    (F(9, 16), z(6), z(6)),
    (F(1, 8), z(8), z(8)),
    (-F(3, 8), z(4), z(4)),
    (F(3, 8), z(4), bar(1 - z(-6))),
    (F(1, 2), z(5) / S2, z(12)),
    (-F(1, 2), z(2), 1 - S2 * z(-17)),
    (-F(1, 4), z(2), z(10)),
    (F(3, 8), z(4), 1 - z(-14)),
    (-F(1, 8), z(4), z(5) / S2),
    (-F(1, 8), z(1) / S2, bar(inv(1 - z(9) / S2))),
    (F(3, 8), z(1) / S2, z(34)),
    (F(3, 8), z(2), bar(1 - S2 * z(-9))),
    (F(1, 16), K.const(F(1, 2)), 1 - z(-8)),
    (-F(3, 8), z(5) / S2, z(26)),
    (F(3, 8), z(1) / S2, z(4)),
    (F(3, 8), z(1) / S2, bar(1 - z(-18))),
    (-F(1, 2), z(6), bar(1 - S2 * z(-9))),
    (-F(1, 4), z(4), z(6)),
]
target4 = ([(1, za(1)), (-1, za(17))] + [(-F(3, 8), z(2 * b)) for b in (1, 3, 7, 9)]
           + [(F(11, 16), z(8)), (F(10, 16), z(16))])
got4, want4, diff4 = check_certificate(tab4, target4, "Table 4 / Thm 4.3 second, Q(zeta_40)")

import cmath


def show_diff(diff, label):
    if not diff:
        return
    print("  residual for", label)
    for k, v in diff.items():
        print("   ", v, "* D(", k.to_complex(), ")")


for d_, lab in ((diff1, "T1"), (diff5, "T5"), (diff2, "T2"), (diff3, "T3"), (diff4, "T4")):
    show_diff(d_, lab)

# Sanity: the residual (if any) evaluated numerically with D.
