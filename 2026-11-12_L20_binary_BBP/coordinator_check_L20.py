# Coordinator's own check (12 Nov 2026) of the formulas lane's Corollary A and Theorem A, from their statements only.
from mpmath import mp, mpf, sqrt, pi, log, zeta, polylog, exp, mpc, im, ceil
mp.dps = 320
chi = {1: 1, 3: 1, 7: 1, 9: 1, 11: -1, 13: -1, 17: -1, 19: -1}          # Kronecker symbol (-20/n)
L = sum(s * zeta(2, mpf(a) / 20) for a, s in chi.items()) / 400        # L(2, chi_-20) via Hurwitz zeta
c20 = [1, 1, -1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, -1, 1, 1, 0]
def c(k):
    r = (k - 1) % 40
    return c20[r] if r < 20 else -c20[r - 20]
S = mpf(0)
for k in range(1, 2300):
    ck = c(k)
    if ck: S += mpf(ck) / (mpf(2) ** ((k + 1) // 2) * k * k)
lhs = L - sqrt(5) / 15 * pi * log(2)
print('L(2,chi_-20) =', mp.nstr(L, 30))
print('Corollary A: |L - (sqrt5/15) pi log2 - (4/3) S| =', mp.nstr(abs(lhs - mpf(4) / 3 * S), 5))
# Theorem A at 120 digits
mp.dps = 130
L = sum(s * zeta(2, mpf(a) / 20) for a, s in chi.items()) / 400
z = lambda a: exp(mpc(0, 1) * pi * a / 20) / sqrt(2)
T = im(polylog(2, z(1))) + im(polylog(2, z(7))) + im(polylog(2, z(9))) - im(polylog(2, z(17)))
print('Theorem A: |sum Im Li2 - ((3 sqrt5/4) L - (pi/4) log 2)| =', mp.nstr(abs(T - (3 * sqrt(5) / 4 * L - pi / 4 * log(2))), 5))
# red control: change c_1 to 2
mp.dps = 60
S2 = S + mpf(1) / 2
print('red control (c_1 + 1): difference', mp.nstr(abs(lhs - mpf(4) / 3 * S2), 5))
