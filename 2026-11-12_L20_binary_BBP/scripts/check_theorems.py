"""Numeric check (mpmath polylog at complex points, a code path independent of the BBP sums) of Theorems A, B1, B2, B3.
z_a = exp(i pi a/20)/sqrt2.  Usage: python -I check_theorems.py digits"""
import sys
import mpmath as mp
P = int(sys.argv[1]); mp.mp.dps = P + 20
z = lambda a: mp.expjpi(mp.mpf(a) / 20) / mp.sqrt(2)
I = {a: mp.im(mp.polylog(2, z(a))) for a in (1, 7, 9, 17)}
chi = {1: 1, 3: 1, 7: 1, 9: 1, 11: -1, 13: -1, 17: -1, 19: -1}
L = mp.fsum(v * mp.zeta(2, mp.mpf(r) / 20) for r, v in chi.items()) / 400
C1, C2 = mp.clsin(2, 2 * mp.pi / 5), mp.clsin(2, 4 * mp.pi / 5)
s5, pl = mp.sqrt(5), mp.pi * mp.log(2)
checks = {
 'A : I1+I7+I9-I17 = (3 sqrt5/4) L - (pi/4) log2': I[1] + I[7] + I[9] - I[17] - (3 * s5 / 4 * L - pl / 4),
 'B1: I7+I17 = (21 Cl2(2pi/5) - Cl2(4pi/5))/16 - (3pi/20) log2': I[7] + I[17] - ((21 * C1 - C2) / 16 - 3 * pl / 20),
 'B2: I1-I17 = (6 sqrt5 L - 11 Cl2(2pi/5) - 10 Cl2(4pi/5))/16 - (pi/40) log2': I[1] - I[17] - ((6 * s5 * L - 11 * C1 - 10 * C2) / 16 - pl / 40),
 'B3: I9-I17 = (6 sqrt5 L - 10 Cl2(2pi/5) + 11 Cl2(4pi/5))/16 - (3pi/40) log2': I[9] - I[17] - ((6 * s5 * L - 10 * C1 + 11 * C2) / 16 - 3 * pl / 40),
 'RED (B1 with 21 -> 20)': I[7] + I[17] - ((20 * C1 - C2) / 16 - 3 * pl / 20),
 'RED (A with 3/4 -> 2/3)': I[1] + I[7] + I[9] - I[17] - (2 * s5 / 3 * L - pl / 4),
}
for k, v in checks.items():
    print(f"{k}: log10|diff| = {float(mp.log10(abs(v))) if v else '-inf'}")
