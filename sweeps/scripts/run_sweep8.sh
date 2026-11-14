#!/bin/zsh
# Sweep 8 (12 Nov 2026 00:1x, coordinator): weight 2, the Q(sqrt-6) family: L(2,chi_-24) with sqrt2, sqrt3, sqrt6 multiples
# and pi times logs of 2, 3, 1+sqrt2, 2+sqrt3; base 2 (periods multiples of 24) and base 3 (periods multiples of 24).
# Controls inside: pi^2 (base 64, period 6k), L(2,chi_-3) (base 3/9 known), Catalan (base 2^12, period 24).
cd "${0:A:h}/.."
T='pi**2;l2**2;l3**2;l2*l3;pi*l2;pi*l3;pi*ls;pi*lt;s2*pi*ls;s3*pi*lt;G;L(-3,2);s3*L(-3,2);L(-8,2);s2*L(-8,2);L(-24,2);s2*L(-24,2);s3*L(-24,2);s6*L(-24,2)'
echo "### (A) base 2, weight 2, periods 24k"
./scripts/guard.sh; nice -n 10 ./venv/bin/python -I scripts/sweep3.py 2 2 "12:24,6:24,12:48,24:48,6:72,12:72" "$T"
echo "### (B) base 3, weight 2, periods 24k"
./scripts/guard.sh; nice -n 10 ./venv/bin/python -I scripts/sweep3.py 3 2 "2:24,6:24,2:48,6:48,12:48,6:72" "$T"
echo "SWEEP8 ALL DONE"
