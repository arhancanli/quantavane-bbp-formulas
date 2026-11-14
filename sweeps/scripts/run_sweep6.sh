#!/bin/zsh
# Sweep 6 (11 Nov 2026 18:2x, coordinator): weight 3 and 4 L-values with natural partners.
# (A) binary weight 3, Q(sqrt2)/Q(sqrt-2) family; (B) binary weight 4: beta(4), L(4,chi_-8), L(4,chi_-3);
# (C) ternary weight 3, Q(sqrt3) family.
cd "${0:A:h}/.."
TA='z3;pi**3;pi**2*l2;l2**3;pi*G;G*l2;L(8,3);s2*L(8,3);pi**2*ls;ls**3;l2**2*ls;l2*ls**2;s2*pi**2*l2;s2*pi**3;pi*L(-8,2);s2*pi*L(-8,2);l2*L(-8,2);s2*l2*L(-8,2);ls*L(-8,2);s2*ls*L(-8,2);pi*L(-3,2);l2*L(-3,2)'
echo "### (A) binary weight 3, Q(sqrt2) family"
./scripts/guard.sh; nice -n 10 ./venv/bin/python -I scripts/sweep3.py 2 3 "4:8,8:16,6:24,12:24,8:32,16:32,12:48,24:48" "$TA"
TB='pi**4;pi**2*l2**2;l2**4;z3*l2;pi*z3;pi**3*l2;pi*l2**3;G*pi*l2;G*l2**2;G*pi**2;L(-4,4);L(-8,4);s2*L(-8,4);L(-3,4);s3*L(-3,4);s2*pi**3*l2;pi**2*L(-8,2);l2**2*L(-8,2);pi*l2*L(-8,2);s2*pi*l2*L(-8,2);s2*pi**2*L(-8,2)'
echo "### (B) binary weight 4"
./scripts/guard.sh; nice -n 10 ./venv/bin/python -I scripts/sweep3.py 2 4 "4:8,8:16,6:24,12:24,8:32,16:32,12:48" "$TB"
TC='z3;pi**3;pi**2*l3;l3**3;s3*pi**3;s3*pi**2*l3;s3*pi*l3**2;L(12,3);s3*L(12,3);pi**2*lt;lt**3;l3**2*lt;l3*lt**2;s3*pi**2*lt;pi*L(-3,2);l3*L(-3,2);s3*l3*L(-3,2);lt*L(-3,2);s3*pi*L(-3,2);pi**2*l2;l2*l3**2'
echo "### (C) ternary weight 3, Q(sqrt3) family"
./scripts/guard.sh; nice -n 10 ./venv/bin/python -I scripts/sweep3.py 3 3 "2:4,6:12,12:24,6:24,24:48,12:48" "$TC"
echo "SWEEP6 ALL DONE"
