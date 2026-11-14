"""Check that every prime factor used by the search is a proven prime (FLINT fmpz.is_prime: a proven primality test,
as opposed to is_probable_prime), and record the factorizations' sizes."""
import json, sys
from flint import fmpz
tot = 0; mx = 0; allok = True
for f in sys.argv[1:]:
    J = json.load(open(f))['summary']
    ps = set()
    for s in J.values():
        for ob in s['obstructions']:
            for p, e in fmpz(ob['N']).factor():
                ps.add(int(p))
    ok = all(fmpz(p).is_prime() for p in ps)
    allok &= ok
    tot += len(ps); mx = max(mx, max(ps))
    print(f"{f}: {len(ps)} distinct prime factors, all proven prime: {ok}, largest {max(ps).bit_length()} bits")
print("ALL PRIME FACTORS PROVEN" if allok else "FAILED")
