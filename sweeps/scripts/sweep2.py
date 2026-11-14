"""Sweep 2: longer ternary lengths for degree 3 (zeta(3) in base 3?).
Targets kept small: zeta(3), pi^2 log 3, log^3 3 (+ odd parts)."""
import sys, time
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import mpmath as mp
from sweep1 import run

T3s = ['z3', 'pi2l3', 'l3_3', 's3pi3', 's3pil3_2', 'pi3', 'pil3_2']

def main():
    ns = [int(x) for x in sys.argv[1].split(',')]
    for n in ns:
        for fn, fd in [(1, 2), (1, 3), (2, 3), (1, 1), (1, 4)]:
            if (n * fn) % fd:
                continue
            m = n * fn // fd
            run(3, 3 ** m, n, T3s)
    print("SWEEP2 DONE", flush=True)

if __name__ == '__main__':
    main()
