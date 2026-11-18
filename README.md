# Binary BBP-type formulas for L(2, χ_d): certificates, search and checks

This repository holds the certificates, the search programs with their outputs, and the checking programs behind
the paper "Binary BBP-type formulas for L(2, χ_d): a Bloch-group criterion, an exhaustive search, and new formulas
for Q(√−2) and Q(√−5)" by Arhan Canli (Quantavane), November 2026.

The paper proves a criterion for when a rational combination of χ_d-twisted Galois orbit sums of the dilogarithm at
base-b points can be certified by five-term relations, reports an exhaustive search in a stated finite range
(base-2, 3 and 5 points with k ≤ 4 in Q(ζ_M), M ≤ 120, and base-2 points with k ≤ 6, M ≤ 168), and proves, by
certificates checked in exact arithmetic:

- binary formulas for L(2, χ_−8) − (3√2/32) π log 2 (base 64) and for L(2, χ_−20) − (√5/15) π log 2 (base 2^20)
  (Theorems 1.1 and 1.2);
- the identities (27), (29), (30) and (37) of Borwein and Broadhurst (1998), for the discriminants −8, −15, −20
  and −24 (Theorems 1.3 and 1.4);
- the value of a dilogarithm combination on the Galois orbit of e^(πi/20)/√2 in terms of D(i/2) and Catalan's
  constant (Theorem 6.3), and formulas for two series of conductor 5 (Proposition 7.1).

A certificate is a list of rational coefficients c_i and pairs (x_i, y_i) of elements of a cyclotomic field such
that the sum of c_i times the five-term relations of the Bloch–Wigner dilogarithm at (x_i, y_i) equals the claimed
combination, class by class, in exact arithmetic.

## What is here

The first three folders keep the names of the lab's working folders.

- `2026-11-11_L8_binary_BBP/`: the field Q(√−2). `certificate.json` has the 8 relations with all 40 arguments
  written exactly; `exceptional_units_q8.json` and `certificate_numeric_keys.json` are the inputs of
  `scripts/verify_cert.py`; `numerics_6000_digits.txt` is a numerical check. `scripts/` also has the
  integer-relation sweep (`sweep3.py`), the S-unit programs and the sparse certificate search.
- `2026-11-12_L20_binary_BBP/`: the field Q(√−5). `certificates/` has the three certificates in Q(ζ_40)
  (`cert_A.json`, `cert_B1.json`, `cert_B2.json`) with the search configurations, `scripts/` the certificate
  search, the exact verifier (`verify_cert_exact.py`) and the numerical programs, `outputs/` their outputs.
  `coordinator_check_L20.py` and `.txt` are a further numerical check written from the statements only.
- `2026-11-12_twisted_bloch_criterion/`: the search of Section 4 and the later certificates.
  `scripts/search_twisted_bloch.py` is the search and `outputs/search_b{2,3,5}_*` its complete output (every orbit,
  relative norm and obstruction vector, and every basis vector with its pair (q, r)). `scripts/check_nu_exact.py`
  recomputes the relative norms without floating point (`outputs/check_nu_exact_*.txt`: all 1320 pairs agree), and
  `scripts/prove_primes.py` proves the prime factors prime (`outputs/prove_primes_output.txt`).
  `certificates/cert_G.json` is the certificate of Theorem 6.3 (250 relations) with its search logs in `outputs/`
  and two red controls (damaged certificates that the verifiers must refuse). `borwein_broadhurst/` has the
  certificates of Theorem 1.4 (`cert_bb15.json`, `cert_bb24.json`, `cert_bb30.json`), their search scripts and
  logs, red controls, and the configurations and logs of the searches for the identity with discriminant −11,
  which gave no certificate. `catalan_row/` has a two-relation certificate in Q(ζ_8) for one row of Table 1.
- `checks/`: the paper's own checks. `third_checker/` is the third exact program, written in a separate session
  from the printed tables alone; `verify_q8_independent.py` is the second program for the certificates in Q(ζ_8)
  and `verify_cert_flint.py` the second verifier (FLINT polynomials) for Theorems 6.3 and 1.4; `cert_bb27.json` is
  the certificate of Theorem 1.3; `prove_cond5.py`, `prove_thm_G.py` and `branch_exact.py` check the exact finite
  steps of the proofs; `rows_q8.py` generates the printed rows of two tables from the exact elements; the
  `numeric_*.py` programs are numerical checks. Each program's output is next to it (`*_output.txt` and the other
  `.txt` files).
- `sweeps/`: the integer-relation sweeps that found the points (Section 9 of the paper): `scripts/` (the sweep
  programs and the shell scripts that ran them) and `logs/` (their logs and four intermediate JSON files).
- `results/`: `run_checks.sh`, which runs the quick checks below, and its output `results.txt`.

`SHA256SUMS` lists every file.

## How to re-run the checks

Commands are run from the repository root. The exact verifiers use the Python 3 standard library only. The
numerical programs need mpmath; the search, `check_nu_exact.py`, `prove_primes.py` and `verify_cert_flint.py` need
python-flint (and SymPy for `check_nu_exact.py`); some search scripts need numpy and scipy.

```sh
sh results/run_checks.sh         # PYTHON=/path/to/python3 to choose the interpreter; about three minutes
```

It runs the third checker on five certificates, the three programs for the exact steps of the proofs (each
compared with its stored output), the first exact verifier on eight certificates and on one red control, the
second program for Q(ζ_8) with a red control, and, when python-flint is installed, the FLINT verifier on four
certificates. The last line of `results/results.txt` must read `ALL PASS`.

One certificate by hand:

```sh
python3 -I 2026-11-12_twisted_bloch_criterion/scripts/verify_cert_exact.py \
    2026-11-12_twisted_bloch_criterion/certificates/cert_G.json
#   CERTIFICATE VERIFIED EXACTLY: sum_i c_i (five-term relation_i) = target combination
```

The search of Section 4 (`scripts/search_twisted_bloch.py b KMAX MFMAX out.json`) and the certificate searches
take longer; their outputs and logs are stored, and we have not run them again from this repository.

## Notes on the files

- Two scripts of `2026-11-12_L20_binary_BBP/scripts/` (`extract40.py`, `span_search_q40.py`) import modules from
  `2026-11-11_L8_binary_BBP/scripts/`; the absolute path of that folder on the lab's machine was replaced by a
  relative one. In the five `sweeps/scripts/run_sweep*.sh`, the line that changed to the working folder now changes
  to `sweeps/`; they call `./venv/bin/python`, a virtual environment with mpmath that is not part of this
  repository. Nothing else was changed.
- `2026-11-12_twisted_bloch_criterion/borwein_broadhurst/scripts/verify_cert_exact.py` imports `cyclo.py`, which
  is in `2026-11-12_twisted_bloch_criterion/scripts/`; put that folder on the module path (as
  `results/run_checks.sh` does).
- Some comments and output lines in `checks/` use the theorem and table numbers of earlier versions of the paper,
  which had two sections fewer before the present Section 5 (for example "Proposition 5.1" in `prove_cond5.py` is
  Proposition 7.1). Each output line also names the field and the number of relations of its certificate.
- The comments of `checks/third_checker/` call its author "the referee": that program was written by an AI agent
  in a separate session, from the printed tables of the paper alone. The word does not mean a journal referee.

## Credit and sources

- J. M. Borwein and D. J. Broadhurst, *Determinations of rational Dedekind-zeta invariants of hyperbolic manifolds
  and Feynman knots and links*, arXiv:hep-th/9811173 (1998), for the identities (27), (29), (30) and (37).
- D. H. Bailey's compendium of BBP-type formulas, whose notation P(s, b, n, A) the paper uses.
- FLINT (through python-flint), mpmath and SymPy, used by the search and by some of the checks.

The paper gives the full list of references.

## Authorship and the role of AI

Author and maintainer: Arhan Canli (Quantavane, https://quantavane.com), who takes responsibility for the results.

The search programs, the certificates, the checking programs and the runs were written and carried out by AI agents
(Claude, by Anthropic) working for Quantavane under Arhan Canli's direction. No mathematician outside Quantavane has
reviewed this repository; the certificates can be checked by anyone with the programs here or with their own.

## How to cite

If you use these certificates or programs, please cite them as below.

```bibtex
@misc{quantavane-bbp-formulas,
  author       = {Canli, Arhan},
  title        = {Binary BBP-type formulas for L(2, chi_d): certificates, search and checks},
  year         = {2026},
  howpublished = {\url{https://github.com/arhancanli/quantavane-bbp-formulas}},
  note         = {Quantavane; CC BY 4.0 for the data and text, MIT licence for the code}
}
```

[CITATION.cff](CITATION.cff) gives the same reference for GitHub's "Cite this repository" button.

## Licence

The certificates, the outputs, the logs and the text are under CC BY 4.0 (`LICENSE`); the code is under the MIT
licence (`LICENSE-CODE`). To cite: see `CITATION.cff`.
