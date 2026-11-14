#!/bin/sh
# Quick exact checks of this repository, about three minutes on one core. Run from the repository root:
#   sh results/run_checks.sh         (PYTHON=/path/to/python3 to choose the interpreter)
# Needs Python 3 with mpmath (one numerical line of checks/prove_cond5.py); python-flint only for step 5.
# Writes results/results.txt. Exit 0 only if every step gave the expected result.
#  1. The third checker (checks/third_checker), written from the printed tables alone: five certificates.
#  2. The exact finite steps of the proofs: checks/prove_cond5.py, prove_thm_G.py, branch_exact.py.
#  3. The first exact verifier on every certificate, and on one red control (a damaged certificate must fail).
#  4. The separate verifier for the two certificates in Q(zeta_8), with a red control.
#  5. If python-flint is installed: the second verifier (FLINT polynomials) on four certificates.
set -u
PY=${PYTHON:-python3}
L8=2026-11-11_L8_binary_BBP; L20=2026-11-12_L20_binary_BBP; TW=2026-11-12_twisted_bloch_criterion
bad=0
expect() {  # expect PATTERN LABEL COMMAND...: the last line of the command's output must contain PATTERN
  pat=$1; label=$2; shift 2
  out=$("$@" 2>&1 | tail -1)
  echo "$label: $out"
  case "$out" in *"$pat"*) ;; *) bad=1; echo "  UNEXPECTED (wanted: $pat)";; esac
}
{
echo "run_checks.sh $(date -u +%Y-%m-%dT%H:%M:%SZ); $($PY --version 2>&1)"

echo "=== 1. third checker"
$PY -I -B checks/third_checker/check_certs.py > results/.tmp 2>&1; cat results/.tmp
[ "$(grep -c 'EXACT MATCH' results/.tmp)" = 5 ] || { bad=1; echo "UNEXPECTED: five EXACT MATCH lines wanted"; }
cmp -s results/.tmp checks/third_checker/check_certs_output.txt && echo "same as checks/third_checker/check_certs_output.txt" || { bad=1; echo "DIFFERS from checks/third_checker/check_certs_output.txt"; }

echo "=== 2. exact steps of the proofs"
for s in prove_cond5 prove_thm_G branch_exact; do
  $PY -I -B checks/$s.py > results/.tmp 2>&1; tail -1 results/.tmp
  cmp -s results/.tmp checks/${s}_output.txt && echo "  same as checks/${s}_output.txt" || { bad=1; echo "  DIFFERS from checks/${s}_output.txt"; }
done

echo "=== 3. first exact verifier"
expect "VERIFIED EXACTLY" "$L8 certificate" $PY -I -B $L8/scripts/verify_cert.py $L8/exceptional_units_q8.json $L8/certificate_numeric_keys.json
for c in cert_A cert_B1 cert_B2; do
  expect "VERIFIED EXACTLY" "$L20 $c" $PY -I -B $L20/scripts/verify_cert_exact.py $L20/certificates/$c.json
done
expect "VERIFIED EXACTLY" "$TW cert_G" $PY -I -B $TW/scripts/verify_cert_exact.py $TW/certificates/cert_G.json
expect "MISMATCH" "$TW cert_G_red (red control)" $PY -I -B $TW/scripts/verify_cert_exact.py $TW/certificates/red_controls/cert_G_red.json
for c in cert_bb15 cert_bb24 cert_bb30; do
  # this copy of the verifier imports cyclo.py, which is in $TW/scripts, so that folder is put on the module path
  expect "VERIFIED EXACTLY" "$TW borwein_broadhurst $c" env PYTHONPATH=$TW/scripts $PY -B $TW/borwein_broadhurst/scripts/verify_cert_exact.py $TW/borwein_broadhurst/certificates/$c.json
done

echo "=== 4. separate verifier in Q(zeta_8)"
expect "VERIFIED EXACTLY" "L8" $PY -I -B checks/verify_q8_independent.py $L8/exceptional_units_q8.json $L8/certificate_numeric_keys.json L8
expect "VERIFIED EXACTLY" "BB27" $PY -I -B checks/verify_q8_independent.py $L8/exceptional_units_q8.json checks/cert_bb27.json BB27
expect "MISMATCH" "BB27 red control" $PY -I -B checks/verify_q8_independent.py $L8/exceptional_units_q8.json checks/cert_bb27_red_control.json BB27

echo "=== 5. second verifier (FLINT)"
if $PY -I -c "import flint" 2>/dev/null; then
  expect "VERIFIED EXACTLY" "cert_G" $PY -I -B checks/verify_cert_flint.py $TW/certificates/cert_G.json
  for c in cert_bb15 cert_bb24 cert_bb30; do
    expect "VERIFIED EXACTLY" "$c" $PY -I -B checks/verify_cert_flint.py $TW/borwein_broadhurst/certificates/$c.json
  done
else
  echo "python-flint is not installed for $PY: skipped"
fi
rm -f results/.tmp
[ $bad = 0 ] && echo "ALL PASS" || echo "SOMETHING FAILED"
} > results/results.txt 2>&1
tail -1 results/results.txt
[ $bad = 0 ]
