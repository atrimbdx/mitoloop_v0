#!/usr/bin/env bash
# Runs every test of tests/. Each test builds synthetic inputs, runs the
# scripts end to end and checks their outputs; it prints PASS <name> or exits
# non-zero.
#
# The .hic route needs hicstraw, java, and JUICER_TOOLS pointing to the
# Juicer Tools jar; without them it is skipped and the .cool route is tested.
#
# Usage
#   bash tests/run_all.sh
#   JUICER_TOOLS=path/to/juicer_tools.jar bash tests/run_all.sh

set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
failed=0
for t in "$here"/test_*.sh; do
    name="$(basename "$t" .sh)"
    start=$SECONDS
    if out="$(bash "$t" 2>&1)"; then
        echo "ok    $name ($((SECONDS - start)) s)"
    else
        echo "FAIL  $name"
        echo "$out" | tail -20 | sed 's/^/      /'
        failed=$((failed + 1))
    fi
done
if [ "$failed" -gt 0 ]; then
    echo "$failed test(s) failed"
    exit 1
fi
echo "all tests passed"
