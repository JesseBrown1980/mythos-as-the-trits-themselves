#!/bin/sh
# Rebuild everything from the star file and check it against the committed receipts.
# Default: replay into replay/ and verify every file against results/SHA256SUMS (results/ is never touched).
# SEAL=1 ./run_all.sh writes a new baseline into results/ instead.
# Needs g++ (C++17) and python3; numpy/scipy/matplotlib only for the figures.
set -e
cd "$(dirname "$0")"
(cd star && sha256sum -c ACER-MYTHOS-BLUE.gguf.sha256)
OUT=replay
[ "${SEAL:-0}" = "1" ] && OUT=results
rm -rf build && mkdir -p build && python3 src/make_gate_includes.py src build
python3 src/star_to_trits.py star/ACER-MYTHOS-BLUE.gguf && mv star.trits build/
for g in tg v1; do
  g++ -O3 -std=c++17 -Wall -Wextra -Ibuild -DGATES_INC="\"repo_gates_$g.inc\"" -o build/net_$g src/star_nested_network.cpp
  rm -rf "$OUT/$g" && mkdir -p "$OUT/$g" && cp build/star.trits "$OUT/$g/" && (cd "$OUT/$g" && ../../build/net_$g > log.txt)
  python3 src/outer_network.py "$OUT/$g/outer_modes.csv" > "$OUT/$g/outer.txt"
done
if [ "$OUT" = "results" ]; then
  (cd results && sha256sum tg/* v1/* > SHA256SUMS)
  echo "SEALED|baseline=results/SHA256SUMS"
else
  (cd replay && sha256sum -c ../results/SHA256SUMS)
  echo "REPLAY|every_file_matches_committed_baseline=1"
fi
