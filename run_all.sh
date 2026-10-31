#!/bin/sh
# Rebuild everything from the star file. Needs g++ (C++17) and python3; numpy/scipy/matplotlib only for the figures.
set -e
cd "$(dirname "$0")"
sha256sum -c star/ACER-MYTHOS-BLUE.gguf.sha256 --quiet 2>/dev/null || (cd star && sha256sum -c ACER-MYTHOS-BLUE.gguf.sha256)
mkdir -p build && python3 src/make_gate_includes.py src build
python3 src/star_to_trits.py star/ACER-MYTHOS-BLUE.gguf && mv star.trits build/
for g in tg v1; do
  g++ -O3 -std=c++17 -Wall -Wextra -Ibuild -DGATES_INC="\"repo_gates_$g.inc\"" -o build/net_$g src/star_nested_network.cpp
  mkdir -p results/$g && cp build/star.trits results/$g/ && (cd results/$g && ../../build/net_$g | tee log.txt)
  python3 src/outer_network.py results/$g/outer_modes.csv | tee results/$g/outer.txt
done
(cd results && sha256sum tg/* v1/* > SHA256SUMS)

if python3 -c "import numpy, scipy, matplotlib" 2>/dev/null; then
  python3 src/read_gguf.py star/ACER-MYTHOS-BLUE.gguf build/star.npz > build/gguf_header.txt
  python3 src/draw_star.py build/star.npz figures/mythos-blue-star.png
  python3 src/draw_net.py
  echo "figures rebuilt"
fi
