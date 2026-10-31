#!/usr/bin/env python3
"""Build the two gate includes from the repository source, changing only:
 - lattice 24x18x12 -> 64x64x3 (the star's own shape: 64x64, slices A,B,C)
 - its main() renamed so this project's main() drives it
 - tg: the Taylor-Green pairing rule from asolaria_taylor_green_regularity.cpp
   (code 40 = (1,1,1,1) kept as a fixed point of every gate)."""
import sys, os
src = open(os.path.join(sys.argv[1], 'asolaria_ternary_fluid_v1.cpp')).read()
a = 'constexpr int NX=24, NY=18, NZ=12'; b = 'int main(int argc,char**argv){'
assert src.count(a) == 1 and src.count(b) == 1
v1 = src.replace(a, 'constexpr int NX=64, NY=64, NZ=3').replace(b, 'int repo_main(int argc,char**argv){')
old = """        for(const auto& entry:cls) {
            const auto& codes=entry.second;"""
new = """        for(const auto& entry:cls) {
            std::vector<int> codes;  // Taylor-Green rule from asolaria_taylor_green_regularity.cpp:
            for(int code:entry.second) if(code!=40) codes.push_back(code);  // 40=(1,1,1,1) pinned"""
assert v1.count(old) == 1
out = sys.argv[2]
open(os.path.join(out, 'repo_gates_v1.inc'), 'w').write(v1)
open(os.path.join(out, 'repo_gates_tg.inc'), 'w').write(v1.replace(old, new))
print('wrote repo_gates_v1.inc, repo_gates_tg.inc')
