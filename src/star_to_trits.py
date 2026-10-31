#!/usr/bin/env python3
"""Read the three 64x64 F32 slices of ACER-MYTHOS-BLUE.gguf and turn every
cell into a trit using integer comparisons of the raw IEEE-754 bit patterns
(positive floats order exactly like their uint32 bits):
    bits == 0                 -> 0   (dark)
    0 < bits <= 0x3F800000    -> 1   (value in (0, 1.0])
    bits >  0x3F800000        -> 2   (value above 1.0)
Writes star.trits: 3*64*64 bytes, order z (A,B,C) then y then x."""
import struct, sys, hashlib
p = sys.argv[1]
b = open(p, 'rb').read()
assert hashlib.sha256(b).hexdigest().startswith('3dd92a42160a38c1')
o = 4
def rd(f):
    global o; v = struct.unpack_from(f, b, o); o += struct.calcsize(f); return v[0]
def s():
    global o; n = rd('<Q'); v = b[o:o+n]; o += n; return v
def skip(t):
    sizes = {0:1,1:1,2:2,3:2,4:4,5:4,6:4,7:1,10:8,11:8,12:8}
    global o
    if t in sizes: o += sizes[t]
    elif t == 8: s()
    elif t == 9:
        et = rd('<I'); n = rd('<Q')
        for _ in range(n): skip(et)
ver = rd('<I'); nt = rd('<Q'); nk = rd('<Q')
for _ in range(nk): s(); skip(rd('<I'))
info = {}
for _ in range(nt):
    name = s().decode(); nd = rd('<I'); dims = [rd('<Q') for _ in range(nd)]; t = rd('<I'); off = rd('<Q')
    info[name] = (dims, t, off)
base = (o + 31) // 32 * 32
out = bytearray(); counts = [0, 0, 0]
for name in ('slice_A', 'slice_B', 'slice_C'):
    dims, t, off = info[name]; assert dims == [64, 64] and t == 0
    words = struct.unpack_from('<4096I', b, base + off)
    for w in words:
        tr = 0 if w == 0 else (1 if w <= 0x3F800000 else 2)
        counts[tr] += 1; out.append(tr)
open('star.trits', 'wb').write(out)
print(f"STAR_TRITS|cells={len(out)}|t0={counts[0]}|t1={counts[1]}|t2={counts[2]}|trit_sum={counts[1]+2*counts[2]}|sha256={hashlib.sha256(out).hexdigest()}")
