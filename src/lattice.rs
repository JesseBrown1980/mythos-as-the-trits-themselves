//! D3Q27 ternary lattice gas on an NX x NY x NZ periodic box. The gate network of
//! Neural-network-inside-a-neural-network-using-trits-to-resolve-Finiteness
//! (asolaria_ternary_fluid_v1.cpp); `pin_uniform` adds the pairing rule of
//! asolaria_taylor_green_regularity.cpp (code 40 = (1,1,1,1) fixed in every gate).

use std::collections::{BTreeMap, BTreeSet};

pub const K: usize = 27;
pub const GATES_PER_CELL: usize = 32;

/// Channel k = (z+1)*9 + (y+1)*3 + (x+1).
pub fn vel(k: usize) -> (i64, i64, i64) {
    (
        (k % 3) as i64 - 1,
        ((k / 3) % 3) as i64 - 1,
        (k / 9) as i64 - 1,
    )
}
pub fn energy(k: usize) -> i64 {
    let (x, y, z) = vel(k);
    x * x + y * y + z * z
}

#[derive(Clone, Copy)]
pub struct Dims {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
}

impl Dims {
    pub fn cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }
    pub fn slots(&self) -> usize {
        self.cells() * K
    }
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        (z * self.ny + y) * self.nx + x
    }
    /// Cell reached from (x, y, z) by velocity offset (dx, dy, dz), periodic.
    pub fn shifted(&self, x: usize, y: usize, z: usize, d: (i64, i64, i64)) -> usize {
        let w = |a: usize, s: i64, n: usize| (a as i64 + s).rem_euclid(n as i64) as usize;
        self.idx(w(x, d.0, self.nx), w(y, d.1, self.ny), w(z, d.2, self.nz))
    }
}

#[derive(Clone)]
pub struct Gate {
    pub ch: [usize; 4],
    pub lut: [u8; 81],
}

fn invariant(ch: &[usize; 4], code: usize) -> [i64; 5] {
    let mut c = code;
    let mut s = [0i64; 5];
    for &k in ch {
        let d = (c % 3) as i64;
        c /= 3;
        let (x, y, z) = vel(k);
        s[0] += d;
        s[1] += d * x;
        s[2] += d * y;
        s[3] += d * z;
        s[4] += d * energy(k);
    }
    s
}

pub fn make_gates(pin_uniform: bool) -> Vec<Gate> {
    let mut pairs: BTreeMap<(i64, i64, i64, i64), Vec<[usize; 2]>> = BTreeMap::new();
    for a in 0..K {
        for b in a + 1..K {
            let (ax, ay, az) = vel(a);
            let (bx, by, bz) = vel(b);
            pairs
                .entry((ax + bx, ay + by, az + bz, energy(a) + energy(b)))
                .or_default()
                .push([a, b]);
        }
    }
    let mut quads: BTreeSet<[usize; 4]> = BTreeSet::new();
    for p in pairs.values() {
        for i in 0..p.len() {
            for j in i + 1..p.len() {
                let mut q = [p[i][0], p[i][1], p[j][0], p[j][1]];
                q.sort_unstable();
                if q[0] != q[1] && q[1] != q[2] && q[2] != q[3] {
                    quads.insert(q);
                }
            }
        }
    }
    let mut out = Vec::new();
    for q in quads {
        let mut lut = [0u8; 81];
        let mut cls: BTreeMap<[i64; 5], Vec<usize>> = BTreeMap::new();
        for (code, slot) in lut.iter_mut().enumerate() {
            *slot = code as u8;
            cls.entry(invariant(&q, code)).or_default().push(code);
        }
        for codes in cls.values() {
            let kept: Vec<usize> = codes
                .iter()
                .copied()
                .filter(|&c| !(pin_uniform && c == 40))
                .collect();
            for pair in kept.chunks_exact(2) {
                lut[pair[0]] = pair[1] as u8;
                lut[pair[1]] = pair[0] as u8;
            }
        }
        for code in 0..81 {
            let dst = usize::from(lut[code]);
            assert_eq!(usize::from(lut[dst]), code, "gate not involutive");
            assert_eq!(
                invariant(&q, code),
                invariant(&q, dst),
                "gate breaks conservation"
            );
        }
        if lut.iter().enumerate().any(|(i, &v)| i != usize::from(v)) {
            out.push(Gate { ch: q, lut });
        }
    }
    out
}

fn apply_one(f: &mut [u8], base: usize, g: &Gate) {
    let code = usize::from(f[base + g.ch[0]])
        + 3 * usize::from(f[base + g.ch[1]])
        + 9 * usize::from(f[base + g.ch[2]])
        + 27 * usize::from(f[base + g.ch[3]]);
    let d = usize::from(g.lut[code]);
    f[base + g.ch[0]] = (d % 3) as u8;
    f[base + g.ch[1]] = ((d / 3) % 3) as u8;
    f[base + g.ch[2]] = ((d / 9) % 3) as u8;
    f[base + g.ch[3]] = ((d / 27) % 3) as u8;
}

fn which(step: usize, cell: usize, j: usize, total: usize) -> usize {
    (step * GATES_PER_CELL + cell * 17 + j * 37) % total
}

pub fn collide(f: &mut [u8], gates: &[Gate], step: usize, cells: usize) {
    for cell in 0..cells {
        for j in 0..GATES_PER_CELL {
            apply_one(f, cell * K, &gates[which(step, cell, j, gates.len())]);
        }
    }
}

pub fn uncollide(f: &mut [u8], gates: &[Gate], step: usize, cells: usize) {
    for cell in 0..cells {
        for j in (0..GATES_PER_CELL).rev() {
            apply_one(f, cell * K, &gates[which(step, cell, j, gates.len())]);
        }
    }
}

/// n_k(x + c_k, t+1) = n_k(x, t)
pub fn stream(d: &Dims, cur: &[u8], next: &mut [u8]) {
    for z in 0..d.nz {
        for y in 0..d.ny {
            for x in 0..d.nx {
                let src = d.idx(x, y, z) * K;
                for k in 0..K {
                    next[d.shifted(x, y, z, vel(k)) * K + k] = cur[src + k];
                }
            }
        }
    }
}

pub fn unstream(d: &Dims, after: &[u8], before: &mut [u8]) {
    for z in 0..d.nz {
        for y in 0..d.ny {
            for x in 0..d.nx {
                let dst = d.idx(x, y, z) * K;
                for k in 0..K {
                    before[dst + k] = after[d.shifted(x, y, z, vel(k)) * K + k];
                }
            }
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct Totals {
    pub mass: i64,
    pub px: i64,
    pub py: i64,
    pub pz: i64,
    pub e2: i64,
}

pub fn totals(f: &[u8]) -> Totals {
    let mut t = Totals::default();
    for cell in f.chunks_exact(K) {
        for (k, &n) in cell.iter().enumerate() {
            let n = i64::from(n);
            let (x, y, z) = vel(k);
            t.mass += n;
            t.px += n * x;
            t.py += n * y;
            t.pz += n * z;
            t.e2 += n * energy(k);
        }
    }
    t
}

pub fn modes(f: &[u8]) -> [i64; K] {
    let mut m = [0i64; K];
    for cell in f.chunks_exact(K) {
        for (k, &n) in cell.iter().enumerate() {
            m[k] += i64::from(n);
        }
    }
    m
}
