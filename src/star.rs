//! The Mythos blue star -> trits, and the 27 shifted copies laid into the lattice.
//! Trit by the raw IEEE-754 bit pattern (positive floats order like their u32 bits):
//! 0 -> 0, 1..=0x3F80_0000 (value in (0, 1.0]) -> 1, above -> 2.

use crate::gguf::Gguf;
use crate::lattice::{vel, Dims, K};

pub const ONE_BITS: u32 = 0x3F80_0000;

pub fn trits(g: &Gguf) -> Vec<u8> {
    let mut out = Vec::new();
    for name in ["slice_A", "slice_B", "slice_C"] {
        let t = g
            .tensors
            .iter()
            .find(|t| t.name == name)
            .expect("slice tensor");
        assert_eq!(t.dims, vec![64, 64]);
        assert_eq!(t.ggml_type, 0);
        out.extend(t.mags.iter().map(|&w| match w {
            0 => 0u8,
            1..=ONE_BITS => 1,
            _ => 2,
        }));
    }
    out
}

/// n_k(x) = star(x + c_k): channel k holds the star shifted by c_k.
pub fn lay_in(d: &Dims, star: &[u8]) -> Vec<u8> {
    let mut f = vec![0u8; d.slots()];
    for z in 0..d.nz {
        for y in 0..d.ny {
            for x in 0..d.nx {
                let cell = d.idx(x, y, z);
                for k in 0..K {
                    f[cell * K + k] = star[d.shifted(x, y, z, vel(k))];
                }
            }
        }
    }
    f
}
