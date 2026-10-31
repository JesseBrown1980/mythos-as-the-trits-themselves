#![deny(clippy::float_arithmetic)]
//! The Rust receipts against the C++/Python history, tick for tick, and the fixed hashes.
use mythos_trits::lattice::make_gates;
use mythos_trits::sha256::hex_digest;
use std::fs;

const H: &str = "history/cxx-python-2026-10-31/results";

fn gate_table_sha(pin: bool) -> String {
    let mut b = Vec::new();
    for g in make_gates(pin) {
        b.extend(g.ch.iter().map(|&k| u8::try_from(k).unwrap()));
        b.extend_from_slice(&g.lut);
    }
    hex_digest(&b)
}

#[test]
fn sha256_vectors() {
    assert_eq!(
        hex_digest(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
}

#[test]
fn gate_tables_match_the_cpp() {
    assert_eq!(
        gate_table_sha(true),
        "a7fca71ac8a85fd3acbe4188e9c246e98e4beeeaa459f3e77d40c0342aa4bb12"
    );
    assert_eq!(
        gate_table_sha(false),
        "bd689114406588ac2f013b4a1cff3def5e6c55ddc084fe27d30335ccede505f5"
    );
}

#[test]
fn star_trits_match_history() {
    let new = fs::read("results/star.trits").unwrap();
    assert_eq!(new, fs::read(format!("{H}/tg/star.trits")).unwrap());
    assert_eq!(
        hex_digest(&new),
        "748ebf67136c3c1c4e344f22dbe2b1ce600cc480900c2c74ca3a6468fccbe0ee"
    );
}

#[test]
fn every_tick_matches_history() {
    for g in ["tg", "v1"] {
        let old = fs::read_to_string(format!("{H}/{g}/outer_modes.csv")).unwrap();
        let old: Vec<String> = old.lines().skip(1).map(str::to_string).collect();
        let new: Vec<String> = fs::read_to_string(format!("results/SERIES-{g}.hbp"))
            .unwrap()
            .lines()
            .map(|l| {
                let get = |k: &str| {
                    l.split('|')
                        .find_map(|kv| kv.strip_prefix(k).and_then(|r| r.strip_prefix('=')))
                        .unwrap()
                        .to_string()
                };
                format!(
                    "{},{},{},{},{}",
                    get("mode"),
                    get("t"),
                    get("modes"),
                    get("Q"),
                    get("aligned")
                )
            })
            .collect();
        assert_eq!(old.len(), 386);
        assert_eq!(old, new, "{g}");
        assert_eq!(
            fs::read(format!("{H}/{g}/rho_frames.bin")).unwrap(),
            fs::read(format!("results/rho-frames-{g}.bin")).unwrap()
        );
    }
}
