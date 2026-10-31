//! mythos-trits: the Mythos blue star as the trits of the nested ternary network.
//! usage: mythos-trits STAR_GGUF OUT_DIR
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic)]

use mythos_trits::gguf;
use mythos_trits::hbp::{index, Hbp};
use mythos_trits::lattice::{
    collide, make_gates, modes, stream, totals, uncollide, unstream, Dims, Gate, K,
};
use mythos_trits::outer;
use mythos_trits::png::{blue, Canvas};
use mythos_trits::sha256::hex_digest;
use mythos_trits::star::{lay_in, trits};
use std::fs;
use std::path::Path;

const D: Dims = Dims {
    nx: 64,
    ny: 64,
    nz: 3,
};
const TMAX: usize = 192;

fn aligned(f: &[u8], star: &[u8]) -> i64 {
    f.chunks_exact(K)
        .zip(star.iter())
        .map(|(cell, &s)| cell.iter().filter(|&&n| n == s).count() as i64)
        .sum()
}

fn rho(f: &[u8]) -> Vec<u8> {
    f.chunks_exact(K)
        .map(|c| c.iter().map(|&n| u32::from(n)).sum::<u32>() as u8)
        .collect()
}

struct GateRun {
    readouts: Vec<[i64; K]>,
    frames: Vec<u8>,
}

/// One gate set: the free run and the gate run, t = 0..=192, with conservation asserted
/// every tick. Writes one SERIES row per tick per mode.
fn run_set(
    tag: &str,
    gates: &[Gate],
    start: &[u8],
    star: &[u8],
    rec: &mut Hbp,
    series: &mut Hbp,
) -> GateRun {
    let t0 = totals(start);
    let mut readouts = Vec::new();
    let mut frames = Vec::new();
    for coll in [false, true] {
        let mode = if coll { "gates" } else { "free" };
        let mut cur = start.to_vec();
        let mut nxt = vec![0u8; D.slots()];
        let mut returns = Vec::new();
        let mut focus_t1 = 0i64;
        for t in 0..=TMAX {
            assert_eq!(totals(&cur), t0, "{tag} {mode}: conservation at t={t}");
            let m = modes(&cur);
            let q: i64 = m.iter().sum();
            let al = aligned(&cur, star);
            if t == 1 {
                focus_t1 = al;
            }
            if t > 0 && cur == start {
                returns.push(t);
            }
            series.row(
                "TICK",
                &[
                    ("set", tag.into()),
                    ("mode", mode.into()),
                    ("t", t.to_string()),
                    ("Q", q.to_string()),
                    ("aligned", al.to_string()),
                    (
                        "modes",
                        m.iter().map(i64::to_string).collect::<Vec<_>>().join(","),
                    ),
                ],
            );
            if coll {
                readouts.push(m);
                if t % 4 == 0 {
                    frames.extend(rho(&cur));
                }
            }
            if t == TMAX {
                break;
            }
            if coll {
                collide(&mut cur, gates, t, D.cells());
            }
            stream(&D, &cur, &mut nxt);
            std::mem::swap(&mut cur, &mut nxt);
        }
        rec.row(
            "RUN",
            &[
                ("set", tag.into()),
                ("mode", mode.into()),
                ("ticks", TMAX.to_string()),
                ("conservation_exact_every_tick", "1".into()),
                ("aligned_at_t1", focus_t1.to_string()),
                ("slots", D.slots().to_string()),
                (
                    "returns_to_start_at",
                    if returns.is_empty() {
                        "none".into()
                    } else {
                        returns
                            .iter()
                            .map(usize::to_string)
                            .collect::<Vec<_>>()
                            .join(",")
                    },
                ),
            ],
        );
    }
    // reversibility, 72 forward and 72 back
    let mut cur = start.to_vec();
    let mut s = vec![0u8; D.slots()];
    for t in 0..72 {
        collide(&mut cur, gates, t, D.cells());
        stream(&D, &cur, &mut s);
        std::mem::swap(&mut cur, &mut s);
    }
    let moved = cur.iter().zip(start.iter()).filter(|(a, b)| a != b).count();
    for t in (0..72).rev() {
        unstream(&D, &cur, &mut s);
        std::mem::swap(&mut cur, &mut s);
        uncollide(&mut cur, gates, t, D.cells());
    }
    assert!(cur == start, "{tag}: reverse failed");
    rec.row(
        "REVERSE",
        &[
            ("set", tag.into()),
            ("steps", "72".into()),
            ("slots_differing_at_t72", moved.to_string()),
            ("every_star_trit_restored", "1".into()),
        ],
    );
    // null control
    let uni = vec![1u8; D.slots()];
    let mut u = uni.clone();
    for t in 0..72 {
        collide(&mut u, gates, t, D.cells());
        stream(&D, &u, &mut s);
        std::mem::swap(&mut u, &mut s);
    }
    rec.row(
        "NULL",
        &[
            ("set", tag.into()),
            ("uniform_unchanged_after_72", u8::from(u == uni).to_string()),
        ],
    );
    // outer network on every gate readout
    let mut ranks = std::collections::BTreeMap::new();
    let (mut cols, mut nl, mut qs) = (true, true, std::collections::BTreeSet::new());
    let mut spread_max = (0i64, 0usize);
    for (t, m) in readouts.iter().enumerate() {
        let j = outer::jacobian(m);
        *ranks.entry(outer::rank(&j)).or_insert(0usize) += 1;
        cols &= outer::column_sums_zero(&j);
        nl &= outer::sum_nl_zero(m);
        qs.insert(m.iter().sum::<i64>());
        let sp = m.iter().max().unwrap_or(&0) - m.iter().min().unwrap_or(&0);
        if sp > spread_max.0 {
            spread_max = (sp, t);
        }
    }
    let last = readouts.last().copied().unwrap_or([0; K]);
    let first_spread = |t: usize| {
        let m = readouts[t];
        m.iter().max().unwrap_or(&0) - m.iter().min().unwrap_or(&0)
    };
    rec.row(
        "OUTER",
        &[
            ("set", tag.into()),
            ("readouts", readouts.len().to_string()),
            (
                "Q_values",
                qs.iter().map(i64::to_string).collect::<Vec<_>>().join(","),
            ),
            (
                "rank_hist",
                ranks
                    .iter()
                    .map(|(r, c)| format!("{r}:{c}"))
                    .collect::<Vec<_>>()
                    .join(","),
            ),
            ("ones_T_J_zero_all", u8::from(cols).to_string()),
            ("sum_NL_zero_all", u8::from(nl).to_string()),
            ("spread_t1", first_spread(1).to_string()),
            ("spread_max", spread_max.0.to_string()),
            ("spread_max_at_t", spread_max.1.to_string()),
            ("spread_t192", first_spread(TMAX).to_string()),
            ("D_mod_P", outer::d_sum().to_string()),
            (
                "drift_law_243_steps_from_t192",
                u8::from(outer::drift_law(&last, 1, 243)).to_string(),
            ),
        ],
    );
    GateRun { readouts, frames }
}

/// The star: slices A, B, C side by side (x4), then all three summed.
fn draw_star(star: &[u8]) -> Vec<u8> {
    let (s, gap) = (4usize, 16usize);
    let w = 4 * 64 * s + 5 * gap;
    let h = 64 * s + 2 * gap;
    let mut c = Canvas::new(w, h, [2, 4, 12]);
    for panel in 0..4 {
        let x0 = gap + panel * (64 * s + gap);
        for y in 0..64 {
            for x in 0..64 {
                let v: i64 = if panel < 3 {
                    i64::from(star[(panel * 64 + y) * 64 + x])
                } else {
                    (0..3).map(|z| i64::from(star[(z * 64 + y) * 64 + x])).sum()
                };
                let max = if panel < 3 { 2 } else { 6 };
                c.rect(x0 + x * s, gap + y * s, s, s, blue(v, max));
            }
        }
    }
    c.png()
}

/// The network: density (A+B+C) at t = 0, 4, 24, 72, 192 and after 72 back, then the
/// 27 outer modes over 193 ticks (rest channel gold).
fn draw_net(tg: &GateRun, restored: &[u8]) -> Vec<u8> {
    let (s, gap) = (3usize, 12usize);
    let tile = 64 * s;
    let w = 6 * tile + 7 * gap;
    let h = tile + 3 * gap + 260;
    let mut c = Canvas::new(w, h, [2, 4, 12]);
    let cells = D.cells();
    let frame = |i: usize| &tg.frames[i * cells..(i + 1) * cells];
    let panels: [&[u8]; 6] = [frame(0), frame(1), frame(6), frame(18), frame(48), restored];
    let sum_z = |f: &[u8], x: usize, y: usize| -> i64 {
        (0..3).map(|z| i64::from(f[D.idx(x, y, z)])).sum()
    };
    let max = panels
        .iter()
        .flat_map(|f| {
            (0..64)
                .flat_map(move |y| (0..64).map(move |x| (x, y)))
                .map(move |(x, y)| sum_z(f, x, y))
        })
        .max()
        .unwrap_or(1);
    for (i, f) in panels.iter().enumerate() {
        let x0 = gap + i * (tile + gap);
        if i == 5 {
            c.rect(x0 - 2, gap - 2, tile + 4, tile + 4, [255, 215, 94]);
        }
        for y in 0..64 {
            for x in 0..64 {
                c.rect(x0 + x * s, gap + y * s, s, s, blue(sum_z(f, x, y), max));
            }
        }
    }
    let (cx, cy, cw, ch) = (gap, tile + 2 * gap, w - 2 * gap, 250usize);
    c.rect(cx, cy, cw, 1, [31, 58, 102]);
    c.rect(cx, cy + ch, cw, 1, [31, 58, 102]);
    let lo = tg
        .readouts
        .iter()
        .flat_map(|m| m.iter())
        .copied()
        .min()
        .unwrap_or(0);
    let hi = tg
        .readouts
        .iter()
        .flat_map(|m| m.iter())
        .copied()
        .max()
        .unwrap_or(1);
    let n = tg.readouts.len() as i64;
    let pt = |t: usize, v: i64| {
        (
            cx as i64 + t as i64 * (cw as i64 - 1) / (n - 1),
            (cy + ch) as i64 - (v - lo) * (ch as i64 - 4) / (hi - lo).max(1),
        )
    };
    for k in (0..K).filter(|&k| k != 13).chain(std::iter::once(13)) {
        let col = if k == 13 {
            [255, 215, 94]
        } else {
            [60, 120, 210]
        };
        for t in 1..tg.readouts.len() {
            let (x0, y0) = pt(t - 1, tg.readouts[t - 1][k]);
            let (x1, y1) = pt(t, tg.readouts[t][k]);
            c.line(x0, y0, x1, y1, col);
        }
    }
    c.png()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: mythos-trits STAR_GGUF OUT_DIR");
        std::process::exit(2);
    }
    let out = Path::new(&args[2]);
    fs::create_dir_all(out).expect("out");
    let raw = fs::read(&args[1]).expect("star");
    let g = gguf::parse(&raw).expect("gguf");
    let star = trits(&g);
    let mut rec = Hbp::new();
    let c = [0, 1, 2].map(|v| star.iter().filter(|&&t| t == v).count());
    let sum: usize = star.iter().map(|&t| usize::from(t)).sum();
    println!(
        "{}",
        rec.row(
            "STAR",
            &[
                ("gguf_bytes", raw.len().to_string()),
                ("gguf_sha256", hex_digest(&raw)),
                ("gguf_version", g.version.to_string()),
                ("tensors", g.tensors.len().to_string()),
                ("bytes_accounted", g.bytes_read_to.to_string()),
                ("cells", star.len().to_string()),
                ("t0", c[0].to_string()),
                ("t1", c[1].to_string()),
                ("t2", c[2].to_string()),
                ("trit_sum", sum.to_string()),
                ("trits_sha256", hex_digest(&star)),
            ],
        )
    );
    let start = lay_in(&D, &star);
    let t = totals(&start);
    println!(
        "{}",
        rec.row(
            "LATTICE",
            &[
                ("dims", "64x64x3".into()),
                ("cells", D.cells().to_string()),
                ("trits", D.slots().to_string()),
                ("star_copies", "27".into()),
                ("mass", t.mass.to_string()),
                (
                    "mass_eq_27x_star",
                    u8::from(t.mass == 27 * sum as i64).to_string()
                ),
                ("px", t.px.to_string()),
                ("py", t.py.to_string()),
                ("pz", t.pz.to_string()),
                ("e2", t.e2.to_string()),
            ],
        )
    );
    let mut files: Vec<(String, usize, String)> = Vec::new();
    let mut tg_run = None;
    for (tag, pin) in [("tg", true), ("v1", false)] {
        let gates = make_gates(pin);
        let mut table = Vec::new();
        for gt in &gates {
            table.extend(gt.ch.iter().map(|&k| k as u8));
            table.extend_from_slice(&gt.lut);
        }
        rec.row(
            "GATES",
            &[
                ("set", tag.into()),
                ("count", gates.len().to_string()),
                (
                    "fix_code40",
                    gates.iter().filter(|x| x.lut[40] == 40).count().to_string(),
                ),
                ("table_sha256", hex_digest(&table)),
            ],
        );
        let mut series = Hbp::new();
        let n0 = rec.len();
        let run = run_set(tag, &gates, &start, &star, &mut rec, &mut series);
        for line in rec.text().lines().skip(n0.saturating_sub(1)) {
            println!("{line}");
        }
        let name = format!("SERIES-{tag}.hbp");
        let text = series.text();
        fs::write(out.join(&name), &text).expect("series");
        files.push((name, text.len(), hex_digest(text.as_bytes())));
        let fname = format!("rho-frames-{tag}.bin");
        fs::write(out.join(&fname), &run.frames).expect("frames");
        files.push((fname, run.frames.len(), hex_digest(&run.frames)));
        if tag == "tg" {
            tg_run = Some(run);
        }
    }
    let tg = tg_run.expect("tg run");
    let star_png = draw_star(&star);
    let net_png = draw_net(&tg, &rho(&start));
    for (name, bytes) in [
        ("mythos-blue-star.png", &star_png),
        ("mythos-star-nested-network.png", &net_png),
    ] {
        fs::write(out.join(name), bytes).expect("png");
        files.push((name.into(), bytes.len(), hex_digest(bytes)));
    }
    fs::write(out.join("star.trits"), &star).expect("trits");
    files.push(("star.trits".into(), star.len(), hex_digest(&star)));
    let text = rec.text();
    fs::write(out.join("RESULT.hbp"), &text).expect("result");
    files.push(("RESULT.hbp".into(), text.len(), hex_digest(text.as_bytes())));
    let idx = index("RESULTINDEX", &files);
    fs::write(out.join("INDEX.hbi"), &idx).expect("index");
    let mut sums = String::new();
    for (f, _, h) in &files {
        sums.push_str(&format!("{h}  {f}\n"));
    }
    sums.push_str(&format!("{}  INDEX.hbi\n", hex_digest(idx.as_bytes())));
    fs::write(out.join("SHA256SUMS"), sums).expect("sums");
    println!("DONE|files={}|json=0", files.len() + 1);
}
