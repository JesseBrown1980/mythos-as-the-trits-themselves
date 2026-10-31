//! The outer network over F_P, P = 1,000,081: 27 shell modes u_n,
//! NL_n = u_n (u_{n+1} - u_{n-1}), J_ii = u_{i+1} - u_{i-1}, J_{i,i+1} = u_i, J_{i,i-1} = -u_i.

pub const P: i64 = 1_000_081;
pub const N: usize = 27;

fn md(a: i64) -> i64 {
    a.rem_euclid(P)
}

fn pow_mod(mut b: i64, mut e: i64) -> i64 {
    let mut r = 1i64;
    b = md(b);
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % P;
        }
        b = b * b % P;
        e >>= 1;
    }
    r
}

pub fn jacobian(u: &[i64; N]) -> Vec<Vec<i64>> {
    let mut j = vec![vec![0i64; N]; N];
    for i in 0..N {
        j[i][i] = md(u[(i + 1) % N] - u[(i + N - 1) % N]);
        j[i][(i + 1) % N] = md(j[i][(i + 1) % N] + u[i]);
        j[i][(i + N - 1) % N] = md(j[i][(i + N - 1) % N] - u[i]);
    }
    j
}

/// Rank mod P by Gauss-Jordan elimination.
pub fn rank(m: &[Vec<i64>]) -> usize {
    let mut a: Vec<Vec<i64>> = m.to_vec();
    let mut r = 0;
    for c in 0..N {
        let Some(p) = (r..N).find(|&i| a[i][c] != 0) else {
            continue;
        };
        a.swap(r, p);
        let inv = pow_mod(a[r][c], P - 2);
        for v in a[r].iter_mut() {
            *v = *v * inv % P;
        }
        let pivot = a[r].clone();
        for (i, row) in a.iter_mut().enumerate() {
            if i != r && row[c] != 0 {
                let f = row[c];
                for (x, &y) in row.iter_mut().zip(pivot.iter()) {
                    *x = md(*x - f * y);
                }
            }
        }
        r += 1;
    }
    r
}

pub fn column_sums_zero(j: &[Vec<i64>]) -> bool {
    (0..N).all(|c| md(j.iter().map(|row| row[c]).sum::<i64>()) == 0)
}

pub fn sum_nl_zero(u: &[i64; N]) -> bool {
    let s: i64 = (0..N)
        .map(|n| md(u[n] * (u[(n + 1) % N] - u[(n + N - 1) % N])))
        .sum();
    md(s) == 0
}

/// D = sum_{n<27} 9^n mod P.
pub fn d_sum() -> i64 {
    (0..N as i64).fold(0, |s, n| md(s + pow_mod(9, n)))
}

/// Run u_{t+1} = u + NL - nu*9^n for `steps`; true if Q(u_t) + t*nu*D == Q(u_0) mod P each step.
pub fn drift_law(u0: &[i64; N], nu: i64, steps: i64) -> bool {
    let d = d_sum();
    let q0 = md(u0.iter().sum());
    let mut u: Vec<i64> = u0.iter().map(|&x| md(x)).collect();
    let mut ok = true;
    for t in 1..=steps {
        let nl: Vec<i64> = (0..N)
            .map(|n| md(u[n] * (u[(n + 1) % N] - u[(n + N - 1) % N])))
            .collect();
        u = (0..N)
            .map(|n| md(u[n] + nl[n] - nu * pow_mod(9, n as i64)))
            .collect();
        let q = md(u.iter().sum());
        ok &= md(q + t * nu % P * d) == q0;
    }
    ok
}
