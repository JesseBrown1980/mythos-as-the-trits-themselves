#!/usr/bin/env python3
"""Outer network: the 27 channel masses M_k (one per star copy) are the 27 shell
modes u_n. Integer-only checks over F_P, P = 1,000,081, as the repository README
specifies:  NL_n = u_n (u_{n+1} - u_{n-1}),  J_ii = u_{i+1}-u_{i-1},
J_{i,i+1} = u_i, J_{i,i-1} = -u_i,  1^T J = 0,  rank 26,
Q(u_t) + t*nu*D == Q(u_0) mod P with D = sum 9^n."""
import csv, sys
P = 1_000_081; N = 27
rows = [r for r in csv.DictReader(open(sys.argv[1])) if r['mode'] == 'gates']
def rank_mod(M):
    M = [row[:] for row in M]; r = 0
    for c in range(N):
        piv = next((i for i in range(r, N) if M[i][c] % P), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; inv = pow(M[r][c], P-2, P)
        M[r] = [x*inv % P for x in M[r]]
        for i in range(N):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(a - f*b) % P for a, b in zip(M[i], M[r])]
        r += 1
    return r
def jac(u):
    J = [[0]*N for _ in range(N)]
    for i in range(N):
        J[i][i] = (u[(i+1) % N] - u[(i-1) % N]) % P
        J[i][(i+1) % N] = (J[i][(i+1) % N] + u[i]) % P
        J[i][(i-1) % N] = (J[i][(i-1) % N] - u[i]) % P
    return J
D = sum(pow(9, n, P) for n in range(N)) % P
assert D == 883_973
ranks = {}; colsum_ok = True; nl_ok = True; Qset = set(); spread = []
for r in rows:
    u = [int(r[f'M{k}']) for k in range(N)]
    Qset.add(sum(u)); spread.append((int(r['tick']), max(u) - min(u)))
    J = jac(u)
    colsum_ok &= all(sum(J[i][j] for i in range(N)) % P == 0 for j in range(N))
    nl_ok &= sum(u[n]*(u[(n+1) % N] - u[(n-1) % N]) for n in range(N)) % P == 0
    rk = rank_mod(J); ranks[rk] = ranks.get(rk, 0) + 1
print(f"OUTER_READOUT|ticks={len(rows)}|Q_values={sorted(Qset)}|rank_hist={ranks}|ones_T_J_zero_all={colsum_ok}|sum_NL_zero_all={nl_ok}")
print("MODE_SPREAD|" + "|".join(f"t{t}={s}" for t, s in spread if t in (0,1,2,4,12,24,48,72,96,144,192)))
# run the outer recurrence itself from the t=192 readout, nu=1, 243 steps
u0 = [int(rows[-1][f'M{k}']) % P for k in range(N)]; u = u0[:]; nu = 1; ok = True
for t in range(1, 244):
    nl = [u[n]*(u[(n+1) % N] - u[(n-1) % N]) for n in range(N)]
    u = [(u[n] + nl[n] - nu*pow(9, n, P)) % P for n in range(N)]
    ok &= (sum(u) + t*nu*D - sum(u0)) % P == 0
print(f"OUTER_RECURRENCE|start=t192_readout|steps=243|drift_law_exact_every_step={ok}|D={D}")
