# t7_quinary.py -- checks for five-var/REPORT.md (27 Oct 2026): split cubic forms with QUINARY blocks H(u1..u5)
# (smooth cubic threefolds Y_H = {H = 0} in P^4), the cyclic cubic fourfolds Z_H = {H + y^3 = 0} in P^5 (sigma: y -> omega y),
# and the count of Hodge classes when every quinary block satisfies (H1') Hg(H^4_prim(Z_H)) = U(V_H, h_H) and
# (H3) Hg(J(Y_H)) = Sp_10, distinct classes being independent.
# (1) Jacobian rings over F_10007 (t1_sk_smooth.hilbert_function): R_H has Hilbert function (1+t)^5, R_{H+y^3} has (1+t)^6;
#     sigma-eigenspaces of H^{3,1}, H^{2,2}_prim, H^{1,3} of Z_H (residues, character omega^(deg_y A + 1)); signature;
#     H^{2,1}, H^{1,2} of Y_H.
# (2) Local Torelli for {Z_H} on the eigenspaces: ranks of (R_H)_3 -> Hom((R_H)_0, (R_H)_3) and -> Hom((R_H)_2, (R_H)_5).
# (3) Block table for n = 5 (Jacobian degree, dimension, character, excess) and its decomposition; totals (2^N+2)/3 and
#     binom(N, N/2) for split forms with blocks in 1, 3, 4, 5 variables.
# (4) Invariant theory: I_n(a,b) = dim (W^a (x) W*^b)^{SL_n} by Schur-Weyl/hook lengths (checked against the weight count
#     t5_quaternary.I_sl5 for n = 5); Sp_10-invariants of V^{(x)r} by oscillating tableaux with <= 5 rows.
# (5) Hodge-class count with quinary blocks: per character and via counts; equal to t5 without quinary blocks;
#     dim H_prim = (2^N+2)/3 and (m,m)-part = binom(N,N/2); specialisation H -> G + u^3 and H -> F + u^3 + v^3 never lowers
#     dim Hdg; the three P^4 in H(x) + H(y) = 0 and h^4 have a Gram matrix of rank 4 (= dim Hdg^4 for H = H').
# (6) Thresholds: first single hypersurface (exhaustive, N <= 62) and first product Z_H^a x E_omega^e with
#     SU(11)-determinant Hodge classes; the Weil-type abelian variety A_H x E_omega^9.
import itertools, json, os, random, sys, types
import flint
from collections import Counter
from fractions import Fraction
from math import comb, factorial, prod

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')
TC = HERE       # t5_quaternary.py is read from the folder of this script (read only)


def _load_t5():
    """t5_quaternary loaded from its source WITHOUT its line t1 = _load('t1_sk_smooth'): t1 has no main guard and rewrites
    ternary-cubics/outputs/t1_sk_smooth.{txt,json} when loaded. t5 uses t1 only in its own jacobian_checks (not called here);
    t5 also loads t3_formula, which has a main guard. Nothing in ternary-cubics/code or /outputs is written."""
    path = os.path.join(TC, 't5_quaternary.py')
    with open(path, encoding='utf-8') as fh:
        src = fh.read()
    line = "t1 = _load('t1_sk_smooth')"
    assert src.count(line) == 1
    mod = types.ModuleType('t5_quaternary')
    mod.__file__ = path
    exec(compile(src.replace(line, 't1 = None'), path, 'exec'), mod.__dict__)
    return mod


t5 = _load_t5()
binom, multinom, nu = t5.binom, t5.multinom, t5.nu
P = 10007


# ---------- helpers copied verbatim from ternary-cubics/code/t1_sk_smooth.py (md5 a7cbf922) ----------
def poly_eval(Pf, v, p):
    tot = 0
    for e, c in Pf.items():
        t = c
        for i, k in enumerate(e):
            if k:
                t = t * pow(v[i], k, p)
        tot += t
    return tot % p


def partials(Pf, nvar):
    D = []
    for i in range(nvar):
        Q = {}
        for e, c in Pf.items():
            if e[i]:
                e2 = list(e); e2[i] -= 1
                Q[tuple(e2)] = (Q.get(tuple(e2), 0) + c * e[i])
        D.append(Q)
    return D


def hilbert_function(Pf, nvar, p, tmax):
    """dim of (F_p[x]/(partials of P))_t for t = 0..tmax, via ranks over F_p (python-flint)."""
    D = partials(Pf, nvar)
    hf = []
    for t in range(tmax + 1):
        mons = [e for e in itertools.product(range(t + 1), repeat=nvar) if sum(e) == t]
        if t < 2:
            hf.append(len(mons)); continue
        idx = {e: i for i, e in enumerate(mons)}
        mus = [e for e in itertools.product(range(t - 1), repeat=nvar) if sum(e) == t - 2]
        M = flint.nmod_mat(len(mus) * len(D), len(mons), p)
        r = 0
        for mu in mus:
            for Q in D:
                for e, c in Q.items():
                    if c % p:
                        M[r, idx[tuple(a + b for a, b in zip(e, mu))]] = c % p
                r += 1
        hf.append(len(mons) - M.rank())
    return hf


def embed(Pf, nvar, offset):
    return {tuple([0] * offset + list(e) + [0] * (nvar - offset - len(e))): c for e, c in Pf.items()}


t1 = types.SimpleNamespace(poly_eval=poly_eval, partials=partials, hilbert_function=hilbert_function, embed=embed)


# ---------- linear algebra mod p ----------
def rref_mod(rows, p):
    piv = []                                    # [(pivot column, fully reduced row)]
    for r in rows:
        v = [x % p for x in r]
        for c, pr in piv:
            if v[c]:
                f = v[c]
                v = [(x - f * y) % p for x, y in zip(v, pr)]
        lead = next((i for i, x in enumerate(v) if x), None)
        if lead is None:
            continue
        inv = pow(v[lead], p - 2, p)
        v = [(x * inv) % p for x in v]
        new = []
        for c, pr in piv:
            if pr[lead]:
                f = pr[lead]
                pr = [(x - f * y) % p for x, y in zip(pr, v)]
            new.append((c, pr))
        piv = new + [(lead, v)]
    return piv


def rank_mod(rows, p):
    return len(rref_mod(rows, p))


def monos(n, d):
    return [e for e in itertools.product(range(d + 1), repeat=n) if sum(e) == d]


class JacRing:
    """(S/J)_d for d <= dmax, S = F_p[x_1..x_n], J = (partials of P); normal forms via RREF of J_d."""

    def __init__(self, Pf, n, p, dmax):
        self.n, self.p, self.data = n, p, {}
        D = t1.partials(Pf, n)
        for d in range(dmax + 1):
            ms = monos(n, d)
            idx = {e: i for i, e in enumerate(ms)}
            rows = []
            if d >= 2:
                for mu in monos(n, d - 2):
                    for Q in D:
                        v = [0] * len(ms)
                        for e, c in Q.items():
                            j = idx[tuple(x + y for x, y in zip(e, mu))]
                            v[j] = (v[j] + c) % p
                        rows.append(v)
            piv = rref_mod(rows, p)
            pc = {c for c, _ in piv}
            self.data[d] = (ms, idx, piv, [i for i in range(len(ms)) if i not in pc])

    def dim(self, d):
        return len(self.data[d][3])

    def basis(self, d):
        ms, idx, piv, b = self.data[d]
        return [ms[i] for i in b]

    def coords(self, poly, d):
        ms, idx, piv, b = self.data[d]
        v = [0] * len(ms)
        for e, c in poly.items():
            v[idx[e]] = (v[idx[e]] + c) % self.p
        for c, pr in piv:
            if v[c]:
                f = v[c]
                v = [(x - f * y) % self.p for x, y in zip(v, pr)]
        return [v[i] for i in b]


# ---------- (1)-(2) Jacobian rings, eigenspaces, local Torelli ----------
def jacobian_checks(p=P, seed=2026, trials=3):
    rnd, rows, ok = random.Random(seed), [], True
    for t in range(trials):
        H = t5.random_form(5, p, rnd)
        hfH = t1.hilbert_function(H, 5, p, 6)
        Z = dict(t1.embed(H, 6, 0))
        Z[(0, 0, 0, 0, 0, 3)] = 1
        hfZ = t1.hilbert_function(Z, 6, p, 7)
        # Z_H: N = 6, H^{4-q,q}_prim = R_{3q-3}, R = R_H (x) C[y]/(y^2); A Omega/(H+y^3)^(q+1) has character omega^(deg_y A + 1)
        eig = {}
        for q in range(5):
            d = 3 * q - 3
            if d < 0:
                continue
            c1 = hfH[d] if 0 <= d <= 6 else 0           # deg_y A = 0: character omega
            c2 = hfH[d - 1] if 1 <= d <= 7 else 0       # A = y * (R_H)_{d-1}: character omega^2
            eig[f'H^({4 - q},{q})'] = (c1, c2)
        # Hodge-Riemann on H^4_prim (weight 4): i^(p-q) int x xbar > 0, i.e. negative on (3,1),(1,3), positive on (2,2)
        sgn = {'H^(3,1)': -1, 'H^(2,2)': 1, 'H^(1,3)': -1}
        sig1 = (sum(v[0] for k, v in eig.items() if sgn.get(k) == 1), sum(v[0] for k, v in eig.items() if sgn.get(k) == -1))
        sig2 = (sum(v[1] for k, v in eig.items() if sgn.get(k) == 1), sum(v[1] for k, v in eig.items() if sgn.get(k) == -1))
        yh = dict(h21=hfH[1], h12=hfH[4])               # Y_H: N = 5, H^{3-q,q} = (R_H)_{3q-2}
        # (2) local Torelli: deformations (R_H)_3 of H; U^{3,1} = (R_H)_0 -> U^{2,2} = (R_H)_3 and
        #     Ubar^{2,2} = y (R_H)_2 -> Ubar^{1,3} = y (R_H)_5, both by multiplication
        R = JacRing(H, 5, p, 5)
        S3 = monos(5, 3)
        rows_a = [R.coords({G: 1}, 3) for G in S3]
        B2 = R.basis(2)
        rows_b = [sum((R.coords({tuple(x + y for x, y in zip(G, A)): 1}, 5) for A in B2), []) for G in S3]
        rk_a, rk_b = rank_mod(rows_a, p), rank_mod(rows_b, p)
        rk_ab = rank_mod([ra + rb for ra, rb in zip(rows_a, rows_b)], p)
        good = (hfH == [1, 5, 10, 10, 5, 1, 0] and hfZ == [1, 6, 15, 20, 15, 6, 1, 0]
                and eig == {'H^(3,1)': (1, 0), 'H^(2,2)': (10, 10), 'H^(1,3)': (0, 1), 'H^(0,4)': (0, 0)}
                and sig1 == (10, 1) and sig2 == (10, 1) and yh == dict(h21=5, h12=5)
                and R.dim(3) == 10 and R.dim(5) == 1 and rk_a == 10 and rk_b == 10 and rk_ab == 10)
        ok &= good
        rows.append(dict(trial=t, hilbert_H=hfH, hilbert_ZH=hfZ, eigen_dims_omega_omega2=eig,
                         signature_omega=sig1, signature_omega2=sig2, Y_H=yh,
                         torelli_rank_31_to_22=rk_a, torelli_rank_22_to_13=rk_b, torelli_rank_total=rk_ab, ok=good))
    return ok, rows


# ---------- (3) block table ----------
def table_checks():
    tab = {n: t5.block_pieces(n) for n in (1, 3, 4, 5)}
    ok = True
    q = {c: sorted(e for a, d, cc, e in tab[5] for _ in range(d) if cc == c) for c in (0, 1, 2)}
    ok &= q[0] == sorted([Fraction(-1, 2)] * 5 + [Fraction(1, 2)] * 5)              # H^3(Y_H): Sp_10, weights -+1/2
    mu = Fraction(-1, 6)                                                               # line of the y-variable
    ubar = sorted([Fraction(0)] * 10 + [Fraction(1)])                                  # Ubar: (2,2) x 10, (1,3) x 1
    u = sorted([Fraction(-1)] + [Fraction(0)] * 10)                                    # U: (3,1) x 1, (2,2) x 10
    ok &= q[1] == sorted(e + mu for e in ubar)                                         # character 1 = Ubar (x) mu
    ok &= q[2] == sorted(e - mu for e in u)                                            # character 2 = U (x) mu-bar
    det1 = sum(q[1])                                                                   # excess of wedge^11 of the chi=1 piece
    ok &= det1 == Fraction(-5, 6)
    for k3, k4, k5, l in itertools.product(range(4), range(3), range(3), range(9)):
        N = 3 * k3 + 4 * k4 + 5 * k5 + l
        if N % 2 or N < 2 or N > 18:
            continue
        dist = {(0, Fraction(0)): 1}
        for n in [3] * k3 + [4] * k4 + [5] * k5 + [1] * l:
            new = {}
            for (c, e), m in dist.items():
                for a, d, cc, ee in tab[n]:
                    key = ((c + cc) % 3, e + ee)
                    new[key] = new.get(key, 0) + m * d
            dist = new
        tot = sum(m for (c, e), m in dist.items() if c == 0)
        ok &= tot == (2 ** N + 2) // 3 and dist.get((0, Fraction(0)), 0) == comb(N, N // 2)
    return ok, {n: [(a, d, c, str(e)) for a, d, c, e in tab[n]] for n in tab}, str(det1)


# ---------- (4) invariant theory ----------
def partitions_of(n, maxlen, maxpart=None):
    maxpart = n if maxpart is None else maxpart
    if n == 0:
        yield ()
        return
    if maxlen == 0:
        return
    for p_ in range(min(n, maxpart), 0, -1):
        for rest in partitions_of(n - p_, maxlen - 1, p_):
            yield (p_,) + rest


def flam(lam):
    n = sum(lam)
    conj = [sum(1 for x in lam if x > j) for j in range(lam[0])] if lam else []
    h = 1
    for i, li in enumerate(lam):
        for j in range(li):
            h *= (li - j - 1) + (conj[j] - i - 1) + 1
    return factorial(n) // h


_ISL = {}


def I_SL(n, a, b):
    """dim (W^{(x)a} (x) W*^{(x)b})^{SL_n} = sum f^lam f^mu over lam |- a, mu |- b (<= n rows) with S_lam = S_mu on SL_n."""
    key = (n, a, b)
    if key not in _ISL:
        tot = 0
        if (a - b) % n == 0:
            for lam in partitions_of(a, n):
                L = list(lam) + [0] * (n - len(lam))
                for m_ in partitions_of(b, n):
                    M = list(m_) + [0] * (n - len(m_))
                    if len({x - y for x, y in zip(L, M)}) == 1:
                        tot += flam(lam) * flam(m_)
        _ISL[key] = tot
    return _ISL[key]


_JSP = {}


def J_sp(r, n=5):
    """dim (V^{(x)r})^{Sp_2n}: oscillating tableaux of length r from 0 to 0 with at most n rows (Sundaram)."""
    if (r, n) in _JSP:
        return _JSP[(r, n)]
    if r % 2:
        return 0
    cur = Counter({(): 1})
    for _ in range(r):
        nxt = Counter()
        for lam, m in cur.items():
            L = list(lam)
            for i in range(min(len(L) + 1, n)):
                if i == len(L):
                    new = L + [1]
                elif i == 0 or L[i - 1] > L[i]:
                    new = L[:]
                    new[i] += 1
                else:
                    continue
                nxt[tuple(new)] += m
            for i in range(len(L)):
                if i == len(L) - 1 or L[i] > L[i + 1]:
                    new = L[:]
                    new[i] -= 1
                    if new[i] == 0:
                        new = new[:i]
                    nxt[tuple(new)] += m
        cur = nxt
    _JSP[(r, n)] = cur.get((), 0)
    return _JSP[(r, n)]


def invariant_checks():
    ok, rows = True, {}
    for a in range(7):
        for b in range(7):
            if a + b <= 8:
                ok &= I_SL(5, a, b) == t5.I_sl5(a, b)                  # independent weight count in t5
    for a in range(13):
        for b in range(13):
            v = I_SL(11, a, b)
            if v:
                rows[f'SL11 I({a},{b})'] = v
            if a + b <= 10:
                ok &= v == (factorial(a) if a == b else 0)            # only contractions below 11 copies
    ok &= I_SL(11, 11, 0) == 1 and I_SL(11, 0, 11) == 1 and I_SL(11, 12, 1) == 11
    sp = {r: J_sp(r) for r in range(0, 15, 2)}
    dfact = {r: prod(range(1, r, 2)) for r in sp}                    # (r-1)!! perfect matchings
    ok &= all(sp[r] == dfact[r] for r in sp if r <= 10) and sp[12] == dfact[12] - 1
    ok &= J_sp(4, 1) == 2 and J_sp(6, 1) == 5                       # SL_2: Catalan numbers
    rows.update({f'Sp10 J({r})': v for r, v in sp.items()})
    return ok, rows


# ---------- (5) Hodge-class counts ----------
EX4 = Fraction(3, 10)      # excess per unit of (a_c - b_c) of an SU(5)-determinant (NOTE Prop C')
EX5 = Fraction(5, 66)      # minus the excess per unit of (a'_d - b'_d) of an SU(11)-determinant (wedge^11 of chi=1: -5/6)


def hodge_chi5(tern, quat, quin, l):
    """Per character. tern: [(label, has_cm)] ('w' = isogenous to E_omega); quat, quin: class labels (equal labels =
    equal block; distinct = independent); l cubes. Returns dim Hdg^m, the parts with SU(5)- and SU(11)-determinants,
    dim H_prim and its (m,m)-part from the block table."""
    k, q, r = len(tern), len(quat), len(quin)
    N = 3 * k + 4 * q + 5 * r + l
    assert N % 2 == 0
    qcl, rcl = sorted(set(quat)), sorted(set(quin))
    total, det4, det5, prim = 1, 0, 0, 0
    for chi in itertools.product(range(3), repeat=k):
        for psi in itertools.product(range(3), repeat=q):
            for phi in itertools.product(range(3), repeat=r):
                for c in itertools.product((1, 2), repeat=l):
                    if (sum(chi) + sum(psi) + sum(phi) + sum(c)) % 3:
                        continue
                    n1, n2, q0 = chi.count(1), chi.count(2), psi.count(0)
                    a, b, l1, l2 = psi.count(1), psi.count(2), c.count(1), c.count(2)
                    K = [s for s in range(k) if chi[s] == 0]
                    mult = 3 ** (n1 + n2) * 6 ** q0
                    prim += mult * 2 ** len(K) * 5 ** (a + b) * 10 ** phi.count(0) * 11 ** (phi.count(1) + phi.count(2))
                    inv, ex, any4, any5 = 1, Fraction(0), False, False
                    for cl in qcl:
                        ac = sum(1 for i in range(q) if quat[i] == cl and psi[i] == 1)
                        bc = sum(1 for i in range(q) if quat[i] == cl and psi[i] == 2)
                        inv *= t5.I(ac, bc)
                        ex += EX4 * (ac - bc)
                        any4 |= ac != bc
                    for cl in rcl:
                        r0 = sum(1 for i in range(r) if quin[i] == cl and phi[i] == 0)
                        ad = sum(1 for i in range(r) if quin[i] == cl and phi[i] == 1)
                        bd = sum(1 for i in range(r) if quin[i] == cl and phi[i] == 2)
                        inv *= I_SL(11, ad, bd) * J_sp(r0)
                        ex -= EX5 * (ad - bd)
                        any5 |= ad != bd
                    if not inv:
                        continue
                    Kw = [s for s in K if tern[s][0] == 'w']
                    x = Fraction(len(Kw), 2) + Fraction(n1 - n2 + a - b + l1 - l2, 6) - ex
                    term = inv * binom(len(Kw), x)
                    groups = {}
                    for s in K:
                        if tern[s][0] != 'w':
                            groups.setdefault(tern[s][0], []).append(s)
                    for ss in groups.values():
                        term *= nu(len(ss), tern[ss[0]][1])
                    total += mult * term
                    det4 += mult * term if any4 else 0
                    det5 += mult * term if any5 else 0
    tab = {n: t5.block_pieces(n) for n in (1, 3, 4, 5)}
    dist = {(0, Fraction(0)): 1}
    for n in [3] * k + [4] * q + [5] * r + [1] * l:
        new = {}
        for (cc, e), m in dist.items():
            for a_, d, c_, e_ in tab[n]:
                key = ((cc + c_) % 3, e + e_)
                new[key] = new.get(key, 0) + m * d
        dist = new
    return dict(hodge=total, det4=det4, det5=det5, dim_prim=prim, hmm_prim=dist.get((0, Fraction(0)), 0), N=N)


def splits(n, parts):
    if parts == 1:
        yield (n,)
        return
    for i in range(n + 1):
        for rest in splits(n - i, parts - 1):
            yield (i,) + rest


def hodge_counts5(kw, others, qsizes, rsizes, l):
    """The same over counts: kw ternary blocks ~ E_omega, others = [(size, has_cm)] further isogeny classes of ternary
    blocks, qsizes / rsizes = sizes of the classes of quaternary / quinary blocks, l cubes."""
    tern_opts = [[((n1, n2), multinom(n1, n2, K), K) for n1, n2, K in splits(kw, 3)]]
    for size, cm in others:
        tern_opts.append([((n1, n2), multinom(n1, n2, K) * nu(K, cm), None) for n1, n2, K in splits(size, 3)])
    quat_opts = [[((q0, a, b), multinom(q0, a, b) * 6 ** q0 * t5.I(a, b)) for q0, a, b in splits(m, 3)] for m in qsizes]
    quin_opts = [[((r0, a, b), multinom(r0, a, b) * J_sp(r0) * I_SL(11, a, b)) for r0, a, b in splits(m, 3)]
                 for m in rsizes]
    quat_opts = [[o for o in opts if o[1]] for opts in quat_opts]
    quin_opts = [[o for o in opts if o[1]] for opts in quin_opts]
    total, det4, det5 = 1, 0, 0
    for l1 in range(l + 1):
        l2 = l - l1
        for topt in itertools.product(*tern_opts):
            n1 = sum(o[0][0] for o in topt)
            n2 = sum(o[0][1] for o in topt)
            f_t = prod(o[1] for o in topt)
            if not f_t:
                continue
            Kw = topt[0][2]
            for qopt in itertools.product(*quat_opts):
                f_q = prod(o[1] for o in qopt)
                a = sum(o[0][1] for o in qopt)
                b = sum(o[0][2] for o in qopt)
                ex4 = sum(EX4 * (o[0][1] - o[0][2]) for o in qopt)
                d4 = any(o[0][1] != o[0][2] for o in qopt)
                for ropt in itertools.product(*quin_opts):
                    f_r = prod(o[1] for o in ropt)
                    a5 = sum(o[0][1] for o in ropt)
                    b5 = sum(o[0][2] for o in ropt)
                    if (n1 + 2 * n2 + a + 2 * b + a5 + 2 * b5 + l1 + 2 * l2) % 3:
                        continue
                    ex = ex4 - sum(EX5 * (o[0][1] - o[0][2]) for o in ropt)
                    x = Fraction(Kw, 2) + Fraction(n1 - n2 + a - b + l1 - l2, 6) - ex
                    term = comb(l, l1) * 3 ** (n1 + n2) * f_t * f_q * f_r * binom(Kw, x)
                    total += term
                    det4 += term if d4 else 0
                    det5 += term if any(o[0][1] != o[0][2] for o in ropt) else 0
    return dict(hodge=total, det4=det4, det5=det5)


def to_counts(tern):
    kw = sum(1 for lab, cm in tern if lab == 'w')
    groups = {}
    for lab, cm in tern:
        if lab != 'w':
            groups.setdefault(lab, [0, cm])[0] += 1
    return kw, [tuple(v) for v in groups.values()]


def sizes(labels):
    return sorted(Counter(labels).values())


W_, CMO = ('w', True), ('i', True)
GEN = [('g%d' % i, False) for i in range(8)]


def name(tern, quat, quin, l):
    t = ''.join('w' if lab == 'w' else ('i' if lab == 'i' else 'E' + lab[1:]) for lab, cm in tern)
    return f"[{t}|{','.join(quat)}|{','.join(quin)}|l={l}]"


CASES = [
    ([], [], ['H'], 1), ([], [], ['H'], 3), ([], [], ['H'], 5), ([W_], [], ['H'], 0), ([GEN[0]], [], ['H'], 0),
    ([], ['G'], ['H'], 1), ([], [], ['H', 'H'], 0), ([], [], ['H', 'K'], 0), ([], [], ['H', 'H'], 2),
    ([], [], ['H', 'K'], 2), ([W_], [], ['H', 'H'], 1), ([], ['G', 'G'], ['H'], 1), ([W_, W_], [], ['H'], 1),
    ([], [], ['H', 'H', 'H'], 1), ([], [], ['H', 'H', 'K'], 1), ([], [], ['H', 'H', 'H', 'H'], 0),
    ([], [], ['H', 'H'], 4), ([GEN[0], GEN[0]], [], ['H', 'H'], 0), ([], ['G'], ['H', 'H'], 2),
]
BIG = [  # counts only
    ([], [], ['H'] * 5, 1), ([], [], ['H'] * 10, 0), ([], [], ['H'] * 10, 2), ([], [], ['H'] * 11, 1),
    ([], [], ['H'] * 11, 3), ([], [], ['H'] * 11, 5), ([W_], [], ['H'] * 11, 2), ([], [], ['H'] * 11 + ['K'] * 11, 0),
    ([], ['G'] * 4, ['H'] * 10, 0),
]


def prop_checks():
    rows, ok = [], True
    for tern, quat, quin, l in CASES:
        r = hodge_chi5(tern, quat, quin, l)
        kw, oth = to_counts(tern)
        rc = hodge_counts5(kw, oth, sizes(quat), sizes(quin), l)
        N = r['N']
        good = (r['hodge'] == rc['hodge'] and r['det5'] == rc['det5'] and r['det4'] == rc['det4']
                and r['dim_prim'] == (2 ** N + 2) // 3 and r['hmm_prim'] == comb(N, N // 2))
        ok &= good
        rows.append(dict(case=name(tern, quat, quin, l), N=N, dim=N - 2, hodge=r['hodge'], det5=r['det5'],
                         counts=rc['hodge'], ok=good))
    for tern, quat, quin, l in BIG:
        kw, oth = to_counts(tern)
        rc = hodge_counts5(kw, oth, sizes(quat), sizes(quin), l)
        N = 3 * len(tern) + 4 * len(quat) + 5 * len(quin) + l
        rows.append(dict(case=name(tern, quat, quin, l), N=N, dim=N - 2, hodge=rc['hodge'], det5=rc['det5'],
                         counts=rc['hodge'], ok=True))
    # without quinary blocks: = t5 (closed formula and counts) on t5's cases and a sweep
    n5 = 0
    for tern, quat, l in t5.CASES:
        N = 3 * len(tern) + 4 * len(quat) + l
        if N % 2 or len(tern) + len(quat) + l > 9:
            continue
        a = t5.hodge_chi(tern, quat, l)['hodge']
        b = hodge_chi5(tern, quat, [], l)['hodge']
        kw, oth = to_counts(tern)
        c = hodge_counts5(kw, oth, sizes(quat), [], l)['hodge']
        ok &= a == b == c
        n5 += 1
    for k in range(4):
        for l in range(9):
            if (3 * k + l) % 2 or 3 * k + l < 4 or 3 * k + l > 12:
                continue
            for labs in itertools.product([W_, CMO, GEN[0], GEN[1]], repeat=k):
                labs = list(labs)
                kw, oth = to_counts(labs)
                ok &= t5.hodge_counts(kw, oth, [], l)['hodge'] == hodge_counts5(kw, oth, [], [], l)['hodge']
                n5 += 1
    return ok, rows, n5


def special_checks():
    """dim Hdg for very general H <= dim Hdg after specialising every class H to G + u^3 (G very general quaternary,
    one class per class) or to F + u^3 + v^3 (F Fermat, or a curve without CM)."""
    ok, rows = True, []
    for tern, quat, quin, l in CASES + BIG:
        kw, oth = to_counts(tern)
        gen = hodge_counts5(kw, oth, sizes(quat), sizes(quin), l)['hodge']
        cls = sorted(set(quin))
        spec_G = hodge_counts5(kw, oth, sizes(quat + ['s' + c for c in quin]), [], l + len(quin))['hodge']
        res = dict(case=name(tern, quat, quin, l), generic=gen, to_G=spec_G)
        good = spec_G >= gen
        if len(tern) + len(quin) <= 6:
            rep_w = list(tern) + [W_] * len(quin)
            rep_g = list(tern) + [('s' + c, False) for c in quin]
            for tag, rep in (('to_Fermat', rep_w), ('to_generic_F', rep_g)):
                kw2, oth2 = to_counts(rep)
                v = hodge_counts5(kw2, oth2, sizes(quat), [], l + 2 * len(quin))['hodge']
                res[tag] = v
                good &= v >= gen
        res['ok'] = good
        ok &= good
        rows.append(res)
    return ok, rows


def eightfold_planes(p=10009, seed=7):
    """H(x) + H(y) = 0 contains Lambda_k = {y = -omega^k x}, k = 0,1,2 (checked mod p, p = 1 mod 3); pairwise disjoint;
    Lambda.Lambda = c_4(N) = [(1+t)^5/(1+3t)]_4 = 11; Gram matrix with h^4 has rank 4."""
    assert p % 3 == 1
    rnd = random.Random(seed)
    w = next(x for x in range(2, p) if pow(x, (p - 1) // 3, p) != 1)
    w = pow(w, (p - 1) // 3, p)
    H = t5.random_form(5, p, rnd)
    okv = True
    for k in range(3):
        for _ in range(20):
            x = [rnd.randrange(p) for _ in range(5)]
            y = [(-pow(w, k, p) * xi) % p for xi in x]
            okv &= (t1.poly_eval(H, x, p) + t1.poly_eval(H, y, p)) % p == 0
    coeff = [sum(comb(5, i) * (-3) ** (j - i) for i in range(j + 1)) for j in range(5)]
    self_int = coeff[4]
    G = [[3, 1, 1, 1], [1, self_int, 0, 0], [1, 0, self_int, 0], [1, 0, 0, self_int]]
    rk = rank_mod(G, 1000003)
    return okv and self_int == 11 and rk == 4, dict(contained=okv, self_intersection=self_int, gram=G, rank=rk)


# ---------- (6) thresholds ----------
def det_threshold_single(Nmax=62):
    """Exhaustive for quinary-determinant classes: they need a class of >= 11 quinary blocks (FFT for SL_11), so
    N >= 55; all configurations with N <= Nmax are listed (ternary: kw ~ E_omega, kg pairwise non-isogenous without CM,
    kc with CM by another field; quaternary classes of size <= 4; quinary classes)."""
    found, n_conf = [], 0
    for rs in ([11], [12], [11, 1]):
        for kw in range(3):
            for kg in range(3):
                for kc in range(2):
                    for qs in ([], [1]):
                        for l in range(8):
                            N = 3 * (kw + kg + kc) + 4 * sum(qs) + 5 * sum(rs) + l
                            if N % 2 or N > Nmax:
                                continue
                            n_conf += 1
                            oth = [(1, False)] * kg + [(1, True)] * kc
                            r = hodge_counts5(kw, oth, qs, rs, l)
                            if r['det5']:
                                found.append(dict(N=N, dim=N - 2, quinary=rs, kw=kw, kg=kg, kc=kc, quaternary=qs, l=l,
                                                  det5=r['det5']))
    Nmin = min(f['N'] for f in found)
    return Nmin, [f for f in found if f['N'] == Nmin], n_conf


def det_threshold_products(amax=12, emax=4):
    """Products Z_H^a x E_omega^e: Hodge classes with an SU(11)-determinant factor. Pieces: H^4_prim(Z_H) = U + Ubar with
    U = (3,1) x 1 + (2,2) x 10 (excess -1, 0), Ubar = (2,2) x 10 + (1,3) x 1; H^*(Z_H) otherwise h^0..h^4 (excess 0);
    H^1(E_omega) = L (1,0) + Lbar (0,1) (excess -+1/2), H^0, H^2 (0). The torus condition det_K g = t^2 is equivalent to
    excess 0 on SU(11)-invariants (one Hodge circle)."""
    res = []
    for a_ in range(amax + 1):
        for e in range(emax + 1):
            cnt = 0
            for nu_, nb in itertools.product(range(a_ + 1), repeat=2):
                if nu_ + nb > a_ or nu_ == nb:
                    continue
                inv = I_SL(11, nu_, nb)
                if not inv:
                    continue
                ex = Fraction(-(nu_ - nb), 11)                     # wedge^11 U has excess -1
                for eL in range(e + 1):
                    for eLb in range(e - eL + 1):
                        if ex - Fraction(eL, 2) + Fraction(eLb, 2) != 0:
                            continue
                        rest = e - eL - eLb
                        cnt += (multinom(nu_, nb, a_ - nu_ - nb) * inv * 5 ** (a_ - nu_ - nb)
                                * multinom(eL, eLb, rest) * 2 ** rest)
            if cnt:
                res.append(dict(a=a_, e=e, dim=4 * a_ + e, det_classes_all_degrees=cnt))
    dmin = min(r['dim'] for r in res)
    return dmin, [r for r in res if r['dim'] == dmin]


def weil_type():
    """A_H = half twist (vGI Lemma 6.2): omega-eigenspace of H^1 has types (1,0) x 1, (0,1) x 10. With E_omega^9 (K acting
    by its CM, (1,0) each) the signature is (10,10); wedge^11 U_{1/2} (x) L^9 has excess 9/2 - 9/2 = 0."""
    sig = (1 + 9, 10 + 0)
    ex = Fraction(-1, 2) * 1 + Fraction(1, 2) * 10 - Fraction(1, 2) * 9
    return sig == (10, 10) and ex == 0, dict(dimension=11 + 9, signature=sig, weil_class_excess=str(ex))


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    lines, out, allok = [], {}, True

    ok, rows = jacobian_checks()
    allok &= ok
    out['jacobian'] = rows
    lines.append(f'(1)+(2) Jacobian rings over F_{P}, 3 random quinary H: Hilbert functions (1+t)^5, (1+t)^6; eigenspaces '
                 f'(omega, omega^2) of H^(3,1), H^(2,2)_prim, H^(1,3) of Z_H = (1,0), (10,10), (0,1); signature (10,1) on '
                 f'each eigenspace; h^(2,1)(Y_H) = h^(1,2) = 5; local Torelli ranks 10, 10, 10 -- {"ok" if ok else "FAIL"}')
    for r in rows:
        lines.append('    ' + json.dumps(r))

    ok, tab, det1 = table_checks()
    allok &= ok
    out['block_table'] = tab
    lines.append(f'(3) block table (a, dim, character, excess); quinary: chi 0 = H^3(Y_H) (weights -+1/2), chi 1 = Ubar (x) mu, '
                 f'chi 2 = U (x) mu-bar, mu = -1/6; excess of wedge^11 of the chi-1 piece = {det1}; totals (2^N+2)/3 and '
                 f'binom(N,N/2) with blocks of 1,3,4,5 variables, N <= 18 -- {"ok" if ok else "FAIL"}')
    for n in (1, 3, 4, 5):
        lines.append(f'    n = {n}: {tab[n]}')

    ok, rows = invariant_checks()
    allok &= ok
    out['invariants'] = rows
    lines.append(f'(4) SL_n invariants by Schur-Weyl (= t5 weight count for n = 5, a+b <= 8); SL_11: only contractions for '
                 f'a+b <= 10, I(11,0) = 1, I(12,1) = 11; Sp_10: (r-1)!! for r <= 10, 10394 for r = 12 -- '
                 f'{"ok" if ok else "FAIL"}')
    lines.append('    ' + ', '.join(f'{k} = {v}' for k, v in rows.items()))

    ok, rows, n5 = prop_checks()
    allok &= ok
    out['hodge_counts'] = rows
    lines.append(f'(5a) Hodge-class count with quinary blocks: per character = counts, dim H_prim = (2^N+2)/3, '
                 f'(m,m)-part = binom(N,N/2); = t5 on {n5} configurations without quinary blocks -- {"ok" if ok else "FAIL"}')
    lines.append('    [ternary | quaternary classes | quinary classes | l]  (equal letters = equal block)')
    for r in rows:
        lines.append(f"    {r['case']:<40} dim {r['dim']:>2}: Hdg {r['hodge']:>14}  SU(11)-det part {r['det5']:>3}"
                     f"  {'ok' if r['ok'] else 'FAIL'}")

    ok, rows = special_checks()
    allok &= ok
    out['specialisation'] = rows
    lines.append(f'(5b) specialisation H -> G + u^3 and H -> F + u^3 + v^3 never lowers dim Hdg ({len(rows)} cases) -- '
                 f'{"ok" if ok else "FAIL"}')
    for r in rows:
        lines.append('    ' + json.dumps(r))

    ok, info = eightfold_planes()
    allok &= ok
    out['eightfold_planes'] = info
    lines.append(f'(5c) H(x) + H(y) = 0: three disjoint P^4 {{y = -omega^k x}}, self-intersection 11, Gram rank with h^4 = '
                 f'{info["rank"]} (= dim Hdg^4 for H = H\') -- {"ok" if ok else "FAIL"}')

    Nmin, mins, nconf = det_threshold_single()
    good = Nmin == 60 and any(m['quinary'] == [11] and m['l'] == 5 and m['kw'] == 0 and m['det5'] == 2 for m in mins)
    allok &= good
    out['det_threshold_single'] = dict(Nmin=Nmin, minimal=mins, configurations=nconf)
    lines.append(f'(6a) SU(11)-determinant Hodge classes in a single hypersurface: first at N = {Nmin} '
                 f'(dimension {Nmin - 2}); {nconf} configurations with N <= 62 searched -- {"ok" if good else "FAIL"}')
    for m in mins:
        lines.append('    ' + json.dumps(m))

    dmin, pm = det_threshold_products()
    good = dmin == 46 and pm == [dict(a=11, e=2, dim=46, det_classes_all_degrees=2)]
    allok &= good
    out['det_threshold_products'] = dict(dmin=dmin, minimal=pm)
    lines.append(f'(6b) products Z_H^a x E_omega^e: first determinant classes in dimension {dmin}: {pm} -- '
                 f'{"ok" if good else "FAIL"}')

    ok, info = weil_type()
    allok &= ok
    out['weil_type'] = info
    lines.append(f'(6c) A_H x E_omega^9: {info} -- {"ok" if ok else "FAIL"}')

    lines.append('ALL CHECKS PASS' if allok else 'SOME CHECK FAILED')
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 't7_quinary.txt'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
    with open(os.path.join(OUT, 't7_quinary.json'), 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(out, fh, indent=1, default=str)
    print('\n'.join(lines))
