# t5_quaternary.py -- checks for NOTE.md Sections 8-10 (Phase 99, 27 Oct 2026): split cubic forms with quaternary
# blocks G (smooth cubic surfaces S_G), the cyclic cubic threefolds Y_G = {G + y^3 = 0}, and Proposition C' (dim Hdg^m
# when every J(Y_G) has Hodge group U(V_G), hypothesis (H1) of Theorem D, and distinct classes are not isogenous).
# (1) Jacobian rings over F_10007 (python-flint, helper from t1_sk_smooth.py): Hilbert function (1+t)^n for random
#     blocks in n = 3, 4 variables and for Y_G; hence the sigma-eigenspaces of H^{2,1}(Y_G), H^{1,2}(Y_G) have
#     dimensions 4 + 1 and 1 + 4 (signature (4,1)).
# (2) The block table: Jacobian degree a, dimension, character (a + n) mod 3, excess (a + n)/3 - n/2.
# (3) I(a, b) = dim (W^{(x)a} (x) W*^{(x)b})^{SL_5} by the Brauer-Klimyk formula, against the first fundamental theorem.
# (4) Proposition C' three ways: closed formula per character (hodge_chi), the same summed over counts (hodge_counts),
#     and an independent count of invariants from explicit weights of the cover SU(5)^{classes} x G_m -> Hg with
#     t = s^5, det g_c = s^15 = t^3 (hodge_weights); and against Proposition C (t3_formula.py) without quaternary blocks.
# (5) Specialisation: replacing each class G by F + u^3 (F Fermat, or a curve without CM) never lowers dim Hdg.
# (6) Where determinant (Weil-type) Hodge classes first occur in a single hypersurface (all configurations, N <= 24).
import importlib.util, itertools, json, os, random, sys
from fractions import Fraction
from math import comb, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


t1 = _load('t1_sk_smooth')
t3 = _load('t3_formula')


def binom(a, b):
    if Fraction(b).denominator != 1:
        return 0
    b = int(b)
    return comb(a, b) if 0 <= b <= a else 0


def multinom(*ks):
    r, s = 1, 0
    for k in ks:
        s += k
        r *= comb(s, k)
    return r


# ---------- (1) Jacobian rings ----------
def random_form(nvar, p, rnd):
    return {e: rnd.randrange(1, p) for e in itertools.product(range(4), repeat=nvar) if sum(e) == 3}


def jacobian_checks(p=10007, seed=99, trials=3):
    rnd, rows, ok = random.Random(seed), [], True
    for n in (3, 4):
        for _ in range(trials):
            hf = t1.hilbert_function(random_form(n, p, rnd), n, p, n + 1)
            good = hf == [comb(n, t) for t in range(n + 1)] + [0]
            ok &= good
            rows.append(dict(block_vars=n, hilbert=hf, ok=good))
    for _ in range(trials):
        G = random_form(4, p, rnd)
        hfG = t1.hilbert_function(G, 4, p, 5)
        Y = dict(t1.embed(G, 5, 0))
        Y[(0, 0, 0, 0, 3)] = 1
        hfY = t1.hilbert_function(Y, 5, p, 6)
        # R_Y = R_G (x) C[y]/(y^2); A.Omega/Y^(p+1) has sigma-character deg_y(A) + 1, so
        # H^{2,1} = R_1 = (R_G)_1 (chi 1) + y (R_G)_0 (chi 2), H^{1,2} = R_4 = (R_G)_4 (chi 1) + y (R_G)_3 (chi 2).
        sig = dict(h21=(hfG[1], hfG[0]), h12=(hfG[4], hfG[3]))
        good = (hfG == [1, 4, 6, 4, 1, 0] and hfY == [1, 5, 10, 10, 5, 1, 0]
                and sig == dict(h21=(4, 1), h12=(1, 4)))
        ok &= good
        rows.append(dict(Y_G=True, hilbert_G=hfG, hilbert_Y=hfY, eigen_dims_chi1_chi2=sig, ok=good))
    return ok, rows


# ---------- (2) block table ----------
def block_pieces(n):
    """[(a, dim, character, excess)] for a smooth cubic block in n variables (R = C[z]/(z^2) for n = 1)."""
    return [(a, comb(n, a), (a + n) % 3, Fraction(a + n, 3) - Fraction(n, 2)) for a in range(n + 1)]


def table_checks():
    tab = {n: block_pieces(n) for n in (1, 3, 4)}
    ok = True
    # quaternary block: character 1 = W (x) mu, character 2 = W-bar (x) mu-bar, W of signature (4,1), mu of excess -1/6
    q1 = sorted(e for a, d, c, e in tab[4] for _ in range(d) if c == 1)
    q2 = sorted(e for a, d, c, e in tab[4] for _ in range(d) if c == 2)
    q0 = [e for a, d, c, e in tab[4] for _ in range(d) if c == 0]
    ok &= q1 == sorted([Fraction(-1, 2) - Fraction(1, 6)] + [Fraction(1, 2) - Fraction(1, 6)] * 4)
    ok &= q2 == sorted([Fraction(1, 2) + Fraction(1, 6)] + [Fraction(-1, 2) + Fraction(1, 6)] * 4)
    ok &= q0 == [0] * 6
    # totals: sum over characters 0 of all blocks = (2^N + 2)/3, its excess-0 part = binom(N, N/2)
    for k3, k4, l in itertools.product(range(5), range(5), range(9)):
        N = 3 * k3 + 4 * k4 + l
        if N % 2 or N == 0 or N > 16:
            continue
        dist = {(0, Fraction(0)): 1}
        for n in [3] * k3 + [4] * k4 + [1] * l:
            new = {}
            for (c, e), m in dist.items():
                for a, d, cc, ee in tab[n]:
                    key = ((c + cc) % 3, e + ee)
                    new[key] = new.get(key, 0) + m * d
            dist = new
        tot = sum(m for (c, e), m in dist.items() if c == 0)
        mid = dist.get((0, Fraction(0)), 0)
        ok &= tot == (2 ** N + 2) // 3 and mid == comb(N, N // 2)
    return ok, {n: [(a, d, c, str(e)) for a, d, c, e in tab[n]] for n in tab}


# ---------- (3) SL_5 invariants ----------
RHO = (4, 3, 2, 1, 0)


def _sign(p):
    s, seen = 1, set()
    for i in range(len(p)):
        if i in seen:
            continue
        j, L = i, 0
        while j not in seen:
            seen.add(j)
            j = p[j]
            L += 1
        s *= -1 if L % 2 == 0 else 1
    return s


def norm5(w):
    return tuple(x - w[4] for x in w)


A5 = {}
for _p in itertools.permutations(range(5)):
    A5[norm5(tuple(RHO[_p[i]] - RHO[i] for i in range(5)))] = _sign(_p)   # w(rho) - rho -> sgn(w)
assert len(A5) == 120
E5 = [tuple(1 if j == i else 0 for j in range(5)) for i in range(5)]


def weights_W(a, b):
    """weight multiset of W^{(x)a} (x) W*^{(x)b} for SL_5, as {weight in Z^5: multiplicity}."""
    dist = {(0,) * 5: 1}
    for sgn in [1] * a + [-1] * b:
        new = {}
        for w, m in dist.items():
            for e in E5:
                v = tuple(x + sgn * y for x, y in zip(w, e))
                new[v] = new.get(v, 0) + m
        dist = new
    return dist


def I_sl5(a, b):
    """dim of SL_5-invariants in W^{(x)a} (x) W*^{(x)b} (Brauer-Klimyk: sum_w sgn(w) m(w rho - rho))."""
    return sum(m * A5.get(norm5(w), 0) for w, m in weights_W(a, b).items())


def fft_checks():
    rows, ok = {}, True
    for a in range(7):
        for b in range(7):
            if a + b > 8:
                continue
            v = I_sl5(a, b)
            rows[f'{a},{b}'] = v
            if a + b <= 4:
                ok &= v == (factorial(a) if a == b else 0)          # only contractions below 5 copies
            if (a - b) % 5:
                ok &= v == 0
    ok &= rows['5,0'] == 1 and rows['0,5'] == 1 and rows['3,3'] == 6 and rows['4,4'] == 24
    ok &= rows['6,1'] == 5 and rows['1,6'] == 5                    # f^{(2,1,1,1,1)} = 5
    return ok, rows


_I = {}


def I(a, b):
    if (a, b) not in _I:
        _I[(a, b)] = I_sl5(a, b) if (a - b) % 5 == 0 else 0
    return _I[(a, b)]


def nu(size, cm):
    if size % 2:
        return 0
    e = size // 2
    return comb(2 * e, e) if cm else comb(2 * e, e) // (e + 1)


# ---------- (4) Proposition C' ----------
def hodge_chi(tern, quat, l):
    """Proposition C' summed over characters. tern: [(label, has_cm)] ('w' = isogenous to E_omega), as in t3;
    quat: [class label] of the quaternary blocks (equal labels = equal G; distinct labels = non-isogenous J(Y_G));
    l: number of terms z^3. Returns dim Hdg^m, the part from determinant terms, dim H_prim and its (m,m)-part."""
    k, q = len(tern), len(quat)
    N = 3 * k + 4 * q + l
    assert N % 2 == 0
    classes = sorted(set(quat))
    total, det, prim, hmm = 1, 0, 0, 0
    for chi in itertools.product(range(3), repeat=k):
        for psi in itertools.product(range(3), repeat=q):
            for c in itertools.product((1, 2), repeat=l):
                if (sum(chi) + sum(psi) + sum(c)) % 3:
                    continue
                n1, n2, q0 = chi.count(1), chi.count(2), psi.count(0)
                a, b, l1, l2 = psi.count(1), psi.count(2), c.count(1), c.count(2)
                K = [s for s in range(k) if chi[s] == 0]
                mult = 3 ** (n1 + n2) * 6 ** q0
                prim += mult * 2 ** len(K) * 5 ** (a + b)
                inv, esum, anydet = 1, 0, False
                for cl in classes:
                    ac = sum(1 for i in range(q) if quat[i] == cl and psi[i] == 1)
                    bc = sum(1 for i in range(q) if quat[i] == cl and psi[i] == 2)
                    inv *= I(ac, bc)
                    esum += Fraction(ac - bc, 5)
                    anydet |= ac != bc
                Kw = [s for s in K if tern[s][0] == 'w']
                x = Fraction(len(Kw), 2) + Fraction(n1 - n2 + a - b + l1 - l2, 6) - Fraction(3, 2) * esum
                term = inv * binom(len(Kw), x)
                groups = {}
                for s in K:
                    if tern[s][0] != 'w':
                        groups.setdefault(tern[s][0], []).append(s)
                for ss in groups.values():
                    term *= nu(len(ss), tern[ss[0]][1])
                total += mult * term
                if anydet:                              # a determinant factor in some class
                    det += mult * term
    # the (m,m)-part of H_prim from the block table (excess 0), as a control
    tab = {n: block_pieces(n) for n in (1, 3, 4)}
    dist = {(0, Fraction(0)): 1}
    for n in [3] * k + [4] * q + [1] * l:
        new = {}
        for (cc, e), m in dist.items():
            for a_, d, c_, e_ in tab[n]:
                key = ((cc + c_) % 3, e + e_)
                new[key] = new.get(key, 0) + m * d
        dist = new
    hmm = dist.get((0, Fraction(0)), 0)
    return dict(hodge=total, det=det, dim_prim=prim, hmm_prim=hmm)


def hodge_counts(kw, others, qsizes, l):
    """The same sum over counts: kw ternary blocks isogenous to E_omega, others = [(size, has_cm)] further isogeny
    classes of ternary blocks, qsizes = [size] of the classes of quaternary blocks, l terms z^3."""
    def splits(n, parts):                      # (n_1, ..., n_parts) with sum n, with multinomial multiplicity
        if parts == 1:
            yield (n,)
            return
        for i in range(n + 1):
            for rest in splits(n - i, parts - 1):
                yield (i,) + rest

    tern_opts = []                              # per ternary class: list of (n1, n2, contribution factor, Kw)
    tern_opts.append([((n1, n2), multinom(n1, n2, K), K) for n1, n2, K in splits(kw, 3)])
    for size, cm in others:
        tern_opts.append([((n1, n2), multinom(n1, n2, K) * nu(K, cm), None) for n1, n2, K in splits(size, 3)])
    quat_opts = [[((q0, a, b), multinom(q0, a, b) * 6 ** q0 * I(a, b)) for q0, a, b in splits(m, 3)] for m in qsizes]
    total, det = 1, 0
    for l1 in range(l + 1):
        l2 = l - l1
        for topt in itertools.product(*tern_opts):
            n1 = sum(o[0][0] for o in topt)
            n2 = sum(o[0][1] for o in topt)
            f_t = 1
            for o in topt:
                f_t *= o[1]
            if not f_t:
                continue
            Kw = topt[0][2]
            for qopt in itertools.product(*quat_opts):
                f_q = 1
                for o in qopt:
                    f_q *= o[1]
                if not f_q:
                    continue
                a = sum(o[0][1] for o in qopt)
                b = sum(o[0][2] for o in qopt)
                if (n1 + 2 * n2 + a + 2 * b + l1 + 2 * l2) % 3:
                    continue
                esum = sum(Fraction(o[0][1] - o[0][2], 5) for o in qopt)
                x = Fraction(Kw, 2) + Fraction(n1 - n2 + a - b + l1 - l2, 6) - Fraction(3, 2) * esum
                term = comb(l, l1) * 3 ** (n1 + n2) * f_t * f_q * binom(Kw, x)
                total += term
                if any(o[0][1] != o[0][2] for o in qopt):
                    det += term
    return dict(hodge=total, det=det)


def hodge_weights(tern, quat, l):
    """Independent count: Hg-invariants of V(chi)/M from explicit weights, for the cover
    prod_c SU(5) x G_m(s) x prod_kappa Hg(E_kappa) -> Hg: W_c has weights e_k s^3, W-bar_c weights -e_k s^-3,
    H^1 of a curve isogenous to E_omega has s^5 and s^-5, the Tate part has s^(10 e_T) (t = s^5 acts on it by its
    Hodge weight t^(2 e_T)), a CM class kappa has its own torus u^(+-1), a class without CM an SL_2 (weights +-1)."""
    k, q = len(tern), len(quat)
    classes = sorted(set(quat))
    tab4 = block_pieces(4)
    mu = {}                                     # excess of the character mu_psi of W (x) mu, W-bar (x) mu-bar
    for psi in (1, 2):
        tot = sum(d * e for a, d, c, e in tab4 if c == psi)
        mu[psi] = (tot - (Fraction(3, 2) if psi == 1 else Fraction(-3, 2))) / 5
    groups = {}
    for s, (lab, cm) in enumerate(tern):
        if lab != 'w':
            groups.setdefault(lab, []).append(s)
    glist = sorted(groups)
    total = 1
    exc = {n: {c: e for a, d, c, e in block_pieces(n) if c} for n in (1, 3)}   # excess of the pieces with chi != 0
    for chi in itertools.product(range(3), repeat=k):
        for psi in itertools.product(range(3), repeat=q):
            for c in itertools.product((1, 2), repeat=l):
                if (sum(chi) + sum(psi) + sum(c)) % 3:
                    continue
                eT = sum(exc[3][x] for x in chi if x) + sum(exc[1][x] for x in c) + sum(mu[x] for x in psi if x)
                s0 = 10 * eT
                assert s0.denominator == 1
                mult = 3 ** (chi.count(1) + chi.count(2)) * 6 ** psi.count(0)
                # factors: per quaternary block with psi != 0 and per ternary block with chi = 0
                dist = {((((0,) * 5),) * len(classes), int(s0), (0,) * len(glist)): 1}
                for i in range(q):
                    if not psi[i]:
                        continue
                    ci = classes.index(quat[i])
                    sg = 1 if psi[i] == 1 else -1
                    new = {}
                    for (su, s, g), m in dist.items():
                        for e in E5:
                            su2 = list(su)
                            su2[ci] = tuple(x + sg * y for x, y in zip(su[ci], e))
                            key = (tuple(su2), s + 3 * sg, g)
                            new[key] = new.get(key, 0) + m
                    dist = new
                for s_ in range(k):
                    if chi[s_]:
                        continue
                    lab = tern[s_][0]
                    new = {}
                    for (su, s, g), m in dist.items():
                        for sg in (1, -1):
                            if lab == 'w':
                                key = (su, s + 5 * sg, g)
                            else:
                                gi = glist.index(lab)
                                g2 = list(g)
                                g2[gi] += sg
                                key = (su, s, tuple(g2))
                            new[key] = new.get(key, 0) + m
                    dist = new
                inv = 0
                for (su, s, g), m in dist.items():
                    if s:
                        continue
                    f = m
                    for w in su:
                        f *= A5.get(norm5(w), 0)
                        if not f:
                            break
                    if not f:
                        continue
                    for gi, lab in enumerate(glist):
                        if tern[groups[lab][0]][1]:         # CM: torus, weight 0 only
                            f *= 1 if g[gi] == 0 else 0
                        else:                               # SL_2: m(0) - m(2), via the sign of the weight
                            f *= 1 if g[gi] == 0 else (-1 if g[gi] == 2 else 0)
                    inv += f
                total += mult * inv
    return total


def to_counts(tern, quat):
    kw = sum(1 for lab, cm in tern if lab == 'w')
    groups = {}
    for lab, cm in tern:
        if lab != 'w':
            groups.setdefault(lab, [0, cm])[0] += 1
    others = [tuple(v) for v in groups.values()]
    qs = {}
    for c in quat:
        qs[c] = qs.get(c, 0) + 1
    return kw, others, sorted(qs.values())


W_ = ('w', True)
GEN = [('g%d' % i, False) for i in range(8)]
CMO = ('i', True)


def name(tern, quat, l):
    t = ''.join('w' if lab == 'w' else ('i' if lab == 'i' else 'E' + lab[1:]) for lab, cm in tern)
    return f"[{t}|{','.join(quat)}|l={l}]"


CASES = [
    ([], ['G'], 1), ([], ['G'], 2), ([], ['G'], 4), ([W_], ['G'], 1), ([GEN[0]], ['G'], 1), ([CMO], ['G'], 1),
    ([], ['G', 'G'], 0), ([], ['G', 'H'], 0), ([], ['G', 'G'], 2), ([W_], ['G', 'G'], 1), ([W_, W_], ['G'], 0),
    ([], ['G', 'G', 'G'], 0), ([], ['G', 'G', 'H'], 0), ([], ['G', 'H', 'K'], 0), ([], ['G', 'G'], 4),
    ([GEN[0], GEN[0]], ['G'], 0), ([GEN[0], GEN[1]], ['G'], 0), ([W_, GEN[0]], ['G'], 1),
    ([], ['G', 'G', 'G', 'G'], 0), ([], ['G', 'G', 'H', 'H'], 0), ([W_, W_], ['G', 'G'], 0),
    ([], ['G'] * 4, 2), ([], ['G'] * 5, 1), ([], ['G'] * 5, 4), ([W_], ['G'] * 5, 1), ([], ['G'] * 5 + ['H'] * 5, 0),
]


def prop_checks():
    rows, ok = [], True
    for tern, quat, l in CASES:
        N = 3 * len(tern) + 4 * len(quat) + l
        if N % 2:
            continue
        r = hodge_chi(tern, quat, l)
        rc = hodge_counts(*to_counts(tern, quat), l)
        heavy = len(quat) + len(tern) + l > 9
        rw = None if heavy else hodge_weights(tern, quat, l)
        good = (r['hodge'] == rc['hodge'] and r['det'] == rc['det'] and (rw is None or rw == r['hodge'])
                and r['dim_prim'] == (2 ** N + 2) // 3 and r['hmm_prim'] == comb(N, N // 2))
        ok &= good
        rows.append(dict(case=name(tern, quat, l), N=N, dim=N - 2, hodge=r['hodge'], det=r['det'],
                         counts=rc['hodge'], weights=rw, ok=good))
    # no quaternary blocks: agreement with Proposition C (t3_formula.hodge_dim) on all t2/t3 cases and more
    t3rows = 0
    for k in range(0, 5):
        for l in range(0, 9):
            N = 3 * k + l
            if N % 2 or N < 4 or N > 12:
                continue
            for labs in itertools.product([W_, CMO, GEN[0], GEN[1]], repeat=k):
                labs = list(labs)
                a = t3.hodge_dim(labs, l)['hodge']
                b = hodge_chi(labs, [], l)['hodge']
                c = hodge_counts(*to_counts(labs, []), l)['hodge']
                ok &= a == b == c
                t3rows += 1
    return ok, rows, t3rows


# ---------- (5) specialisation ----------
def special_checks():
    ok, rows = True, []
    for tern, quat, l in CASES:
        N = 3 * len(tern) + 4 * len(quat) + l
        if N % 2 or len(quat) + len(tern) + l > 11:
            continue
        gen = hodge_counts(*to_counts(tern, quat), l)['hodge']
        cls = sorted(set(quat))
        for mode in ('fermat', 'generic F'):
            rep = {c: (W_ if mode == 'fermat' else ('s%d' % i, False)) for i, c in enumerate(cls)}
            tern2 = list(tern) + [rep[c] for c in quat]
            spec = hodge_counts(*to_counts(tern2, []), l + len(quat))['hodge']
            good = spec >= gen
            ok &= good
            rows.append(dict(case=name(tern, quat, l), mode=mode, generic=gen, special=spec, ok=good))
    return ok, rows


# ---------- (6) first determinant classes ----------
def partitions(n, maxpart=None):
    maxpart = n if maxpart is None else maxpart
    if n == 0:
        yield []
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield [p] + rest


def det_threshold(Nmax=24):
    found = []
    for q in range(0, 7):
        for qs in partitions(q):
            for kw in range(0, 9):
                for kg in range(0, 9):
                    for l in range(0, 25):
                        N = 3 * (kw + kg) + 4 * q + l
                        if N % 2 or N > Nmax or N < 4:
                            continue
                        r = hodge_counts(kw, [(1, False)] * kg, qs, l)
                        if r['det']:
                            found.append(dict(N=N, dim=N - 2, quaternary_classes=qs, kw=kw, kg=kg, l=l,
                                              det=r['det'], hodge=r['hodge']))
    Nmin = min(f['N'] for f in found) if found else None
    return Nmin, [f for f in found if f['N'] == Nmin], len(found)


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    lines, out, allok = [], {}, True

    ok, rows = jacobian_checks()
    allok &= ok
    out['jacobian'] = rows
    lines.append(f'(1) Jacobian rings over F_10007: Hilbert functions (1+t)^n for n = 3, 4 and Y_G; sigma-eigenspace '
                 f'dimensions (chi 1, chi 2) of H^(2,1), H^(1,2): (4,1), (1,4) -- {"ok" if ok else "FAIL"}')
    for r in rows:
        lines.append('    ' + json.dumps(r))

    ok, tab = table_checks()
    allok &= ok
    out['block_table'] = tab
    lines.append(f'(2) block table (a, dim, character, excess); character-1 part of G = W (x) mu with W of signature '
                 f'(4,1), mu of excess -1/6; totals (2^N+2)/3 and binom(N,N/2) for N <= 16 -- {"ok" if ok else "FAIL"}')
    for n in (1, 3, 4):
        lines.append(f'    n = {n}: {tab[n]}')

    ok, rows = fft_checks()
    allok &= ok
    out['I_sl5'] = rows
    lines.append(f'(3) I(a,b) = dim (W^a (x) W*^b)^SL5 (Brauer-Klimyk): a! for a = b, 0 unless a = b mod 5, '
                 f'I(5,0) = 1, I(6,1) = 5 -- {"ok" if ok else "FAIL"}')
    lines.append('    ' + ', '.join(f'I({k}) = {v}' for k, v in rows.items() if v))

    ok, rows, nt3 = prop_checks()
    det_of = {r['case']: r['det'] for r in rows}
    ok &= det_of['[|G,G,G,G,G|l=4]'] == 2 and det_of['[|G,G,G,G,G,H,H,H,H,H|l=0]'] == 2
    ok &= all(r['det'] == 0 for r in rows if r['case'].count('G') <= 4)   # Theorem D: at most 4 copies
    allok &= ok
    out['prop_Cprime'] = rows
    lines.append(f"(4) Proposition C': closed formula = sum over counts = explicit-weight count (where run), "
                 f"dim H_prim = (2^N+2)/3, (m,m)-part = binom(N,N/2); and = Proposition C on {nt3} configurations "
                 f"without quaternary blocks -- {'ok' if ok else 'FAIL'}")
    lines.append('    [ternary blocks | quaternary classes | l]  (w: E_omega, i: CM by Z[i], E#: no CM; '
                 'equal letters = equal G)')
    for r in rows:
        lines.append(f"    {r['case']:<34} dim {r['dim']:>2}: Hdg {r['hodge']:>7}  determinant part {r['det']:>3}"
                     f"  (counts {r['counts']}, weights {r['weights']}) {'ok' if r['ok'] else 'FAIL'}")

    ok, rows = special_checks()
    allok &= ok
    out['specialisation'] = rows
    lines.append(f'(5) specialisation G -> F + u^3 never lowers dim Hdg ({len(rows)} comparisons) -- '
                 f'{"ok" if ok else "FAIL"}')
    for r in rows:
        lines.append(f"    {r['case']:<34} {r['mode']:<10} generic {r['generic']:>7} <= special {r['special']:>7}")

    Nmin, mins, nfound = det_threshold()
    out['det_threshold'] = dict(Nmin=Nmin, minimal=mins, configurations_with_det=nfound)
    good = Nmin == 24 and any(m['quaternary_classes'] == [5] and m['l'] == 4 and m['kw'] == 0 and m['kg'] == 0
                              and m['det'] == 2 for m in mins)
    allok &= good
    lines.append(f'(6) determinant (Weil-type) Hodge classes in a single hypersurface, all configurations with '
                 f'N <= 24 (quaternary classes, kw blocks ~ E_omega, kg pairwise non-isogenous blocks without CM, '
                 f'l cubes): first at N = {Nmin} (dimension {Nmin - 2 if Nmin else None}); '
                 f'{nfound} configurations have them -- {"ok" if good else "FAIL"}')
    for m in mins:
        lines.append(f"    classes {m['quaternary_classes']}, kw {m['kw']}, kg {m['kg']}, l {m['l']}: "
                     f"dim Hdg {m['hodge']}, determinant part {m['det']}")

    lines.append('ALL CHECKS PASS' if allok else 'SOME CHECK FAILED')
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 't5_quaternary.txt'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
    with open(os.path.join(OUT, 't5_quaternary.json'), 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(out, fh, indent=1, default=str)
    print('\n'.join(lines))
