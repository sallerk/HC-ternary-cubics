# t7_eckardt.py -- checks for NOTE.md Section 11 (Phase 101, 27 Oct 2026): local Torelli for the cyclic cubic
# threefolds Y_G = {G + y^3 = 0}, Eckardt blocks G = C(x1,x2,x3) + w^2 L(x1,x2,x3), and dim Hdg^m with Eckardt blocks.
# (1) G generic, over F_10007: the pairing (R_G)_3 x (R_G)_1 -> (R_G)_4 = F_p has rank 4, so the differential of the
#     period map of {Y_G} to the 4-ball (multiplication by B in (R_G)_3, from W^{2,1} = (R_G)_1 to W^{1,2} = (R_G)_4)
#     is bijective.
# (2) G = C + w^2 L with random C, L: Hilbert function 1,4,6,4,1,0 (smooth); two singular controls (L tangent to C at a
#     point; C nodal) have R_5 != 0; the tau-graded Hilbert functions (tau: w -> -w); the residue characters (tau acts on
#     the class of A Omega/(G + y^3)^(p+1) by det(tau) (A o tau) = -(A o tau)) give W = W^- (+) W^+ of signatures (3,1)
#     and (1,0); the pairing (R_G)_3^+ x (R_G)_1^+ -> (R_G)_4 has rank 3 (local Torelli for the Eckardt family on W^-);
#     the odd pairing has rank 1; the Fermat surface is of this form (sympy).
# (3) dim Hdg^m with Eckardt blocks (Proposition C''), two ways: a closed formula per character (each Eckardt factor
#     W = W^- (+) L, SL_4 invariants I_4 by Brauer-Klimyk) and explicit weights of the cover
#     SU(5)^{(H1) classes} x SU(4)^{Eckardt classes} x G_m(s) -> Hg with t = s^10 (W_c: s^6, W^-_e: s^5, L: s^10, Tate
#     part s^(20 e_T)). With (H1) blocks only, both agree with t5; generic <= Eckardt <= Fermat specialisation.
# (4) Where Eckardt determinant classes first occur in a single hypersurface: all configurations with N <= 20 that
#     contain four blocks of one Eckardt class (fewer cannot carry them, as I_4(a, b) = 0 unless a = b mod 4).
import importlib.util, itertools, json, os, random, sys
from fractions import Fraction
from math import comb, factorial, prod

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')
P = 10007


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


t1 = _load('t1_sk_smooth')
t5 = _load('t5_quaternary')


# ---------- linear algebra over F_p ----------
def rref(rows, ncols, p):
    M = [[x % p for x in r] for r in rows]
    piv, r = [], 0
    for c in range(ncols):
        k = next((i for i in range(r, len(M)) if M[i][c]), None)
        if k is None:
            continue
        M[r], M[k] = M[k], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == len(M):
            break
    return M[:r], piv


def rank(rows, ncols, p):
    return len(rref(rows, ncols, p)[1]) if rows else 0


def kernel(rows, ncols, p):
    """basis of {v : r.v = 0 for every row r}."""
    R, piv = rref(rows, ncols, p)
    out = []
    for f in (c for c in range(ncols) if c not in piv):
        v = [0] * ncols
        v[f] = 1
        for r, c in zip(R, piv):
            v[c] = -r[f] % p
        out.append(v)
    return out


# ---------- (1), (2) Jacobian rings ----------
def mons(n, d):
    return [e for e in itertools.product(range(d + 1), repeat=n) if sum(e) == d]


def jac_rows(G, n, d, p, parity=None, wv=3):
    """J_d = sum_i S_{d-2} dG/du_i in the monomial basis of S_d. With parity = 0 or 1, only the part of that parity in
    u_wv (for G even in u_wv: dG/du_i has parity [i == wv])."""
    cols = [e for e in mons(n, d) if parity is None or e[wv] % 2 == parity]
    idx = {e: i for i, e in enumerate(cols)}
    rows = []
    if d >= 2:
        D = t1.partials(G, n)
        for m in mons(n, d - 2):
            for i, Q in enumerate(D):
                if parity is not None and (m[wv] + (i == wv)) % 2 != parity:
                    continue
                v = [0] * len(cols)
                for e, c in Q.items():
                    j = idx[tuple(a + b for a, b in zip(e, m))]
                    v[j] = (v[j] + c) % p
                rows.append(v)
    return cols, rows


def hilbert(G, n, p, dmax, parity=None):
    out = []
    for d in range(dmax + 1):
        cols, rows = jac_rows(G, n, d, p, parity)
        out.append(len(cols) - rank(rows, len(cols), p))
    return out


def pairing_rank(G, n, p, lin, cub):
    """rank of (ell(a b)) for a in lin, b in cub, where ell spans the annihilator of J_4 (R_4 = F_p)."""
    cols, rows = jac_rows(G, n, 4, p)
    ker = kernel(rows, len(cols), p)
    assert len(ker) == 1
    ell = dict(zip(cols, ker[0]))
    M = [[ell[tuple(x + y for x, y in zip(a, b))] for b in cub] for a in lin]
    return rank(M, len(cub), p)


def eckardt(C, L):
    """G(x1, x2, x3, w) = C(x) + w^2 L(x), for a ternary cubic C and a linear form L (dicts over exponents)."""
    G = {e + (0,): c for e, c in C.items()}
    for e, c in L.items():
        G[e + (2,)] = G.get(e + (2,), 0) + c
    return G


def generic_checks(p=P, seed=101, trials=3):
    rnd, ok, rows = random.Random(seed), True, []
    for _ in range(trials):
        G = {e: rnd.randrange(1, p) for e in mons(4, 3)}
        hf = hilbert(G, 4, p, 5)
        r = pairing_rank(G, 4, p, mons(4, 1), mons(4, 3))
        good = hf == [1, 4, 6, 4, 1, 0] and r == 4
        ok &= good
        rows.append(dict(hilbert=hf, pairing_rank=r, ok=good))
    return ok, rows


def eckardt_checks(p=P, seed=102, trials=3):
    rnd, ok, rows = random.Random(seed), True, []
    lin_even = [e for e in mons(4, 1) if e[3] % 2 == 0]
    lin_odd = [e for e in mons(4, 1) if e[3] % 2]
    cub_even = [e for e in mons(4, 3) if e[3] % 2 == 0]
    cub_odd = [e for e in mons(4, 3) if e[3] % 2]
    for _ in range(trials):
        C = t1.random_ternary(p, rnd)
        L = {e: rnd.randrange(1, p) for e in mons(3, 1)}
        G = eckardt(C, L)
        hf = hilbert(G, 4, p, 5)
        he = hilbert(G, 4, p, 4, parity=0)
        ho = hilbert(G, 4, p, 4, parity=1)
        # W^{2,1} = (R_G)_1, W^{1,2} = (R_G)_4 (omega-part); tau = -(parity sign): W^- <-> w-even, W^+ <-> w-odd
        wm, wp = (he[1], he[4]), (ho[1], ho[4])
        r_even = pairing_rank(G, 4, p, lin_even, cub_even)
        r_odd = pairing_rank(G, 4, p, lin_odd, cub_odd)
        good = (hf == [1, 4, 6, 4, 1, 0] and he == [1, 3, 4, 3, 1] and ho == [0, 1, 2, 1, 0]
                and wm == (3, 1) and wp == (1, 0) and r_even == 3 and r_odd == 1 and len(cub_even) == 13)
        ok &= good
        rows.append(dict(hilbert=hf, hilbert_w_even=he, hilbert_w_odd=ho, W_minus=wm, W_plus=wp,
                         pairing_rank_even=r_even, pairing_rank_odd=r_odd, ok=good))
    # singular controls: (a) L tangent to C at (1:0:0); (b) C nodal at (0:0:1)
    while True:
        C = {e: rnd.randrange(1, p) for e in mons(3, 3)}
        C[(3, 0, 0)] = 0
        if t1.smooth_plane_cubic(C, p):
            break
    L = {(0, 1, 0): C[(2, 1, 0)], (0, 0, 1): C[(2, 0, 1)]}
    hf_a = hilbert(eckardt(C, L), 4, p, 5)
    Cn = {(1, 1, 1): 1, (3, 0, 0): 1, (0, 3, 0): 1}
    hf_b = hilbert(eckardt(Cn, {e: rnd.randrange(1, p) for e in mons(3, 1)}), 4, p, 5)
    ctrl = hf_a[5] != 0 and hf_b[5] != 0
    ok &= ctrl
    rows.append(dict(control_tangent=hf_a, control_nodal=hf_b, singular_detected=ctrl))
    # the Fermat surface: u3 = (s + w)/2, u4 = (s - w)/2 gives u1^3 + u2^3 + s^3/4 + (3/4) s w^2
    u1, u2, s, w = sp.symbols('u1 u2 s w')
    fer = sp.expand(u1**3 + u2**3 + ((s + w) / 2)**3 + ((s - w) / 2)**3
                    - (u1**3 + u2**3 + s**3 / 4 + sp.Rational(3, 4) * s * w**2)) == 0
    ok &= fer
    rows.append(dict(fermat_is_C_plus_w2L=fer))
    return ok, rows


# ---------- (3) SL_4 invariants and Proposition C'' ----------
RHO4 = (3, 2, 1, 0)


def norm4(w):
    return tuple(x - w[3] for x in w)


A4 = {}
for _p in itertools.permutations(range(4)):
    A4[norm4(tuple(RHO4[_p[i]] - RHO4[i] for i in range(4)))] = t5._sign(_p)   # w(rho) - rho -> sgn(w)
assert len(A4) == 24
E4 = [tuple(1 if j == i else 0 for j in range(4)) for i in range(4)]


def I4_raw(a, b):
    """dim (W^{(x)a} (x) W*^{(x)b})^{SL_4} (Brauer-Klimyk)."""
    dist = {(0,) * 4: 1}
    for sg in [1] * a + [-1] * b:
        new = {}
        for w, m in dist.items():
            for e in E4:
                v = tuple(x + sg * y for x, y in zip(w, e))
                new[v] = new.get(v, 0) + m
        dist = new
    return sum(m * A4.get(norm4(w), 0) for w, m in dist.items())


_I4 = {}


def I4(a, b):
    if (a, b) not in _I4:
        _I4[(a, b)] = I4_raw(a, b) if (a - b) % 4 == 0 else 0
    return _I4[(a, b)]


def fft4_checks():
    rows, ok = {}, True
    for a in range(7):
        for b in range(7):
            if a + b > 8:
                continue
            v = I4_raw(a, b)
            rows[f'{a},{b}'] = v
            if a + b <= 3:
                ok &= v == (factorial(a) if a == b else 0)
            if (a - b) % 4:
                ok &= v == 0
    ok &= rows['4,0'] == 1 and rows['0,4'] == 1 and rows['4,4'] == 24 and rows['5,1'] == 4 and rows['1,5'] == 4
    return ok, rows


def _setup(tern, quat):
    gen = sorted({lab for lab, kd in quat if kd == 'gen'})
    eck = sorted({lab for lab, kd in quat if kd == 'eck'})
    assert not set(gen) & set(eck)
    groups = {}
    for s, (lab, cm) in enumerate(tern):
        if lab != 'w':
            groups.setdefault(lab, []).append(s)
    return gen, eck, groups


def _characters(k, q, l):
    for chi in itertools.product(range(3), repeat=k):
        for psi in itertools.product(range(3), repeat=q):
            for c in itertools.product((1, 2), repeat=l):
                if (sum(chi) + sum(psi) + sum(c)) % 3 == 0:
                    yield chi, psi, c


def hodge_E_chi(tern, quat, l):
    """Proposition C'' summed over characters. tern: [(label, has_cm)] as in t3/t5; quat: [(label, kind)] with kind
    'gen' ((H1)) or 'eck' (Eckardt with Hg(B_G) = U(W^-)); equal labels = equal G, distinct labels = non-isogenous;
    l terms z^3. Returns dim Hdg^m and the parts with a determinant factor of an (H1) or an Eckardt class."""
    k, q = len(tern), len(quat)
    assert (3 * k + 4 * q + l) % 2 == 0
    gen, eck, groups = _setup(tern, quat)
    total, det_g, det_e = 1, 0, 0
    for chi, psi, c in _characters(k, q, l):
        n1, n2, q0 = chi.count(1), chi.count(2), psi.count(0)
        a, b, l1, l2 = psi.count(1), psi.count(2), c.count(1), c.count(2)
        eT = Fraction(-n1 + n2 - a + b - l1 + l2, 6)          # excess of lambda'; Hg acts on it by t^(2 eT)
        Kw = sum(1 for s in range(k) if chi[s] == 0 and tern[s][0] == 'w')
        base = 3 ** (n1 + n2) * 6 ** q0
        for ss in groups.values():
            base *= t5.nu(sum(1 for s in ss if chi[s] == 0), tern[ss[0]][1])
        if not base:
            continue
        inv_g, tg, anyg = 1, Fraction(0), False
        for cl in gen:
            ac = sum(1 for i in range(q) if quat[i][0] == cl and psi[i] == 1)
            bc = sum(1 for i in range(q) if quat[i][0] == cl and psi[i] == 2)
            inv_g *= t5.I(ac, bc)
            tg += Fraction(3 * (ac - bc), 5)                    # a determinant of W_c has t-weight 3
            anyg |= ac != bc
        if not inv_g:
            continue
        opts = []
        for cl in eck:
            ae = sum(1 for i in range(q) if quat[i][0] == cl and psi[i] == 1)
            be = sum(1 for i in range(q) if quat[i][0] == cl and psi[i] == 2)
            o = []
            for x1 in range(ae + 1):          # x1 of the psi = 1 factors use L, x2 of the psi = 2 factors use L-bar
                for x2 in range(be + 1):
                    f = comb(ae, x1) * comb(be, x2) * I4(ae - x1, be - x2)
                    if f:                     # a determinant of W^- has t-weight 2, L has t-weight 1
                        o.append((f, Fraction(ae - x1 - be + x2, 2) + x1 - x2, ae - x1 != be - x2))
            opts.append(o)
        for combo in itertools.product(*opts):
            te = sum((o[1] for o in combo), Fraction(0))
            j = Fraction(Kw, 2) - eT - tg / 2 - te / 2          # number of E_omega factors in K that take L
            term = base * inv_g * prod(o[0] for o in combo) * t5.binom(Kw, j)
            total += term
            det_g += term if anyg else 0
            det_e += term if any(o[2] for o in combo) else 0
    return dict(hodge=total, det_gen=det_g, det_eck=det_e)


def hodge_E_weights(tern, quat, l):
    """Independent count: Hg-invariants from explicit weights of the cover
    prod_{(H1) classes} SU(5) x prod_{Eckardt classes} SU(4) x G_m(s) x prod_kappa Hg(E_kappa) -> Hg, t = s^10:
    W_c has e_k s^6, W^-_e has f_k s^5, the line L = W^+_e and H^{1,0} of a curve ~ E_omega have s^10 (conjugates:
    negatives), the Tate part s^(20 e_T); a CM class kappa has its own torus, a class without CM an SL_2."""
    k, q = len(tern), len(quat)
    gen, eck, groups = _setup(tern, quat)
    glist = sorted(groups)
    total = 1
    for chi, psi, c in _characters(k, q, l):
        n1, n2, q0 = chi.count(1), chi.count(2), psi.count(0)
        eT = Fraction(-n1 + n2 - psi.count(1) + psi.count(2) - c.count(1) + c.count(2), 6)
        s0 = 20 * eT
        assert s0.denominator == 1
        mult = 3 ** (n1 + n2) * 6 ** q0
        dist = {((((0,) * 5),) * len(gen), (((0,) * 4),) * len(eck), int(s0), (0,) * len(glist)): 1}
        for i in range(q):
            if not psi[i]:
                continue
            sg = 1 if psi[i] == 1 else -1
            lab, kd = quat[i]
            new = {}
            for (u5, u4, s, g), m in dist.items():
                if kd == 'gen':
                    ci = gen.index(lab)
                    for e in t5.E5:
                        v = list(u5)
                        v[ci] = tuple(x + sg * y for x, y in zip(u5[ci], e))
                        key = (tuple(v), u4, s + 6 * sg, g)
                        new[key] = new.get(key, 0) + m
                else:
                    ci = eck.index(lab)
                    for e in E4:
                        v = list(u4)
                        v[ci] = tuple(x + sg * y for x, y in zip(u4[ci], e))
                        key = (u5, tuple(v), s + 5 * sg, g)
                        new[key] = new.get(key, 0) + m
                    key = (u5, u4, s + 10 * sg, g)
                    new[key] = new.get(key, 0) + m
            dist = new
        for s_ in range(k):
            if chi[s_]:
                continue
            lab = tern[s_][0]
            new = {}
            for (u5, u4, s, g), m in dist.items():
                for sg in (1, -1):
                    if lab == 'w':
                        key = (u5, u4, s + 10 * sg, g)
                    else:
                        g2 = list(g)
                        g2[glist.index(lab)] += sg
                        key = (u5, u4, s, tuple(g2))
                    new[key] = new.get(key, 0) + m
            dist = new
        inv = 0
        for (u5, u4, s, g), m in dist.items():
            if s:
                continue
            f = m
            for w in u5:
                f *= t5.A5.get(t5.norm5(w), 0)
            for w in u4:
                f *= A4.get(norm4(w), 0)
            if not f:
                continue
            for gi, lab in enumerate(glist):
                if tern[groups[lab][0]][1]:                     # CM: torus, weight 0 only
                    f *= 1 if g[gi] == 0 else 0
                else:                                           # SL_2: m(0) - m(2)
                    f *= 1 if g[gi] == 0 else (-1 if g[gi] == 2 else 0)
            inv += f
        total += mult * inv
    return total


def kinds(quat, kind_of):
    return [(c, kind_of(c)) for c in quat]


def prop_checks():
    rows, ok = [], True
    for tern, quat, l in t5.CASES:
        N = 3 * len(tern) + 4 * len(quat) + l
        if N % 2:
            continue
        nb = len(tern) + len(quat) + l
        name = t5.name(tern, quat, l)
        # regression: only (H1) blocks
        rg = hodge_E_chi(tern, kinds(quat, lambda c: 'gen'), l)
        r5 = t5.hodge_chi(tern, quat, l)
        good = rg['hodge'] == r5['hodge'] and rg['det_gen'] == r5['det'] and rg['det_eck'] == 0
        if nb <= 9:
            good &= hodge_E_weights(tern, kinds(quat, lambda c: 'gen'), l) == r5['hodge']
        row = dict(case=name, N=N, dim=N - 2, generic=r5['hodge'], regression_ok=good)
        if len(quat) <= 5 and nb <= 10:
            gen = t5.hodge_counts(*t5.to_counts(tern, quat), l)['hodge']
            re = hodge_E_chi(tern, kinds(quat, lambda c: 'eck'), l)
            we = hodge_E_weights(tern, kinds(quat, lambda c: 'eck'), l) if nb <= 9 else None
            tern2 = list(tern) + [t5.W_] * len(quat)
            fer = t5.hodge_counts(*t5.to_counts(tern2, []), l + len(quat))['hodge']
            g2 = gen <= re['hodge'] <= fer and (we is None or we == re['hodge'])
            cls = sorted(set(quat))
            mixed = None
            if len(cls) >= 2:                   # first class Eckardt, the others (H1)
                mixed = hodge_E_chi(tern, kinds(quat, lambda c: 'eck' if c == cls[0] else 'gen'), l)['hodge']
                g2 &= gen <= mixed <= re['hodge']
            good &= g2
            row.update(eckardt=re['hodge'], eckardt_weights=we, eckardt_det=re['det_eck'], mixed=mixed, fermat=fer)
        row['ok'] = good
        ok &= good
        rows.append(row)
    return ok, rows


# ---------- (4) first Eckardt determinant classes ----------
def eck_det_threshold(Nmax=20):
    """All configurations with N <= Nmax containing four blocks of one Eckardt class e0, with the remaining <= 4
    variables in: cubes, one more quaternary block (class e0, another Eckardt class, or (H1)), or one ternary block
    (~ E_omega, CM by Z[i], or without CM) and one cube."""
    rest = [('l=0', [], [], 0), ('l=2', [], [], 2), ('l=4', [], [], 4),
            ('+e0', [], [('e0', 'eck')], 0), ('+e1', [], [('e1', 'eck')], 0), ('+g', [], [('g', 'gen')], 0),
            ('+F_w,l=1', [t5.W_], [], 1), ('+F_i,l=1', [t5.CMO], [], 1), ('+F_gen,l=1', [t5.GEN[0]], [], 1)]
    found = []
    for lab, tern, extra, l in rest:
        quat = [('e0', 'eck')] * 4 + extra
        N = 3 * len(tern) + 4 * len(quat) + l
        if N % 2 or N > Nmax:
            continue
        r = hodge_E_chi(tern, quat, l)
        found.append(dict(config='4 e0 ' + lab, N=N, dim=N - 2, hodge=r['hodge'], det_eck=r['det_eck']))
    with_det = [f for f in found if f['det_eck']]
    Nmin = min(f['N'] for f in with_det) if with_det else None
    return Nmin, found


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    lines, out, allok = [], {}, True

    ok, rows = generic_checks()
    allok &= ok
    out['generic'] = rows
    lines.append(f'(1) G generic over F_{P}: Hilbert function 1,4,6,4,1,0 and the pairing (R_G)_3 x (R_G)_1 -> (R_G)_4 '
                 f'of rank 4 (local Torelli for Y_G on the 4-ball) -- {"ok" if ok else "FAIL"}')
    for r in rows:
        lines.append('    ' + json.dumps(r))

    ok, rows = eckardt_checks()
    allok &= ok
    out['eckardt'] = rows
    lines.append(f'(2) G = C + w^2 L over F_{P}: smooth; R_G by w-parity 1,3,4,3,1 (even) and 0,1,2,1,0 (odd); '
                 f'W^- (3,1), W^+ (1,0); pairing ranks 3 (even, local Torelli on the Eckardt family) and 1 (odd); '
                 f'singular controls detected; Fermat surface of this form -- {"ok" if ok else "FAIL"}')
    for r in rows:
        lines.append('    ' + json.dumps(r))

    ok, rows4 = fft4_checks()
    allok &= ok
    out['I_sl4'] = rows4
    lines.append(f'(3a) I_4(a,b) = dim (W^a (x) W*^b)^SL4 (Brauer-Klimyk): a! for a = b, a + b <= 3; 0 unless a = b '
                 f'mod 4; I_4(4,0) = 1, I_4(4,4) = 24, I_4(5,1) = 4 -- {"ok" if ok else "FAIL"}')
    lines.append('    ' + ', '.join(f'I4({k}) = {v}' for k, v in rows4.items() if v))

    ok, rows = prop_checks()
    val = {r['case']: r for r in rows}
    ok &= val['[|G|l=2]']['generic'] == 13 and val['[|G|l=2]']['eckardt'] == 15
    allok &= ok
    out['prop_Csecond'] = rows
    lines.append(f"(3b) Proposition C'': (H1) blocks only = t5 (formula and weights); all quaternary classes Eckardt: "
                 f"formula = weights (where run), generic <= Eckardt <= Fermat, generic <= mixed <= Eckardt; "
                 f"G + z1^3 + z2^3: 13 generic, 15 Eckardt -- {'ok' if ok else 'FAIL'}")
    lines.append('    [ternary blocks | quaternary classes | l]: generic / mixed (first class Eckardt) / Eckardt '
                 '(Eckardt det part) / Fermat')
    for r in rows:
        if 'eckardt' in r:
            lines.append(f"    {r['case']:<34} dim {r['dim']:>2}: {r['generic']:>7} / {str(r['mixed']):>7} / "
                         f"{r['eckardt']:>7} ({r['eckardt_det']}) / {r['fermat']:>7}  weights {r['eckardt_weights']}"
                         f"  {'ok' if r['ok'] else 'FAIL'}")
        else:
            lines.append(f"    {r['case']:<34} dim {r['dim']:>2}: {r['generic']:>9} (regression only)"
                         f"  {'ok' if r['ok'] else 'FAIL'}")

    Nmin, found = eck_det_threshold()
    out['eck_det_threshold'] = dict(Nmin=Nmin, configurations=found)
    good = Nmin == 18 and any(f['config'] == '4 e0 l=2' and f['det_eck'] == 2 for f in found)
    allok &= good
    lines.append(f'(4) Eckardt determinant classes in a single hypersurface (configurations with four blocks of one '
                 f'Eckardt class, N <= 20): first at N = {Nmin} (dimension {Nmin - 2 if Nmin else None}) -- '
                 f'{"ok" if good else "FAIL"}')
    for f in found:
        lines.append(f"    {f['config']:<14} N = {f['N']:>2}: dim Hdg {f['hodge']:>7}, Eckardt determinant part "
                     f"{f['det_eck']}")

    lines.append('ALL CHECKS PASS' if allok else 'SOME CHECK FAILED')
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 't7_eckardt.txt'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
    with open(os.path.join(OUT, 't7_eckardt.json'), 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(out, fh, indent=1, default=str)
    print('\n'.join(lines))
