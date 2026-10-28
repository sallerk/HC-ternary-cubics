# t6_joins.py -- NOTE.md Section 6 (Phase 99, 27 Oct 2026): joins of planes in the cubic tenfolds F_1+F_2+F_3+F_4 = 0.
# (1) Symbolic (SymPy, exact in Q(omega)[lambda], omega^2 + omega + 1 = 0): for the Hesse cubic
#     F = x^3 + y^3 + z^3 - 3 lambda x y z, the 9 flexes f = g f_0 and tangent vectors v(f) = g v_0, g = phi^i tau^j
#     (phi: (x,y,z) -> (y,z,x), tau: (x,y,z) -> (x, omega y, omega^2 z), f_0 = (0,1,-1), v_0 = (1,-lambda,0)):
#     F(a f + t v(f)) = (1 - lambda^3) t^3; the flexes are distinct; for f != g the tangent lines meet in
#     r = x f + v(f) = x' g + kappa v(g) with kappa = omega^K(f,g), K(f,g) in Z/3 independent of lambda.
# (2) Planes of F_a(x_a) + F_b(x_b) = 0 of rank 2 (Koike): S = <(f,0), (0,g), (v(f), -rho v(g))>, rho = omega^k mu_a/mu_b,
#     mu^3 = 1 - lambda^3: 9 * 9 * 3 = 243 planes. Joins S (+) S' over the three pairings of {1,2,3,4}: 3 * 243^2 = 177147
#     linear P^5's in the tenfold.
# (3) Intersection dimensions: the model dim(A n B) = e + dim(L_A n kappa^-1 L_B) (e = number of blocks with the same flex,
#     L = the plane's equations in the coordinates c_i of v(f_i)) against direct ranks over F_p with explicit Hesse
#     forms, for two primes; and the cubic vanishes on every join tested.
# (4) Gram matrices: linear m-spaces meeting in P^r in a cubic 2m-fold meet with multiplicity (1 - (-2)^(r+1))/3, so with
#     h^m (degree 3) the primitive Gram matrix is 3G - 1 = -((-2)^dim(A n B)). Its rank over F_q is a lower bound for the
#     rank over Q. Planes in the fourfold: dim Hdg^2_prim = 18 for non-isogenous curves (Koike Thm 1.7; Prop C).
#     Joins in the tenfold: dim Hdg^5_prim = 486 for pairwise non-isogenous curves not isogenous to E_omega (Prop C).
import itertools, json, os, random, sys
import sympy as sp
import flint

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')
PAIRINGS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))


# ---------- (1) symbolic ----------
lam, a_, t_, w = sp.symbols('lambda a t w')


def red(e):
    return sp.rem(sp.expand(e), w ** 2 + w + 1, w)


def Fh(v):
    return v[0] ** 3 + v[1] ** 3 + v[2] ** 3 - 3 * lam * v[0] * v[1] * v[2]


def phi(v):
    return (v[1], v[2], v[0])


def tau(v):
    return (v[0], w * v[1], w ** 2 * v[2])


def gamma(i, j, v):
    for _ in range(j):
        v = tau(v)
    for _ in range(i):
        v = phi(v)
    return tuple(red(c) for c in v)


FLEX = [(i, j) for i in range(3) for j in range(3)]
f0, v0 = (0, 1, -1), (1, -lam, 0)


def symbolic_checks():
    fs = {ij: gamma(*ij, f0) for ij in FLEX}
    vs = {ij: gamma(*ij, v0) for ij in FLEX}
    ok = True
    for ij in FLEX:
        ok &= red(Fh(tuple(a_ * x + t_ * y for x, y in zip(fs[ij], vs[ij]))) - (1 - lam ** 3) * t_ ** 3) == 0
    for p, q in itertools.combinations(FLEX, 2):
        cr = sp.Matrix(fs[p]).cross(sp.Matrix(fs[q]))
        ok &= any(red(c) != 0 for c in cr)                       # distinct points of P^2
    K = {}
    for p in FLEX:
        for q in FLEX:
            if p == q:
                continue
            ok &= red(sp.Matrix([fs[p], fs[q], vs[p]]).det()) != 0  # q is not on the tangent line at p
            sols = [e for e in range(3)
                    if red(sp.Matrix([fs[p], fs[q], [x - w ** e * y for x, y in zip(vs[p], vs[q])]]).det()) == 0]
            ok &= len(sols) == 1
            K[(p, q)] = sols[0] if sols else None
    return ok, fs, vs, K


# ---------- (2)-(3) over F_p ----------
def find_prime(start):
    p = start
    while not (sp.isprime(p) and p % 3 == 1):
        p += 1
    return p


def setup_field(p, seed):
    rnd = random.Random(seed)
    while True:
        om = pow(rnd.randrange(2, p - 1), (p - 1) // 3, p)
        if om != 1:
            break
    lams, mus = [], []
    while len(lams) < 4:
        L = rnd.randrange(2, p - 1)
        c = (1 - pow(L, 3, p)) % p
        if c == 0 or pow(c, (p - 1) // 3, p) != 1:
            continue
        m = int(sp.ntheory.residue_ntheory.nthroot_mod(c, 3, p))
        assert pow(m, 3, p) == c
        lams.append(L)
        mus.append(m)
    return om, lams, mus


def numeric_flexes(p, om, L):
    """f_ij, v_ij over F_p for the Hesse parameter L (same group elements as the symbolic part)."""
    def ev(c):
        e = sp.Poly(sp.expand(c), lam, w)
        return sum(int(co) * pow(L, i, p) * pow(om, j, p) for (i, j), co in e.terms()) % p
    fs, vs = {}, {}
    for ij in FLEX:
        fs[ij] = [ev(c) for c in gamma(*ij, f0)]
        vs[ij] = [ev(c) for c in gamma(*ij, v0)]
    return fs, vs


class Field:
    def __init__(self, p, seed):
        self.p = p
        self.om, self.lams, self.mus = setup_field(p, seed)
        self.fs, self.vs = [], []
        for L in self.lams:
            f, v = numeric_flexes(p, self.om, L)
            self.fs.append(f)
            self.vs.append(v)
        self.inv = lambda x: pow(x, p - 2, p)

    def rho(self, a, b, k):
        return self.mus[a] * self.inv(self.mus[b]) * pow(self.om, k, self.p) % self.p

    def join(self, J):
        """basis (6 vectors in F_p^12) of the join J = (pairing index, {block: flex}, {pair: label})."""
        pi, flex, lab = J
        p, rows = self.p, []
        for i in range(4):
            v = [0] * 12
            v[3 * i:3 * i + 3] = self.fs[i][flex[i]]
            rows.append(v)
        for (a, b) in PAIRINGS[pi]:
            r = self.rho(a, b, lab[(a, b)])
            v = [0] * 12
            v[3 * a:3 * a + 3] = self.vs[a][flex[a]]
            v[3 * b:3 * b + 3] = [(-r * x) % p for x in self.vs[b][flex[b]]]
            rows.append(v)
        return rows

    def h_vanishes(self, rows, rnd, trials=4):
        p = self.p
        for _ in range(trials):
            cs = [rnd.randrange(p) for _ in rows]
            x = [sum(c * r[j] for c, r in zip(cs, rows)) % p for j in range(12)]
            val = 0
            for i in range(4):
                X, Y, Z = x[3 * i:3 * i + 3]
                val += X ** 3 + Y ** 3 + Z ** 3 - 3 * self.lams[i] * X * Y * Z
            if val % p:
                return False
        return True

    def direct_dim(self, A, B, nrows=6, ncols=12):
        M = flint.nmod_mat(2 * nrows, ncols, self.p)
        for r, row in enumerate(A + B):
            for c, x in enumerate(row):
                if x:
                    M[r, c] = x
        return 2 * nrows - M.rank()


def random_join(rnd):
    pi = rnd.randrange(3)
    flex = {i: FLEX[rnd.randrange(9)] for i in range(4)}
    lab = {pr: rnd.randrange(3) for pr in PAIRINGS[pi]}
    return (pi, flex, lab)


class Model:
    """dim(A n B) = e + dim(L_A n kappa^-1 L_B), the second term by a 4 x 4 rank over F_p (the mu's cancel)."""
    def __init__(self, K, fld):
        self.K, self.f, self.memo = K, fld, {}

    def kexp(self, f, g):
        return 0 if f == g else self.K[(f, g)]

    def dim(self, A, B):
        (pa, fa, la), (pb, fb, lb) = A, B
        e = sum(1 for i in range(4) if fa[i] == fb[i])
        Ks = tuple(self.kexp(fa[i], fb[i]) for i in range(4))
        key = (pa, pb, tuple(sorted(la.items())), tuple(sorted(lb.items())), Ks)
        if key not in self.memo:
            f, p = self.f, self.f.p
            rows = []
            for (a, b) in PAIRINGS[pa]:
                v = [0] * 4
                v[a], v[b] = 1, (-f.rho(a, b, la[(a, b)])) % p
                rows.append(v)
            for (a, b) in PAIRINGS[pb]:             # B: d_b = -sigma d_a with d_i = kappa_i c_i
                v = [0] * 4
                v[a] = pow(f.om, (-Ks[a]) % 3, p)
                v[b] = (-f.rho(a, b, lb[(a, b)]) * pow(f.om, (-Ks[b]) % 3, p)) % p
                rows.append(v)
            M = flint.nmod_mat(4, 4, p)
            for r, row in enumerate(rows):
                for c, x in enumerate(row):
                    M[r, c] = x
            self.memo[key] = 4 - M.rank()
        return e + self.memo[key]


def plane_model_dim(K, P, Q):
    """planes in F_1 + F_2: P = (f1, f2, k); dim(P n Q) = e + [k_P = k_Q + K_1 - K_2 mod 3]."""
    (f1, f2, k), (g1, g2, l) = P, Q
    e = (f1 == g1) + (f2 == g2)
    K1 = 0 if f1 == g1 else K[(f1, g1)]
    K2 = 0 if f2 == g2 else K[(f2, g2)]
    return e + (1 if (k - l - K1 + K2) % 3 == 0 else 0)


def rank_mod(rowsfun, n, q=1000000007):
    M = flint.nmod_mat(n, n, q)
    for i in range(n):
        for j in range(n):
            M[i, j] = rowsfun(i, j) % q
    return M.rank()


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    lines, out, allok = [], {}, True

    ok, fs, vs, K = symbolic_checks()
    allok &= ok
    out['kappa_exponents'] = {f'{p}->{q}': v for (p, q), v in K.items()}
    hist = {e: sum(1 for v in K.values() if v == e) for e in range(3)}
    lines.append(f'(1) Hesse cubic, exact in Q(omega)[lambda]: F(a f + t v(f)) = (1 - lambda^3) t^3 at the 9 flexes, '
                 f'9 distinct flexes, kappa(f,g) = omega^K unique and independent of lambda for all 72 ordered pairs '
                 f'-- {"ok" if ok else "FAIL"}; K = 0, 1, 2 occur {hist[0]}, {hist[1]}, {hist[2]} times')
    for i in range(3):
        lines.append('    K((0,0) -> (i,j)), i = %d: %s' % (i, [K[((0, 0), (i, j))] if (i, j) != (0, 0) else '-'
                                                             for j in range(3)]))

    # (3) model vs direct, two primes
    fields = [Field(find_prime(10 ** 9), 7), Field(find_prime(2 * 10 ** 9 + 11), 8)]
    rnd = random.Random(2026)
    planes = [(f1, f2, k) for f1 in FLEX for f2 in FLEX for k in range(3)]
    for fld in fields:
        mdl = Model(K, fld)
        # planes in F_1 + F_2 (blocks 0, 1 of the field), all 243^2 pairs
        def pbasis(P):
            f1, f2, k = P
            r = fld.rho(0, 1, k)
            return [list(fld.fs[0][f1]) + [0, 0, 0], [0, 0, 0] + list(fld.fs[1][f2]),
                    list(fld.vs[0][f1]) + [(-r * x) % fld.p for x in fld.vs[1][f2]]]
        pb = [pbasis(P) for P in planes]
        bad = sum(1 for i in range(243) for j in range(243)
                  if fld.direct_dim(pb[i], pb[j], 3, 6) != plane_model_dim(K, planes[i], planes[j]))
        joins = [random_join(rnd) for _ in range(400)]
        jb = [fld.join(J) for J in joins]
        inX = all(fld.h_vanishes(b, rnd) for b in jb)
        tested, badj = 0, 0
        for i in range(len(joins)):
            for j in range(i, len(joins)):
                if (i * 7 + j) % 5 == 0 or i == j:
                    tested += 1
                    badj += fld.direct_dim(jb[i], jb[j]) != mdl.dim(joins[i], joins[j])
        # pairs that share flexes (the interesting cases): perturb one join
        near = 0
        for _ in range(3000):
            A = random_join(rnd)
            pi = rnd.randrange(3)
            flex = dict(A[1])
            for i in range(4):
                if rnd.random() < 0.5:
                    flex[i] = FLEX[rnd.randrange(9)]
            B = (pi, flex, {pr: rnd.randrange(3) for pr in PAIRINGS[pi]})
            near += fld.direct_dim(fld.join(A), fld.join(B)) != mdl.dim(A, B)
        good = bad == 0 and badj == 0 and near == 0 and inX
        allok &= good
        lines.append(f'(3) p = {fld.p}: planes in F_1+F_2: model = direct for all 243^2 pairs (mismatches {bad}); '
                     f'joins: cubic vanishes on 400 random joins: {inX}; model = direct on {tested} random pairs '
                     f'(mismatches {badj}) and on 3000 pairs sharing flexes (mismatches {near}) -- '
                     f'{"ok" if good else "FAIL"}')

    # (4) Gram ranks (model, exact intersection dimensions)
    mdl = Model(K, fields[0])
    q = 1000000007
    rk_planes = rank_mod(lambda i, j: -(-2) ** plane_model_dim(K, planes[i], planes[j]), 243, q)
    ok4 = rk_planes <= 18
    lines.append(f'(4a) fourfold F_1 + F_2: rank of the primitive Gram matrix of the 243 planes = {rk_planes} '
                 f'(dim Hdg^2_prim = 18 for non-isogenous curves; at most 18 required) -- {"ok" if ok4 else "FAIL"}')
    out['rank_planes_prim'] = rk_planes
    allok &= ok4
    rnd = random.Random(99)
    res = {}
    for label, pis, n in (('pairing 12|34', (0,), 900), ('pairings 12|34, 13|24', (0, 1), 1300),
                          ('all three pairings', (0, 1, 2), 1300), ('all three pairings', (0, 1, 2), 2000)):
        js = []
        while len(js) < n:
            J = random_join(rnd)
            if J[0] in pis:
                js.append(J)
        rk = rank_mod(lambda i, j: -(-2) ** mdl.dim(js[i], js[j]), n, q)
        res[f'{label}, n={n}'] = rk
        bound = 324 if len(pis) == 1 else 486
        okr = rk <= bound
        allok &= okr
        lines.append(f'(4b) tenfold, {label}: {n} random joins, rank of the primitive Gram matrix mod {q} = {rk} '
                     f'(Hodge bound {bound}) -- {"ok" if okr else "FAIL"}')
    out['rank_joins_prim'] = res
    full = res.get('all three pairings, n=2000') == 486
    lines.append(f'(4c) the joins and h^5 span a subspace of Hdg^5 of dimension 1 + {res.get("all three pairings, n=2000")}'
                 f'{" = 487 = dim Hdg^5 for pairwise non-isogenous curves not isogenous to E_omega" if full else ""}')
    out['span_is_Hdg5'] = full
    lines.append('ALL CHECKS PASS' if allok else 'SOME CHECK FAILED')
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 't6_joins.txt'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
    with open(os.path.join(OUT, 't6_joins.json'), 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(out, fh, indent=1, default=str)
    print('\n'.join(lines))
