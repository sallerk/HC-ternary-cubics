# t1_sk_smooth.py -- checks for Sections 1 and 3 of the ternary-cubics NOTE.md, for degree d = 3 with a ternary
# block split off.  Consistency checks only; the proofs use no computation.
# (1) phi((x:y),(x':y')) = (y'x : eps y x'), eps = -1, maps X1 x X2 into W = {f + g = 0}, symbolically, for
#     f = (generic ternary cubic) + z^3 and g = generic ternary cubic (symbolic coefficients).
# (2) Fulton's degree identity (Prop. 4.4 in Fulton's Intersection theory, as in the quartic note's check) for
#     dim X2 = 2 (a cubic surface) and 1 <= dim X1 <= 10: deg(psi) = 3.
# (3) Fibre count over F_p, p = 1 mod 3, by brute force over all of X1(F_p) x X2(F_p): every point (x:x') of
#     W(F_p) with f(x) != 0 has exactly #{y in F_p : y^3 = -f(x)} preimages; the image lies in W.
#     Cases: W = {F1(x) + F2(x') = 0} in P^5 (X1, X2 cubic surfaces) and W = {F1(x) + z^3 + F2(x') = 0} in P^6.
# (4) Smoothness (Lemma 1.1): for split cubics with random ternary blocks over F_p (p = 10007), the Jacobian ring
#     R = F_p[x]/(dh) has Hilbert function binom(N, t) in degrees t = 0..N and R_{N+1} = 0, so the partials have no
#     common zero but 0 (p != 3).  Control: a block with a nodal plane cubic gives R_{N+1} != 0.
import itertools, json, os, random
import sympy as sp
import flint

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')
out, lines = {}, []

# ---------- (1) symbolic identity ----------
x = sp.symbols('x0:3'); z, y = sp.symbols('z y'); xp = sp.symbols('a0:3'); yp = sp.Symbol('yp')
mons3 = [e for e in itertools.product(range(4), repeat=3) if sum(e) == 3]
cf = sp.symbols('c0:10'); ce = sp.symbols('e0:10')
tern = lambda X, C: sum(c * X[0] ** e[0] * X[1] ** e[1] * X[2] ** e[2] for c, e in zip(C, mons3))
f = lambda X, Z: tern(X, cf) + Z ** 3
g = lambda X: tern(X, ce)
eps = -1
img = f(tuple(yp * t for t in x), yp * z) + g(tuple(eps * y * t for t in xp))
ok_hom = sp.expand(img - (yp ** 3 * f(x, z) + eps ** 3 * y ** 3 * g(xp))) == 0
# on X1 x X2: f(x, z) = -y^3 and g(x') = -yp^3
ok_land = sp.expand(yp ** 3 * (-y ** 3) + eps ** 3 * y ** 3 * (-yp ** 3)) == 0
out['phi_lands_in_W'] = bool(ok_hom and ok_land)
lines.append('(1) phi maps X1 x X2 to W (symbolic, f = ternary + z^3, g = ternary, d = 3, eps = -1): %s'
             % out['phi_lands_in_W'])

# ---------- (2) Fulton degree identity, s = dim X2 = 2 ----------
A, B = sp.symbols('A B')
fulton = {}
s = 2
for r in range(1, 11):
    n = r + s
    ser = sp.series(sp.series((1 + A + B) ** n / ((1 + A) * (1 + B)), A, 0, r).removeO(), B, 0, s).removeO()
    coeff = sp.Poly(sp.expand(ser), A, B).coeff_monomial(A ** (r - 1) * B ** (s - 1))
    fulton[r] = int(sp.binomial(n, r) - coeff)        # deg(psi) = 3 * this
out['fulton_deg_psi_over_3_s2'] = fulton
lines.append('(2) Fulton: deg(psi)/3 for dim X2 = 2, dim X1 = 1..10: %s  (all 1: %s)'
             % (list(fulton.values()), all(v == 1 for v in fulton.values())))


# ---------- helpers over F_p ----------
def proj_points(n, p):
    """Normalised points of P^n(F_p) (first non-zero coordinate 1)."""
    pts = []
    for lead in range(n + 1):
        for rest in itertools.product(range(p), repeat=n - lead):
            pts.append((0,) * lead + (1,) + rest)
    return pts


def normalise(v, p):
    for c in v:
        if c % p:
            inv = pow(c, p - 2, p)
            return tuple((t * inv) % p for t in v)
    return None


def poly_eval(P, v, p):
    tot = 0
    for e, c in P.items():
        t = c
        for i, k in enumerate(e):
            if k:
                t = t * pow(v[i], k, p)
        tot += t
    return tot % p


def random_ternary(p, rnd):
    while True:
        P = {e: rnd.randrange(p) for e in mons3}
        if smooth_plane_cubic(P, p):
            return P


def partials(P, nvar):
    D = []
    for i in range(nvar):
        Q = {}
        for e, c in P.items():
            if e[i]:
                e2 = list(e); e2[i] -= 1
                Q[tuple(e2)] = (Q.get(tuple(e2), 0) + c * e[i])
        D.append(Q)
    return D


def hilbert_function(P, nvar, p, tmax):
    """dim of (F_p[x]/(partials of P))_t for t = 0..tmax, via ranks over F_p (python-flint)."""
    D = partials(P, nvar)
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


def smooth_plane_cubic(P, p):
    return hilbert_function(P, 3, p, 4) == [1, 3, 3, 1, 0]


def embed(P, nvar, offset):
    return {tuple([0] * offset + list(e) + [0] * (nvar - offset - len(e))): c for e, c in P.items()}


def add(*Ps):
    S = {}
    for P in Ps:
        for e, c in P.items():
            S[e] = S.get(e, 0) + c
    return S


# ---------- (3) fibre counts ----------
def fibre_check(p, F1, F2, with_z, seed):
    a = 4 if with_z else 3                      # number of variables of f
    fdict = add(embed(F1, a, 0), {(0, 0, 0, 3): 1} if with_z else {})
    X1 = [v for v in proj_points(a, p) if (poly_eval(fdict, v[:a], p) + pow(v[a], 3, p)) % p == 0]
    X2 = [v for v in proj_points(3, p) if (poly_eval(F2, v[:3], p) + pow(v[3], 3, p)) % p == 0]
    fib = {}
    for u in X1:
        for w in X2:
            yy, yyp = u[a], w[3]
            if yy == 0 and yyp == 0:
                continue                            # Z = {y = y' = 0}
            img = normalise(tuple((yyp * t) % p for t in u[:a]) + tuple((-yy * t) % p for t in w[:3]), p)
            assert img is not None
            assert (poly_eval(fdict, img[:a], p) + poly_eval(F2, img[a:], p)) % p == 0, 'image not in W'
            fib[img] = fib.get(img, 0) + 1
    cubes = {}
    for t in range(p):
        cubes[pow(t, 3, p)] = cubes.get(pow(t, 3, p), 0) + 1
    W = [v for v in proj_points(a + 2, p) if (poly_eval(fdict, v[:a], p) + poly_eval(F2, v[a:], p)) % p == 0]
    nU = bad = 0
    for v in W:
        fx = poly_eval(fdict, v[:a], p)
        if fx == 0:
            continue
        nU += 1
        want = cubes.get((-fx) % p, 0)
        if fib.get(v, 0) != want:
            bad += 1
    return dict(p=p, dimW=a + 1, X1=len(X1), X2=len(X2), W=len(W), U=nU, mismatches=bad)


rnd = random.Random(20261025)
fib_results = []
for p in (7, 13):
    F1 = random_ternary(p, rnd); F2 = random_ternary(p, rnd)
    for with_z in (False, True):
        if with_z and p > 7:
            continue
        r = fibre_check(p, F1, F2, with_z, 0)
        fib_results.append(r)
        lines.append('(3) p=%d dim W=%d: #X1(F_p)=%d #X2(F_p)=%d #W(F_p)=%d, points with f(x)!=0: %d, mismatches: %d'
                     % (r['p'], r['dimW'], r['X1'], r['X2'], r['W'], r['U'], r['mismatches']))
out['fibre_counts'] = fib_results

# ---------- (4) smoothness via the Jacobian ring ----------
P = 10007
sm = []
for (k, l) in ((2, 0), (2, 1), (1, 3)):
    N = 3 * k + l
    blocks = [random_ternary(P, rnd) for _ in range(k)]
    h = add(*[embed(B, N, 3 * i) for i, B in enumerate(blocks)],
            *[{tuple(3 if j == 3 * k + i else 0 for j in range(N)): 1} for i in range(l)])
    hf = hilbert_function(h, N, P, N + 1)
    want = [sp.binomial(N, t) for t in range(N + 1)] + [0]
    sm.append(dict(k=k, l=l, N=N, hilbert=hf, expected=[int(w) for w in want], ok=hf == want))
    lines.append('(4) k=%d l=%d N=%d: Hilbert function of R_h %s, expected binom(N,t) and R_{N+1}=0: %s'
                 % (k, l, N, hf, hf == want))
# control: nodal block y^2 z - x^3 - x^2 z
nodal = {(0, 2, 1): 1, (3, 0, 0): P - 1, (2, 0, 1): P - 1}
N = 6
h = add(embed(nodal, N, 0), embed(random_ternary(P, rnd), N, 3))
hf = hilbert_function(h, N, P, N + 1)
sm.append(dict(control='nodal block', N=N, hilbert=hf, R_N1_nonzero=hf[-1] != 0))
lines.append('(4) control, nodal block + smooth block, N=6: Hilbert function %s, R_7 != 0: %s' % (hf, hf[-1] != 0))
out['smoothness'] = sm

os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, 't1_sk_smooth.json'), 'w', encoding='utf-8') as fh:
    json.dump(out, fh, indent=1, default=str)
with open(os.path.join(OUT, 't1_sk_smooth.txt'), 'w', encoding='utf-8', newline='\n') as fh:
    fh.write('\n'.join(lines) + '\n')
print('\n'.join(lines))
