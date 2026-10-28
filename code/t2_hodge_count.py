# t2_hodge_count.py -- dim Hdg^m(X) for X = {sum_s F_s(x_s, y_s, z_s) + sum_j z_j^3 = 0} in P^(N-1),
# N = 3k + l = 2m + 2 (k ternary blocks, l Fermat terms).  Consistency checks only; the proofs use no computation.
#
# (1) DELSARTE count (independent of Proposition C: no elliptic curves, no Shioda-Katsura isomorphism).
#     Blocks 'F' = x^3 + y^3 + z^3 (E_F = E_omega, j = 0) and 'J' = y^2 z + x^3 + x z^2 (j = 1728, CM by Z[i]),
#     and l terms z^3.  For each block: a monomial basis of its Jacobian ring, the characters of its diagonal
#     symmetry group (including the factor of the residue form Omega), multiplicity one checked.  The Galois group
#     of Q(mu_D) (D = 12 if a J block occurs, else 3) acts on characters by multiplication with t in (Z/D)^*.
#     dim Hdg^m = 1 + #{basis tensors of degree m+1 all of whose Galois conjugates have degree m+1}.
# (2) MODEL count: the decomposition of Proposition C (V(chi) = M (x) lambda_chi (x) (x)_{s in K} H^1(E_s)),
#     with the Hodge group invariants found by weight counting (sign patterns for the tori, m_0 - m_2 for SL_2),
#     not by the closed formula; the weight of lambda_chi is read off from the Fermat character that is left
#     when the blocks are split off (alpha = (chi'_j, chi_s for s not in K)).
# (3) Fermat cubic (k = 0) via Shioda's criterion; (4) h^{m,m}_prim = binom(N, m+1), dim H_prim = (2^N + 2)/3.
import itertools, json, math, os
from fractions import Fraction
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')

# ------------------------------------------------------------------ Delsarte character count
BLOCKS = {
    'F': [(3, 0, 0), (0, 3, 0), (0, 0, 3)],      # x^3 + y^3 + z^3
    'J': [(0, 2, 1), (3, 0, 0), (1, 0, 2)],      # y^2 z + x^3 + x z^2
}


def monomials(deg, nvar):
    return [e for e in itertools.product(range(deg + 1), repeat=nvar) if sum(e) == deg]


def jacobian_basis(block):
    """Monomial basis of C[x,y,z]/(dF/dx, dF/dy, dF/dz), F = sum of the monomials of the block."""
    X = sympy.symbols('x y z')
    F = sum(X[0] ** a * X[1] ** b * X[2] ** c for (a, b, c) in block)
    partials = [sympy.diff(F, v) for v in X]
    basis, hilb = [], []
    for deg in range(0, 6):
        mons = monomials(deg, 3)
        gens = []
        if deg >= 2:
            for e in monomials(deg - 2, 3):
                mon = X[0] ** e[0] * X[1] ** e[1] * X[2] ** e[2]
                gens += [sympy.expand(mon * P) for P in partials]

        def vec(expr):
            P = sympy.Poly(expr, *X)
            return [P.coeff_monomial(X[0] ** a * X[1] ** b * X[2] ** c) for (a, b, c) in mons]
        cur = sympy.Matrix([vec(g) for g in gens]) if gens else sympy.zeros(0, len(mons))
        r = cur.rank() if gens else 0
        chosen = []
        for mono in mons:
            row = sympy.Matrix([[1 if mm == mono else 0 for mm in mons]])
            test = row if cur.shape[0] == 0 else cur.col_join(row)
            rt = test.rank()
            if rt > r:
                chosen.append(mono); cur = test; r = rt
        hilb.append(len(chosen))
        basis += chosen
    assert hilb == [1, 3, 3, 1, 0, 0], hilb
    return basis


def block_group(block, D):
    return [g for g in itertools.product(range(D), repeat=3)
            if all(sum(i * a for i, a in zip(mono, g)) % D == 0 for mono in block)]


_JB = {}


def block_info(name, D):
    """(basis, characters): character of the class of the monomial e = (i, j, k) under g = (a, b, c) is
    (i+1) a + (j+1) b + (k+1) c mod D (the +1 is the factor of Omega)."""
    if name == 'z':
        G = [(a,) for a in range(D) if (3 * a) % D == 0]
        B = [(0,), (1,)]
    else:
        G = block_group(BLOCKS[name], D)
        if name not in _JB:
            _JB[name] = jacobian_basis(BLOCKS[name])
        B = _JB[name]
    chars = {}
    for e in B:
        chi = tuple(sum((i + 1) * a for i, a in zip(e, g)) % D for g in G)
        assert chi not in chars, ('multiplicity > 1', name, e, chars[chi])
        chars[chi] = e
    return B, G, chars


def delsarte_count(blocks, l):
    k = len(blocks); N = 3 * k + l
    assert N % 2 == 0
    m = N // 2 - 1
    D = 12 if 'J' in blocks else 3
    names = list(blocks) + ['z'] * l
    units = [t for t in range(1, D) if math.gcd(t, D) == 1]
    # per factor: list of (degree, degrees of the Galois conjugates for t in units)
    factors = []
    for nm in names:
        B, G, chars = block_info(nm, D)
        rows = []
        for x in B:
            conj = []
            for t in units:
                chi = tuple((t * (sum((i + 1) * a for i, a in zip(x, g)) % D)) % D for g in G)
                assert chi in chars, 'Galois conjugate character not realised'
                conj.append(sum(chars[chi]))
            rows.append((sum(x), tuple(conj)))
        factors.append(rows)
    hodge = hmm = total = 0
    for e in itertools.product(*factors):
        total += 1
        if sum(r[0] for r in e) != m + 1:
            continue
        hmm += 1
        if all(sum(r[1][u] for r in e) == m + 1 for u in range(len(units))):
            hodge += 1
    return dict(hodge=1 + hodge, prim=hodge, hmm_prim=hmm, dimR=total, D=D)


# ------------------------------------------------------------------ model count (Proposition C, by weights)
def sl2_invariants(n):
    """dim of SL_2-invariants in (C^2)^{(x) n}: m_0 - m_2 (weight multiplicities)."""
    w = [sum(s) for s in itertools.product((1, -1), repeat=n)]
    return w.count(0) - w.count(2)


def model_count(types, l):
    """types[s] = ('w',) (E_s isogenous to E_omega), ('cm', label) (CM, not ~ E_omega), ('gen', label) (no CM)."""
    k = len(types); N = 3 * k + l
    assert N % 2 == 0
    m = N // 2 - 1
    total = hmm = prim = 0
    for chi in itertools.product(range(3), repeat=k):
        for chip in itertools.product((1, 2), repeat=l):
            if (sum(chi) + sum(chip)) % 3:
                continue
            n12 = sum(1 for c in chi if c)
            K = [s for s in range(k) if chi[s] == 0]
            alpha = list(chip) + [c for c in chi if c]          # Fermat character left over
            L, A = len(alpha), sum(alpha)
            assert A % 3 == 0
            q = Fraction(A, 3) - 1; p = (L - 2) - q              # Hodge type of the Fermat piece V(alpha)
            cl = p - q                                           # weight of lambda_chi: z^{p-q}
            mult = 3 ** n12
            prim += mult * 2 ** len(K)
            # Hodge types: sign +1 = type (1,0), -1 = type (0,1) for each H^1(E_s), s in K; (m,m) <=> weight 0
            hmm += mult * sum(1 for sg in itertools.product((1, -1), repeat=len(K)) if cl + sum(sg) == 0)
            # Hodge group: Hg(E_omega) x prod over isogeny classes (Imai)
            gen = {}
            tori = []
            for s in K:
                t = types[s]
                if t[0] == 'gen':
                    gen.setdefault(t[1], []).append(s)
                else:
                    tori.append(s)
            fg = 1
            for ss in gen.values():
                fg *= sl2_invariants(len(ss))
            if fg == 0:
                continue
            cnt = 0
            for sg in itertools.product((1, -1), repeat=len(tori)):
                sig = dict(zip(tori, sg))
                if cl + sum(sig[s] for s in tori if types[s][0] == 'w') != 0:
                    continue
                labels = set(types[s][1] for s in tori if types[s][0] == 'cm')
                if all(sum(sig[s] for s in tori if types[s][0] == 'cm' and types[s][1] == lab) == 0
                       for lab in labels):
                    cnt += 1
            total += mult * fg * cnt
    return dict(hodge=1 + total, prim=total, hmm_prim=hmm, dim_prim=prim)


def fermat_shioda(N):
    """Fermat cubic in N variables: alpha in {1,2}^N, sum = 0 mod 3; V(alpha) in H^{n-q,q} with q = |alpha|/3 - 1;
    Hodge iff |t alpha| = 3N/2 for t = 1, 2."""
    cnt = hmm = 0
    for a in itertools.product((1, 2), repeat=N):
        if sum(a) % 3:
            continue
        if 2 * sum(a) == 3 * N:
            hmm += 1
            if 2 * sum((2 * x) % 3 for x in a) == 3 * N:
                cnt += 1
    return dict(hodge=1 + cnt, prim=cnt, hmm_prim=hmm)


def main():
    T = {'F': ('w',), 'J': ('cm', 'i')}
    out, lines, allok = {}, [], True
    lines.append('(1)+(4) DELSARTE character count (independent) vs MODEL count of Proposition C; blocks F = x^3+y^3+z^3,')
    lines.append('        J = y^2z+x^3+xz^2, plus l terms z^3; N = 3k+l even, 4 <= N <= 16')
    for N in range(4, 17, 2):
        for k in range(0, N // 3 + 1):
            l = N - 3 * k
            for combo in itertools.combinations_with_replacement('FJ', k):
                tag = '%s|l=%d' % (''.join(combo) or '-', l)
                dl = delsarte_count(list(combo), l)
                md = model_count([T[b] for b in combo], l)
                m = N // 2 - 1
                hmm_want = math.comb(N, m + 1)
                prim_want = (2 ** N + 2) // 3
                ok = (dl['hodge'] == md['hodge'] and dl['hmm_prim'] == md['hmm_prim'] == hmm_want
                      and dl['dimR'] == 8 ** k * 2 ** l and md['dim_prim'] == prim_want)
                allok &= ok
                out[tag] = dict(N=N, k=k, l=l, delsarte=dl, model=md, agree=ok)
                lines.append('   N=%2d dim=%2d %-12s Delsarte Hdg=%-6d model Hdg=%-6d h^{m,m}_prim=%-6d dim H_prim=%-6d %s'
                             % (N, N - 2, tag, dl['hodge'], md['hodge'], dl['hmm_prim'], md['dim_prim'],
                                'AGREE' if ok else 'DISAGREE'))
    lines.append('(3) Fermat cubic via Shioda criterion vs Delsarte (k = 0) and model (all blocks F)')
    for N in range(4, 17, 2):
        fs = fermat_shioda(N)
        dl = out['-|l=%d' % N]['delsarte']
        k = N // 3
        md_list = []
        for kk in range(0, k + 1):
            if (N - 3 * kk) >= 0:
                md_list.append(model_count([('w',)] * kk, N - 3 * kk)['hodge'])
        ok = fs['hodge'] == dl['hodge'] == 1 + math.comb(N, N // 2) and all(v == fs['hodge'] for v in md_list)
        allok &= ok
        out['fermat_N%d' % N] = dict(shioda=fs, model_all_F=md_list, agree=ok)
        lines.append('   N=%2d  Shioda Hdg=%-6d Delsarte=%-6d 1+binom(N,N/2)=%-6d model (k F-blocks, k=0..%d)=%s %s'
                     % (N, fs['hodge'], dl['hodge'], 1 + math.comb(N, N // 2), k, md_list, 'AGREE' if ok else 'DISAGREE'))
    lines.append('ALL AGREE' if allok else 'SOME DISAGREE')
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 't2_hodge_count.json'), 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1)
    with open(os.path.join(OUT, 't2_hodge_count.txt'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
