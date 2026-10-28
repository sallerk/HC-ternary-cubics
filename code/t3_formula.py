# t3_formula.py -- evaluates the closed formula of Proposition C (NOTE.md Section 5) for dim Hdg^m(X),
#   X = {F_1(x_1,y_1,z_1) + ... + F_k(x_k,y_k,z_k) + z_1^3 + ... + z_l^3 = 0} in P^(N-1), N = 3k + l = 2m + 2,
# and compares it with the Delsarte counts of t2_hodge_count.py and with values from the literature.
# Written from the statement of the proposition, independently of the model count in t2_hodge_count.py.
#
# dim Hdg^m(X) = 1 + sum_chi 3^(n1+n2) binom(|K_w|, |K_w|/2 + j(chi) - |K|/2) prod_kappa nu(K_kappa),
# chi = (chi_1..chi_k; chi'_1..chi'_l) in (Z/3)^k x {1,2}^l with sum = 0 mod 3; n_c = #{s : chi_s = c};
# K = {s : chi_s = 0}; l2 = #{j : chi'_j = 2}; j(chi) = (m + 1 - n1 - 2 n2 - l2)/3 (term 0 unless an integer);
# K_w = {s in K : E_s isogenous to E_omega (j = 0)}; K \ K_w = disjoint union of the K_kappa by isogeny class;
# nu = 0 for odd |K_kappa|, otherwise binom(2e, e) (CM) or binom(2e, e)/(e+1) (no CM), e = |K_kappa|/2.
# Input per block: (isogeny class label, has CM); the label 'w' is the class of E_omega.
import itertools, json, os
from fractions import Fraction
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')


def binom(a, b):
    if b.denominator != 1:
        return 0
    b = int(b)
    return comb(a, b) if 0 <= b <= a else 0


def hodge_dim(blocks, l):
    k = len(blocks); N = 3 * k + l
    assert N % 2 == 0
    m = N // 2 - 1
    total, hmm, prim = 1, 0, 0
    for chi in itertools.product(range(3), repeat=k):
        for chip in itertools.product((1, 2), repeat=l):
            if (sum(chi) + sum(chip)) % 3:
                continue
            n1, n2, l2 = chi.count(1), chi.count(2), chip.count(2)
            K = [s for s in range(k) if chi[s] == 0]
            j = Fraction(m + 1 - n1 - 2 * n2 - l2, 3)
            mult = 3 ** (n1 + n2)
            prim += mult * 2 ** len(K)                                  # dim V(chi)
            hmm += mult * binom(len(K), j)                              # dim of its (m,m)-part
            Kw = [s for s in K if blocks[s][0] == 'w']
            term = binom(len(Kw), Fraction(len(Kw), 2) + j - Fraction(len(K), 2))
            classes = {}
            for s in K:
                if blocks[s][0] != 'w':
                    classes.setdefault(blocks[s][0], []).append(s)
            for ss in classes.values():
                if len(ss) % 2:
                    term = 0
                    break
                e = len(ss) // 2
                term *= comb(2 * e, e) if blocks[ss[0]][1] else comb(2 * e, e) // (e + 1)
            total += mult * term
    return dict(hodge=total, hmm_prim=hmm, dim_prim=prim)


if __name__ == '__main__':
    F = ('w', True)     # x^3 + y^3 + z^3: E_F = E_omega (j = 0)
    J = ('i', True)     # y^2 z + x^3 + x z^2: j = 1728, CM by Z[i], not isogenous to E_omega
    G = ('g', False)    # a curve without CM
    H = ('h', False)    # a second curve without CM, not isogenous to G
    I = ('k', False)    # a third
    L = ('n', False)    # a fourth
    out, lines, ok = {}, [], True
    with open(os.path.join(OUT, 't2_hodge_count.json'), encoding='utf-8') as fh:
        dels = json.load(fh)
    lines.append('Proposition C vs the Delsarte counts of t2_hodge_count.py (F = x^3+y^3+z^3, J = y^2z+x^3+xz^2, l terms z^3)')
    for tag, rec in dels.items():
        if not tag.endswith(tuple('0123456789')) or '|' not in tag:
            continue
        name, l = tag.split('|l=')
        l = int(l)
        bl = [] if name == '-' else [F if c == 'F' else J for c in name]
        r = hodge_dim(bl, l)
        N = 3 * len(bl) + l
        want = rec['delsarte']['hodge']
        good = (r['hodge'] == want and r['dim_prim'] == (2 ** N + 2) // 3
                and r['hmm_prim'] == comb(N, N // 2))
        ok &= good
        out[tag] = dict(r, delsarte=want, agree=good)
        lines.append('   %-12s dim=%2d  formula %6d  Delsarte %6d  h^{m,m}_prim %6d  dim H_prim %6d  %s'
                     % (tag, N - 2, r['hodge'], want, r['hmm_prim'], r['dim_prim'], 'agree' if good else 'DISAGREE'))
    lines.append('Literature controls (cubic fourfolds F_1(x) + F_2(y) = 0)')
    controls = [
        ('E_1, E_2 not isogenous', [G, H], 0, 19, 'Koike Thm 1.7: transcendental rank 4, so 23 - 4 = 19'),
        ('E_1 ~ E_2 without CM', [G, G], 0, 20, '19 + rk Hom(E_1, E_2) = 20'),
        ('E_1 ~ E_2 with CM', [J, J], 0, 21, 'Beauville (Koike Prop 1.6): maximal, h^{2,2} = 21'),
        ('E_1 = E_2 = E_omega', [F, F], 0, 21, 'Fermat cubic fourfold: h^{2,2} = 21'),
        ('E_1 ~ E_omega, E_2 no CM', [F, G], 0, 19, '19 + 0'),
    ]
    for name, bl, l, want, src in controls:
        r = hodge_dim(bl, l)
        good = r['hodge'] == want
        ok &= good
        out['control: ' + name] = dict(r, want=want, source=src, agree=good)
        lines.append('   %-28s formula %4d  expected %4d  (%s)  %s' % (name, r['hodge'], want, src,
                                                                        'agree' if good else 'DISAGREE'))
    lines.append('Blocks without CM (G, H, I, L pairwise non-isogenous, no CM); the very general members')
    for name, bl, l in [('GH|l=0', [G, H], 0), ('GG|l=0', [G, G], 0), ('G|l=5', [G], 5), ('GH|l=2', [G, H], 2),
                        ('GG|l=2', [G, G], 2), ('GHI|l=1', [G, H, I], 1), ('GGG|l=1', [G, G, G], 1),
                        ('GHIL|l=0', [G, H, I, L], 0), ('GGHH|l=0', [G, G, H, H], 0), ('GGGG|l=0', [G] * 4, 0),
                        ('GHI|l=3', [G, H, I], 3), ('G|l=9', [G], 9)]:
        r = hodge_dim(bl, l)
        out[name] = r
        lines.append('   %-12s dim=%2d  formula %6d  h^{m,m}_prim %6d' % (name, 3 * len(bl) + l - 2, r['hodge'],
                                                                         r['hmm_prim']))
    lines.append('ALL AGREE' if ok else 'SOME DISAGREE')
    with open(os.path.join(OUT, 't3_formula.txt'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
    with open(os.path.join(OUT, 't3_formula.json'), 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1)
    print('\n'.join(lines))
