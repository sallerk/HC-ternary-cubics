# t4_tenfold_range.py -- Corollary 5.2(ii) of the note: for the cubic tenfolds F_1 + F_2 + F_3 + F_4 = 0 (k = 4, l = 0),
# evaluates the formula of Proposition C (t3_formula.hodge_dim) over all configurations of the four curves E_s:
# a set partition of {1,2,3,4} into isogeny classes, each class labelled 'w' (isogenous to E_omega; at most one class),
# 'cm' (complex multiplication, not isogenous to E_omega) or 'gen' (no complex multiplication).
# Checks 487 <= dim Hdg^5 <= 925, with 487 exactly when the four classes are singletons and none is labelled 'w',
# and that the value depends only on the configuration up to permuting the curves.
import importlib.util, itertools, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')
_spec = importlib.util.spec_from_file_location('t3_formula', os.path.join(HERE, 't3_formula.py'))
t3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(t3)


def set_partitions(s):
    if not s:
        yield []
        return
    first, rest = s[0], s[1:]
    for p in set_partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    parts = list(set_partitions([0, 1, 2, 3]))
    configs, bad = [], 0
    for part in parts:
        for kinds in itertools.product(('w', 'cm', 'gen'), repeat=len(part)):
            if kinds.count('w') > 1:
                continue                                    # only one class is the class of E_omega
            blocks = [None] * 4
            for ci, (cls, kd) in enumerate(zip(part, kinds)):
                for s in cls:
                    blocks[s] = ('w', True) if kd == 'w' else ('c%d' % ci, kd == 'cm')
            d = t3.hodge_dim(blocks, 0)['hodge']
            generic = len(part) == 4 and 'w' not in kinds
            if (d == 487) != generic or not 487 <= d <= 925:
                bad += 1
            orbit = tuple(sorted((len(c), kd) for c, kd in zip(part, kinds)))   # up to permuting the curves
            configs.append(dict(partition=[[s + 1 for s in c] for c in part], kinds=list(kinds),
                                orbit=[list(o) for o in orbit], hodge=d))
    orbits = {}
    for c in configs:
        orbits.setdefault(tuple(map(tuple, c['orbit'])), set()).add(c['hodge'])
    invariant = all(len(v) == 1 for v in orbits.values())
    values = sorted({c['hodge'] for c in configs})
    lines = ['Cubic tenfolds F_1+F_2+F_3+F_4 = 0: dim Hdg^5 from Proposition C over all configurations of E_1..E_4',
             'set partitions of {1,2,3,4}: %d' % len(parts),
             'configurations (classes labelled w / cm / gen, at most one w): %d' % len(configs),
             'configurations up to permuting the curves: %d' % len(orbits),
             'value depends only on the configuration up to permutation: %s' % invariant,
             'min %d, max %d' % (values[0], values[-1]),
             'distinct values (%d): %s' % (len(values), ' '.join(map(str, values))),
             'configurations violating "487 <= dim <= 925, and 487 iff four singleton classes, none w": %d' % bad,
             'value 925 attained by: %s' % '; '.join(str(o) for o, v in sorted(orbits.items()) if 925 in v),
             'per configuration up to permutation, (class size, label) ... : dim Hdg^5']
    for o, v in sorted(orbits.items()):
        lines.append('   %-52s %d' % (' '.join('%d%s' % (n, kd) for n, kd in o), min(v)))
    ok = (len(parts) == 15 and invariant and bad == 0 and values[0] == 487 and values[-1] == 925)
    lines.append('ALL CHECKS PASS' if ok else 'SOME CHECKS FAIL')
    with open(os.path.join(OUT, 't4_tenfold_range.txt'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
    with open(os.path.join(OUT, 't4_tenfold_range.json'), 'w', encoding='utf-8') as fh:
        json.dump(dict(configurations=configs, values=values, ok=ok), fh, indent=1)
    print('\n'.join(lines))
