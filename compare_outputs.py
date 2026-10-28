"""compare_outputs.py -- compare the outputs of a rerun with the reference outputs.

Usage: python3 -I compare_outputs.py REFERENCE_FOLDER RERUN_FOLDER

Every .txt and .json file of REFERENCE_FOLDER is compared with the file of the same name in RERUN_FOLDER, line by
line. Line ends are ignored (the reference outputs were written on Windows, some with CRLF line ends); nothing else
is masked, as the scripts print no running times. Files of RERUN_FOLDER without a reference are listed. The last
line is "ALL EQUAL" if every reference file has an equal counterpart, and "NOT ALL EQUAL" otherwise.
"""
import os
import sys


def lines(path):
    with open(path, 'rb') as fh:
        return fh.read().decode('utf-8').splitlines()


def main(ref, new):
    names = sorted(n for n in os.listdir(ref) if n.endswith(('.txt', '.json')))
    ok = True
    for name in names:
        p = os.path.join(new, name)
        if not os.path.isfile(p):
            print(f'{name}: MISSING in the rerun')
            ok = False
            continue
        with open(os.path.join(ref, name), 'rb') as fa, open(p, 'rb') as fb:
            same_bytes = fa.read() == fb.read()
        a, b = lines(os.path.join(ref, name)), lines(p)
        if a == b:
            print(f'{name}: equal' + (' (byte for byte)' if same_bytes else ' (up to line ends)'))
            continue
        ok = False
        k = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
        print(f'{name}: DIFFERS from line {k + 1} ({len(a)} reference lines, {len(b)} rerun lines)')
        print(f'  reference: {a[k] if k < len(a) else "<end of file>"}')
        print(f'  rerun:     {b[k] if k < len(b) else "<end of file>"}')
    extra = sorted(n for n in os.listdir(new) if n.endswith(('.txt', '.json')) and n not in names)
    for name in extra:
        print(f'{name}: in the rerun only')
    print('ALL EQUAL' if ok else 'NOT ALL EQUAL')
    return 0 if ok else 1


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
