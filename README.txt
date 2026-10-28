Supporting code and data for the note
"Cubic surfaces, Eckardt points and the Hodge conjecture for sums of cubic forms"
================================================================================

Status: accompanies the note (K. Saller, 2026), DOI 10.5281/zenodo.23746290; not peer reviewed.

Setting. X in P^(N-1) is a smooth cubic hypersurface F_1 + ... + F_k + w_1^3 + ... + w_l^3 = 0 whose blocks F_s are
cubic forms in disjoint sets of variables: ternary blocks (each defines a smooth plane cubic curve E_s), quaternary
blocks G (smooth cubic surfaces S_G; Y_G = {G + y^3 = 0} is a cyclic cubic threefold, J_G its intermediate
Jacobian), among them the Eckardt blocks G = C(x_1, x_2, x_3) + w^2 L(x_1, x_2, x_3), and quinary blocks H
(Y_H = {H = 0} is a cubic threefold and Z_H = {H + y^3 = 0} a cyclic cubic fourfold). E_omega is the elliptic curve
with j = 0; the scripts write F for the block x^3 + y^3 + z^3 (E_F = E_omega) and J for y^2 z + x^3 + x z^2
(j(E_J) = 1728).

Role of the computations. Section 12 of the note lists the items (C1)-(C8); Section 3 below says which script makes
each of them. Except for Lemma 6.2(a), Proposition 6.3, the numerical values in Section 11 and in Theorem D(b), and
Remarks 11.2 and 11.3, the proofs use no computation; the other items are safeguards. The arithmetic is exact: over
Q and Q(omega)[lambda] (SymPy), over finite fields (python-flint and plain Python), and with integers and fractions
for the counts. Where the note says "random", the choices come from seeded generators (fixed seeds in the scripts).

Not checked by any script: the proofs; the cited results (among them those of Shioda and Katsura, Allcock, Carlson
and Toledo, Achter, Schoen, Markman and the author's earlier notes); and the hypotheses on Hodge groups under which
the dimensions of Section 11 are counted, which the scripts assume.

Labels. The headers of the scripts name statements of the working notes NOTE.md and of a working report
five-var/REPORT.md, neither of which is part of this repository (the report was never saved), in the numbering of
the day each script was written, and steps of the working log ("Phase 99", "Phase 101"); some outputs repeat these
labels. Section 4 translates them into the numbering of the note.

Licence: the code (*.py, *.sh) is under the MIT licence; the outputs and this text are under CC BY 4.0 (see the
LICENSE files of the repository). No third-party files are included.


1. Requirements
---------------
Python 3 with SymPy (t1, t2, t6 and t7_eckardt, and t5, which loads t1) and python-flint (t1, t6 and t7_quinary, and
t5 and t7_eckardt, which load t1); t3 and t4 use the standard library only, as does compare_outputs.py. Tested
versions: Python 3.12.9, SymPy 1.13.1 and python-flint 0.9.0, on Windows 10; these are the versions that Section 12
of the note names.


2. Layout and how to run
------------------------
code/      t1_sk_smooth.py, t2_hodge_count.py, t3_formula.py, t4_tenfold_range.py, t5_quaternary.py, t6_joins.py,
           t7_eckardt.py and t7_quinary.py
outputs/   the reference outputs: for each script NAME.txt, the report that it writes (and prints), and NAME.json,
           the same results in machine-readable form
run_all.sh, compare_outputs.py

Each script writes NAME.txt and NAME.json into ../outputs, relative to its own folder, and loads the scripts it uses
from its own folder: t3 reads t2_hodge_count.json; t4 loads t3; t5 loads t1 and t3; t7_eckardt loads t1 and t5; and
t7_quinary loads t5 without the line of t5 that loads t1. t1 has no main guard, so every load of t1 runs its checks
again and rewrites t1_sk_smooth.txt and t1_sk_smooth.json, with the same content (its random choices are seeded); in
a full run this happens once in t5 and twice in t7_eckardt (directly and through t5). run_all.sh therefore copies
code/ to OUT/code and runs the scripts there, so that the reference outputs in outputs/ are not touched.

From the top folder of this repository:
    sh run_all.sh rerun
    python3 -I compare_outputs.py outputs rerun/outputs
run_all.sh runs every script once, in the order t1, t2, t3, t4, t5, t6, t7_eckardt, t7_quinary (t2 must run before
t3), each as python3 -I -B SCRIPT, with its screen output in OUT/logs/NAME.log, and prints a line "NAME exit N" for
each run. The shell script calls python3; set the environment variable PYTHON to use another interpreter name (for
example PYTHON=python). compare_outputs.py compares each .txt and .json file with the reference, line by line, and
ignores the line ends (the .json files of t1-t4 have CRLF line ends, from Windows); the scripts print no running
times, so nothing else is masked. Its expected last line is "ALL EQUAL". Running time: about a minute and a half for
the whole run on the machine of Section 5.


3. The scripts and the items of Section 12 of the note
------------------------------------------------------
Statement numbers are those of the note; (1), (2), ... are the checks named in the header and in the output of each
script.

(C1)  t1, checks (1) and (3): Proposition 3.1, the map phi, symbolically, and the fibre counts over F_7 and F_13 by
      enumerating X1(F_p) x X2(F_p).
(C2)  t1, check (4): Lemma 2.1, the Hilbert functions of the Jacobian rings over F_10007, with the nodal control.
      t1 also makes check (2), the degree identity deg(psi) = 3 of Fulton's Proposition 4.4 for dim X2 = 2 and
      1 <= dim X1 <= 10, which the note does not list.
(C3)  t2: Proposition 5.1 for the blocks F and J and terms z^3 with N even, 4 <= N <= 16 (76 cases): (1) the
      Delsarte count from the characters of the diagonal symmetry group, (2) the second implementation, which counts
      invariants by weights, (3) the Fermat cubic by Shioda's criterion, (4) h^{m,m}_prim and dim H^{2m}_prim; the
      values printed in (C3) are in t2_hodge_count.txt. t3: Proposition 5.1 evaluated from its statement, against
      the Delsarte counts of t2; it also evaluates the formula for some blocks without complex multiplication and
      compares it, for five cubic fourfolds F_1(x) + F_2(y) = 0, with values from [Koi21] (the note's reference)
      and of the Fermat cubic fourfold.
(C4)  t4: Corollary 5.2(ii): the 227 configurations of the four curves (38 up to permuting them) and the 18 values.
(C5)  t5: (1) the Jacobian rings over F_10007 and the eigenspaces of sigma^* on H^{2,1}(Y_G) and H^{1,2}(Y_G)
      (Lemma 7.1); (2) the block table; (3) the numbers I_5(a, b) by the Brauer-Klimyk rule, against Lemma 2.8;
      (4) Proposition 11.1 without Eckardt blocks, three ways, on 23 configurations (the explicit weights for at
      most nine blocks), among them the values of Theorem D(b), and against Proposition 5.1 on 467 configurations
      without quaternary blocks; (5) the 46 specialisations G -> F + u^3; (6) the search over all configurations
      with N <= 24 for the first part of Remark 11.2.
(C6)  t6: (1) Lemma 6.2(a), the table K(f, g), symbolically, for all 72 ordered pairs of distinct flexes; (2) the
      243 planes and the joins; (3) Lemma 6.2(b) against direct ranks over F_p for p = 1000000009 and 2000000011,
      and the cubic vanishing on 400 random joins; (4) the ranks of the primitive Gram matrices modulo 10^9 + 7
      (Lemma 6.1 and Proposition 6.3): 18, 324, 486, 486 and 486.
(C7)  t7_eckardt: (1) the rank 4 of the pairing for generic G (Proposition 7.4(a)); (2) Eckardt blocks: smoothness,
      the Hilbert functions of Lemma 7.6(a), the Hodge numbers of W^- and W^+, the ranks 3 and 1 of the even and odd
      pairings (Proposition 7.7), the two singular controls, and the Fermat surface in the form of Lemma 7.5;
      (3a) the numbers I_4(a, b); (3b) Proposition 11.1 with Eckardt blocks, two ways, on 22 configurations, against
      t5 for blocks satisfying (H1), and the inequalities between the columns of Table 1; (4) the second part of
      Remark 11.2.
(C8)  t7_quinary: (1)-(2) the Jacobian rings, the eigenspaces of H^4(Z_H) and the ranks of the multiplication maps
      of Section 8; (3) the pieces of a quinary block; (4) the invariants of SL_n by Schur-Weyl duality, against t5
      for n = 5, and of Sp_10; (5a) the dimensions of Section 11 with quinary blocks, against t5 on 233
      configurations without quinary blocks; (5b) the 28 specialisations; (6a)-(6c) the thresholds of Remark 11.3,
      (6a) for N <= 62. Check (5c), on the three P^4 in H(x) + H(y) = 0, is not used by the note.


4. Labels in the scripts
------------------------
The headers name statements of NOTE.md, in the numbering of the day each script was written. In the numbering of
the note:
  t1      "Sections 1 and 3 of the ternary-cubics NOTE.md"     Lemma 2.1 and Proposition 3.1
          "Lemma 1.1"                                           Lemma 2.1
          "the quartic note's check"                            the same check for the note on sums of binary
                                                                quartics ([Sal26] in the note)
  t2, t3, t4  "Proposition C"                                   Proposition 5.1
  t5      "NOTE.md Sections 8-10"                               Sections 7, 9 and 11
          "Proposition C'"                                      Proposition 11.1 without Eckardt blocks
          "hypothesis (H1) of Theorem D"                        hypothesis (H1) of Theorem C
  t6      "NOTE.md Section 6"                                   Section 6 (Lemmas 6.1, 6.2 and Proposition 6.3)
          "Koike Thm 1.7", "Prop C"                             [Koi21, Theorem 1.7], Proposition 5.1
  t7_eckardt  "NOTE.md Section 11"                              Section 7 (Proposition 7.4 to Proposition 7.7)
                                                                and Section 11
          "Proposition C''"                                     Proposition 11.1 with Eckardt blocks
  t7_quinary  "five-var/REPORT.md"                              Section 8 and the quinary cases of Section 11
          "(H1')", "(H3)" and the independence of distinct      condition (Q) of Section 8
          classes
"Phase 99" and "Phase 101" are steps of the working log, also not part of this repository.


5. Reference outputs and the rerun of 28 October 2026
-----------------------------------------------------
Reference outputs. The files in outputs/ were made in the working tree, each by the last version of its script:
those of t2, t3 and t4 on 25 October 2026, those of t5, t6, t7_eckardt and t7_quinary on 27 October; those of t1
were last rewritten on 28 October (t1 rewrites them whenever it runs or is loaded, Section 2). The platform and the
versions of these runs are not recorded in them.

Changes against the working tree. t7_quinary.py and its outputs were in a separate folder of the working tree, and
the script read t5_quaternary.py from the folder of the other scripts there; it now reads it from its own folder
(line 26). All other scripts and outputs are unchanged. run_all.sh, compare_outputs.py, the README files, the
licences and .gitattributes are new.

The rerun. When this repository was assembled, on 28 October 2026, run_all.sh was run from its files on the author's
machine (Windows 10, the sh of Git for Windows, and the versions listed in Section 1), with PYTHON=python and a
folder outside the repository as OUT. It took 87 seconds: t1 about 16, t5, t6 and t7_eckardt between 18 and 26
each (t5 and t7_eckardt run the checks of t1 again when they load it), and the other scripts a few seconds or less.
All 8 runs ended with exit status 0, and compare_outputs.py found all 16 files equal to the reference outputs, byte
for byte ("ALL EQUAL"). The reference outputs were made on the same machine, so the comparison shows only that two
runs there agree.
