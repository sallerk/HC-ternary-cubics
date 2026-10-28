# HC-ternary-cubics

Code and data supporting the note *Cubic surfaces, Eckardt points and the Hodge conjecture for sums of cubic forms*
(K. Saller, 2026), DOI: [10.5281/zenodo.23746290](https://doi.org/10.5281/zenodo.23746290).

**Subject.** The note studies smooth cubic hypersurfaces F₁ + … + F_k + w₁³ + … + w_l³ = 0 whose blocks F_s are
cubic forms in disjoint sets of at most five variables. For ternary blocks, which define smooth plane cubic curves,
the Hodge conjecture holds for every such hypersurface and for the products of such hypersurfaces and of curves whose
Jacobians are isogenous to products of elliptic curves, by the inductive structure of Shioda and Katsura (this case
was stated before; see the introduction of the note). The same structure gives a formula for the dimension of the
space of Hodge classes. For blocks G in four variables, which define smooth cubic surfaces, the Hodge classes are
governed by the intermediate Jacobian of the cubic threefold G + y³ = 0, an abelian fivefold with multiplication by
Q(√−3) of signature (4,1). For very general blocks the note proves the Hodge conjecture when each block occurs at
most four times, and for any multiplicities if the Weil classes of certain abelian eightfolds are algebraic. Very
general blocks whose cubic surface has an Eckardt point may occur any number of times, by work of Schoen and a
theorem of Markman (a preprint), and very general blocks in five variables at most ten times.

**Role of the computations.** Section 12 of the note lists eight computations, (C1)–(C8). Except for Lemma 6.2(a),
Proposition 6.3, the numerical values in Section 11 and in Theorem D(b), and Remarks 11.2 and 11.3, the proofs use
no computation; the other items are safeguards, such as independent recounts of the dimension formulas and checks
over finite fields. The arithmetic is exact. [`README.txt`](README.txt) gives the requirements, the commands, the script
behind each item, how the labels in the scripts translate into the numbering of the note, and the result of the
rerun made when this repository was assembled.

**Status:** not peer reviewed.

**Disclosure:** the arguments, the computations and the code were produced with the assistance of an AI system
(Claude, by Anthropic).

## Licences

- **Code** (`*.py`, `*.sh`): MIT License, see [`LICENSE`](LICENSE).
- **Everything else** (README files and computed outputs): Creative Commons Attribution 4.0 International
  (CC BY 4.0), see [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).
- No third-party files are included.
