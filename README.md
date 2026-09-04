# Five-point monodromy relations for AdS string integrals

At five points and weight 2 the monodromy relation expands into **25 scalar relations**
among **100 distinct Z-integrals**. This repository derives those relations symbolically
and verifies all 25 of them numerically.

At the kinematic point `{s12, s23, s34, s45, s15} = {1/5, 3/10, 2/5, 3/5, 7/20}`:

| | |
|---|---|
| relations checked | 25 / 25 |
| Z-integrals evaluated | 100 (no failures) |
| worst residual | `2.97e-8` |
| best residual | `1.81e-10` |
| max Z | `535.09` |
| worst residual relative to max Z | `5.5e-11` |

The residuals are quadrature error, not physics.

## Reproduce the check

Needs nothing but Mathematica — no NCAlgebra, no PolyLogTools, no HPL, no ginsh.

Open `numericalChecks5pt.nb` and evaluate it, or headless:

```bash
wolframscript -code '
  Get["zNumerics.m"]; monData = Get["monEqs5_W2.m"];
  sMon = {1/5, 3/10, 2/5, 3/5, 7/20};
  zList = Union@Cases[monData["equations"][[All,2]], _Zc, Infinity];
  zVals = Association[(# -> ZNum[Sequence @@ (List @@ #), sMon]["value"]) & /@ zList];
  Print[Max[Table[Abs[numS[monData["equations"][[k,2]] /. zVals, sMon]], {k, 25}]]]'
```

About 4-5 minutes: 100 integrals, each split into five sector-decomposed pieces. The
notebook adds an optional appendix that re-derives the convergence domain symbolically,
which takes a few minutes more.

## Files

| file | what it is | needs |
|---|---|---|
| `numericalChecks5pt.nb` | the numerical verification — the 100 integrals and the 25 residuals | nothing |
| `zNumerics.m` | all the numerical machinery, self-contained | nothing |
| `monEqs5_W2.m` | the derived relations: `<\|"Wmax", "rho", "equations", "integrals"\|>` | nothing |
| `monodromyRelations5pt.nb` | symbolic derivation — regenerates `monEqs5_W2.m` | NCAlgebra, PolyLogTools, HPL |
| `drinfeldAssociator.m` | the Drinfeld associator as a truncated series | NCAlgebra, PolyLogTools, HPL |

`monEqs5_W2.m` is committed, so the derivation notebook is optional: it is there to show
where the relations come from, not to be run before the check.

## The relations

`monodromyRelations5pt.nb` implements

```
Z(tau_1|rho) e^{i pi E_12} + Z(I_n|rho)
  + Sum_{j=3}^{n-1} Z(tau_j|rho) prodL_{k=3}^{j} C_{2->k} = 0
```

with chambers `tau_1 = (2 1 3 ... n)`, `tau_2 = I_n`, `tau_j = (1 3 ... j 2 j+1 ... n)`,
and crossing factors

```
C_{2->k} = Phi(e_2k, tauhat_k(e^<_{k-1})) e^{-i pi E_2k} Phi(tauhat_{k-1}(e^<_{k-1}), e_2k).
```

At `n = 5` the five independent braid generators are `x1 = e13, x2 = e23, x3 = e34,
y1 = e12, y2 = e24`; all braid relations reduce to the five commutators `[y_j, x_i]`, so
`{x-word} ** {y-word}` is a basis and every expression has a normal form. Expanding
the relation in that basis and truncating at total degree `Wmax` gives one scalar
equation per basis word: 6 at weight <= 1 and 19 more at weight 2.

`Wmax = 2` takes seconds, `Wmax = 3` about two minutes.

## How the numerics work

Two things make the check cheap enough to run:

**No ginsh at run time.** PolyLogTools evaluates MPLs by forking an external `ginsh`
process — about 10 ms per call, which inside `NIntegrate` (10^4–10^6 samples) means hours
per integral. But at `Wmax = 2` the chamber alphabets are `{0, Y, 1}` and `{0, 1}`, so
every MPL that appears has weight <= 2 and an elementary closed form. `zNumerics.m`
substitutes those directly: 9.8 ms -> 7.8 micros per integrand call, a factor of ~1250. The
closed forms were checked against ginsh over all weight <= 2 words at six random rational
points; the discrepancy was exactly zero.

**Sector decomposition at the hard corner.** `X = Y U` maps the simplex to the unit
square and factorizes three of the four corners. The remaining corner `U = Y = 1` goes
like `a^(s23-1) b^s24 (a+b)^(s34-1)` — integrable but not factorizable, and adaptive
quadrature stalls there. Splitting `a < b` from `b < a` makes each sector a product of
pure powers. `ZNum` therefore integrates five pieces per Z.

**The kinematic point matters.** The relation mixes the four chambers `12345`, `13245`,
`13425` and `21345`, so the point must lie in the convergence domain of all four at once:

```
s12, s23, s34, s45, s15 > 0        and        s13, s24, s25 > -1.
```

A point chosen for one chamber generally fails this. `{1/5, 3/10, 3/10, 6/5, 6/5}`
satisfies the first group but has `s25 = -11/5`, and the two chambers it misses return
`1e7` and `1e33`. `convExponents` in `zNumerics.m` reads the facet exponents off the
sector-decomposed integrand, so the domain is checked rather than assumed. At the point
used here the worst facet exponent in the four chambers is `-4/5`, `-7/10`, `-13/20` and
`-4/5` — all comfortably above `-1`. The appendix of `numericalChecks5pt.nb` computes
them; it takes a few minutes and nothing else depends on it.

## License

MIT — see `LICENSE`.
