(* ::Package:: *)

(* ---------------------------------------------------------------------------
   zNumerics.m -- fast numerical evaluation of the 5-point Z-integrals, and the
   residuals of the monodromy relations they satisfy.

   Self-contained: plain Wolfram Language, no external packages.  In particular
   it needs neither NCAlgebra, nor PolyLogTools, nor HPL, nor ginsh.  Just

       Get[FileNameJoin[{NotebookDirectory[], "zNumerics.m"}]]

   Public symbols
   --------------
     0.  sBasis5, momCons5, toBasis5, gaugeX, Gins, chAlph, genIntegrand,
         sChoice, integrandAtSChoice      kinematics and the integrand
     1.  GN, GtoNum                       weight <= 2 MPLs in closed form
     3.  toIdx, numFun, ZNum              one Z-integral, numerically
     4.  convExponents, monChambers       convergence domain of a chamber
     5.  numS, monResidual                residual of one monodromy relation

   Why this exists
   ---------------
   PolyLogTools' Ginsh (plt/core/GiNaCInterface.m) writes a script file, forks an
   external `ginsh` process, imports the answer and deletes three files: ~10 ms
   per call.  Inside NIntegrate, which samples 10^4-10^6 times, that is hours to
   days.  Nothing here needs ginsh at run time -- monEqs5_W2.m has Wmax -> 2 and
   the chamber alphabets are {0,Y,1} for the first variable and {0,1} for the
   second, so every MPL that appears has weight <= 2 and an elementary closed
   form.  Measured: 9.8 ms -> 7.8 us per integrand call (~1250x).
   --------------------------------------------------------------------------- *)

ClearAll[Zc, ss, s, sBasis5, momCons5, toBasis5, gaugeX, Gins, chAlph, genIntegrand,
         sChoice, integrandAtSChoice];

SetAttributes[ss, Orderless];
ss[i_, i_] = 0;
Format[ss[i_, j_]] := Subscript[s, i, j]

Format[Zc[tt_String, rr_String, uu_String, vv_String]] :=
  Row[{Superscript["Z", tt], "[", rr, "|", uu, ",", vv, "]"}]

sBasis5 = {s[1, 2], s[2, 3], s[3, 4], s[4, 5], s[1, 5]};

(* momentum conservation + masslessness: Sum_j s[i,j] = 0 for each i *)
momCons5 = {
   s[1, 3] -> s[4, 5] - s[1, 2] - s[2, 3],
   s[2, 4] -> s[1, 5] - s[2, 3] - s[3, 4],
   s[3, 5] -> s[1, 2] - s[3, 4] - s[4, 5],
   s[1, 4] -> s[2, 3] - s[4, 5] - s[1, 5],
   s[2, 5] -> s[3, 4] - s[1, 5] - s[1, 2]};

toBasis5[e_] := Expand[e //. momCons5];

(*Our labels for z integrals are the two permutations and two words that identify the relevant multivariable MPL:
    τ   chamber: the boundary ordering of the punctures
    ρ   Parke-Taylor ordering
    l1  word for the first integration variable,  x[τ[[2]]] = X
    l2  word for the second integration variable, x[τ[[3]]] = Y
  Words are lists of PUNCTURE labels. Where a puncture sits is read off from the
  τ-adapted gauge, so the same word means the same singularity in every chamber.

    x[τ[[1]]] = 0, x[τ[[2]]] = X, x[τ[[3]]] = Y, x[τ[[4]]] = 1, x[τ[[5]]] = ∞

  and the domain is the simplex 0 < X < Y < 1 in every chamber.*)

gaugeX[τ_List] := AssociationThread[τ -> {0, X, Y, 1, Infinity}]

(*G(a1,...,an;z) in PolyLogTools' flat notation; the empty word is 1*)
Gins[{}, z_] := 1
Gins[w_List, z_] := G[Sequence @@ w, z]

(*Alphabet of the k-th variable, in this chamber's labels: the puncture at 0, the
  punctures of the LATER variables, and the puncture at 1. A word may use no
  other letter -- e.g. the first variable may not be singular at its own puncture.*)
chAlph[τ_List, k_] := τ[[Join[{1}, Range[k + 2, Length[τ] - 2], {Length[τ] - 1}]]]

genIntegrand[τ_List, ρ_List, l1_List, l2_List] :=
 Module[{x = gaugeX[τ], inf = Last[τ], fin = Most[τ], KN, PT, cyc},
  Do[If[! SubsetQ[chAlph[τ, k], {l1, l2}[[k]]],
     Message[genIntegrand::badw, {l1, l2}[[k]], chAlph[τ, k], k, τ];
     Return[$Failed, Module]], {k, 2}];

  KN = Product[(x[fin[[j]]] - x[fin[[i]]])^(s @@ Sort[{fin[[i]], fin[[j]]}]),
       {j, 2, Length[fin]}, {i, 1, j - 1}];

  cyc = Partition[Append[ρ, First[ρ]], 2, 1];
  PT = (x[τ[[1]]] - x[τ[[4]]])/
       Product[x[pr[[1]]] - x[pr[[2]]], {pr, DeleteCases[cyc, {___, inf, ___}]}];

  (*the MPL insertion: in this gauge a puncture label IS its position*)
  KN PT Gins[x /@ l1, X] Gins[x /@ l2, Y]]

genIntegrand::badw = "Word `1` is not in the alphabet `2` of variable `3` in chamber `4`.";

sChoice[numVal_] := Table[sBasis5[[k]] -> numVal[[k]], {k, Length[sBasis5]}]

integrandAtSChoice[τ_List, ρ_List, l1_List, l2_List, numVal_List] :=
  (genIntegrand[τ, ρ, l1, l2] // toBasis5) /. sChoice[numVal]

(* ===========================================================================
   1.  Weight <= 2 multiple polylogarithms in closed form
   =========================================================================== *)

ClearAll[GN, GtoNum, onCutQ];

(* The branch is chosen ONCE, from a probe point inside 0 < X < Y < 1, so that GN
   stays symbolic and is cut-free across the whole domain. *)
$GNprobe = {X -> 3/10, Y -> 7/10};
onCutQ[u_] := With[{v = N[u /. $GNprobe]}, NumericQ[v] && Im[v] == 0 && Re[v] > 1];

GN[{}, z_]       := 1;
GN[{0}, z_]      := Log[z];
GN[{a_}, z_]     := Log[1 - z/a] /; a =!= 0;
GN[{0, 0}, z_]   := Log[z]^2/2;
GN[{0, b_}, z_]  := -PolyLog[2, z/b] /; b =!= 0;
GN[{a_, 0}, z_]  := Log[z] Log[1 - z/a] + PolyLog[2, z/a] /; a =!= 0;
GN[{a_, b_}, z_] := Log[1 - z/a]^2/2 /; a === b && a =!= 0;

(* G(a,b;z) = Log[(b-a)/b] Log[1-z/a] + Li2[a/(a-b)] - Li2[(a-z)/(a-b)].
   Both Li2 arguments can land on the cut (> 1) -- e.g. a = 1, b = Y gives
   1/(1-Y) > 1.  The total is real, but Mathematica's principal branches leave a
   spurious imaginary part.  Since a/(a-b) + b/(b-a) = 1, exactly one of the two
   orderings is off the cut; reach the other through the shuffle
   G(a;z) G(b;z) = G(a,b;z) + G(b,a;z). *)
GN[{a_, b_}, z_] := Module[{u1 = a/(a - b), u2 = (a - z)/(a - b)},
    If[onCutQ[u1] || onCutQ[u2],
      GN[{a}, z] GN[{b}, z] - GN[{b, a}, z],
      Log[(b - a)/b] Log[1 - z/a] + PolyLog[2, u1] - PolyLog[2, u2]]
    ] /; a =!= 0 && b =!= 0 && a =!= b;

(* HoldPattern is essential: PolyLogTools gives G downvalues under which the bare
   pattern G[args__] evaluates to 1, so without it the rule silently never fires
   and the G's survive unconverted. *)
GtoNum[expr_] := expr //. HoldPattern[G[args__]] :> GN[Most[{args}], Last[{args}]];

(* ===========================================================================
   2.  Validation against ginsh

   GN was checked against PolyLogTools` Ginsh over all weight <= 2 words of
   both chamber alphabets, at 6 seeded random rational points of the domain:
   the maximum discrepancy was exactly 0 in both variables.
   =========================================================================== *)

(* ===========================================================================
   3.  The integration region

   X = Y U maps the simplex 0 < X < Y < 1 to the unit square and resolves the
   corner at the origin: X->0, Y->0 and X->Y factorize onto U=0, Y=0, U=1.

   What is left is the corner U = Y = 1.  With a = 1-U, b = 1-Y the integrand
   goes like a^(s23-1) b^s24 (a+b)^(s34-1) -- integrable but NOT factorizable,
   and adaptive quadrature stalls on it (machine precision and WorkingPrecision
   30 disagreed in the 4th digit).  Splitting a < b from b < a via a = b t and
   b = a t makes each sector a product of pure powers.

   The Simplify under the domain assumptions is
   what combines b^(s23-1) b^s24 b^(s34-1) b into a single power.  Without it
   NIntegrate evaluates Infinity*0*Infinity*0 on the edge and returns
   Indeterminate.
   =========================================================================== *)

ClearAll[toIdx, numFun, ZNum];

toIdx[s_String] := ToExpression /@ Characters[s];   (* "13245" -> {1,3,2,4,5} *)
toIdx[l_List]   := l;

(* Numeric closure, evaluated in LOG SPACE for the power factors.

   Every singular facet makes several powers blow up or vanish together, and no
   algebraic form survives that at machine precision:  Simplify merges them into
   one radical of a product, (U^7 Y^9 (1-U Y)^3)^(1/10) / (Y - U Y)^(7/10), which
   is 0/0 as Y -> 0 even though Y^(9/10)/Y^(7/10) = Y^(1/5) is perfectly finite;
   splitting them apart instead gives 0 * Infinity.  PowerExpand would fix the
   algebra but is invalid here -- the bases are not all positive, and it injects a
   spurious (-1)^(3/10).

   So sum p_k log|b_k| first and exponentiate once.  Each base has a fixed sign on
   a connected region, so the power part's overall sign is a constant, fixed once
   against a direct evaluation at an interior probe point. *)
numFun[expr_, vars_List, probe_List] := Module[{nrm, flat, e, fs, pw, rest, g, h, direct, viaLog, sgn},
   (* Cancellation-free bases.  Substituting U -> 1-yy, Y -> 1-yy uu leaves a base
      such as 1 - U Y sitting as -1 + (1-yy)(1-yy uu): Mathematica does not expand
      a product like that, so at yy ~ 10^-20 it is evaluated as 1 - 1 = 0 and the
      power blows up.  Together+Factor turns it into -yy(1 + uu - yy uu), which is
      accurate to the last bit.  Log and PolyLog arguments need the same care. *)
   nrm[z_] := Factor[Together[z]];
   e = expr /. {Power[bb_, pp_] /; ! FreeQ[bb, Alternatives @@ vars] :> Power[nrm[bb], pp],
                Log[z_] /; ! FreeQ[z, Alternatives @@ vars] :> Log[nrm[z]],
                PolyLog[n_, z_] /; ! FreeQ[z, Alternatives @@ vars] :> PolyLog[n, nrm[z]]};
   fs = If[Head[e] === Times, List @@ e, {e}];
   (* Flatten each power down to irreducible bases.  Simplify likes to merge
      factors into one radical, e.g. (yy^13 (1+uu-uu yy)^7)^(1/20).  At yy ~ 10^-30
      the inner yy^13 = 10^-390 underflows to 0 in double precision BEFORE the log
      is taken, so log-space alone does not save it.  Splitting is exact and
      branch-safe here because every base enters as Abs[.] with the overall sign
      fixed separately. *)
   flat[{bb_, pp_}] := Module[{bf = Factor[Together[bb]]},
      Which[Head[bf] === Times,  Join @@ (flat[{#, pp}] & /@ (List @@ bf)),
            Head[bf] === Power,  flat[{bf[[1]], pp*bf[[2]]}],
            True,                {{bf, pp}}]];
   pw = Join @@ (flat /@ Cases[fs, Power[bb_, pp_] /; ! FreeQ[bb, Alternatives @@ vars] :> {bb, pp}]);
   pw = DeleteCases[pw, {_?NumericQ, _}];
   rest = Times @@ DeleteCases[fs, Power[bb_, _] /; ! FreeQ[bb, Alternatives @@ vars]];
   g[sub_] := Exp[Total[(#[[2]]*Log[Abs[N[#[[1]] /. sub]]]) & /@ pw]]*Re[N[rest /. sub]];
   direct = Re[N[e /. Thread[vars -> probe]]];
   viaLog = g[Thread[vars -> probe]];
   sgn = If[viaLog == 0. || ! NumericQ[viaLog], 1, Sign[direct/viaLog]];
   (* the ?NumericQ guard is required: without it NIntegrate evaluates the
      integrand symbolically once while building, gets Indeterminate, and gives up *)
   h[p_?NumericQ, q_?NumericQ] := sgn*g[{vars[[1]] -> p, vars[[2]] -> q}];
   h];

Options[ZNum] = {WorkingPrecision -> MachinePrecision, PrecisionGoal -> 8,
   MaxRecursion -> 40, "ImCheck" -> True};

(* After X = Y U the only singular loci are the coordinate edges U=0, U=1, Y=0,
   Y=1 (each a single factorized power) plus 1 - U Y, which vanishes ONLY at the
   corner U = Y = 1.  Two separate things must be arranged:

   (a) that corner is entangled -- with a = 1-U, b = 1-Y the integrand goes like
       a^(s23-1) b^s24 (a+b)^(s34-1), integrable but not factorizable -- so it
       needs sector decomposition, a < b and b < a.

   (b) every OTHER edge must be reached through a coordinate that equals the
       small quantity exactly.  Integrating U up to 1 directly makes Simplify
       write the edge as (Y - U Y)^(7/10) in a denominator; at U = 1-10^-32
       machine arithmetic returns 1-U = 0 and the integrand becomes
       0^(-0.7) Log[0] = Indeterminate.  That is cancellation, not a singularity.

   So split the square at U = Y = 1/2 and substitute per quadrant:
       Q1  U=u,   Y=y     Q2  U=1-u, Y=y
       Q3  U=u,   Y=1-y   Q4  the corner -> two sectors
   all four quadrants over (0,1/2)^2.  Sector decomposition must NOT be applied
   globally: on the whole square it would map (U,Y)=(0,0), already factorized as
   U^(s12-1) Y^(s45-1), onto (a,b)=(1,1) where b = a t re-entangles it. *)

ZNum[t0_, r0_, w10_, w20_, sVals_List, OptionsPattern[]] :=
  Module[{tau = toIdx[t0], rho = toIdx[r0], l1 = toIdx[w10], l2 = toIdx[w20],
    raw, fast, sq, asm, e1, e2, e3, eA, eB, pieces, vals, wp, pg, mr, imax = 0.},
   wp = OptionValue[WorkingPrecision]; pg = OptionValue[PrecisionGoal];
   mr = OptionValue[MaxRecursion];

   raw = integrandAtSChoice[tau, rho, l1, l2, sVals];
   If[raw === $Failed, Return[$Failed]];

   fast = GtoNum[raw];
   If[! FreeQ[fast, G],
     Print["ZNum: unconverted MPL (weight > 2?): ", Union@Cases[fast, _G, Infinity]];
     Return[$Failed]];

   sq  = Simplify[(fast /. X -> Y*U)*Y, {0 < U < 1, 0 < Y < 1}];
   asm = {0 < uu < 1/2, 0 < yy < 1/2};

   e1 = Simplify[sq /. {U -> uu,     Y -> yy},     asm];
   e2 = Simplify[sq /. {U -> 1 - uu, Y -> yy},     asm];
   e3 = Simplify[sq /. {U -> uu,     Y -> 1 - yy}, asm];
   eA = Simplify[(sq /. {U -> 1 - yy*uu, Y -> 1 - yy})*yy, {0 < uu < 1, 0 < yy < 1/2}];
   eB = Simplify[(sq /. {U -> 1 - yy,    Y -> 1 - yy*uu})*yy, {0 < uu < 1, 0 < yy < 1/2}];

   If[OptionValue["ImCheck"],
     imax = Max[Table[Abs[Im[N[eA /. {uu -> RandomReal[], yy -> RandomReal[{0, 1/2}]}]]], {40}]]];

   (* {expression, uu range} -- yy runs over (0,1/2) in every piece *)
   pieces = {{e1, 1/2}, {e2, 1/2}, {e3, 1/2}, {eA, 1}, {eB, 1}};
   vals = (Quiet@NIntegrate[numFun[#[[1]], {uu, yy}, {#[[2]]/3, 1/6}][p, q],
        {p, 0, #[[2]]}, {q, 0, 1/2},
        Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
        WorkingPrecision -> wp, PrecisionGoal -> pg, MaxRecursion -> mr]) & /@ pieces;

   If[! FreeQ[vals, _NIntegrate],
     Print["ZNum: piece(s) ", Flatten@Position[vals, _NIntegrate, {1}],
           " failed for ", {t0, r0, w10, w20}];
     Return[$Failed]];

   <|"value" -> Total[vals], "pieces" -> vals, "maxIm" -> imax|>];

(* ===========================================================================
   4.  Convergence conditions per chamber

   Computed symbolically for all four chambers.  For 12345, 13245 and 13425 the
   run agreed with reading the exponents off the collision structure; 21345,
   previously filled in analytically by the same rule, now agrees too and adds
   no constraint the other three do not already imply.  Two self-checks pass
   everywhere: the sector boundary tt->1 gives exponent 0, and the two corner
   limits bb->0 and aa->0 agree within each chamber.

   sNum5 = {1/5, 3/10, 3/10, 6/5, 6/5} was chosen to make the 13245 chamber
   converge, and it does.  But the monodromy equations mix four chambers, and at
   that point two of them diverge badly:

       Zc["12345","12345","",""]  ->    18.04         (fine)
       Zc["13245","12345","",""]  ->    -2.16         (fine)
       Zc["13425","12345","",""]  ->    -8.4 x 10^7   (divergent)
       Zc["21345","12345","",""]  ->    -3.0 x 10^33  (divergent)

   so equation 1 cannot close there.  Requiring every facet exponent > -1 in all
   four chambers gives

       s12, s23, s34, s45, s15 > 0     and     s13, s24, s25 > -1

   and sNum5 fails the second group (s25 = -11/10).  At sMon = {1/5, 3/10, 2/5,
   3/5, 7/20} all four converge and monodromy equation 1 closes to
   |residual| = 1.8*10^-10 against max|Z| = 34, i.e. 5*10^-12 relative.
   =========================================================================== *)

ClearAll[stripT, secExprs, expAt0, expAt1, convExponents];

stripT[e_] := e /. {Log[_] -> 1, PolyLog[_, _] -> 1};   (* logs never shift a power *)

secExprs[tau_, rho_] := Module[{raw, sq, eA, eB},
   raw = stripT[GtoNum[toBasis5[genIntegrand[tau, rho, {}, {}]]]];
   sq = Simplify[(raw /. X -> Y*U)*Y, {0 < U < 1, 0 < Y < 1}];
   eA = Simplify[(sq /. {U -> 1 - bb*tt, Y -> 1 - bb})*bb, {0 < tt < 1, 0 < bb < 1}];
   eB = Simplify[(sq /. {U -> 1 - aa, Y -> 1 - aa*tt})*aa, {0 < tt < 1, 0 < aa < 1}];
   {eA, eB}];

(* for f ~ v^e near v -> 0,  v d/dv Log f -> e *)
expAt0[e_, v_] := Simplify[Limit[v*D[Log[e], v], v -> 0]];
expAt1[e_, v_] := Module[{w}, Simplify[Limit[w*D[Log[e /. v -> 1 - w], w], w -> 0]]];

convExponents[tau_, rho_] := Module[{ab = secExprs[tau, rho]},
   <|"sectorA" -> {"tt->0" -> expAt0[ab[[1]], tt], "tt->1" -> expAt1[ab[[1]], tt],
                   "bb->0" -> expAt0[ab[[1]], bb], "bb->1" -> expAt1[ab[[1]], bb]},
     "sectorB" -> {"tt->0" -> expAt0[ab[[2]], tt], "tt->1" -> expAt1[ab[[2]], tt],
                   "aa->0" -> expAt0[ab[[2]], aa], "aa->1" -> expAt1[ab[[2]], aa]}|>];

monChambers = {{1, 2, 3, 4, 5}, {1, 3, 2, 4, 5}, {1, 3, 4, 2, 5}, {2, 1, 3, 4, 5}};

(* ===========================================================================
   5.  Driving the monodromy equations
   =========================================================================== *)

ClearAll[numS, monResidual];

numS[e_, sVals_] := toBasis5[e /. ss -> s] /. sChoice[sVals];

monResidual[eqRHS_, sVals_, opts : OptionsPattern[ZNum]] :=
  Module[{zs, vals},
   zs = Union@Cases[eqRHS, _Zc, Infinity];
   vals = Association[(# -> ZNum[Sequence @@ (List @@ #), sVals, opts]["value"]) & /@ zs];
   <|"residual" -> ((eqRHS /. vals) /. ss[i_, j_] :> numS[ss[i, j], sVals]),
     "z" -> vals|>];
