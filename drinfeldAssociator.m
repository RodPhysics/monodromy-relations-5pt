(* ::Package:: *)

(* ---------------------------------------------------------------------------
   drinfeldAssociator.m -- the Drinfeld associator as a truncated series in two
   non-commuting letters.

   Requires NCAlgebra for the non-commutative products, PolyLogTools for G and
   HPL for HPL, all three of which the notebook loads before this file.

   Convention, matching the draft:

       Phi(e0, e1) = 1 + Sum_w e_w G_reg(w; 1),

   the sum running over all words w of length 1..weight in the two letters, with
   the multiple polylogarithms regularised at z = 1 (the divergent HPL pieces are
   dropped).  The series stops at weight 4, where RegHPLat1's expansion stops --
   hence Wmax <= 4 everywhere downstream.
   --------------------------------------------------------------------------- *)

Needs["NCAlgebra`"];

(* all NC words of length 2..maxLength in {a,b}, preceded by the two letters *)
GenerateNCFactor[a_, b_, maxLength_] := Module[{words},
  words = Flatten[Table[Tuples[{a, b}, n], {n, 2, maxLength}], 1];
  Flatten[Join[{a}, {b}, NonCommutativeMultiply @@ # & /@ words]]]

(* the same words as 0/1 sequences, i.e. as MPL indices *)
GenerateNCWords[maxLength_] := Flatten[Table[Tuples[{0, 1}, n], {n, 1, maxLength}], 1]

(* G(w; z) regularised at z -> 1 *)
RegHPLat1[w__] := Normal[Series[G[w, x], {x, 1, 4}]] /. x -> 1 /. HPL[a__] -> 0

Drinfeld[a_, b_, weight_] := 1 + Sum[
   GenerateNCFactor[a, b, weight][[k]] (RegHPLat1 @@ GenerateNCWords[weight][[k]]),
   {k, 1, Length[GenerateNCWords[weight]]}]
