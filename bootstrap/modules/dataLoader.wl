(* ::Package:: *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
(* :!CodeAnalysis::Disable::UppercaseVariable:: *)

BeginPackage["dataLoader`"]

(* A basis is an association that depends on a sector association defined like so
 *
 * sectors = <| "labels"  -> {...}, "index" -> AssociationThread[Range[n] -> labels],
 *              "weights" -> {...}, "S" -> ..., "T" -> ... |>
 * 
 * basis   = <| "id"    -> "...", "sectors" -> <ref>, "invariant" -> z,
 *              "slots" -> {{a,b}, ...}, "index" -> AssociationThread[slots -> Range[n]],
 *              "degeneracy" -> {...}, "degPos" -> ..., "nondegPos" -> ...|>
 *)

MinimalModularData::uasage = "Get the modular data of a minimal model";
FoldedMinimalModularData::usage = "Obtain parameters for folded minimal models";

Begin["`Private`"]

(* ADE modular invariants *)
kacLabels[p_, pp_] := Join[
    Flatten[Table[{r, s}, {r, 1, Floor[(pp - 1)/2]}, {s, 1, p - 1}], 1],
    If[EvenQ[pp], Table[{pp/2, s}, {s, 1, Floor[(p - 1)/2]}], {}]];

(* {b, ...} -> Sum_b (Sum_{a in b} e_a) (Sum_{a in b} e_a)^T *)
blocks[n_, bs_] := Total[Outer[Times, #, #] & /@
    Map[Table[Boole[MemberQ[#, a]], {a, n - 1}] &, bs]];

ade[n_, "A"] := IdentityMatrix[n - 1];
ade[n_, "D"] /; EvenQ[n] && n >= 6 := Table[
    If[Mod[n, 4] == 2,
        Boole[OddQ[a]] (Boole[a == b] + Boole[a + b == n]),   (* extension   *)
        If[OddQ[a], Boole[a == b], Boole[a + b == n]]],       (* automorphism *)
    {a, n - 1}, {b, n - 1}];
ade[12, "E"] := blocks[12, {{1, 7}, {4, 8}, {5, 11}}];
ade[18, "E"] := blocks[18, {{1, 17}, {5, 13}, {7, 11}, {3, 9, 15}}] - blocks[18, {{3, 15}}];
ade[30, "E"] := blocks[30, {{1, 11, 19, 29}, {7, 13, 17, 23}}];

MinimalModularInvariant[p_, pp_, invariant_] := Module[
    {
        labels = kacLabels[p, pp], mr, ms, lift, z
    },
        (* the D/E factor sits on the label with even Coxeter number *)
        {mr, ms} = If[EvenQ[pp],
            {ade[pp, invariant], IdentityMatrix[p - 1]},
            {IdentityMatrix[pp - 1], ade[p, invariant]}];

        (* sum over both lifts of each label, halving the (r,s) ~ (pp-r,p-s) double count *)
        lift[{r_, s_}] := {{r, s}, {pp - r, p - s}};
        z[{r_, s_}, {a_, b_}] := mr[[r, a]] ms[[s, b]];
        Table[Sum[z[x, y], {x, lift[f]}, {y, lift[g]}]/2, {f, labels}, {g, labels}]
    ]

MinimalModularData[p_:5, pp_:4, invariant_: "A"] := Module[
    {
        sectors = <| |>, c, h, s, t
    },
        c = 1 - 6 (p - pp)^2 / (p pp);

        (* Labels for each of the virasoro irreps *)
        sectors["labels"] = kacLabels[p, pp];
        sectors["index"] = AssociationThread[Range@Length[sectors["labels"]] -> sectors["labels"]];

        (* Conformal Weights *)
        h[r_, s_] := h[r,s] = ((p r - pp s)^2 - (p - pp)^2) / (4 p pp);
        sectors["weights "] = hv @@ # & /@ sectors["labels"];

        (* Modular Data *)
        s[a_, b_] := 2 Sqrt[2 / (p pp)] (-1)^(1 + a[[2]] b[[1]] + a[[1]] b[[2]]) Sin[a[[1]] b[[1]] Pi p / pp] Sin[a[[2]] b[[2]] Pi pp / p];
        t[a_, e_:1] := Exp[2 Pi I e (h[a[[1]], a[[2]]] - c/24)];
        sectors["S"] = Outer[s, sectors["labels"], sectors["labels"], 1];
        sectors["T"] = DiagonalMatrix[t /@ sectors["labels"]];

        (* The actual basis of ishibashi states with all the relevant data  *)
        <| 
            "id" -> StringJoin["M(", ToString[p], ", ", ToString[pp], ")"],
            "sectors" -> sectors,
            "invariant" -> MinimalModularInvariant[p, pp, invariant]
        |>
]

FoldedMinimalFixedPointModularData[p_:5, pp_:4] := Module[
    {
        cv, hv, sv, tv, H, i, ho, t, sMaker, s, nH, HIndex, svMat, tv2, twist
    },
        cv = 1 - 6 (p - pp)^2 / (p pp);
        hv[r_, s_] := hv[r,s] = ((p r - pp s));
        sv[a_, b_] := sv[a,b] = 2 Sqrt[2 / (p pp)] (-1)^(1 + a[[2]] b[[1]] + a[[1]] b[[2]]) Sin[a[[1]] b[[1]] Pi p / pp] Sin[a[[2]] b[[2]] Pi pp / p];
        tv[a_, e_:1] := tv[a,e] = Exp[2 Pi I e (hv[a[[1]], a[[2]]] - cv/24)];
    ]

End[]
EndPackage[]

(* :!CodeAnalysis::EndBlock:: *)
